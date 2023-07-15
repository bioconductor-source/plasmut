importance_samples <- function(y, n,
                               a, b,
                               S,
                               prior.weight,
                               model,
                               dna_source){
    a.n <- a + y
    b.n <- n - y + b
    x <- rbinom(S, 1, prior.weight)
    g.samples <- rep(NA, S)
    g.samples[x==1] <- rbeta(sum(x), a, b)
    g.samples[x==0] <- rbeta(S-sum(x), a.n, b.n)
    ## probability of CTC
    posterior.samples <- rbeta(S, a.n, b.n)
    prior.samples <- rbeta(S, a, b)
    ##
    ## probability of theta with respect to prior and posterior
    ## p(theta | model )
    ##
    g.prob.posterior <- dbeta(g.samples, a.n, b.n)
    g.prob.prior <- dbeta(g.samples, a, b)
    ## the density of g is the mixture density
    lik <- dbinom(y, n, g.samples)
    ## probability of theta with respect to g
    g.density <- prior.weight * g.prob.prior +
        (1-prior.weight) * g.prob.posterior
    result <- lik * g.prob.prior / g.density
    ml <- mean(result)
    samples <- tibble(posterior=posterior.samples,
                      prior=prior.samples,
                      g=g.samples,
                      model=model,
                      dna_source=dna_source)
    densities <- tibble(posterior=g.prob.posterior,
                        prior=g.prob.prior,
                        g=g.density,
                        x=x,
                        lik=lik,
                        model=model,
                        dna_source=dna_source)
    result <- list(samples=samples,
                   densities=densities,
                   marglik=ml)
    result
}

#' Estimate the marginal likelihood of observing somatic mutations from CTCs present in buffy coat
#'
#'
#' p(y_w | theta_w, n_w, model_S) x p(theta_w| Model_S)
#' theta_w | model_S  ~ beta(1, 999)  ## sequencing error or CTC
#' @param dat tibble containing elements `y`and `n`.  The `y` and `n` vectors should be named
#' @param params a list with named elements that must include the following:
#' `gamma`: weight given to heavier-tailed prior distribution
#' `a`: prior expectation for number of somatic variants observed in the WBC sequencing data (either by error or from a CTC)
#' `b`: prior expectation for number of WBCs not containing the variant
#'
#' @examples
#' ## asssume we observed y_w=500 and n_w = 1000
#' ## this is clearly germline, but to go through the machinery we would obtain the probability of observing 500 someatic
#' ## variants (y_w) under a model where the prior for theta_w has most of the mass near zero
#' ##    - the importance sampler does not really make sense here as it samples from a distribution that is fatter tailed than the
#' ##      posterior -> it
#'
wbc_somatic <- function(dat, params){
    S <- params$montecarlo.samples
    prior.weight <- params$prior.weight
    params <- params$ctc
    dat <- filter(dat, analyte=="buffy coat")
    y <- dat$y
    n <- dat$n
    a <- params$a
    b <- params$b
    result <- importance_samples(y, n,
                                 a, b,
                                 S,
                                 prior.weight,
                                 model="Somatic",
                                 dna_source="CTC")
}

#' Estimate the marginal likelihood that variants identified in cell-free DNA are derived from tumor cells (ctDNA-derived)
#'
plasma_somatic <- function(dat, params){
    prior.weight <- params$prior.weight
    S <- params$montecarlo.samples
    params <- params$ctdna
    dat <- filter(dat, analyte=="plasma")
    y <- dat$y
    n <- dat$n
    a <- params$a
    b <- params$b
    result <- importance_samples(y, n,
                                 a, b,
                                 S,
                                 prior.weight,
                                 model="Somatic",
                                 dna_source="ctDNA")
    return(result)
}

#' Estimate the marginal likelihood that mutations in buffy coat and cfDNA reflect CH or correspond to germline mutations.
#'
#' If germline, the allele frequency should be 50%.  The prior should be diffuse enough to handle CHIP mutations (potentially << 50%).
#'
model_w <- function(dat, params){
    ##    - We use g to propose values of theta and evaluate the likelihood of the observed data in plasma and WBCs given this value
    ##    - If we use a beta (1, 10) for example, this should be diffuse enough to capture germline
    prior.weight <- params$prior.weight
    S <- params$montecarlo.samples
    params <- params$chip
    y <- setNames(dat$y, dat$analyte)
    n <- setNames(dat$n, dat$analyte)
    a <- params$a
    b <- params$b
    ##
    ## With single theta and conjugate (beta) prior, posterior is
    ## beta(a + y_w + y_p, b + (n.w + n.p)-(y.w+y.p))
    ##
    a.n <- a + sum(y)
    b.n <- sum(n)- sum(y) + b
    g.samples <- rep(NA, S)
    x <- rbinom(S, 1, prior.weight)
    g.samples[x == 1] <- rbeta(sum(x), a, b)
    g.samples[x == 0] <- rbeta(S-sum(x), a.n, b.n)
    prior.samples <- rbeta(S, a, b)
    posterior.samples <- rbeta(S, a.n, b.n)
    ##
    ##
    ##
    g.prob.posterior <- dbeta(g.samples, a.n, b.n)
    g.prob.prior <- dbeta(g.samples, a, b)
    g.density <- prior.weight * g.prob.prior +
        (1-prior.weight)*g.prob.posterior
    not.finite <- !is.finite(g.density)
    if(any(not.finite)) browser("Not finite")
    loglik.p <- dbinom(y["plasma"], n["plasma"], g.samples, log=TRUE)
    loglik.w <- dbinom(y["buffy coat"], n["buffy coat"], g.samples, log=TRUE)
    lik <- exp(loglik.p + loglik.w)
    result <- lik * g.prob.prior/g.density
    ml <- mean(result)
    samples <- tibble(posterior=posterior.samples,
                      prior=prior.samples,
                      g=g.samples,
                      model="WBC",
                      dna_source="WBC")
    densities <- tibble(posterior=g.prob.posterior,
                        prior=g.prob.prior,
                        g=g.density,
                        x=x,
                        lik=lik,
                        model="WBC",
                        dna_source="WBC")
    result <- list(samples=samples,
                   densities=densities,
                   marglik=ml)
    result
}

#' Importance sampler to estimate marginal likelihoods and Bayes factors
#'
#' @examples
#' param.list <- list(ctc=list(a=1, b=9999),
#'                    ctdna=list(a=1, b=9),
#'                    chip=list(a=1, b=9),
#'                    montecarlo.samples=50e3,
#'                    prior.weight=0.1)
#' dat <- tibble(y=c(4, 1),
#'               n=c(1000, 1000),
#'               analyte=c("plasma", "buffy coat"),
#'               mutation="mutA",
#'               sample_id="id1")
#' importance_sampler(dat, param.list)
#' @export
importance_sampler <- function(dat, params, save_montecarlo=TRUE){
    ##
    ## Probability we see a somatic alteration when
    ## sequencing WBCs
    A <- wbc_somatic(dat, params)
    ## Probability we see the somatic alteration in
    ## cfDNA
    B <- plasma_somatic(dat, params)
    ml.tumor <- A$marglik * B$marglik
    ##ns.params <- nonsomatic_params(1, 10, gamma=0.2)
    model.not.tumor <- model_w(dat, params)
    ml.not.tumor <- model.not.tumor$marglik
    bf <- ml.tumor/ml.not.tumor
    stats <- tibble(ctc=A$marglik,
                    ctdna=B$marglik,
                    chip=model.not.tumor$marglik,
                    bayesfactor=bf)
    if(!save_montecarlo) return(stats)
    densities <- bind_rows(A$densities,
                           B$densities,
                           model.not.tumor$densities)
    samples <- bind_rows(A$samples,
                         B$samples,
                         model.not.tumor$samples)
    list(data=dat,
         densities=densities,
         samples=samples,
         bayesfactor=stats,
         params=params)
}
