#' Compute the likelihood of a given plasma mutation being somatic in origin
#' @import tidyverse
#' @import rjags
#' @import ggmcmc
#' @param data is a list with 4 named items (np = total distinct reads plasma, yp = mutant reads plasma, nw = total distinct reads wbc, yw = mutant reads wbc) 
#' @param rho_ctc_parameter fraction of cells in buffy coat layer expected to be circulating tumor cells, default 0. If CTC suspected, 0.001 is a good starting point. 
#' @param a_w_somatic shape beta parameter for informing wbc maf in somatic CTC model
#' @param b_w_somatic rate beta parameter for informing wbc maf in somatic CTC model
#' @param a_p_somatic shape beta parameter for informating plasma maf in somatic CTC and No CTC model
#' @param b_p_somatic rate beta parameter for informating plasma maf in somatic CTC and No CTC model
#' @param a_theta shape beta parameter for informing plasma and wbc maf in hematopoietic model
#' @param b_theta rate beta parameter for informing plasma and wbc maf in hematopoietic model
#' @param prior.odds prior odds of a mutation being somatic P(somatic) / P(hematopoietic), default to 1
#' @param chains number of monte carlo chains, default 2
#' @param adapt number of iterations for burn in of each monte carlo chain, default 2000
#' @param iter number of iterations for each chain, default 1e5
#' @param nthin thinning value, take every nth value of chain, default 10
#' @param samples # of samples for Monte Carlo integration
#' @param prior.odds P(mutation being somatic) / P(mutation being hematopoietic)
#' @return list with two named vectors - bf (bayes factor) and p (probbability of tumor-specific) for all mutations supplied into function
#' @examples 
#' compute_p_tumor(list(np=100, yp=5, nw=100, yw=10))
#'@export

compute_p_tumor <- function(data, rho_ctc_parameter=0, a_w_somatic=1, b_w_somatic=10, a_p_somatic=1, b_p_somatic=10, a_theta=1, b_theta=10, prior.odds=1, chains=2, adapt=2000, iter=1e5, nthin=10){
  
  data$a_p <- a_p_somatic
  data$b_p <- b_p_somatic
  data$a_theta <- a_theta
  data$b_theta <- b_theta
  data$a_w <- a_w_somatic
  data$b_w <- b_w_somatic
  data$rho <- rho_ctc_parameter
  
  if (data$rho == 0){
    somatic_model <- system.file(file.path("JAGS", "somtic-model-no-tail.jag"), package = "plasmut", mustWork = TRUE)
  } else{
    somatic_model <- system.file(file.path("JAGS", "somatic-model.jag"), package = "plasmut", mustWork = TRUE)
  }
  
  hemato_model <- system.file(file.path("JAGS", "hemato-model.jag"), package="plasmut", mustWork = TRUE)

  #number of mutations to estimate probabilities for
  data$N <- length(data$np)

  fit.s <- rjags::jags.model(somatic_model, data=data, n.chains=chains, n.adapt=adapt)

  samples.s <- rjags::coda.samples(fit.s,
  	                         variable.names="out",
  	                          n.iter=iter,
  	                          thin=nthin)

  fit.h <- rjags::jags.model(hemato_model, data=data, n.chains=chains, n.adapt=adapt)

  samples.h <- rjags::coda.samples(fit.h,
  	                         variable.names="out",
  	                          n.iter=iter,
  	                          thin=nthin)

  somatic <- ggmcmc::ggs(samples.s) %>% filter(grepl("out", Parameter)) %>% group_by(Chain, Parameter) %>% summarize(v=mean(value)) %>% ungroup() %>% group_by(Parameter) %>% summarize(lik=mean(v))


  hemato <- ggmcmc::ggs(samples.h) %>% filter(grepl("out", Parameter)) %>% group_by(Chain, Parameter) %>% summarize(v=mean(value)) %>% ungroup() %>% group_by(Parameter) %>% summarize(lik=mean(v))

  bf <- somatic$lik / hemato$lik

  posterior.odds <- bf * prior.odds

  posterior.odds <- unlist(lapply(posterior.odds, function(x) ifelse(!is.finite(x), 1e100, x)))

  p <- posterior.odds / (1 + posterior.odds)

  return(list(bf=bf, p=p))
}

