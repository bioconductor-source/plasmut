#' Strip a column of a dataframe into a vector
#' @import tidyverse
#' @param data a data frame.
#' @param col.idx the column name or column number
#' @return a vector of the selected column
strip.column <- function(data, col.idx=1){

    x <- data %>% as.list()

    return(x[[col.idx]])
}

#' Compute the likelihood of a given plasma mutation being somatic in origin
#'
#' @param y.p # of distinct mutant reads in plasma
#' @param n.p # of distinct reads in plasma
#' @param y.w # of distinct mutant reads in wbc
#' @param n.w # of distinct reads in wbc
#' @param thetas samples from a beta distribution
#' @param p.theta density of beta computed at thetas (i.e. dbeta(thetas, a, b))
#' @return the likelihood of the mutation being somatic in origin
lik.somatic <- function(y.p, n.p, y.w, n.w, thetas, p.theta){

    mean(dbinom(y.p, size=n.p, thetas) * p.theta) * dbinom(y.w, size=n.w, prob=0)
}

#' Compute the likelihood of a given plasma mutation being CHIP in origin
#' 
#' @param y.p # of distinct mutant reads in plasma
#' @param n.p # of distinct reads in plasma
#' @param y.w # of distinct mutant reads in wbc
#' @param n.w # of distinct reads in wbc
#' @param thetas samples from a beta distribution
#' @param p.theta density of beta computed at thetas (i.e. dbeta(thetas, a, b))
#' @return the likelihood of the mutation being hematopoietic in origin
lik.chip <- function(y.p, n.p, y.w, n.w, thetas, p.theta){

    mean(dbinom(y.p, size=n.p, thetas) * dbinom(y.w, size=n.w, thetas) * p.theta)
}


#' Compute the likelihood of a given plasma mutation being somatic in origin
#'
#' @param pdmr vector of distinct mutant reads in plasma for different mutations
#' @param pdr vector of distinct reads in plasma for different mutations
#' @param wdmr vector of distinct mutant reads in wbc for different mutations
#' @param wdr vector of distinct reads in wbc for different mutations
#' @param data data frame with the pdmr, pdr, wdmr, and wdr information as column and rows correspond to different mutations
#' @param cols vector of column names to use (cols[1] refers to pdmr, etc.)
#' @param samples # of samples for Monte Carlo integration
#' @param prior.odds P(mutation being somatic) / P(mutation being hematopoietic)
#' @return the posterior odds of the mutation being somatic
#'@export
posterior.odds <- function(pdmr, pdr, wdmr, wdr, data=NULL, cols=NULL, samples=1e6, prior.odds=1){

    if (!is.null(data)){

        if (!is.null(cols)){

            pdmr <- strip.column(data, cols[1])
            pdr <- strip.column(data, cols[2])
            wdmr <- strip.column(data, cols[3])
            wdr <- strip.column(data, cols[4])

        } else{
            pdmr <- strip.column(data, 1)
            pdr <- strip.column(data, 2)
            wdmr <- strip.column(data, 3)
            wdr <- strip.column(data, 4)
        }
    }

    a <- 0.5
    b <- 0.5
    theta <- rbeta(samples, a, b)
    p.theta <- dbeta(theta, a, b)

    lik.SOMATIC <- unlist(lapply(seq_along(pdmr), function(i) lik.somatic(pdmr[i], pdr[i], wdmr[i], wdr[i], theta, p.theta)) )

    lik.CHIP <- unlist(lapply(seq_along(pdmr), function(i) lik.chip(pdmr[i], pdr[i], wdmr[i], wdr[i], theta, p.theta)))

    bayes.factor <- lik.SOMATIC / lik.CHIP

    posterior.odds <-  bayes.factor * prior.odds

    posterior.odds <- ifelse(!is.finite(posterior.odds), 1e100, posterior.odds)

    return(posterior.odds)

}

#' Posterior Probability of mutation being somatic (tumor specific) computed from Posterior Odds
#' @param po a vector of posterior odds
#' @return a vector of posterior probabilities
#' @examples
#' posterior.probability(1)
#' @export
posterior.probability <- function(po){

    return ( po / (1 + po) )

}



