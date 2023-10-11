#' This data is an example dataset to show how to use the package
#'
#'
#' A cohort of metastatic colorectal cancer patients whose plasma and buffy 
#' coat were sequenced as part of the CAIRO5 trial. 
#' The cohort and analyses are described here: 
#' \url{https://pubmed.ncbi.nlm.nih.gov/36534496/}
#'
#'
#'
#' @docType data
#' @name crcseq
#' @return An example DNA sequencing dataset of matched plasma 
#' and wbc colorectal cancer samples
#' crcseq
NULL



#' This 
#'
#'
#' Processed matched sequencing output from the CAIRO5 clinical trial.
#' The data in this sheet is presented in this paper: 
#' The cohort and analyses are described here: 
#' \url{https://pubmed.ncbi.nlm.nih.gov/36534496/}
#' 
#' The data in this csv file (inst/extdata/cairo5-matched-sequencing.csv)
#' is processed according to inst/script/crcseq/R to create crcseq.rda 
#' that is fed through plasmut to estimate the probability
#' of tumor specific alterations
#'
#' @docType data
#' @name cairo5-matched-sequencing
#' @return Number of distinct and mutant reads in white blood cells
#' and cell-free DNA tissues for each patient's mutations (each row)
#' cairo5-matched-sequencing
NULL











