## plasmut: stratifying mutations observed in cell-free DNA and estimating the probability of tumor-specific alterations

Adith S. Arun, Robert B. Scharpf
Oct 6, 2022

### Overview

The package accomplishes three things: a) Bayesian estimate of the probability that a mutation is tumor-specific, specifically when zero white blood cell mutant reads are observed, b) classification of white blood cell mutants as germline or hematopoietic, and c) visualization of mutant classes on cell-free DNA MAF and white blood cell DNA MAF logarithmic axes. 

For each mutation, we require plasma total distinct mutant reads, plasma total distinct reads, white blood cell distinct mutant reads, and white blood cell distinct reads.

### Installation

``` r
## install.packages("devtools")
devtools::install_github("cancer-genomics/plasmut")
```


### Usage

Please see the vignette in `doc/plasmut.html`
