## plasmut: Weighing the evidence of tumor-specific mutations in cell-free DNA

Adith S. Arun and Robert B. Scharpf  

### Overview

Mutation-based approaches for detection of cancer from cell-free DNA (cfDNA) using liquid biopsies have the potential to track a patient’s response to treatment, enabling effective and timely decisions on therapy. However, mutations arising from clonal hematopoeisis (CH) are common and tumor biopsies for definitive identification of the origin of these mutations is not always available. Sequencing of matched cells from buffy coat and the absence of mutations in these cells has been used as a test to rule-out CH, but uneven sequencing depths between matched cell-free DNA and buffy coat samples and the potential for contamination of buffy coat with circulating tumor cells (CTCs) are not captured by rule-based analyses. This package estimates Bayes factors that weigh the evidence of competing CH- and tumor-origin models of cfDNA mutations detected in cfDNA, requiring only the allele frequencies of high quality alignments available from standard mutation callers.

### Installation

``` r
## install.packages("devtools")
devtools::install_github("cancer-genomics/plasmut")
```
This package is also available through Bioconductor. 


### Usage

Please see the vignette in `doc/plasmut.html`. 

We illustrate this approach on a dataset of cfDNA and matched buffy coat sequencing for patients with metastatic colorectal cancer [van ’t Erve et al. 2023)](https://pubmed.ncbi.nlm.nih.gov/36534496/). Below, we select four mutations and run the importance sampler for these candidate mutations independently.

``` r
params <- list(ctdna = list(a = 1, b = 9), 
               ctc = list(a = 1, b = 10^3), 
               chip = list(a= 1, b = 9), 
               montecarlo.samples = 50e3, 
               prior.weight = 0.1)
muts <- unite(crcseq, "uid", c(patient, gene), remove = FALSE) %>% 
        group_by(uid) %>% nest()
#Each element of the data column contains a table with the variant and total allele counts in plasma and buffy coat. 
#Run the importance sampler
muts$IS <- muts$data %>% map(importance_sampler, params)
fun <- function(x){
    result <- x$data %>%
        select(-position) %>%
        mutate(bayes_factor = x$IS$bayesfactor$bayesfactor)
    return(result)
}
bf <- apply(muts, 1, fun) 
bf %>% do.call(rbind, .) %>%
    as_tibble() %>%
    select(patient, gene, aa, bayes_factor) %>%
    rename(log_bf=bayes_factor) %>%
    distinct()
```
