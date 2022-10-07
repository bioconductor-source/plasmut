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

Please see the vignette in `doc/plasmut.html` or the below documentation.  

An example workflow of the plasmut package is shown using data from the CAIRO5 study (van't Erve and Medina et. al 2022). For a group of patients with cancer, both cell-free DNA and white blood cell DNA were sequenced. The distinct and total reads for each mutation present in cell-free DNA are displayed in the data. Corresponding white blood cell total and distinct reads are displayed as well. 

From this information, we wish to assign labels to each of the mutations. Specifically, for mutations seen in cell-free DNA and not in white blood cells (distinct reads wbc = 0), we estimate the probability that the mutation arises from the tumor. This is done in a Bayesian framework and approximates integrals as Monte Carlo sums. Please see the paper for a workthrough of the math underlying this. Additionally, for mutations seen in white blood cells and cell-free DNA, we classify them as either hematopoietic or germline variants in an unsupervised manner.

``` r
library(plasmut)
library(tidyverse)
library(readxl)
```

Read in the data and convert columns to their natural types. 
``` r
extdir <- system.file("extdata", package="plasmut")
fpath <- file.path(extdir, "cairo5-matched-sequencing.csv")
df <- read_csv(fpath) %>% type.convert(as.is=TRUE)
```

Call the plasmut function posterior.odds() and posterior.probability() that computes the probability of a mutation being tumor derived. 
``` r
prob.tumor.specific <- posterior.odds(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads")) %>% posterior.probability()
```

Call the plasmut function stratify_wbc_vars that identifies which white blood cell mutants are germline versus hematopoietic. 
``` r
wbc.var.type <- stratify_wbc_vars(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads"))
```

Compute mutant allele fractions (MAFs) and clean data.
``` r
df <- df %>% mutate(prob.tumor.specific=prob.tumor.specific, wbc.var.type=wbc.var.type) %>% mutate(wbcmaf=`wbc distinct mutant reads`/`wbc distinct reads`, cfdnamaf=`cfDNA distinct mutant reads`/`cfDNA distinct reads`) %>% mutate(cfdnamaf=ifelse(is.nan(cfdnamaf), 0, cfdnamaf))
df <- df %>% filter(g_include == "yes")
```

We estimate the probability of a mutation being somatic (tumor-derived) using only information from cell-free DNA and white blood cells with posterior.odds() and posterior.probability(). We can assess the accuracy of this method by sequencing the tumor and observing whether the mutation is found there. We sequenced 42 cases where the tumor displayed the same variant as the cell-free DNA. Our estimation procedure assigned a probability of tumor of 1 to 41 of these cases and a probability of 0.987 to the other, which was consistent with the tumor sequencing results. 

``` r
df %>% select(cfdnamaf, wbcmaf, prob.tumor.specific, type) %>% filter(type=="Tumor-confirmed") %>% mutate(prob.tumor.specific = round(prob.tumor.specific, 3)) %>% group_by(prob.tumor.specific) %>% summarize(category=unique(type), n=n())
```

Finally, we can visualize the results of this approach. The final object here is a ggplot object that can be further customized. 
``` r
cfvar_plot <- geom_cfvar(df)
cfvar_plot + theme(axis.title=element_text(size=16), axis.text=element_text(size=13), legend.text=element_text(size=12), legend.title=element_text(size=14))
```

