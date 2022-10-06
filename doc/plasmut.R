## ----setup, message=FALSE-----------------------------------------------------
library(plasmut)

library(tidyverse)
library(readxl)

## ----load data----------------------------------------------------------------
extdir <- system.file("extdata", package="plasmut")
fpath <- file.path(extdir, "cairo5-matched-sequencing.csv")

df <- read_csv(fpath) %>% type.convert(as.is=TRUE)

## ----compute probability of tumor-derived mutation----------------------------
prob.tumor.specific <- posterior.odds(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads")) %>% posterior.probability()

## ----stratify wbc variants----------------------------------------------------
wbc.var.type <- stratify_wbc_vars(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads"))

## ----data cleaning and organization-------------------------------------------
df <- df %>% mutate(prob.tumor.specific=prob.tumor.specific, wbc.var.type=wbc.var.type) %>% mutate(wbcmaf=`wbc distinct mutant reads`/`wbc distinct reads`, cfdnamaf=`cfDNA distinct mutant reads`/`cfDNA distinct reads`) %>% mutate(cfdnamaf=ifelse(is.nan(cfdnamaf), 0, cfdnamaf))

df <- df %>% filter(g_include == "yes")

## ----estimation accuracy------------------------------------------------------
df %>% select(cfdnamaf, wbcmaf, prob.tumor.specific, type) %>% filter(type=="Tumor-confirmed") %>% mutate(prob.tumor.specific = round(prob.tumor.specific, 3)) %>% group_by(prob.tumor.specific) %>% summarize(category=unique(type), n=n())

## ----visualization, fig.align="center", fig.width=7, fig.height=7-------------
cfvar_plot <- geom_cfvar(df)

cfvar_plot

