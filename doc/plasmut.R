## ----setup, message=FALSE-----------------------------------------------------
library(plasmut)

library(tidyverse)
library(readxl)

## ----load data----------------------------------------------------------------
extdir <- system.file("extdata", package="plasmut")
fpath <- file.path(extdir, "cairo5-matched-sequencing.csv")

df <- read_csv(fpath, show_col_types=FALSE) %>% type.convert(as.is=TRUE)

df

## ----compute bayes factors----------------------------------------------------
bayesfactors <- bayes.factor(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads")) 

## ----compute probability of tumor specific alteration-------------------------
posteriorodds <- posterior.odds(bayesfactors)
prob <- posterior.probability(posteriorodds)

## ----stratify wbc variants----------------------------------------------------
wbc.var.type <- stratify_wbc_vars(data=df, cols=c("cfDNA distinct mutant reads","cfDNA distinct reads","wbc distinct mutant reads","wbc distinct reads"))

## ----data cleaning and organization-------------------------------------------
df <- df %>% mutate(prob.tumor.specific=prob, wbc.var.type=wbc.var.type, bayes.factor=bayesfactors) %>% mutate(wbcmaf=`wbc distinct mutant reads`/`wbc distinct reads`, cfdnamaf=`cfDNA distinct mutant reads`/`cfDNA distinct reads`) %>% mutate(cfdnamaf=ifelse(is.nan(cfdnamaf), 0, cfdnamaf))

df <- df %>% filter(g_include == "yes")

df %>% relocate(Gene, Patient, prob.tumor.specific, bayes.factor, cfdnamaf, wbcmaf) %>% head()

## ----estimation accuracy------------------------------------------------------
df %>% select(cfdnamaf, wbcmaf, prob.tumor.specific, type) %>% filter(type=="Tumor-confirmed") %>% mutate(prob.tumor.specific = round(prob.tumor.specific, 3)) %>% group_by(prob.tumor.specific) %>% summarize(category=unique(type), n=n())

## ----visualization, fig.align="center", fig.width=9, fig.height=9, warning=FALSE----
cluster.annotate <- df %>% select(wbc.var.type, cfdnamaf, wbcmaf) %>% group_by(wbc.var.type) %>% summarize(wbclevel=mean(wbcmaf)) %>% arrange(desc(wbclevel)) %>% filter(!is.na(wbc.var.type))

df <- df %>% mutate(status = wbc.var.type) %>% mutate(status = ifelse(is.na(status), "Likely tumor-specific", wbc.var.type)) %>% mutate(status = ifelse(status==cluster.annotate$wbc.var.type[1], "Germline variants", status)) %>% mutate(status = ifelse(status==cluster.annotate$wbc.var.type[2], "Hematopoietic variants", status)) %>% mutate(prob=status) %>% mutate(prob=ifelse(status=="Likely tumor-specific", prob.tumor.specific, prob)) %>% mutate(prob=ifelse(status=="Germline variants", 2, prob), prob=ifelse(status=="Hematopoietic variants", 0, prob)) %>% mutate(prob=as.numeric(prob))

df <- df %>% mutate(shapestatus=ifelse(type=="Tumor-confirmed", "Tumor confirmed", "Tumor not sequenced"))

cfvar_plot <- ggplot(df, aes(x=wbcmaf, y=cfdnamaf)) +
    geom_jitter(aes(color = prob, fill = status, shape=shapestatus),
                fill = "black", stat = "identity", lwd = 3, width = 0.1, alpha = 0.7) +
    scale_color_gradientn(colors = c("red3", "orange", "blue", "black", "black", "orange"), 
                          na.value = "transparent",
                          breaks=c(0.6, 0.8, 1, 2, 0),
                          labels=c(0.6, 0.8, 1,"Germline variants", "Hematopoietic variants"),
                          limits=c(0 , 3)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),
                       breaks=c(0.000084, 0.001, 0.01,  0.1, 0.5),
                       trans='log10',
                       oob = scales::squish_infinite) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       breaks=c(0.000084, 0.001, 0.01,  0.1, 1),
                       trans='log10',
                       oob = scales::squish_infinite) +
    scale_shape_manual(values = c(17, 16)) + 
    labs(x="White blood cell MAF (%)",
          y="cfDNA MAF (%)")+
    theme_bw() + 
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.tag = element_text()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.text.x = element_text(vjust=1.2, hjust= 0.3)) +
    theme(legend.direction = "vertical", legend.box = "horizontal",
          legend.background = element_rect(fill = "transparent"),
          legend.position = c(0.47, 0.88),
          legend.key = element_rect(size = 5), #legend.box.background = element_rect(colour = "black"),
          legend.key.size = unit(0.4, "cm")) +
    guides(color=guide_legend(title="Probability\nTumor-Specific", nrow = 3),
           fill = guide_legend(title=""),
           shape = guide_legend(title="")) +
    geom_vline(xintercept=c(0.00015), linetype="dashed", size = 0.2) +
    geom_vline(xintercept=c(0.00009), linetype="dashed", size = 0.2, color = NA) +
    geom_vline(xintercept=c(0.25), linetype="dashed", size = 0.2) +
    geom_abline(slope=1, intercept=0, color="gray", alpha = 0.4, size = 0.4)

cfvar_plot + theme(axis.title=element_text(size=16), axis.text=element_text(size=13), legend.text=element_text(size=12), legend.title=element_text(size=14))

## ----lower probability tumor specific mutations-------------------------------
df %>% filter(wbcmaf == 0) %>% arrange(prob) %>% select(prob, `cfDNA distinct mutant reads`, `cfDNA distinct reads`, `wbc distinct mutant reads`, `wbc distinct reads`, type, Gene, Patient) %>% slice(1:10)

