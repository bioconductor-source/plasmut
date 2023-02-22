## ----setup, message=FALSE-----------------------------------------------------
library(plasmut)
library(tidyverse)
library(rjags)
library(ggmcmc)

## ----load data----------------------------------------------------------------
extdir <- system.file("extdata", package="plasmut")
fpath <- file.path(extdir, "cairo5-matched-sequencing.csv")

df <- read_csv(fpath, show_col_types=FALSE) %>% type.convert(as.is=TRUE)

## ----sample columns of data frame to show for paper---------------------------
df %>% select(Gene, Patient, `cfDNA distinct mutant reads`, `cfDNA distinct reads`, `wbc distinct reads`, `wbc distinct mutant reads`) %>% slice(1:4)

## ----get relevant data--------------------------------------------------------
data <- df %>% select(8:11) %>% magrittr::set_colnames(c("np", "yp", "nw", "yw")) %>% mutate(yw=ifelse(is.na(yw), 0, yw), nw = ifelse(is.na(nw), 1e8, nw)) %>% as.list()

## ----compute prob tumor, warning=FALSE----------------------------------------
model_results <- compute_p_tumor(data)

## ----stratify wbc variants----------------------------------------------------
cols <- c("cfDNA distinct mutant reads",
          "cfDNA distinct reads",
          "wbc distinct mutant reads",
          "wbc distinct reads")

wbc.var.type <- stratify_wbc_vars(data=df,
                                  cols=cols)

## ----data cleaning and organization-------------------------------------------
df <- df %>%
    mutate(prob.tumor.specific=model_results$p, bayes.factor=model_results$bf, 
           wbc.var.type=wbc.var.type) %>%
    mutate(wbcmaf=`wbc distinct mutant reads`/`wbc distinct reads`,
           cfdnamaf=`cfDNA distinct mutant reads`/`cfDNA distinct reads`) %>%
    mutate(cfdnamaf=ifelse(is.nan(cfdnamaf), 0, cfdnamaf)) 

df <- df %>% filter(g_include == "yes")

df %>% relocate(Gene, Patient, prob.tumor.specific, bayes.factor, cfdnamaf, wbcmaf) %>% head()

## ----estimation accuracy------------------------------------------------------
df %>% select(cfdnamaf, wbcmaf, prob.tumor.specific, type) %>%
    filter(type=="Tumor-confirmed") %>%
    ## @ again, this is an odds and not a probability
    mutate(prob.tumor.specific = round(prob.tumor.specific, 3)) %>%
    group_by(prob.tumor.specific) %>%
    summarize(category=unique(type), n=n())

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

## ----bayes factor relationship to cfDNA, fig.align="center", fig.width=5, fig.height=5, warning=FALSE, message=FALSE, echo=FALSE, results='hide'----
cor.val <- df %>% select(prob.tumor.specific, bayes.factor, cfdnamaf, wbcmaf, Gene, Patient) %>% mutate(logbf=log(bayes.factor)) %>% filter(!is.infinite(logbf)) %>% select(cfdnamaf, logbf) %>% cor(method="spearman") %>% as_tibble() %>% slice(1) %>% pull(logbf) %>% round(2)

corr <- paste0("Corr: ", cor.val)

plt <- df %>% select(prob.tumor.specific, bayes.factor, cfdnamaf, wbcmaf, Gene, Patient) %>% mutate(logbf=log(bayes.factor)) %>% filter(!is.infinite(logbf)) %>% ggplot(aes(x=logbf, y=cfdnamaf)) + geom_point(size=1.5) + theme_bw() + ylab("cfDNA MAF") + xlab("log(Bayes Factor)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + geom_smooth(method="loess", se = F) 

#+ geom_text(x=400, y=0.22, label=corr, size=4)

plt

## ----tile plot, fig.height=6, fig.width=10, fig.align="center"----------------
plt.df <- df %>% select(status, Gene, Patient, bayes.factor) %>% mutate(border=ifelse(bayes.factor > 18.01, "> 18.01", "Else")) %>% select(-bayes.factor) %>% mutate(alt.source=ifelse(status=="Likely tumor-specific", "Tumor", status)) %>% select(-status)

ordered.levels <- plt.df %>% group_by(Gene) %>% summarize(n=n()) %>% arrange(desc(n)) %>% pull(Gene)

gene.numbers <- plt.df %>% select(Gene) %>% distinct() %>% mutate(gene.num=1:n())

plt.df <- plt.df %>% left_join(gene.numbers, by=c("Gene"="Gene"))

fig <- plt.df %>% ggplot(aes(y=factor(Gene, levels=ordered.levels), x=factor(Patient), fill=alt.source, color=border)) + geom_tile(size=1.25, alpha=0.7, width=0.85, height=0.75) + theme_minimal() + ylab("Gene") + xlab("Patient") + theme(panel.grid.major=element_blank()) + theme(axis.text.x=element_blank()) + scale_color_manual(values=c("firebrick", "darkblue"), name="Bayes Factor") + scale_fill_manual(values=c("mediumorchid", "cornflowerblue", "sienna"), name="Alteration\nSource")

fig

