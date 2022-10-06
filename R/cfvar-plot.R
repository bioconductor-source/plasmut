#' visualize probability of mutation being somatic against maf of wbc and plasma
#' @import tidyverse
#' @import scales
#' @param df dataframe with somatic probabilities computed by mutchips, wbc var type categorized by mutchips, and wbcmaf and cfdna maf computed as specificed in vignette
#' @return a ggplot2 plot that can be further customized
#' @export

geom_cfvar <- function(df){

  cluster.annotate <- df %>% select(wbc.var.type, cfdnamaf, wbcmaf) %>% group_by(wbc.var.type) %>% summarize(wbclevel=mean(wbcmaf)) %>% arrange(desc(wbclevel)) %>% filter(!is.na(wbc.var.type))

  df <- df %>% mutate(status = wbc.var.type) %>% mutate(status = ifelse(is.na(status), "Likely tumor-specific", wbc.var.type)) %>% mutate(status = ifelse(status==cluster.annotate$wbc.var.type[1], "Germline variants", status)) %>% mutate(status = ifelse(status==cluster.annotate$wbc.var.type[2], "Hematopoietic variants", status)) %>% mutate(prob=status) %>% mutate(prob=ifelse(status=="Likely tumor-specific", prob.tumor.specific, prob)) %>% mutate(prob=ifelse(status=="Germline variants", 2, prob), prob=ifelse(status=="Hematopoietic variants", 0, prob)) %>% mutate(prob=as.numeric(prob))

  cfvar_plot <- ggplot(df, aes(x=wbcmaf, y=cfdnamaf)) +
    geom_jitter(aes(color = prob, fill = status),
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
           fill = guide_legend(title="")) +
    geom_vline(xintercept=c(0.00015), linetype="dashed", size = 0.2) +
    geom_vline(xintercept=c(0.00009), linetype="dashed", size = 0.2, color = NA) +
    geom_vline(xintercept=c(0.25), linetype="dashed", size = 0.2) +
    geom_abline(slope=1, intercept=0, color="gray", alpha = 0.4, size = 0.4)

  return(cfvar_plot)

}
