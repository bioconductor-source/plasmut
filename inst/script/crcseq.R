library(tidyverse)
library(magrittr)
extdir <- system.file("extdata", package="plasmut")
fpath <- file.path(extdir, "cairo5-matched-sequencing.csv")
df <- read_csv(fpath, show_col_types=FALSE) %>%
    type.convert(as.is=TRUE)
cfdna <- df %>%
    select(Gene, Patient, `cfDNA distinct mutant reads`,
           `cfDNA distinct reads`,
           `AAChange`,
           `Nucleotide Position (hg19)`) %>%
    set_colnames(c("gene", "patient", "y", "n",
                   "aa", "position")) %>%
    mutate(analyte="plasma")
buffy <- df %>%
    select(Gene, Patient, `wbc distinct reads`,
           `wbc distinct mutant reads`,
           `AAChange`,
           `Nucleotide Position (hg19)`) %>%
    set_colnames(c("gene", "patient", "n", "y",
                   "aa", "position")) %>%
    mutate(analyte="buffy coat")
crc.data <- bind_rows(cfdna, buffy) %>%
    arrange(patient, aa)
p12 <- crc.data %>%
    filter(patient == 12,
           gene %in% c("APC", "HRAS"))
p13 <- crc.data %>%
    filter(patient == 13, gene=="FGFR2")
p157 <- crc.data %>%
    filter(patient == 157, aa=="M237K")
dat <- bind_rows(p12, p13, p157) %>%
    arrange(patient, gene) %>%
    select(patient, gene, aa, position, analyte, y, n)
crcseq <- dat

output_dir <- "../../data"
if (!dir.exists(output_dir)){
    dir.create(output_dir)
}
save(crcseq, file=file.path(output_dir, "crcseq.rda"))
