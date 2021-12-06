args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})


sqanti_class <- args[1] # sqanti _classification_processed
prefix <- args[2]
outpath <- args[3]
NE_thresh <- args[4]
NFLR_thresh <- args[5]
min_sample_perc <- args[6]


x.gene = read.table(sqanti_class, header=TRUE, sep="\t")

########### filtering ###########

x.gene$total_samples <- x.gene %>% dplyr::select(starts_with("NFLR.")) %>% ncol()

pass = x.gene %>% dplyr::select(starts_with("NFLR.")) %>% .[] > NFLR_thresh

x.pass <- x.gene %>%
  mutate(
    n_pass = rowSums(pass),
    perc_pass = round((n_pass/total_samples)*100,2)
  ) %>%
  dplyr::filter(NE >= NE_thresh, perc_pass >= min_sample_perc)


write.table(x.pass, str_c(outpath, "/", prefix, "_classification_filtered.txt"), col.names = TRUE, row.names=FALSE, sep="\t", quote = FALSE)
