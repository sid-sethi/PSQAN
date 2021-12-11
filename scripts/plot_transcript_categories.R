args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
})


sqanti_class <- args[1] # sqanti _classification_processed or _classification_filtered
prefix <- args[2]
outpath <- args[3]

# ggplot plotting theme parameters
mytheme <- theme_bw()  +
  theme(plot.title = element_text(lineheight=.4, size=12, hjust=0.5),
        legend.key.size = unit(0.40,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour="BLACK",size=0.4),
        #panel.border = element_blank(),
        #axis.line = element_line(color = "black", size=0.4),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12,angle = 90),
        panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )


x.gene = read.table(sqanti_class, header=TRUE, sep="\t")

if(nrow(x.gene) == 0){
  cat("Input Sqanti file has 0 rows\nExiting...............\n\n")
  quit()
}

x.gene$Isoform_class <- forcats::fct_recode(x.gene$Isoform_class, "Coding known\n(complete match)" = "Coding known (complete match)", "Coding known\n(alternate 3'/5' end)" = "Coding known (alternate 3/5 end)")


n_samples <- x.gene %>%
  dplyr:::select(starts_with("FL")) %>% 
  ncol()

########### plots ###########

if(n_samples > 1){
  
  ## Number of unique transcripts vs. transcript category ##
  
  tc.count <- x.gene %>%
    dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
    pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
    dplyr::filter(NFLR >0) %>%
    group_by(Isoform_class, Sample) %>%
    dplyr::summarise(count = n(), .groups = "keep") %>%
    group_by(Isoform_class) %>%
    summarise(count_mean = mean(count), count_sd = sd(count), .groups = "keep")
  
  
  p1 <- ggplot(data=tc.count, aes(x=Isoform_class, y = count_mean, fill=Isoform_class)) +
    geom_bar(color="black", size=0.3, width=0.7, stat="identity", alpha=0.8) +
    geom_errorbar(aes(ymin = count_mean-count_sd, ymax = count_mean+count_sd), width = 0.2) +
    scale_fill_brewer(palette = "Dark2") +
    guides(fill="none") +
    mytheme +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      axis.text.x = element_text(size=10, angle=35, hjust=1)
    ) +
    scale_x_discrete(name = "Transcript category") +
    scale_y_continuous(name = "Number of unique transcripts")
  
  
  png(str_c(outpath, "/", prefix, "_tc_count_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()
  
  
  ## expression of transcript categories ##
  
  ## per sample ##
  
  tc.nflr <- x.gene %>%
    dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
    pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
    dplyr::filter(NFLR >0)
  
  
  p2 <- ggplot(data=tc.nflr, aes(x=Isoform_class, y = log2(NFLR), fill=Isoform_class)) +
    geom_boxplot(lwd = 0.55, alpha=0.8, outlier.shape = NA) +
    geom_jitter(position=position_jitter(0.22), shape=21, size=1, alpha=0.5) +
    scale_fill_brewer(palette = "Dark2") +
    guides(fill="none") +
    mytheme +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      axis.text.x = element_text(size=10, angle=35, hjust=1)
    ) +
    scale_x_discrete(name = "Transcript category") +
    scale_y_continuous(name = "Transcript expression\nlog2(NFLR)")
  
  
  png(str_c(outpath, "/", prefix, "_tc_NFLR_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p2)
  dev.off()
  
  
  
}else{
  
  x.gene$NFLR_mean <- x.gene$NFLR
  
}



png(str_c(outpath, "/", prefix, "_tc_count_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

ggplot(data=x.gene, aes(x=Isoform_class, fill=Isoform_class)) +
  geom_bar(color="black", size=0.3, width=0.7, stat="count", alpha=0.8) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill="none") +
  mytheme +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size=0.4),
    axis.text.x = element_text(size=10, angle=35, hjust=1)
  ) +
  scale_x_discrete(name = "Transcript category") +
  scale_y_continuous(name = "Number of unique transcripts")

dev.off()




## expression of transcript categories ##

## per transcript ##


png(str_c(outpath, "/", prefix, "_tc_NFLR_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

ggplot(data=x.gene, aes(x=Isoform_class, y = log2(NFLR_mean), fill=Isoform_class)) +
  geom_boxplot(lwd = 0.55, alpha=0.8, outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.22), shape=21, size=1, alpha=0.5) +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill="none") +
  mytheme +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size=0.4),
    axis.text.x = element_text(size=10, angle=35, hjust=1)
  ) +
  scale_x_discrete(name = "Transcript category") +
  scale_y_continuous(name = "Transcript expression\nacross samples; log2(NFLR)")

dev.off()
