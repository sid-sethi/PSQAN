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
  print(str_c(Sys.time(), " - ", sqanti_class))
  cat("Input Sqanti file has 0 rows\nExiting...............\n\n")
  quit()
}


n_samples <- x.gene %>%
  dplyr:::select(starts_with("FL")) %>% 
  ncol()

x.gene$Isoform_class <- forcats::fct_recode(x.gene$Isoform_class, "Coding known\n(complete match)" = "Coding known (complete match)", "Coding known\n(alternate 3'/5' end)" = "Coding known (alternate 3/5 end)")



#### Transcripts (ranked) vs. NFLR ######

if(n_samples > 1){
  
  ## with error bars #####
  x.ranked <- x.gene %>%
    dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
    pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
    group_by(isoform, Isoform_class) %>%
    summarise(NFLR_mean = mean(NFLR), NFLR_sd = sd(NFLR), .groups = "keep") %>%
    arrange(desc(NFLR_mean)) %>%
    tibble::rowid_to_column(., "isoform_index")
  
  
  p1 <- ggplot(data=x.ranked, aes(x=isoform_index, fill = Isoform_class)) +
    geom_bar(aes(y = NFLR_mean), color=NA, size=0.3, width=0.8, stat="identity") +
    geom_errorbar(aes(ymin = NFLR_mean-NFLR_sd, ymax = NFLR_mean+NFLR_sd), width = 0.2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(name = "Transcript expression\nacross samples; NFLR") +
    scale_x_continuous(name = "Transcripts ranked by expression (NFLR)") +
    guides(fill=guide_legend(title="Transcript category")) +
    mytheme +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      legend.position = c(0.80, 0.80)
    )
  
  png(str_c(outpath, "/", prefix, "_transcriptsRanked_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()
  
  
}else{
  
  x.gene$NFLR_mean <- x.gene$NFLR
  
}




#### without error bars ########

x.ranked <- x.gene %>%
  arrange(desc(NFLR_mean)) %>%
  tibble::rowid_to_column(., "isoform_index")
  #mutate(isoform_index = 1:nrow(x.gene))


png(str_c(outpath, "/", prefix, "_transcriptsRanked_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

ggplot(data=x.ranked, aes(x=isoform_index, fill = Isoform_class)) +
  geom_bar(aes(y = NFLR_mean), color=NA, size=0.3, width=0.8, stat="identity") +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(name = "Transcript expression\nacross samples; NFLR") +
  scale_x_continuous(name = "Transcripts ranked by expression (NFLR)") +
  guides(fill=guide_legend(title="Transcript category")) +
  mytheme +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size=0.4),
    legend.position = c(0.80, 0.80)
  )
dev.off()
