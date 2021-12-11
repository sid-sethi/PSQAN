args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
})


sqanti_class <- args[1] # sqanti _classification_processed
prefix <- args[2]
outpath <- args[3]


x.gene = read.table(sqanti_class, header=TRUE, sep="\t")


n_samples <- x.gene %>%
  dplyr:::select(starts_with("FL")) %>% 
  ncol()


counts.long <- x.gene %>%
  dplyr::select(isoform, starts_with("NFLR")) %>%
  pivot_longer(cols = c(-isoform), names_to = "sample", values_to = "read_count")

# if multiple samples, remove primer tags
counts.long$sample <- counts.long$sample %>% str_replace("NFLR.*\\.", "") %>% str_replace("_3p", "")



################# No. of unique transcripts detected at each NFLR threshold ###############

if(n_samples > 1){
  threshs <- x.gene$NFLR_mean %>% round(2) %>% sort() %>% unique()
}else{
  threshs <- x.gene$NFLR %>% round(2) %>% sort() %>% unique()
}



res <- data.frame()

for(i in threshs){

  # calculating number of unique transcripts detected
  transcripts_per_sample <- counts.long %>%
    dplyr::filter(read_count >=i) %>%
    group_by(sample) %>%
    summarise(n_transcripts = n(), .groups = "keep") %>%
    mutate(NFLR_threshold = i)

  res <- rbind(res, transcripts_per_sample)
  rm(transcripts_per_sample)

}


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



if(n_samples > 1){
  
  
  # plotting: transcripts detected vs. mean NFLR (mean of NFLR across all samples)
  
  df1 <- res %>%
    dplyr::filter(sample %in% "NFLR_mean")
  
  sd <- res %>%
    dplyr::filter(!sample %in% "NFLR_mean") %>%
    group_by(NFLR_threshold) %>%
    summarise(sd = sd(n_transcripts), .groups = "keep")
  
  df1 <- left_join(df1, sd, by = "NFLR_threshold")
  write.table(df1, str_c(outpath, "/", prefix, "_NFLR_curve_main.txt"), quote=FALSE, sep="\t", row.names=FALSE)
  

  
  p1 <- ggplot(data=df1, aes(x=factor(NFLR_threshold), y = n_transcripts)) +
    geom_line(group=1) +
    geom_ribbon(aes(factor(NFLR_threshold), ymax = n_transcripts + sd, ymin = n_transcripts - sd), alpha = 0.2, fill = "blue", group=1) +
    scale_x_discrete(name = expression(paste("Normalised expression (", NFLR[T], ")"))) +
    scale_y_continuous(name = "No. of transcripts detected") +
    mytheme +
    guides(fill="none") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      axis.text.x = element_text(size=5, angle=90, hjust=1, vjust=0.5)
    )
    
    
  png(str_c(outpath, "/", prefix, "_NFLR_curve_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()
  
  
  # plotting: transcripts detected vs. NFLR per sample
  
  df2 <- res %>%
    dplyr::filter(!sample %in% "NFLR_mean")
  
  
  p2 <- ggplot(data=df2, aes(x=factor(NFLR_threshold), y = n_transcripts, group=sample)) +
    geom_line(aes(color=sample), alpha=0.5) +
    scale_x_discrete(name = expression(paste("Normalised expression (", NFLR[i], ")"))) +
    scale_y_continuous(name = "No. of transcripts detected") +
    mytheme +
    guides(color=guide_legend(title="Sample")) +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      legend.title = element_text(size=8),
      legend.text = element_text(size = 7),
      axis.text.x = element_text(size=5, angle=90, hjust=1, vjust=0.5)
    )
  
  png(str_c(outpath, "/", prefix, "_NFLR_curve_perSample.png"), bg="transparent",units="in",width = 5.25, height= 3.75 ,res=600)
  plot(p2)
  dev.off()
  
  
}else{
  
  # plotting: transcripts detected vs. NFLR
  
  p3 <- ggplot(data=res, aes(x=factor(NFLR_threshold), y = n_transcripts)) +
    geom_line(group=1) +
    scale_x_discrete(name = expression(paste("Normalised expression (", NFLR[T], ")"))) +
    scale_y_continuous(name = "No. of transcripts detected") +
    mytheme +
    guides(fill="none") +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size=0.4),
      axis.text.x = element_text(size=5, angle=90, hjust=1, vjust=0.5)
    )
  
  png(str_c(outpath, "/", prefix, "_NFLR_curve_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p3)
  dev.off()
  
}

