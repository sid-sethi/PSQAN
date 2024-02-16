options(warn=-1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

option_list = list(
  make_option(c("-a", "--abundance"), type="character", default=NULL, help="input NORMALISED abundance file: gene_normalised_abundance.txt [REQUIRED] [default= %default]", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=getwd(), help="output file directory path [default= cwd]", metavar="character"),
  make_option(c("-l", "--log_file"), type="character", default=NULL, help=" optional log file to redirect error messages [OPTIONAL] [default= %default]", metavar="character"),
  make_option(c("-u", "--utils_path"), type="character", default=".", help="full path of functions file - utils.R [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$abundance)){
  message("\n", Sys.time(), ": *--- Missing required options ---*")
  print_help(opt_parser)
  stop("\n", call.=FALSE)
}

if(! is.null(opt$log_file)){
  log <- file(opt$log_file, open = "wt")
  sink(log, append = FALSE, type = "message")
}

# source functions
source(str_c(opt$utils_path, "/utils.R"))

x.gene = read.table(opt$abundance, header=TRUE, sep="\t")

tryCatch(
  expr = {
    
    # validate number of rows
    if(! nrow(x.gene) >= 1){
      stop(Sys.time(), ": Found ", nrow(x.gene), " rows, expected >= 1 rows\n")
    }

    # validate NFLR column
    NFLR_col <- stringr::str_starts(colnames(x.gene), "NFLR") %>% any()
    if(! NFLR_col){
      stop(Sys.time(), ": Did not find any 'NFLR' column(s) in abundance file\n")
    } else{
      message(Sys.time(), ": File = ", opt$abundance)
      message(Sys.time(), ': *--file succesfully validated !!')
    }

  },
  error = function(e) {
    message(Sys.time(), ": File = ", opt$abundance)
    message(Sys.time(), ": *--Error in ", opt$abundance_file_type, " file validation:\n", e$message, "\n")
    stop("*----Stopping further execution due to an input file validation error")
  }
)


# calculating number of unique transcripts at every NFLR value
res <- .calculate_transcripts_detected_numbers(x.gene)


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



if(res$n_samples > 1){

  message(Sys.time(), ": Multiple samples detected - ", res$n_samples - 1, " samples")
  message(Sys.time(), ': plotting - transcripts detected vs. NFLR_mean (mean NFLR across all samples)')
  
  df1 <- res$data %>%
    dplyr::filter(sample %in% "NFLR_mean")
  
  sd <- res$data %>%
    dplyr::filter(!sample %in% "NFLR_mean") %>%
    group_by(NFLR_threshold) %>%
    summarise(sd = sd(n_transcripts), .groups = "keep")
  
  df1 <- left_join(df1, sd, by = "NFLR_threshold")
  write.table(df1, str_c(opt$outpath, "/NFLR_curve_values.txt"), quote=FALSE, sep="\t", row.names=FALSE)
  
  
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
    
    
  png(str_c(opt$outpath, "/NFLR_curve_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()
  
  
  message(Sys.time(), ': plotting - transcripts detected vs. NFLR per sample')

  df2 <- res$data %>%
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
  
  png(str_c(opt$outpath, "/NFLR_curve_perSample.png"), bg="transparent",units="in",width = 5.25, height= 3.75 ,res=600)
  plot(p2)
  dev.off()
  
  
}else{
  
  message(Sys.time(), ": Single sample detected - ", res$n_samples - 1, " sample")
  message(Sys.time(), ': plotting - transcripts detected vs. NFLR')
  
  p3 <- ggplot(data=res$data, aes(x=factor(NFLR_threshold), y = n_transcripts)) +
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
  
  png(str_c(opt$outpath, "/NFLR_curve_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p3)
  dev.off()
  
}


message(Sys.time(), ': All done !')
sink()