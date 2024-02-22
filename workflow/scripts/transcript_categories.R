options(warn=-1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(optparse)
})

option_list = list(
  make_option(c("-a", "--abundance"), type="character", default=NULL, help="input NORMALISED abundance file: gene_normalised_abundance.txt [REQUIRED] [default= %default]", metavar="character"),
  make_option(c("--abundance_file_type"), type="character", default="SQANTI", help=" abundance file type, accepted choices: ['SQANTI', 'TALON'] [default= %default]", metavar="character"),
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

if (opt$abundance_file_type != "SQANTI" && opt$abundance_file_type != "TALON"){
  message("\n", Sys.time(), ": *--- Incorrect --abundance_file_type value ---*")
  message(Sys.time(), ": *--- accepted choices: ['SQANTI', 'TALON'] ---*")
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


if(opt$abundance_file_type == "SQANTI"){
  x.gene$Isoform_class <- forcats::fct_recode(
      x.gene$Isoform_class,
      "Coding known\n(complete match)" = "Coding known (complete match)",
      "Coding known\n(alternate 3'/5' end)" = "Coding known (alternate 3/5 end)"
  )
}

n_samples <- x.gene %>%
  dplyr:::select(starts_with("NFLR")) %>% 
  ncol()


########### plots ###########
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
  
  message(Sys.time(), ": Multiple samples detected - ", n_samples - 1, " samples")

  res <- .quantify_transcript_categories(x.gene)  

  message(Sys.time(), ': plot per sample - number of unique transcripts in each transcript category')
  p1 <- ggplot(data=res$tc.count, aes(x=Isoform_class, y = count_mean, fill=Isoform_class)) +
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
  
  png(str_c(opt$outpath, "/tc_count_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()
  
  
  ## expression of transcript categories ##
  message(Sys.time(), ': plot per sample - expression of each transcript category')
  
  p2 <- ggplot(data=res$tc.nflr, aes(x=Isoform_class, y = log2(NFLR), fill=Isoform_class)) +
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
  
  png(str_c(opt$outpath, "/tc_NFLR_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p2)
  dev.off() 
  
} else {
  
  message(Sys.time(), ": Single sample detected - ", n_samples - 1, " sample")
  x.gene$NFLR_mean <- x.gene %>% dplyr::select(starts_with("NFLR")) %>% pull(1)
  
}


message(Sys.time(), ': plot across sample - number of unique transcripts in each transcript category')
png(str_c(opt$outpath, "/tc_count_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

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
message(Sys.time(), ': plot across sample - expression of each transcript category')
## per transcript ##

png(str_c(opt$outpath, "/tc_NFLR_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

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


message(Sys.time(), ': All done !')
sink()