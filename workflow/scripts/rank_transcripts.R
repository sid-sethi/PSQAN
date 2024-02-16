options(warn=-1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(optparse)
})

option_list = list(
  make_option(c("-a", "--abundance"), type="character", default=NULL, help="input abundance file: gene_normalised_abundance.txt or filteres_transcripts.txt [REQUIRED] [default= %default]", metavar="character"),
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


#### Transcripts (ranked) vs. NFLR ######
res <- .rank_transcripts(x.gene)

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

  message(Sys.time(), ': plotting ranked transcripts with sd across samples as error bars')

  p1 <- ggplot(data=res$x.ranked.mean_sd, aes(x=isoform_index, fill = Isoform_class)) +
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

  png(str_c(opt$outpath, "/transcriptsRanked_perSample.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)
  plot(p1)
  dev.off()

}


message(Sys.time(), ': plotting ranked transcripts without error bars')

png(str_c(opt$outpath, "/transcriptsRanked_main.png"), bg="transparent",units="in",width = 4.25, height= 3.75 ,res=600)

ggplot(data=res$x.ranked, aes(x=isoform_index, fill = Isoform_class)) +
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


message(Sys.time(), ': All done !')
sink()