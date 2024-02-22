options(warn=-1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})


option_list = list(
  make_option(c("-a", "--abundance"), type="character", default=NULL, help="input NORMALISED abundance file: gene_normalised_abundance.txt [REQUIRED] [default= %default]", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=getwd(), help="output file directory path [default= cwd]", metavar="character"),
  make_option(c("--NFLR_perSample_thresh"), type="double", default=0.3, help="NFLR threshold applied per sample [default= %default]", metavar="double"),
  make_option(c("--NFLR_mean_thresh"), type="double", default=0.3, help="NFLR mean threshold [default= %default]", metavar="double"),
  make_option(c("--min_sample_perc"), type="double", default=0, help="minimum number of samples (%) which should pass NFLR threshold [default= %default]", metavar="double"),
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


########### filtering ###########
x.pass <- .filter_transcripts(x.gene, opt$NFLR_perSample_thresh, opt$NFLR_mean_thresh, opt$min_sample_perc)


# validate number of retained transcripts
if(! nrow(x.pass) >= 1){
  message(Sys.time(), ": WARNING - NO transcripts were retained after filtering on provided read count thresholds")
  message(Sys.time(), ": WARNING - run PSQAN again with lenient thresholds")
} else {
  message(Sys.time(), ": Found ", nrow(x.pass), " transcripts after filtering!")
}

write.table(x.pass, str_c(opt$outpath, "/filtered_transcripts.txt"), col.names = TRUE, row.names=FALSE, sep="\t", quote = FALSE)



message(Sys.time(), ': All done !')
sink()