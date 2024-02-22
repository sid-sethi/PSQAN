options(warn=-1)

suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(optparse)
  library(pointblank)
})


option_list = list(
  make_option(c("-a", "--abundance"), type="character", default=NULL, help="input abundance file. SQANTI: *classification.txt; TALON: *talon_read_annot.tsv [REQUIRED] [default= %default]", metavar="character"),
  make_option(c("--abundance_file_type"), type="character", default="SQANTI", help=" abundance file type, accepted choices: ['SQANTI', 'TALON'] [default= %default]", metavar="character"),
  make_option(c("-g", "--gene_id"), type="character", default=NULL, help="Ensembl gene id of gene of interest [REQUIRED] [default= %default]", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=getwd(), help="output file directory path [default= cwd]", metavar="character"),
  make_option(c("-l", "--log_file"), type="character", default=NULL, help=" optional log file to redirect error messages [OPTIONAL] [default= %default]", metavar="character"),
  make_option(c("-p", "--perc_As"), type="double", default=80, help="maximum value of percent As allowed downstream of transcript termination site. For TALON input, this number would be divided by 100 to convert into fraction [default= %default]", metavar="double"),
  make_option(c("-u", "--utils_path"), type="character", default=".", help="full path of functions file - utils.R [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$abundance) || is.null(opt$gene_id)){
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

x = read.table(opt$abundance, header=TRUE, sep="\t")


if(opt$abundance_file_type == "TALON"){

  opt$perc_As <- round(opt$perc_As/100, 2)
  
  # validate required columns
  req_cols <- c("annot_gene_id", "fraction_As", "annot_transcript_id", "gene_novelty", "transcript_novelty", "ISM_subtype", "dataset")

  tryCatch(
    expr = {
      
      # validate number of rows
      if(! nrow(x) >= 1){
        stop(Sys.time(), ": Found ", nrow(x), " rows, expected >= 1 rows\n")
      }
      
      x %>% 
        pointblank::expect_col_exists(columns = req_cols) %>%
        pointblank::expect_col_is_character(columns = c("annot_gene_id", "annot_transcript_id", "gene_novelty", "transcript_novelty", "ISM_subtype", "dataset"))

      # validate fraction_As column type
      perc_col_type <- class(x$fraction_As)
      if(perc_col_type != "numeric" && perc_col_type != "integer"){
        stop(Sys.time(), ": Column type of 'fraction_As' should be either numeric or integer; instead found column type as ", perc_col_type, "\n")
      }

      # validate gene_id
      gene_val <- opt$gene_id %in% x$annot_gene_id
      if(! gene_val){
        message(Sys.time(), ": Please check your Gene ID and/or ", opt$abundance_file_type, " results file")
        message(Sys.time(), ": Make sure your Gene ID is present in ", opt$abundance_file_type, " results")
        stop(Sys.time(), ": Did not find ", opt$gene_id, " in annot_gene_id column")
      } else{
        message(Sys.time(), ": File = ", opt$abundance)
        message(Sys.time(), ": *-- ", opt$abundance_file_type, " file succesfully validated !!")
      }

    },
    error = function(e) {
      message(Sys.time(), ": File = ", opt$abundance)
      message(Sys.time(), ": *--Error in ", opt$abundance_file_type, " file validation:\n", e$message, "\n")
      stop("*----Stopping further execution due to an input file validation error")
    }
  )


  # select data for gene of interest
  x.gene <- x %>% 
    dplyr::filter(annot_gene_id == opt$gene_id, fraction_As <= opt$perc_As)

  # process Talon file
  x.gene <- .process_talon_file(x.gene)

  if(! nrow(x.gene) >= 1){
    message(Sys.time(), ": File = ", opt$abundance)
    message(Sys.time(), ": *--- No transcripts detected for gene ", opt$gene_id, " after intial filtering - fraction_As <= ", opt$perc_As)
    stop("*----Stopping further execution")
  }

  # validate FL column
  FL_col <- stringr::str_starts(colnames(x.gene), "FL") %>% any()
  if(! FL_col){
    message(Sys.time(), ": File = ", opt$abundance)
    message(Sys.time(), ": *--Error in ", opt$abundance_file_type, " processed data frame")
    message(Sys.time(), ": Did not find any FL column after processing TALON file")
    stop("*----Stopping further execution due to an internal data frame validation error")
  }

}


if(opt$abundance_file_type == "SQANTI"){

  # validate required columns
  req_cols <- c("associated_gene", "perc_A_downstream_TTS", "RTS_stage")

  tryCatch(
    expr = {
      
      # validate number of rows
      if(! nrow(x) >= 1){
        stop(Sys.time(), ": Found ", nrow(x), " rows, expected >= 1 rows\n")
      }
      
      x %>% 
        pointblank::expect_col_exists(columns = req_cols) %>%
        pointblank::expect_col_is_character(columns = "associated_gene") %>% 
        pointblank::expect_col_is_logical(columns = "RTS_stage")

      # validate perc_A_downstream_TTS column type
      perc_col_type <- class(x$perc_A_downstream_TTS)
      if(perc_col_type != "numeric" && perc_col_type != "integer"){
        stop(Sys.time(), ": Column type of 'perc_A_downstream_TTS' should be either numeric or integer; instead found column type as ", perc_col_type, "\n")
      }

      # validate FL column
      FL_col <- stringr::str_starts(colnames(x), "FL") %>% any()
      if(! FL_col){
        stop(Sys.time(), ": Did not find any 'FL' column(s)\n")
      }

      # validate gene_id
      gene_val <- opt$gene_id %in% x$associated_gene
      if(! gene_val){
        message(Sys.time(), ": Please check your Gene ID and/or ", opt$abundance_file_type, " results file")
        message(Sys.time(), ": Make sure your Gene ID is present in ", opt$abundance_file_type, " results")
        stop(Sys.time(), ": Did not find ", opt$gene_id, " in associated_gene column")
      } else{
        message(Sys.time(), ": File = ", opt$abundance)
        message(Sys.time(), ": *-- ", opt$abundance_file_type, " file succesfully validated !!")
      }

    },
    error = function(e) {
      message(Sys.time(), ": File = ", opt$abundance)
      message(Sys.time(), ": *--Error in ", opt$abundance_file_type,  " file validation:\n", e$message, "\n")
      stop("*----Stopping further execution due to an input file validation error")
    }
  )


  # select data for gene of interest
  x.gene <- x %>%
    dplyr::filter(associated_gene == opt$gene_id, perc_A_downstream_TTS <= opt$perc_As, RTS_stage == "FALSE")

  if(! nrow(x.gene) >= 1){
    message(Sys.time(), ": File = ", opt$abundance)
    message(Sys.time(), ": *--- No transcripts detected for gene ", opt$gene_id, " after intial filtering - perc_A_downstream_TTS <= ", opt$perc_As, ", RTS_stage == FALSE")
    stop("*----Stopping further execution")
  }


  #####################################################################
  ################ Adding custom isoform categorisation ###############
  #####################################################################

  message(Sys.time(), ': Adding PSQAN isoform categorisation')
  x.gene <- .add_psqan_isoform_categorisation(x.gene)
  

}



###########################################################
################## Normalisation ##########################
###########################################################

message(Sys.time(), ': Performing transcript abundence normalisation....')
x.gene <- .add_normalised_counts(x.gene, opt$abundance_file_type)


write.table(x.gene, str_c(opt$outpath, "/gene_normalised_abundance.txt"), col.names = TRUE, row.names=FALSE, sep="\t", quote = FALSE)


message(Sys.time(), ': All done !')
sink()