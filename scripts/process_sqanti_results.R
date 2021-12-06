args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(ORFik)
  library(Biostrings)
})


sqanti_class <- args[1] # sqanti classification.txt file
sqanti_fasta <- args[2] # sqanti corrected.fasta file
prefix <- args[3]
outpath <- args[4]
gene_id <- args[5] # multiple genes can be input separated by ";"


x = read.table(sqanti_class, header=TRUE, sep="\t") %>%
  dplyr::filter(perc_A_downstream_TTS <= 80, RTS_stage %in% "FALSE")


# select data for gene of interest
if(!is.na(gene_id)) {
  if(!gene_id == ""){
    ids <- str_split(gene_id, ";") %>% unlist()
    # select gene data
    x.gene <- x %>% dplyr::filter(associated_gene %in% ids)
  }else{
    x.gene <- x
  }
}else{
  x.gene <- x
}



# # adding mean FL counts
 x.gene$FL_mean = x.gene %>% dplyr::select(starts_with("FL.")) %>% rowMeans()


#####################################################################
################ Adding custom isoform categorisation ###############
#####################################################################


x.gene <- x.gene %>%
  dplyr::mutate(
    Isoform_class = case_when(

      coding == "non_coding" &
        structural_category != "full-splice_match" ~ "Non-coding novel",

      coding == "non_coding" &
        structural_category == "full-splice_match" ~ "Non-coding known",

      coding == "coding" &
        predicted_NMD == "TRUE" &
        structural_category != "full-splice_match"~ "NMD novel",

      coding == "coding" &
        predicted_NMD == "TRUE" &
        structural_category == "full-splice_match"~ "NMD known",

      predicted_NMD == "FALSE" &
        structural_category != "full-splice_match" ~ "Coding novel",

      predicted_NMD == "FALSE" &
        structural_category == "full-splice_match" &
        subcategory == "reference_match" ~ "Coding known (complete match)",

      predicted_NMD == "FALSE" &
        structural_category == "full-splice_match" &
        subcategory != "reference_match" ~ "Coding known (alternate 3'/5' end)"

    )
  )



###########################################################
################ ORF prediction using ORFik ###############
###########################################################

seqs <- Biostrings::readDNAStringSet(sqanti_fasta, format = "fasta")

# filter non-gene related transcripts
seqs <- seqs[x.gene$isoform]

orfs.gr <- ORFik::findORFs(seqs, longestORF = TRUE, startCodon = "ATG") %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(seqnames = names(seqs)[as.integer(.$names)]) %>%
  group_by(seqnames) %>%
  slice_max(width, n=1, with_ties=FALSE) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)


orfs.cds <- getSeq(seqs, orfs.gr)
names(orfs.cds) <- seqnames(orfs.gr)
orfs.aa <- Biostrings::translate(orfs.cds)


cds.df <- orfs.cds %>%
  as.data.frame() %>%
  rownames_to_column("isoform") %>%
  plyr::rename(c("x" = "ORFik_cds"))


aa.df <- orfs.aa %>%
  as.data.frame() %>%
  rownames_to_column("isoform") %>%
  plyr::rename(c("x" = "ORFik_aa"))


x.gene <- x.gene %>%
  dplyr::left_join(cds.df, by = "isoform") %>%
  dplyr::left_join(aa.df, by = "isoform")




###########################################################
################## Normalisation ##########################
###########################################################

# normalise FLR per transcript per sample

x.norm <- x.gene %>%
  dplyr:::select(starts_with("FL.")) %>%
  mutate_at(vars(starts_with("FL.")), ~./sum(.)*100) %>%
  dplyr::rename_with(~str_replace(., "FL", "NFLR"))

x.gene <- bind_cols(x.gene, x.norm)

# normalised expression per transcript across all samples; NE = mean(NFLR)

x.gene$NE = x.gene %>% dplyr::select(starts_with("NFLR.")) %>% rowMeans()


write.table(x.gene, str_c(outpath, "/", prefix, "_classification_processed.txt"), col.names = TRUE, row.names=FALSE, sep="\t", quote = FALSE)
