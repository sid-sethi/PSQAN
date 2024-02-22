
#' Process Talon data to a gerneralise format
#' @param x data frame containing talon data
#' @return a data frame
.process_talon_file <- function(x){
    
    x <- x %>%
        dplyr::group_by(annot_gene_id, annot_transcript_id, gene_novelty, transcript_novelty, ISM_subtype, dataset) %>% 
        summarise(abundance = n(), .groups = "keep") %>%
        ungroup() %>%
        tidyr::pivot_wider(id_cols = annot_gene_id:ISM_subtype, names_from = dataset, values_from = abundance, names_prefix = "FL.") %>%
        dplyr::rename(
            "associated_gene" = annot_gene_id,
            "isoform" = annot_transcript_id,
            "Isoform_class" = transcript_novelty
        )

    return(x)

}


#' Add custom PSQAN Isoform categorisation
#' @param x.gene data frame containing sqanti data
#' @return a data frame
.add_psqan_isoform_categorisation <- function(x.gene){

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

            #is.na(predicted_NMD) &
            #  structural_category != "full-splice_match" ~ "Coding novel",

            predicted_NMD == "FALSE" &
                structural_category == "full-splice_match" &
                subcategory == "reference_match" ~ "Coding known (complete match)",

            predicted_NMD == "FALSE" &
                structural_category == "full-splice_match" &
                subcategory %in% c("alternative_3end", "alternative_3end5end", "alternative_5end") ~ "Coding known (alternate 3'/5' end)",

            TRUE ~ "Other"

            )
        )

    return(x.gene)

}


#' normalise full length read counts
#' @param x.gene data frame containing FL count data
#' @param file_type abundance file type: TALON or SQANTI
#' @return a data frame
.add_normalised_counts <- function(x.gene, file_type) {

    n_samples <- x.gene %>%
        dplyr:::select(starts_with("FL")) %>%
        ncol()

    # normalise FLR per transcript per sample
    if(n_samples > 1) {

        if(file_type == "TALON"){
            exact_samples <- n_samples
        } else {
            exact_samples <- n_samples -1
        }

        # replace NA abundance as 0. NAs are found in TALON file
        x.gene <- x.gene %>% 
            mutate_at(vars(starts_with("FL.")), ~replace_na(., 0))

        message(Sys.time(), ": Detected ", exact_samples, " samples in ", file_type , " file")
        
        # # adding mean FL counts
        x.gene$FL_mean = x.gene %>% dplyr::select(starts_with("FL.")) %>% rowMeans()

        x.norm <- x.gene %>%
            dplyr:::select(starts_with("FL.")) %>%
            mutate_at(vars(starts_with("FL.")), ~./sum(.)*100) %>%
            dplyr::rename_with(~str_replace(., "FL", "NFLR"))

        x.norm[is.na(x.norm)] <- 0

        # normalised expression per transcript across all samples; NFLR(T) = mean(NFLR)
        x.norm$NFLR_mean = rowMeans(x.norm)

    }else{

        # replace NA abundance as 0. NAs are found in TALON file
        x.gene <- x.gene %>% 
            mutate_at(vars(starts_with("FL")), ~replace_na(., 0))
        
        message(Sys.time(), ": Detected ", n_samples, " sample in ", file_type, " file")

        x.norm <- x.gene %>%
            dplyr:::select(starts_with("FL")) %>%
            mutate_at(vars(starts_with("FL")), ~./sum(.)*100) %>%
            dplyr::rename_with(~str_replace(., "FL", "NFLR"))

    }


    x.gene <- bind_cols(x.gene, x.norm)

    return(x.gene)

}



#' calculating number of unique transcripts at every NFLR value
#' @param x.gene data frame containing normalised count data
#' @return a list containing a data frame and number of samples
.calculate_transcripts_detected_numbers <- function(x.gene){

    n_samples <- x.gene %>%
        dplyr:::select(starts_with("NFLR")) %>% 
        ncol()

    counts.long <- x.gene %>%
        dplyr::select(isoform, starts_with("NFLR")) %>%
        pivot_longer(cols = c(-isoform), names_to = "sample", values_to = "read_count")

    # if multiple samples, remove primer tags
    counts.long$sample <- counts.long$sample %>% str_replace("NFLR.*\\.", "") %>% str_replace("_3p", "")

    
    ################# No. of unique transcripts detected at each NFLR threshold ###############
    message(Sys.time(), ': Extracting NFLR values')
    if(n_samples > 1){
        threshs <- x.gene$NFLR_mean %>% round(2) %>% sort() %>% unique()
    }else{
        threshs <- x.gene %>% select(starts_with("NFLR")) %>% pull(var = 1) %>% round(2) %>% sort() %>% unique()
    }


    message(Sys.time(), ': calculating number of unique transcripts at every NFLR value')

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

    result <- list()
    result$data <- res
    result$n_samples <- n_samples

    return(result) 

}



#' quantify transcript categories metrics for multisample
#' @param x.gene data frame containing normalised count data
#' @return a list containing two data frames
.quantify_transcript_categories <- function(x.gene){

    ## Number of unique transcripts vs. transcript category ##
    tc.count <- x.gene %>%
        dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
        pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
        dplyr::filter(NFLR >0) %>%
        group_by(Isoform_class, Sample) %>%
        dplyr::summarise(count = n(), .groups = "keep") %>%
        group_by(Isoform_class) %>%
        summarise(count_mean = mean(count), count_sd = sd(count), .groups = "keep")

    ## expression of transcript categories ##
    ## per sample ##

    tc.nflr <- x.gene %>%
        dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
        pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
        dplyr::filter(NFLR >0)
    
    result <- list()
    result$tc.count <- tc.count
    result$tc.nflr <- tc.nflr

    return(result)

}



#' filter transcripts based on provided NFLR thresholds
#' @param x.gene data frame containing normalised count data
#' @param NFLR_perSample_thresh minimum normalised expression per sample
#' @param NFLR_mean_thresh minimum normalised expression averaged across all samples
#' @param min_sample_perc minimum number of samples (%) which should pass NFLR_perSample threshold
#' @return a data frame
.filter_transcripts <- function(x.gene, NFLR_perSample_thresh, NFLR_mean_thresh, min_sample_perc){

    n_samples <- x.gene %>%
        dplyr:::select(starts_with("NFLR")) %>%
        ncol()

    if(n_samples > 1) {

        message(Sys.time(), ": Multiple samples detected - ", n_samples - 1, " samples")

        x.gene$total_samples <- x.gene %>% dplyr::select(starts_with("NFLR.")) %>% ncol()

        pass = x.gene %>% dplyr::select(starts_with("NFLR.")) %>% .[] >= NFLR_perSample_thresh

        x.pass <- x.gene %>%
            mutate(
                n_pass = rowSums(pass),
                perc_pass = round((n_pass/total_samples)*100,2)
            ) %>%
            dplyr::filter(NFLR_mean >= NFLR_mean_thresh, perc_pass >= min_sample_perc)

    } else {

        message(Sys.time(), ": Single sample detected - ", n_samples, " sample")

        x.pass <- x.gene %>%
            mutate(
                total_samples = n_samples,
                n_pass = NA,
                perc_pass = NA
            ) %>%
            dplyr::filter_at(vars(starts_with("NFLR")), any_vars(. >= NFLR_perSample_thresh))

    }

    return(x.pass)


}



#' rank transcripts based on their expression
#' @param x.gene data frame containing normalised count data
#' @return a list of three elements
.rank_transcripts <- function(x.gene){

    result <- list()
    
    result$n_samples <- x.gene %>%
        dplyr:::select(starts_with("NFLR")) %>%
        ncol()

    if(result$n_samples > 1){

        message(Sys.time(), ": Multiple samples detected - ", result$n_samples - 1, " samples")
        
        result$x.ranked.mean_sd <- x.gene %>%
            dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
            pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
            group_by(isoform, Isoform_class) %>%
            summarise(NFLR_mean = mean(NFLR), NFLR_sd = sd(NFLR), .groups = "keep") %>%
            arrange(desc(NFLR_mean)) %>%
            tibble::rowid_to_column(., "isoform_index")
 
    } else {

        message(Sys.time(), ": Single sample detected - ", result$n_samples, " sample")
        x.gene$NFLR_mean <- x.gene %>% dplyr::select(starts_with("NFLR")) %>% pull(var = 1)
        result$x.ranked.mean_sd <- NULL

    }

    result$x.ranked <- x.gene %>%
        arrange(desc(NFLR_mean)) %>%
        tibble::rowid_to_column(., "isoform_index")
        #mutate(isoform_index = 1:nrow(x.gene))

    
    return(result)

}