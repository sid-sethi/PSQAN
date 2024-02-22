test_that("SQANTI multi-sample", {
  
  file <- "../test_data/sqanti_classification_multiSample.txt"
  gene_id <- "ENSG00000TEST"
  perc_As <- 80
  abundance_file_type <- "SQANTI"
  min_exp_perSample <- 0
  min_exp_mean <- 0.3
  min_sample_perc <- 8
  
  x <- read.table(file, header=TRUE, sep="\t")
  
  x.gene <- x %>%
    dplyr::filter(associated_gene == gene_id, perc_A_downstream_TTS <= perc_As, RTS_stage == "FALSE")
  
  x.gene <- .add_psqan_isoform_categorisation(x.gene)
  x.gene <- .add_normalised_counts(x.gene, abundance_file_type)

  req_cols <- c("associated_gene", "isoform", "Isoform_class", "FL_mean", "NFLR_mean")
  expect_contains(colnames(x.gene), req_cols)
  expect_true(stringr::str_starts(colnames(x.gene), "FL.") %>% any())
  expect_true(stringr::str_starts(colnames(x.gene), "NFLR.") %>% any())
  expect_equal(nrow(x.gene), 326)

  res <- .calculate_transcripts_detected_numbers(x.gene)
  expect_is(res, "list")
  expect_contains(colnames(res$data), c("sample", "n_transcripts", "NFLR_threshold"))
  expect_equal(nrow(res$data), 1583)
  expect_equal(res$n_samples, 39)
  expect_true("NFLR_mean" %in% res$data$sample)

  df1 <- res$data %>%
    dplyr::filter(sample %in% "NFLR_mean")
  
  sd <- res$data %>%
    dplyr::filter(!sample %in% "NFLR_mean") %>%
    group_by(NFLR_threshold) %>%
    summarise(sd = sd(n_transcripts), .groups = "keep")
  
  df1 <- left_join(df1, sd, by = "NFLR_threshold")
  expect_contains(colnames(df1), c("sample", "n_transcripts", "NFLR_threshold", "sd"))
  expect_equal(nrow(df1), 41)


  res <- .quantify_transcript_categories(x.gene) 
  expect_is(res, "list")
  expect_contains(colnames(res$tc.count), c("Isoform_class", "count_mean", "count_sd"))
  expect_equal(nrow(res$tc.count), 5)
  expect_contains(colnames(res$tc.nflr), c("Isoform_class", "isoform", "Sample", "NFLR"))
  expect_equal(nrow(res$tc.nflr), 2365)


  x.pass <- .filter_transcripts(x.gene, min_exp_perSample, min_exp_mean, min_sample_perc)
  expect_contains(colnames(x.pass), c("total_samples", "n_pass", "perc_pass"))
  expect_equal(nrow(x.pass), 22)

  res <- .rank_transcripts(x.pass)
  expect_is(res, "list")
  expect_contains(colnames(res$x.ranked.mean_sd), c("isoform_index", "Isoform_class", "isoform", "NFLR_mean", "NFLR_sd"))
  expect_equal(nrow(res$x.ranked.mean_sd), 22)
  expect_equal(res$n_samples, 39)
  expect_contains(colnames(res$x.ranked), c("isoform_index"))
  expect_equal(nrow(res$x.ranked), 22)

})



test_that("SQANTI single-sample", {
  
  file <- "../test_data/sqanti_classification_singleSample.txt"
  gene_id <- "ENSG00000TEST"
  perc_As <- 80
  abundance_file_type <- "SQANTI"
  min_exp_perSample <- 0.3
  min_exp_mean <- 0
  min_sample_perc <- 0
  
  x <- read.table(file, header=TRUE, sep="\t")
  
  x.gene <- x %>%
    dplyr::filter(associated_gene == gene_id, perc_A_downstream_TTS <= perc_As, RTS_stage == "FALSE")
  
  x.gene <- .add_psqan_isoform_categorisation(x.gene)
  x.gene <- .add_normalised_counts(x.gene, abundance_file_type)

  req_cols <- c("associated_gene", "isoform", "Isoform_class")
  expect_contains(colnames(x.gene), req_cols)
  expect_true(stringr::str_starts(colnames(x.gene), "FL") %>% any())
  expect_true(stringr::str_starts(colnames(x.gene), "NFLR") %>% any())
  expect_equal(nrow(x.gene), 326)

  res <- .calculate_transcripts_detected_numbers(x.gene)
  expect_is(res, "list")
  expect_contains(colnames(res$data), c("sample", "n_transcripts", "NFLR_threshold"))
  expect_equal(nrow(res$data), 31)
  expect_equal(res$n_samples, 1)


  x.pass <- .filter_transcripts(x.gene, min_exp_perSample, min_exp_mean, min_sample_perc)
  expect_contains(colnames(x.pass), c("total_samples", "n_pass", "perc_pass"))
  expect_equal(nrow(x.pass), 22)
  
  res <- .rank_transcripts(x.pass)
  expect_is(res, "list")
  expect_equal(res$n_samples, 1)
  expect_contains(colnames(res$x.ranked), c("isoform_index", "NFLR_mean"))
  expect_equal(nrow(res$x.ranked), 22)

})