test_that("TALON mult-sample", {
  
  file <- "../test_data/example_talon_read_annot_multiSample.tsv"
  gene_id <- "SIRV1.1"
  perc_As <- 0.80
  abundance_file_type <- "TALON"
  min_exp_perSample <- 0
  min_exp_mean <- 0.3
  min_sample_perc <- 0
  
  x <- read.table(file, header=TRUE, sep="\t")
  
  x.gene <- x %>% 
    dplyr::filter(annot_gene_id == gene_id, fraction_As <= perc_As)

  x.gene <- .process_talon_file(x.gene)
  x.gene <- .add_normalised_counts(x.gene, abundance_file_type)
  
  req_cols <- c("associated_gene", "isoform", "Isoform_class", "FL_mean", "NFLR_mean")
  expect_contains(colnames(x.gene), req_cols)
  expect_true(stringr::str_starts(colnames(x.gene), "FL.") %>% any())
  expect_true(stringr::str_starts(colnames(x.gene), "NFLR.") %>% any())
  expect_equal(nrow(x.gene), 30)

  res <- .calculate_transcripts_detected_numbers(x.gene)
  expect_is(res, "list")
  expect_contains(colnames(res$data), c("sample", "n_transcripts", "NFLR_threshold"))
  expect_equal(nrow(res$data), 47)
  expect_equal(res$n_samples, 3)
  expect_true("NFLR_mean" %in% res$data$sample)

  df1 <- res$data %>%
    dplyr::filter(sample %in% "NFLR_mean")
  
  sd <- res$data %>%
    dplyr::filter(!sample %in% "NFLR_mean") %>%
    group_by(NFLR_threshold) %>%
    summarise(sd = sd(n_transcripts), .groups = "keep")
  
  df1 <- left_join(df1, sd, by = "NFLR_threshold")
  expect_contains(colnames(df1), c("sample", "n_transcripts", "NFLR_threshold", "sd"))
  expect_equal(nrow(df1), 16)

  res <- .quantify_transcript_categories(x.gene) 
  expect_is(res, "list")
  expect_contains(colnames(res$tc.count), c("Isoform_class", "count_mean", "count_sd"))
  expect_equal(nrow(res$tc.count), 5)
  expect_contains(colnames(res$tc.nflr), c("Isoform_class", "isoform", "Sample", "NFLR"))
  expect_equal(nrow(res$tc.nflr), 44)
  

  x.pass <- .filter_transcripts(x.gene, min_exp_perSample, min_exp_mean, min_sample_perc)
  expect_contains(colnames(x.pass), c("total_samples", "n_pass", "perc_pass"))
  expect_equal(nrow(x.pass), 9)

  res <- .rank_transcripts(x.pass)
  expect_is(res, "list")
  expect_contains(colnames(res$x.ranked.mean_sd), c("isoform_index", "Isoform_class", "isoform", "NFLR_mean", "NFLR_sd"))
  expect_equal(nrow(res$x.ranked.mean_sd), 9)
  expect_equal(res$n_samples, 3)
  expect_contains(colnames(res$x.ranked), c("isoform_index"))
  expect_equal(nrow(res$x.ranked), 9)

})



test_that("TALON single-sample", {
  
  file <- "../test_data/example_talon_read_annot_singleSample.tsv"
  gene_id <- "SIRV1.1"
  perc_As <- 0.80
  abundance_file_type <- "TALON"
  min_exp_perSample <- 0.3
  min_exp_mean <- 0
  min_sample_perc <- 0
  
  x <- read.table(file, header=TRUE, sep="\t")
  
  x.gene <- x %>% 
    dplyr::filter(annot_gene_id == gene_id, fraction_As <= perc_As)

  x.gene <- .process_talon_file(x.gene)
  x.gene <- .add_normalised_counts(x.gene, abundance_file_type)
  
  req_cols <- c("associated_gene", "isoform", "Isoform_class")
  
  expect_contains(colnames(x.gene), req_cols)
  expect_true(stringr::str_starts(colnames(x.gene), "FL.") %>% any())
  expect_true(stringr::str_starts(colnames(x.gene), "NFLR.") %>% any())
  expect_equal(nrow(x.gene), 19)

  res <- .calculate_transcripts_detected_numbers(x.gene)
  expect_is(res, "list")
  expect_contains(colnames(res$data), c("sample", "n_transcripts", "NFLR_threshold"))
  expect_equal(nrow(res$data), 10)
  expect_equal(res$n_samples, 1)


  x.pass <- .filter_transcripts(x.gene, min_exp_perSample, min_exp_mean, min_sample_perc)
  expect_contains(colnames(x.pass), c("total_samples", "n_pass", "perc_pass"))
  expect_equal(nrow(x.pass), 10)

  res <- .rank_transcripts(x.pass)
  expect_is(res, "list")
  expect_equal(res$n_samples, 1)
  expect_contains(colnames(res$x.ranked), c("isoform_index", "NFLR_mean"))
  expect_equal(nrow(res$x.ranked), 10)
  
})
