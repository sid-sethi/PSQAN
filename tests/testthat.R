options(warn=-1)
suppressPackageStartupMessages({
    library(tidyverse)
    library(testthat)
})

# load functions
source("./workflow/scripts/utils.R")

testthat::test_dir("./tests/testthat/", reporter=c("progress", "minimal", "location"))