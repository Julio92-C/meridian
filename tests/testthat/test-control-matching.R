# Tests for metadata.control_match_cols: each sample is compared only against
# the negative controls that share its value in the first listed metadata
# column that has any (e.g. its sequencing run, else its batch).
suppressPackageStartupMessages({ library(dplyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
R_DIR <- file.path("..", "..", "R")
source(file.path(R_DIR, "02_clean_data.R"))

logged <- character(0)
pipeline_log <<- function(cfg, msg) logged <<- c(logged, msg)

meta <- data.frame(
  sample  = c("s1", "s2", "s3", "s4", "c1", "c2", "c3"),
  seq_run = c("R1", "R2", "R4", "R9", "R1", "R2", "R7"),
  batch   = c("B1", "B1", "B2", "B3", "B1", "B1", "B2")
)
controls <- c("c1", "c2", "c3")

test_that("without match columns every sample uses every control (old behaviour)", {
  sets <- control_sets_for_samples(meta, "sample", c("s1", "s3"), controls, NULL)
  expect_equal(sets, list(s1 = controls, s3 = controls))
})

test_that("a sample uses its own run's controls, else its batch's, else none", {
  logged <<- character(0)
  sets <- control_sets_for_samples(meta, "sample", c("s1", "s2", "s3", "s4"), controls,
                                   c("seq_run", "batch"))
  expect_equal(sets$s1, "c1")                 # R1 has its own control
  expect_equal(sets$s2, "c2")                 # R2 has its own control
  expect_equal(sets$s3, "c3")                 # R4 has none -> batch B2's control
  expect_equal(sets$s4, character(0))         # neither R9 nor B3 has a control
  expect_true(any(grepl("s4", logged)))       # unmatched samples are logged
})

test_that("a column with no metadata row keeps every control and is not logged", {
  logged <<- character(0)
  sets <- control_sets_for_samples(meta, "sample", c("s1", "s1_EXCLUSIVE_species"), controls,
                                   c("seq_run", "batch"))
  expect_equal(sets$s1_EXCLUSIVE_species, controls)
  expect_false(any(grepl("EXCLUSIVE", logged)))
})

test_that("a match column missing from the metadata is an error", {
  expect_error(control_sets_for_samples(meta, "sample", "s1", controls, "lane"),
               "control_match_cols.*lane")
})

test_that("a control missing from the metadata is an error when matching", {
  expect_error(control_sets_for_samples(meta, "sample", "s1", c(controls, "c9"), "seq_run"),
               "c9")
})

test_that("control_max_by_sample takes the row max over each sample's own controls", {
  long <- data.frame(taxid = c("t1", "t1", "t1", "t2"), sample = c("s1", "s2", "s4", "s1"),
                     count = c(10, 10, 10, 1),
                     c1 = c(5, 5, 5, 2), c2 = c(50, 50, 50, 0), c3 = c(0, 0, 0, 0))
  sets <- list(s1 = "c1", s2 = "c2", s4 = character(0))
  expect_equal(control_max_by_sample(long, sets), c(5, 50, -Inf, 2))
})
