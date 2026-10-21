# Stage safety net: run_pipeline.R appends caught stage errors to
# <reports>/stage_errors.csv; R/14 stage_status() must report those stages as
# "failed" even when they wrote partial artifacts.
source(file.path("..", "..", "R", "00_setup.R"))
source(file.path("..", "..", "R", "14_manifest.R"))

test_that("read_stage_errors is empty without a file and reads written rows", {
  root <- file.path(tempdir(), "stage_err_test"); unlink(root, recursive = TRUE)
  cfg <- list(project_root = root, outputs = list(log_file = "Reports/pipeline.log"),
              stages = list(network = TRUE, taxonomy = TRUE))
  expect_equal(nrow(read_stage_errors(cfg)), 0L)
  dir.create(file.path(root, "Reports"), recursive = TRUE)
  utils::write.table(data.frame(stage = "network", error = "boom, with comma"),
                     stage_errors_file(cfg), sep = ",", row.names = FALSE,
                     qmethod = "double")
  errs <- read_stage_errors(cfg)
  expect_equal(errs$stage, "network"); expect_equal(errs$error, "boom, with comma")
  # failed beats a recorded duration (partial run); other stages unaffected
  expect_equal(stage_status(cfg, "network",  list(network = 3)), "failed")
  expect_equal(stage_status(cfg, "taxonomy", list(taxonomy = 3)), "complete")
})

test_that("upstream_failed flags GE/network stages when normalisation failed", {
  expect_equal(upstream_failed("resistome", "normalisation"), "normalisation")
  expect_equal(upstream_failed("virulome", data.frame(stage = "normalisation", error = "x")),
               "normalisation")
  for (s in c("mobilome", "network"))
    expect_equal(upstream_failed(s, c("taxonomy", "normalisation")), "normalisation")
  expect_length(upstream_failed("resistome", character(0)), 0)
  expect_length(upstream_failed("resistome", "taxonomy"), 0)
  expect_length(upstream_failed("taxonomy", "normalisation"), 0)
})
