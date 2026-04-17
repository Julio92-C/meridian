# 00_setup.R — load libraries, resolve paths, helper utilities.
# Sourced first by run_pipeline.R; all other R/ modules assume these are loaded.

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(paletteer)
})

#' Resolve a path relative to the project_root defined in the config.
resolve_path <- function(cfg, key) {
  path <- cfg$inputs[[key]]
  if (is.null(path)) path <- cfg$outputs[[key]]
  file.path(cfg$project_root, path)
}

#' Log a step to both stdout and the run log file.
pipeline_log <- function(cfg, msg) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line  <- sprintf("[%s] %s", stamp, msg)
  message(line)
  log_path <- file.path(cfg$project_root, cfg$outputs$log_file)
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  cat(line, "\n", file = log_path, append = TRUE)
}

#' Read and validate the study config.
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  stopifnot(
    !is.null(cfg$project_root),
    dir.exists(cfg$project_root),
    !is.null(cfg$metadata$file),
    !is.null(cfg$metadata$sample_id_col)
  )
  cfg
}
