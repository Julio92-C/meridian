#!/usr/bin/env Rscript
# run_pipeline.R — master runner.
# Usage:
#   Rscript run_pipeline.R projects/chicken_batch1/config.yaml

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_pipeline.R <config.yaml>")
cfg_path <- args[[1]]

repo_root <- dirname(normalizePath(sys.frame(1)$ofile %||% sub("--file=", "",
  grep("--file=", commandArgs(), value = TRUE))))
if (length(repo_root) == 0 || !nzchar(repo_root)) {
  repo_root <- getwd()
}

# Source every R/NN_*.R module in order.
`%||%` <- function(a, b) if (is.null(a)) b else a
mods <- sort(list.files(file.path(repo_root, "R"),
                        pattern = "^\\d+_.*\\.R$",
                        full.names = TRUE))
invisible(lapply(mods, source))

cfg <- load_config(cfg_path)
pipeline_log(cfg, sprintf("Starting pipeline for study '%s'", cfg$study$id))

inputs <- load_inputs(cfg)

cleaned <- if (isTRUE(cfg$stages$clean_data))
  clean_data(inputs, cfg) else NULL

if (isTRUE(cfg$stages$normalisation) && !is.null(cleaned))
  normalise_data(cleaned, cfg)

if (isTRUE(cfg$stages$taxonomy) && !is.null(cleaned))
  run_taxonomy(cleaned, cfg)

if (isTRUE(cfg$stages$relative_abundance) && !is.null(cleaned))
  run_relative_abundance(cleaned, cfg)

if (isTRUE(cfg$stages$alpha_diversity) && !is.null(cleaned))
  run_alpha_diversity(cleaned, cfg)

if (isTRUE(cfg$stages$beta_diversity) && !is.null(cleaned))
  run_beta_diversity(cleaned, cfg)

if (isTRUE(cfg$stages$differential_abundance) && !is.null(cleaned))
  run_aldex(cleaned, cfg)

if (isTRUE(cfg$stages$resistome) && !is.null(cleaned))
  run_resistome(cleaned, cfg)

if (isTRUE(cfg$stages$virulome) && !is.null(cleaned))
  run_virulome(cleaned, cfg)

if (isTRUE(cfg$stages$mobilome) && !is.null(cleaned))
  run_mobilome(cleaned, cfg)

if (isTRUE(cfg$stages$network) && !is.null(cleaned))
  run_network(cleaned, cfg)

if (isTRUE(cfg$stages$report)) {
  report_src <- file.path(repo_root, "templates", "report.qmd")
  report_out <- file.path(cfg$project_root, cfg$outputs$report_html)
  if (file.exists(report_src) && nzchar(Sys.which("quarto"))) {
    pipeline_log(cfg, "Rendering Quarto HTML report")
    system2("quarto", c("render", report_src,
                        "--to", "html",
                        "-P", paste0("config=", cfg_path),
                        "--output", basename(report_out)))
    out_tmp <- file.path(dirname(report_src), basename(report_out))
    if (file.exists(out_tmp)) file.rename(out_tmp, report_out)
  } else {
    pipeline_log(cfg, "Quarto not found or report.qmd missing — skipping HTML")
  }
}

pipeline_log(cfg, "Pipeline finished")
