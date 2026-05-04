#!/usr/bin/env Rscript
# run_pipeline.R — master runner.
# Usage:
#   Rscript run_pipeline.R projects/chicken_batch1/config.yaml

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_pipeline.R <config.yaml>")
cfg_path <- args[[1]]

ofile    <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
script   <- ofile %||% (if (length(file_arg)) file_arg[[1]] else NULL)
repo_root <- if (!is.null(script)) dirname(normalizePath(script)) else getwd()

# Source every R/NN_*.R module in order.
mods <- sort(list.files(file.path(repo_root, "R"),
                        pattern = "^\\d+_.*\\.R$",
                        full.names = TRUE))
invisible(lapply(mods, source))

cfg <- load_config(cfg_path)
pipeline_log(cfg, sprintf("Starting pipeline for study '%s'", cfg$study$id))

# Build the ordered list of stages that will actually run, so we can show
# "[i/n  pct%] stage_name" progress in the log as each one starts/finishes.
enabled_names <- c(
  "load_inputs",
  if (isTRUE(cfg$stages$clean_data))             "clean_data",
  if (isTRUE(cfg$stages$normalisation))          "normalisation",
  if (isTRUE(cfg$stages$taxonomy))               "taxonomy",
  if (isTRUE(cfg$stages$relative_abundance))     "relative_abundance",
  if (isTRUE(cfg$stages$alpha_diversity))        "alpha_diversity",
  if (isTRUE(cfg$stages$beta_diversity))         "beta_diversity",
  if (isTRUE(cfg$stages$differential_abundance)) "differential_abundance",
  if (isTRUE(cfg$stages$resistome))              "resistome",
  if (isTRUE(cfg$stages$virulome))               "virulome",
  if (isTRUE(cfg$stages$mobilome))               "mobilome",
  if (isTRUE(cfg$stages$network))                "network",
  if (isTRUE(cfg$stages$report))                 "report"
)
total_stages <- length(enabled_names)
stage_idx    <- 0L
pipeline_log(cfg, sprintf("Plan: %d stage(s) — %s",
                          total_stages, paste(enabled_names, collapse = ", ")))

run_stage <- function(name, fn) {
  if (!(name %in% enabled_names)) return(invisible(NULL))
  stage_idx <<- stage_idx + 1L
  pct <- round(100 * (stage_idx - 1) / total_stages)
  pipeline_log(cfg, sprintf("[%d/%d %3d%%] %s — start",
                            stage_idx, total_stages, pct, name))
  t0 <- Sys.time()
  result <- fn()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  pipeline_log(cfg, sprintf("[%d/%d %3d%%] %s — done in %.1fs",
                            stage_idx, total_stages, pct, name, elapsed))
  result
}

inputs  <- run_stage("load_inputs", function() load_inputs(cfg))
cleaned <- run_stage("clean_data",  function() clean_data(inputs, cfg))

if (!is.null(cleaned)) {
  run_stage("normalisation",          function() normalise_data(cleaned, cfg))
  run_stage("taxonomy",               function() run_taxonomy(cleaned, cfg))
  run_stage("relative_abundance",     function() run_relative_abundance(cleaned, cfg))
  run_stage("alpha_diversity",        function() run_alpha_diversity(cleaned, cfg))
  run_stage("beta_diversity",         function() run_beta_diversity(cleaned, cfg))
  run_stage("differential_abundance", function() run_aldex(cleaned, cfg))
  run_stage("resistome",              function() run_resistome(cleaned, cfg))
  run_stage("virulome",               function() run_virulome(cleaned, cfg))
  run_stage("mobilome",               function() run_mobilome(cleaned, cfg))
  run_stage("network",                function() run_network(cleaned, cfg))
}

run_stage("report", function() {
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
})

pipeline_log(cfg, sprintf("Pipeline finished — %d/%d stages complete",
                          stage_idx, total_stages))
