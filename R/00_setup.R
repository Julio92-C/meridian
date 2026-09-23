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
  # Runtime override: MERIDIAN_PROJECT_ROOT wins over the config's host path,
  # so the SAME config.yaml runs unchanged natively and in-container. With the
  # env var unset, native behaviour is identical (backward-compatible).
  env_root <- Sys.getenv("MERIDIAN_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) cfg$project_root <- env_root
  stopifnot(
    !is.null(cfg$project_root),
    dir.exists(cfg$project_root),
    !is.null(cfg$metadata$file),
    !is.null(cfg$metadata$sample_id_col)
  )
  cfg
}

# Allowed values for cfg$stats$padjust_method. Matches `stats::p.adjust`
# methods plus "none" (no correction). BH is the field standard for
# omics-style multiple testing; BY is the more conservative variant for
# arbitrary-dependence FDR control.
PADJUST_METHODS_ALLOWED <- c(
  "BH", "fdr",            # Benjamini-Hochberg (PRDS dependence)
  "BY",                   # Benjamini-Yekutieli (arbitrary dependence)
  "bonferroni",           # FWER, most conservative
  "holm", "hochberg",     # FWER step-down / step-up
  "hommel",
  "none"                  # explicit opt-out
)

#' Resolve the configured p-value adjustment method. Defaults to "BH" when
#' unset so legacy configs continue to behave the same. Errors fast on an
#' unknown method name to avoid silent typos.
padjust_method <- function(cfg) {
  m <- cfg$stats$padjust_method %||% "BH"
  if (!m %in% PADJUST_METHODS_ALLOWED) {
    stop(sprintf(
      "cfg$stats$padjust_method = '%s' is not allowed (choose: %s)",
      m, paste(PADJUST_METHODS_ALLOWED, collapse = ", ")
    ))
  }
  m
}

#' Convenience wrapper: adjust a vector of raw p-values using the
#' configured method. NA-tolerant. Use this at every call site instead of
#' a bare `stats::p.adjust(p, method = "BH")` so the cfg knob actually
#' threads through.
padjust_p <- function(p, cfg) {
  m <- padjust_method(cfg)
  if (m == "none") return(p)
  stats::p.adjust(p, method = m)
}

#' Cross-module stats accumulator for families that span more than one
#' R/ module (e.g. resistome / virulome / mobilome each run their own
#' alpha-diversity KW per metric). Each module appends a row with
#' `padj_summary_record`; a finaliser reads the accumulated CSV, groups by
#' `family`, family-adjusts within each group with the configured method,
#' and writes a `<datasets>/stats/<surface>_padj_summary.csv`.
.padj_summary_path <- function(cfg, surface) {
  file.path(cfg$project_root, cfg$outputs$datasets_dir, "stats",
            sprintf("%s_raw.csv", surface))
}

#' Clear any prior accumulator for `surface`. Call once at the top of a
#' run (run_pipeline.R does this) so re-runs don't append to stale rows.
padj_summary_reset <- function(cfg, surface) {
  p <- .padj_summary_path(cfg, surface)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(p)) file.remove(p)
  invisible(p)
}

#' Append one (surface, family, key, raw_p) row to the accumulator. Safe
#' to call repeatedly; creates the CSV with a header on first call.
padj_summary_record <- function(cfg, surface, family, key, raw_p) {
  p <- .padj_summary_path(cfg, surface)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(family = as.character(family),
                    key    = as.character(key),
                    raw_p  = as.numeric(raw_p),
                    stringsAsFactors = FALSE)
  if (!file.exists(p)) {
    utils::write.table(row, p, sep = ",", row.names = FALSE,
                       col.names = TRUE)
  } else {
    utils::write.table(row, p, sep = ",", row.names = FALSE,
                       col.names = FALSE, append = TRUE)
  }
  invisible(p)
}

#' Read the accumulator, family-adjust within each `family` group, write
#' a summary CSV. Returns the summary tibble (invisibly) so callers can
#' surface it directly if needed. No-op when the accumulator is missing
#' (e.g. surface wasn't exercised this run).
padj_summary_finalise <- function(cfg, surface) {
  in_path  <- .padj_summary_path(cfg, surface)
  if (!file.exists(in_path)) return(invisible(NULL))
  raw      <- utils::read.csv(in_path, stringsAsFactors = FALSE)
  if (nrow(raw) == 0) return(invisible(NULL))
  method   <- padjust_method(cfg)
  raw$p_adj <- NA_real_
  for (f in unique(raw$family)) {
    idx <- raw$family == f
    raw$p_adj[idx] <- padjust_p(raw$raw_p[idx], cfg)
  }
  raw$method <- method
  out_path <- sub("_raw\\.csv$", "_padj_summary.csv", in_path)
  readr::write_csv(raw, out_path)
  invisible(raw)
}

# `%||%` is needed by padjust_method but the rest of the pipeline only
# defines it inside individual modules. Define it here too so 00_setup is
# self-contained.
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}
