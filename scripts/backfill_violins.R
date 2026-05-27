#!/usr/bin/env Rscript
# Re-run only the four stages that produce KW-annotated violin PNGs so the
# upper-right annotation tweak in R/06 + utils_ge_profile.R takes effect
# without re-running the full pipeline (DAA alone is 2.5 minutes).
#
# Usage:
#   Rscript scripts/backfill_violins.R projects/test_real_data/config.yaml

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/backfill_violins.R <config.yaml>")
cfg_path <- args[[1]]

repo_root <- dirname(dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))[[1]]
)))
mods  <- sort(list.files(file.path(repo_root, "R"),
                          pattern = "^\\d+_.*\\.R$",  full.names = TRUE))
utils <- sort(list.files(file.path(repo_root, "R"),
                          pattern = "^utils_.*\\.R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))
cfg <- load_config(cfg_path)

# R/06 only uses cleaned$noncontaminants when no precomputed diversity CSV
# exists; R/09/10/11 always read genetable_normdata.csv from disk. So the
# `cleaned` arg is unused in this backfill path.
inputs  <- load_inputs(cfg)
cleaned <- clean_data(inputs, cfg)

run_alpha_diversity(cleaned, cfg)
run_resistome(cleaned, cfg)
run_virulome(cleaned, cfg)
run_mobilome(cleaned, cfg)
cat("\nviolin backfill complete.\n")
