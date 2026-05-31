#!/usr/bin/env Rscript
# reemit_manifests.R — re-emit manifest.json for one or more existing runs
# without rerunning the full pipeline. Useful when iterating on R/13_manifest.R
# or after pipeline output paths change.
#
# Usage:
#   Rscript scripts/reemit_manifests.R                   # both reference configs
#   Rscript scripts/reemit_manifests.R projects/foo/config.yaml [more...]
#
# The writer falls back to filesystem evidence when stage_times is empty
# (see stage_status() in R/13_manifest.R), so re-emitting works as long as
# the per-stage primary artifacts from the previous run are still on disk.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- c("projects/chicken_batch1/config.yaml",
            "projects/chicken_batch2/config.yaml")
}

mods  <- sort(list.files("R", pattern = "^\\d+_.*\\.R$", full.names = TRUE))
utils <- sort(list.files("R", pattern = "^utils_.*\\.R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))

for (cfg_path in args) {
  cat("---", cfg_path, "---\n")
  cfg <- load_config(cfg_path)
  out <- write_manifest(cfg, stage_times = list(), pipeline_repo_root = ".")
  cat("wrote:", out, "\n\n")
}
