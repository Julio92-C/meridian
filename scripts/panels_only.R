# Re-run only the panels stage against an already-populated <project>/test_run/.
# Skips every upstream stage — relies on the .rds siblings written by
# save_panel_ggplot() and on the per-stage figure scanners in R/14_manifest.R
# (which build_stages_index() walks). Cuts a ~4 minute full-pipeline iteration
# down to ~45s while tuning composite theme overrides.
#
# Usage:
#   Rscript scripts/panels_only.R [config.yaml]
#   (default: projects/chicken_batch1/config.yaml)

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[[1]] else "projects/chicken_batch1/config.yaml"
stopifnot(file.exists(cfg_path))

`%||%` <- function(a, b) if (is.null(a)) b else a

mods <- sort(list.files("R", pattern = "^\\d+_.*\\.R$", full.names = TRUE))
utils <- sort(list.files("R", pattern = "^utils_.*\\.R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))

cfg <- load_config(cfg_path)
t0 <- Sys.time()
run_panels(cfg)
cat(sprintf("\npanels_only: %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
