# One-off harness to re-run only R/04_taxonomy.R against an already-cleaned
# test_run/ state. Runs load_inputs + clean_data (cheap) and then
# run_taxonomy in isolation so the Venn / heatmap / etc. can be re-rendered
# quickly while iterating on cfg$taxonomy$* knobs (e.g. heatmap_fontsize_*).
#
# Usage:
#   Rscript scripts/taxonomy_only.R [config.yaml]
#   (default: projects/chicken_batch1/config.yaml)

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[[1]] else "projects/chicken_batch1/config.yaml"
stopifnot(file.exists(cfg_path))

`%||%` <- function(a, b) if (is.null(a)) b else a

mods  <- sort(list.files("R", pattern = "^\\d+_.*\\.R$", full.names = TRUE))
utils <- sort(list.files("R", pattern = "^utils_.*\\.R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))

cfg <- load_config(cfg_path)

t0      <- Sys.time()
inputs  <- load_inputs(cfg)
t_load  <- Sys.time()
cleaned <- clean_data(inputs, cfg)
t_clean <- Sys.time()
run_taxonomy(cleaned, cfg)
t_tax   <- Sys.time()

dur <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
cat(sprintf(
  "\ntaxonomy_only timing:\n  load_inputs %6.1fs\n  clean_data  %6.1fs\n  taxonomy    %6.1fs\n  TOTAL       %6.1fs\n",
  dur(t0, t_load), dur(t_load, t_clean),
  dur(t_clean, t_tax), dur(t0, t_tax)
))
