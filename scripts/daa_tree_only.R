# One-off harness: re-render only the DAA tree (R/08_daa_tree.R) from the
# ALDEx2 taxa tables already in <datasets_dir>/aldex2/taxa/. Runs load_inputs +
# clean_data for the lineage, then run_daa_tree.
#
# Usage:
#   Rscript scripts/daa_tree_only.R [config.yaml]

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[[1]] else "projects/chicken_batch1/config.yaml"
stopifnot(file.exists(cfg_path))

`%||%` <- function(a, b) if (is.null(a)) b else a

mods  <- sort(list.files("R", pattern = "^[0-9]+_.*[.]R$", full.names = TRUE))
utils <- sort(list.files("R", pattern = "^utils_.*[.]R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))

cfg     <- load_config(cfg_path)
inputs  <- load_inputs(cfg)
cleaned <- clean_data(inputs, cfg)
run_daa_tree(cleaned, cfg)
