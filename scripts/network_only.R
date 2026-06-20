# One-off harness to re-run just R/12_network.R on a given config.
# Mirrors scripts/alpha_only.R: sources every R/NN_*.R + R/utils_*.R,
# then runs load_inputs -> clean_data -> run_network in isolation.
#
# Usage:
#   Rscript scripts/network_only.R [config.yaml]
#   (default: projects/hospital_microbiome/config.yaml)

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[[1]] else "projects/hospital_microbiome/config.yaml"
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
run_network(cleaned, cfg)
t_net   <- Sys.time()

dur <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
cat(sprintf(
  "\nnetwork_only timing:\n  load_inputs %6.1fs\n  clean_data  %6.1fs\n  network     %6.1fs\n  TOTAL       %6.1fs\n",
  dur(t0, t_load), dur(t_load, t_clean), dur(t_clean, t_net), dur(t0, t_net)
))
