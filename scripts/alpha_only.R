# One-off harness to time just R/06_alpha_diversity.R on a given config.
# Runs load_inputs + clean_data (cheap) and then run_alpha_diversity in
# isolation so a single stage's perf can be measured without sitting
# through the rest of the pipeline.
#
# Usage:
#   Rscript scripts/alpha_only.R [config.yaml]
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
run_alpha_diversity(cleaned, cfg)
t_alpha <- Sys.time()

dur <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
cat(sprintf(
  "\nalpha_only timing:\n  load_inputs    %6.1fs\n  clean_data     %6.1fs\n  alpha_diversity %6.1fs\n  TOTAL          %6.1fs\n",
  dur(t0, t_load), dur(t_load, t_clean),
  dur(t_clean, t_alpha), dur(t0, t_alpha)
))
