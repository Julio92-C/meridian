#!/usr/bin/env Rscript
# Re-render the sankey PNG from a saved sankey_long.csv. Lets you iterate
# on .network_render_sankey_png without re-running the full pipeline.
# Usage: Rscript scripts/rerender_sankey.R <config.yaml>

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/rerender_sankey.R <config.yaml>")
cfg_path <- args[[1]]

for (f in sort(list.files("R", pattern = "^\\d+_.*\\.R$|^utils_.*\\.R$",
                          full.names = TRUE))) source(f)

cfg <- load_config(cfg_path)
csv <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                 "network/sankey_long.csv")
df  <- readr::read_csv(csv, show_col_types = FALSE)
cat(sprintf("Loaded %s: %d rows, %d samples, %d genes\n",
            csv, nrow(df), dplyr::n_distinct(df$sample),
            dplyr::n_distinct(df$GENE)))

fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "network")
scfg    <- (cfg$network %||% list())$sankey %||% list()
.network_render_sankey_png(df, file.path(fig_dir, "sankey_overall.png"),
                           scfg, cfg)
cat("Rendered: ", file.path(fig_dir, "sankey_overall.png"), "\n", sep = "")
