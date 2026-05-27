#!/usr/bin/env Rscript
# Backfill the new pipeline PNG outputs (volcano per pair, degree distribution)
# against an existing test_run/ without re-executing the full pipeline.
# Reads the already-on-disk TSVs and topology CSV and writes the new PNGs to
# the same fig dirs the next R/08 / R/12 invocation would use.
#
# Usage:
#   Rscript scripts/backfill_new_pngs.R projects/test_real_data/config.yaml

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/backfill_new_pngs.R <config.yaml>")
cfg_path <- args[[1]]

repo_root <- dirname(dirname(normalizePath(
  sub("^--file=", "",
      grep("^--file=", commandArgs(), value = TRUE))
)))
source(file.path(repo_root, "R", "00_setup.R"))
cfg <- load_config(cfg_path)

`%||%` <- function(a, b) if (is.null(a)) b else a

ds_root  <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
fig_root <- file.path(cfg$project_root, cfg$outputs$figures_dir)

# ---- 1. Per-pair volcano PNGs (mirrors R/08 §5b) -----------------------------
alpha <- cfg$stats$alpha %||% 0.05
effect_thresh <- cfg$differential_abundance$effect_threshold %||% 0.8

for (level in c("taxa", "gene")) {
  fig_dir <- file.path(fig_root, "differential_abundance", level)
  ds_dir  <- file.path(ds_root,  "aldex2",                 level)
  if (!dir.exists(ds_dir)) next
  pair_tsvs <- list.files(ds_dir, pattern = "^[^.]+_vs_[^.]+\\.tsv$",
                           full.names = TRUE)
  for (pf in pair_tsvs) {
    cn <- tools::file_path_sans_ext(basename(pf))
    df <- readr::read_tsv(pf, show_col_types = FALSE)
    qcol <- intersect(c("wi.eBH", "we.eBH"), colnames(df))[1]
    if (is.na(qcol) || !"effect" %in% colnames(df)) next
    df$neglogq <- -log10(pmax(df[[qcol]], 1e-300))
    df$sig     <- df[[qcol]] < alpha
    p <- ggplot2::ggplot(df,
          ggplot2::aes(x = .data$effect, y = .data$neglogq,
                       color = .data$sig)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::geom_vline(xintercept = c(-effect_thresh, effect_thresh),
                          linetype = "dashed", colour = "grey60") +
      ggplot2::geom_hline(yintercept = -log10(alpha),
                          linetype = "dashed", colour = "red") +
      ggplot2::scale_color_manual(values = c(`TRUE` = "#d62728",
                                              `FALSE` = "#7f7f7f"),
                                   labels = c("ns", sprintf("q < %.2f", alpha))) +
      ggplot2::labs(
        title = sprintf("%s — effect vs -log10(q)", cn),
        x = "ALDEx2 effect (CLR diff)",
        y = sprintf("-log10(%s)", qcol),
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     text = ggplot2::element_text(size = 12))
    ggplot2::ggsave(file.path(fig_dir, paste0("volcano_", cn, ".png")),
                    p, width = 7, height = 5, dpi = 300)
    cat(sprintf("wrote %s/volcano_%s.png\n", level, cn))
  }
}

# ---- 2. Network degree distribution PNG (mirrors R/12) -----------------------
topo_csv <- file.path(ds_root, "network", "topology.csv")
if (file.exists(topo_csv)) {
  topo <- readr::read_csv(topo_csv, show_col_types = FALSE)
  top_n <- (cfg$network$degree_distribution_top_n %||% 50)
  topo <- topo[order(-topo$degree), , drop = FALSE]
  topo <- utils::head(topo, top_n)
  dd <- ggplot2::ggplot(topo,
        ggplot2::aes(x = stats::reorder(.data$node, -.data$degree),
                     y = .data$degree, fill = .data$kind)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "Degree",
                  title = sprintf("Top %d nodes by degree", nrow(topo)),
                  fill = "Node kind") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      text        = ggplot2::element_text(size = 12)
    )
  ggplot2::ggsave(file.path(fig_root, "network", "degree_distribution.png"),
                  dd, width = max(8, 0.18 * nrow(topo) + 4),
                  height = 5.5, dpi = 200, bg = "white")
  cat("wrote network/degree_distribution.png\n")
}
