# 07_beta_diversity.R — Bray-Curtis → PCoA + PERMANOVA (adonis2, 9999 perms).
# Matches the Chicken batch 1 paper methodology.

run_beta_diversity <- function(cleaned, cfg) {
  pipeline_log(cfg, "Beta diversity (PCoA + PERMANOVA)")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "beta_diversity")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid   <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]

  counts <- cleaned$noncontaminants |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)
  mat <- as.matrix(counts[, -1])
  rownames(mat) <- counts$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  bc <- vegan::vegdist(mat, method = "bray")
  pcoa <- ape::pcoa(bc)
  scores <- as.data.frame(pcoa$vectors[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$sample <- rownames(scores)
  scores <- dplyr::inner_join(scores, meta, by = c("sample" = sid))

  meta_for_adonis <- meta[match(rownames(mat), meta[[sid]]), , drop = FALSE]
  permanova <- vegan::adonis2(
    reformulate(group, "bc"),
    data = meta_for_adonis,
    permutations = cfg$stats$permanova_permutations %||% 9999
  )
  capture.output(permanova,
                 file = file.path(fig_dir, "permanova.txt"))
  pipeline_log(cfg, sprintf("PERMANOVA %s R2 = %.3f, p = %.4g",
                            group, permanova$R2[1], permanova$`Pr(>F)`[1]))

  pct <- round(pcoa$values$Relative_eig[1:2] * 100, 1)
  p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(level = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", pct[1]),
      y = sprintf("PC2 (%.1f%%)", pct[2]),
      title = sprintf("Bray–Curtis PCoA — %s", group)
    )
  ggplot2::ggsave(file.path(fig_dir, "pcoa.png"), p,
                  width = 6, height = 5, dpi = 300)

  invisible(list(scores = scores, permanova = permanova))
}
