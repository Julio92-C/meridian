# 11_mobilome.R — MGE profile analysis. Uses PlasmidFinder rows from ABRicate.

run_mobilome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Mobilome (PlasmidFinder)")
  .run_ge_block(cleaned, cfg, db = "plasmidfinder", label = "MGE",
                out_subdir = "mobilome")
}

# Shared implementation for resistome / virulome / mobilome — the three
# blocks differ only by which ABRicate DATABASE column value they filter on
# and by the figure subdirectory.
.run_ge_block <- function(cleaned, cfg, db, label, out_subdir) {
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, out_subdir)
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, out_subdir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid   <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]

  df <- cleaned$abri_kraken2
  df <- df[tolower(df$DATABASE) == db, , drop = FALSE]
  if (nrow(df) == 0) {
    pipeline_log(cfg, sprintf("No %s hits — skipping %s", db, out_subdir))
    return(invisible(NULL))
  }
  df <- dplyr::inner_join(df, meta, by = c("sample" = sid))

  # Richness + Shannon per sample on GE counts ---------------------------
  gene_mat <- df |>
    dplyr::group_by(sample, GENE) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = GENE, values_from = n, values_fill = 0)
  mat <- as.matrix(gene_mat[, -1])
  rownames(mat) <- gene_mat$sample

  alpha <- data.frame(
    sample   = rownames(mat),
    richness = rowSums(mat > 0),
    shannon  = vegan::diversity(mat, "shannon")
  )
  alpha <- dplyr::inner_join(alpha, meta, by = c("sample" = sid))
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  kw <- kruskal.test(reformulate(group, "shannon"), data = alpha)
  pipeline_log(cfg, sprintf("%s Shannon ~ %s KW p = %.4g", label, group, kw$p.value))

  # Beta + PERMANOVA on GE counts ----------------------------------------
  bc <- vegan::vegdist(mat, method = "bray")
  sample_meta <- meta[match(rownames(mat), meta[[sid]]), , drop = FALSE]
  permanova <- vegan::adonis2(
    reformulate(group, "bc"),
    data = sample_meta,
    permutations = cfg$stats$permanova_permutations %||% 9999
  )
  capture.output(permanova, file = file.path(fig_dir, "permanova.txt"))

  pcoa <- ape::pcoa(bc)
  scores <- as.data.frame(pcoa$vectors[, 1:2]); colnames(scores) <- c("PC1","PC2")
  scores$sample <- rownames(scores)
  scores <- dplyr::inner_join(scores, meta, by = c("sample" = sid))
  p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) + ggplot2::stat_ellipse() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = sprintf("%s profile — PCoA", label))
  ggplot2::ggsave(file.path(fig_dir, "pcoa.png"), p, width = 6, height = 5, dpi = 300)

  # Venn of genes shared across groups -----------------------------------
  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    gene_sets <- split(df$GENE, df[[group]])
    gene_sets <- lapply(gene_sets, unique)
    VennDiagram::venn.diagram(
      x = gene_sets,
      filename = file.path(fig_dir, "venn.png"),
      fill = paletteer::paletteer_d("ggsci::default_nejm")[seq_along(gene_sets)],
      main = sprintf("%s shared across %s groups", label, group)
    )
  }

  # pHeatmap -------------------------------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && ncol(mat) > 1) {
    pheatmap::pheatmap(
      log1p(t(mat)),
      filename = file.path(fig_dir, "heatmap.png"),
      fontsize_row = 6, fontsize_col = 7
    )
  }

  invisible(list(alpha = alpha, permanova = permanova, scores = scores))
}
