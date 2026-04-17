# 04_taxonomy.R — Venn diagrams, pHeatmap and per-level taxonomy summaries.
# Paper methodology: taxa richness + Shannon diversity underpin these plots;
# Venn diagrams are reported across treatment groups.

run_taxonomy <- function(cleaned, cfg) {
  pipeline_log(cfg, "Taxonomy (Venn + heatmap)")
  df <- cleaned$noncontaminants
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "taxonomy")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  group <- cfg$metadata$group_cols[[1]]
  sid   <- cfg$metadata$sample_id_col
  meta  <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )

  joined <- dplyr::inner_join(df, meta, by = c("sample" = sid))

  # Venn of taxa shared across groups -------------------------------------
  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    taxa_sets <- split(joined$name, joined[[group]])
    taxa_sets <- lapply(taxa_sets, unique)
    venn_path <- file.path(fig_dir, "taxa_venn.png")
    VennDiagram::venn.diagram(
      x = taxa_sets, filename = venn_path,
      fill = paletteer::paletteer_d("ggsci::default_nejm")[seq_along(taxa_sets)],
      main = "Shared taxa across treatments"
    )
  }

  # pHeatmap of top-N taxa -----------------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    top_taxa <- joined |>
      dplyr::group_by(name) |>
      dplyr::summarise(total = sum(count, na.rm = TRUE)) |>
      dplyr::slice_max(total, n = 50) |>
      dplyr::pull(name)

    mat <- joined |>
      dplyr::filter(name %in% top_taxa) |>
      dplyr::group_by(sample, name) |>
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0) |>
      tibble::column_to_rownames("name") |>
      as.matrix()

    pheatmap::pheatmap(
      log1p(mat),
      filename = file.path(fig_dir, "taxa_heatmap.png"),
      show_rownames = TRUE, show_colnames = TRUE,
      fontsize_row = 6, fontsize_col = 7
    )
  }

  invisible(NULL)
}
