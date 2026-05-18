# 04_taxonomy.R — Venn diagram + annotated heatmap of taxa across treatment
# groups. Mirrors the methodology used in kraken2_VennDiagram.R and
# taxaNorm_pheatmap.R from the reference chicken_batch1 scripts, but driven
# entirely by config (no hardcoded sample IDs, treatments, or colour maps).
#
# Steps
#   1. Apply shared taxa-name cleanup (clean_taxa_names from utils_taxa.R):
#      drop clade-group noise, keep species-level only, strip NCBI brackets,
#      abbreviate genus, truncate long names.
#   2. Join with metadata so we know which treatment each sample belongs to.
#   3. Build a (taxon x group) presence matrix by summing counts within
#      each group.
#   4. Venn diagram of the presence matrix (skipped for <2 or >5 groups,
#      since VennDiagram only supports 2-5 sets).
#   5. Heatmap of relative abundance for the top-N taxa (by mean % across
#      samples), with Treatment column annotation and Core / Accessory /
#      Unique row annotation derived from across-group occupancy.

# Build a named colour vector for a categorical variable.
# Falls back to a paletteer palette if the user has not supplied an explicit
# mapping under cfg$taxonomy[[key]].
group_palette <- function(levels, cfg, key = "treatment_colors") {
  user_map <- cfg$taxonomy[[key]]
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- tryCatch(
    as.character(paletteer::paletteer_d("ggsci::default_nejm")),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(length(levels), palette = "Dark 3")
  }
  pal <- rep_len(pal, length(levels))
  names(pal) <- levels
  pal
}

run_taxonomy <- function(cleaned, cfg) {
  pipeline_log(cfg, "Taxonomy (Venn + heatmap)")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "taxonomy")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  df <- clean_taxa_names(cleaned$noncontaminants, cfg)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Taxonomy: no rows left after clean_taxa_names — skipping")
    return(invisible(NULL))
  }

  group <- cfg$metadata$group_cols[[1]]
  sid   <- cfg$metadata$sample_id_col
  meta  <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  joined <- dplyr::inner_join(df, meta, by = c("sample" = sid))
  if (!group %in% colnames(joined)) {
    stop(sprintf("Taxonomy: group column '%s' not found in metadata", group))
  }

  # Per-(taxon, sample) count threshold. Evaluated row-by-row in long form,
  # so a taxon can be kept in samples where it crosses the threshold and
  # dropped in samples where it does not. Mirrors the `filter(count > 5)`
  # step in GT Bracken_VennDiagram.R; floor configurable via
  # cfg$taxonomy$min_count_per_sample (default 5).
  min_count <- cfg$taxonomy$min_count_per_sample %||% 5
  n_before  <- nrow(joined)
  joined    <- dplyr::filter(joined, count >= min_count)
  pipeline_log(cfg, sprintf(
    "Taxonomy: %d/%d (taxon, sample) rows kept with count >= %d",
    nrow(joined), n_before, min_count
  ))
  if (nrow(joined) == 0) {
    pipeline_log(cfg, "Taxonomy: nothing left after count filter — skipping")
    return(invisible(NULL))
  }

  # ---- Group-level (taxon x group) presence matrix -----------------------
  group_mat <- joined |>
    dplyr::group_by(name, .data[[group]]) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(group),
      values_from = count,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("name") |>
    as.matrix()
  presence_mat <- (group_mat > 0) + 0L

  # ---- Venn diagram (group-level) ---------------------------------------
  n_grp <- ncol(presence_mat)
  if (n_grp < 2) {
    pipeline_log(cfg, "Taxonomy: <2 groups — Venn skipped")
  } else if (n_grp > 5) {
    pipeline_log(cfg, sprintf(
      "Taxonomy: %d groups — Venn supports 2-5 only, skipped", n_grp
    ))
  } else if (requireNamespace("VennDiagram", quietly = TRUE)) {
    sets <- lapply(seq_len(n_grp), function(j) {
      rownames(presence_mat)[presence_mat[, j] == 1]
    })
    names(sets) <- colnames(presence_mat)

    # Optional inspection aid: dump the per-group species list to the log.
    # Off by default — set cfg$taxonomy$log_venn_sets: true to re-enable.
    if (isTRUE(cfg$taxonomy$log_venn_sets)) {
      species_lines <- unlist(lapply(names(sets), function(g) {
        c(sprintf("  %s (%d species):", g, length(sets[[g]])),
          paste0("    - ", sets[[g]]))
      }))
      pipeline_log(cfg, sprintf(
        "Taxonomy Venn sets:\n%s",
        paste(species_lines, collapse = "\n")
      ))
    }

    fill_pal <- group_palette(names(sets), cfg)
    venn_path <- file.path(fig_dir, "taxa_venn.png")
    VennDiagram::venn.diagram(
      x = sets,
      filename = venn_path,
      imagetype = "png",
      height = 2400,
      width = 2400,
      resolution = 300,
      fill = unname(fill_pal),
      alpha = 0.5,
      cat.cex = 1.3,
      cat.fontface = "bold",
      cex = 1.5,
      main = "Shared species across treatment groups",
      main.cex = 1.4,
      margin = 0.08,
      disable.logging = TRUE
    )
  }

  # ---- Heatmap (relative abundance, top-N by mean %) --------------------
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    pipeline_log(cfg, "Taxonomy: pheatmap not available — heatmap skipped")
    return(invisible(NULL))
  }

  top_n <- cfg$taxonomy$top_n %||% 50

  sample_pct <- joined |>
    dplyr::group_by(sample) |>
    dplyr::mutate(pct = count / sum(count, na.rm = TRUE) * 100) |>
    dplyr::ungroup()

  top_taxa <- sample_pct |>
    dplyr::group_by(name) |>
    dplyr::summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop") |>
    dplyr::slice_max(mean_pct, n = top_n) |>
    dplyr::pull(name)

  heat_mat <- sample_pct |>
    dplyr::filter(name %in% top_taxa) |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(pct = sum(pct, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample, values_from = pct, values_fill = 0) |>
    tibble::column_to_rownames("name") |>
    as.matrix()
  heat_mat <- heat_mat[order(rowMeans(heat_mat), decreasing = TRUE), , drop = FALSE]

  ann_col <- joined |>
    dplyr::select(sample, dplyr::all_of(group)) |>
    dplyr::distinct() |>
    tibble::column_to_rownames("sample")
  ann_col <- ann_col[colnames(heat_mat), , drop = FALSE]

  # Order samples by treatment so all replicates of a group sit together,
  # then compute the between-group gap positions for pheatmap.
  ann_col  <- ann_col[order(ann_col[[group]]), , drop = FALSE]
  heat_mat <- heat_mat[, rownames(ann_col), drop = FALSE]
  grp_run  <- as.character(ann_col[[group]])
  grp_lvls <- unique(grp_run)
  grp_sizes <- as.integer(table(grp_run)[grp_lvls])
  gaps_col <- if (length(grp_sizes) > 1) head(cumsum(grp_sizes), -1L) else NULL

  ann_colors <- list()
  ann_colors[[group]] <- group_palette(grp_lvls, cfg)

  # Size the canvas to the matrix so labels stay readable as top_n / sample
  # count change. Inches; pheatmap renders at 72 dpi by default so these
  # produce roughly 900x900 to 1200x1200 PNGs at typical settings.
  hm_width  <- max(10, 0.45 * ncol(heat_mat) + 4)
  hm_height <- max(8,  0.22 * nrow(heat_mat) + 2)

  pheatmap::pheatmap(
    heat_mat,
    filename = file.path(fig_dir, "taxa_heatmap.png"),
    width  = hm_width,
    height = hm_height,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    gaps_col = gaps_col,
    color = grDevices::colorRampPalette(c("#15b379", "yellow", "#f2615a"))(100),
    border_color = NA,
    fontsize_row = 9,
    fontsize_col = 10,
    fontsize = 10
  )

  invisible(NULL)
}
