# 05_relative_abundance.R — stacked bar plot of relative abundance per sample
# (top-N taxa, remainder collapsed into "Others"), plus a horizontal total-
# count bar for the same taxa. Mirrors relativeAbundance.R and
# bracken_relativeAbundance.R from the reference chicken_batch1 scripts but
# parameterised by config — no hardcoded sample IDs, treatments, or palettes.
#
# Steps
#   1. Apply shared taxa-name cleanup (clean_taxa_names from utils_taxa.R):
#      drop clade-group noise, keep species-level only, strip NCBI brackets,
#      abbreviate genus, truncate long names.
#   2. Apply a per-(taxon, sample) count floor.
#   3. Optionally join metadata so the stacked bar can be facetted by group.
#   4. Compute per-sample relative abundance %.
#   5. Identify the top-N taxa by total count; collapse everything else into
#      "Others" so the palette only has to cover top_n + 1 levels.
#   6. Stacked bar plot (PNG + plotly HTML) and horizontal total-count plot
#      (PNG, log10 scale with median line and count labels).

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build a palette sized to exactly `n` levels. Caller picks the palette name
# via cfg$relative_abundance$palette; falls back to ggsci::default_igv (matches
# GT bracken_relativeAbundance.R) and finally to hcl.colors if paletteer fails.
relab_palette <- function(n, cfg) {
  name <- cfg$relative_abundance$palette %||% "ggsci::default_igv"
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name, n = n)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    # paletteer_d with n > length(palette) errors on some palettes; cycle a
    # broader default rather than failing the stage.
    pal <- tryCatch(
      as.character(paletteer::paletteer_d(name)),
      error = function(e) NULL
    )
  }
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

run_relative_abundance <- function(cleaned, cfg) {
  pipeline_log(cfg, "Relative abundance")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "relative_abundance")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- 1. Clean taxa names ------------------------------------------------
  df <- clean_taxa_names(cleaned$noncontaminants, cfg)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Relative abundance: no rows after clean_taxa_names — skipping")
    return(invisible(NULL))
  }

  # ---- 2. Per-(taxon, sample) count floor ---------------------------------
  ra_cfg    <- cfg$relative_abundance %||% list()
  min_count <- ra_cfg$min_count_per_sample %||%
               cfg$filters$min_count_per_sample %||% 5
  n_before  <- nrow(df)
  df <- dplyr::filter(df, !is.na(name), count > min_count)
  pipeline_log(cfg, sprintf(
    "Relative abundance: %d/%d (taxon, sample) rows kept with count > %d",
    nrow(df), n_before, min_count
  ))
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Relative abundance: nothing left after count filter — skipping")
    return(invisible(NULL))
  }

  # ---- 3. Optional metadata join for facetting ----------------------------
  facet_by <- ra_cfg$facet_by %||% cfg$metadata$group_cols[[1]]
  if (!is.null(facet_by) && nzchar(facet_by)) {
    sid  <- cfg$metadata$sample_id_col
    meta <- readr::read_csv(
      file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
    )
    if (facet_by %in% colnames(meta)) {
      df <- dplyr::inner_join(df, meta[, c(sid, facet_by)],
                              by = c("sample" = sid))
    } else {
      pipeline_log(cfg, sprintf(
        "Relative abundance: facet column '%s' not in metadata — single panel",
        facet_by
      ))
      facet_by <- NULL
    }
  }

  # ---- 4. Per-sample % ----------------------------------------------------
  df_pct <- df |>
    dplyr::group_by(sample) |>
    dplyr::mutate(percentage = count / sum(count, na.rm = TRUE) * 100) |>
    dplyr::ungroup()

  # ---- 5. Top-N + collapse rest into "Others" -----------------------------
  top_n         <- ra_cfg$top_n         %||% 50
  others_label  <- ra_cfg$others_label  %||% "Others"

  top_taxa <- df_pct |>
    dplyr::group_by(name) |>
    dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::slice_max(total, n = top_n, with_ties = FALSE) |>
    dplyr::pull(name)

  df_pct <- df_pct |>
    dplyr::mutate(name = ifelse(name %in% top_taxa, name, others_label)) |>
    dplyr::group_by(sample, name, dplyr::across(dplyr::any_of(facet_by))) |>
    dplyr::summarise(
      count      = sum(count, na.rm = TRUE),
      percentage = sum(percentage, na.rm = TRUE),
      .groups    = "drop"
    )

  # Order taxa globally by total % so the legend matches the largest stacks.
  taxa_levels <- df_pct |>
    dplyr::group_by(name) |>
    dplyr::summarise(total = sum(percentage, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::pull(name)
  # "Others" pinned to the bottom of the stack regardless of size.
  taxa_levels <- c(setdiff(taxa_levels, others_label),
                   intersect(others_label, taxa_levels))
  df_pct$name <- factor(df_pct$name, levels = taxa_levels)

  n_levels <- length(taxa_levels)
  palette  <- relab_palette(n_levels, cfg)
  pipeline_log(cfg, sprintf(
    "Relative abundance: %d top taxa + %s (palette sized to %d)",
    min(top_n, length(top_taxa)),
    if (others_label %in% taxa_levels) sprintf("'%s'", others_label) else "no overflow",
    n_levels
  ))

  # ---- 6a. Stacked bar plot ----------------------------------------------
  show_legend <- ra_cfg$show_legend %||% TRUE
  legend_pos  <- if (isTRUE(show_legend)) "top" else "none"

  p <- ggplot2::ggplot(df_pct,
                       ggplot2::aes(x = factor(sample), y = percentage, fill = name)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)", fill = "Taxa") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = legend_pos,
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      text            = ggplot2::element_text(size = 12)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 10))

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                                 scales = "free_x", nrow = 1)
  }

  # Size canvas to sample count so labels stay readable across studies.
  n_samples <- length(unique(df_pct$sample))
  pw <- max(10, 0.5 * n_samples + 3)
  ph <- if (isTRUE(show_legend)) max(7, 6 + ceiling(n_levels / 10) * 0.4) else 6

  ggplot2::ggsave(file.path(fig_dir, "stacked_bar.png"), p,
                  width = pw, height = ph, dpi = 300)

  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    # selfcontained=TRUE still extracts plotly assets into <name>_files/ during
    # render; on re-run dir.create() warns that the libdir already exists.
    # Warnings are benign housekeeping — suppress to keep the log clean.
    suppressWarnings(
      htmlwidgets::saveWidget(
        plotly::ggplotly(p),
        file = file.path(fig_dir, "stacked_bar.html"),
        selfcontained = TRUE
      )
    )
  }

  # ---- 6b. Horizontal total-count plot ------------------------------------
  min_for_species <- ra_cfg$min_count_for_species %||%
                     cfg$filters$min_count_for_species %||% 2000

  df_species <- df_pct |>
    dplyr::filter(name != others_label) |>
    dplyr::group_by(name) |>
    dplyr::summarise(totCount = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(totCount >= min_for_species) |>
    dplyr::arrange(dplyr::desc(totCount))

  if (nrow(df_species) > 0) {
    name_colors <- setNames(
      palette[match(as.character(df_species$name), taxa_levels)],
      df_species$name
    )

    sp <- ggplot2::ggplot(
        df_species,
        ggplot2::aes(x = stats::reorder(name, totCount), y = totCount, fill = name)
      ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_hline(yintercept = stats::median(df_species$totCount),
                          linetype = "dashed", colour = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(round(totCount / 1000, 1), "K")),
        hjust = -0.05, size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = name_colors) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Taxon", y = "Total count", fill = "Taxon") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 12)) +
      ggplot2::coord_flip()

    sph <- max(5, 0.25 * nrow(df_species) + 2)
    ggplot2::ggsave(file.path(fig_dir, "species_count.png"), sp,
                    width = 9, height = sph, dpi = 300)
  } else {
    pipeline_log(cfg, sprintf(
      "Relative abundance: no taxa with total count >= %d — species_count.png skipped",
      min_for_species
    ))
  }

  invisible(df_pct)
}
