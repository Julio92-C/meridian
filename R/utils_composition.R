# utils_composition.R — Shared streamgraph composition helper.
#
# Companion to utils_prevalence.R. Renders the per-treatment-group
# smoothed stacked area mirroring templates/panels_ref/relativate_abundace.png
# (relative abundance %, x = treatment level, y = 0–100, fill = category).
#
# Used by:
#   R/05_relative_abundance.R   (taxa species + per-rank rollups)
#   R/09_resistome.R            (drug class)
#   R/10_virulome.R             (VF functional category)
#   R/11_mobilome.R             (replicon family)
#
# Aggregation: long-form (sample, category, value, group) is reduced to a
# per-group composition by:
#   1. converting value to % of within-sample total,
#   2. averaging % per (group, category) across samples in that group.
# The resulting per-group means are stacked + smoothed via ggstream when
# available; falls back to geom_area with stat = "identity" when ggstream
# is missing. Top-N categories surface in the legend; the rest are
# collapsed into "Others" so the palette stays readable.

`%||%` <- function(a, b) if (is.null(a)) b else a

build_stream_composition <- function(df,
                                      category_col,
                                      group_col,
                                      value_col      = "count",
                                      sample_col     = "sample",
                                      palette        = NULL,
                                      top_n          = 12,
                                      others_label   = "Others",
                                      title          = NULL,
                                      stream_type    = "proportional") {
  for (pkg in c("ggplot2", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("build_stream_composition: required package '", pkg, "' not available")
    }
  }
  needed <- c(category_col, group_col, value_col, sample_col)
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) {
    stop("build_stream_composition: missing column(s) ",
         paste(miss, collapse = ", "))
  }

  df <- df[!is.na(df[[category_col]]) &
             nzchar(as.character(df[[category_col]])) &
             !is.na(df[[group_col]]), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  # 1. Within-sample % so groups with different total counts contribute
  #    comparably.
  df_pct <- df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::mutate(.sample_total = sum(.data[[value_col]], na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(.pct = ifelse(.data$.sample_total > 0,
                                 100 * .data[[value_col]] /
                                       .data$.sample_total, 0))

  # 2. Top-N globally by mean %, collapse the rest into Others. Computing
  #    the cap on the unaggregated frame so a category that's dominant in
  #    just one group still surfaces.
  global_means <- df_pct |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(mean_pct = mean(.data$.pct, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$mean_pct))
  top_cats <- utils::head(as.character(global_means[[category_col]]),
                          max(1L, as.integer(top_n)))
  df_pct[[category_col]] <- ifelse(
    as.character(df_pct[[category_col]]) %in% top_cats,
    as.character(df_pct[[category_col]]), others_label
  )
  df_pct <- df_pct |>
    dplyr::group_by(.data[[sample_col]], .data[[category_col]],
                     .data[[group_col]]) |>
    dplyr::summarise(.pct = sum(.data$.pct, na.rm = TRUE),
                     .groups = "drop")

  # 3. Expand to (sample, category) so absent categories contribute zero
  #    to the per-group mean. Without this the stack tops out below 100%
  #    because sparse categories' means are computed over present samples
  #    only and don't compose with present-category means.
  sample_group <- df_pct |>
    dplyr::distinct(.data[[sample_col]], .data[[group_col]])
  cats <- unique(c(top_cats, others_label))
  cats <- intersect(cats, unique(as.character(df_pct[[category_col]])))
  full_grid <- tidyr::expand_grid(
    sample_group,
    !!category_col := cats
  )
  df_pct <- dplyr::left_join(full_grid, df_pct,
                              by = c(sample_col, group_col, category_col)) |>
    dplyr::mutate(.pct = tidyr::replace_na(.data$.pct, 0))

  # 4. Per-(group, category) mean %.
  per_group <- df_pct |>
    dplyr::group_by(.data[[group_col]], .data[[category_col]]) |>
    dplyr::summarise(pct = mean(.data$.pct, na.rm = TRUE),
                     .groups = "drop")

  # Category order for the stack (largest at bottom, Others pinned bottom).
  cat_levels <- c(top_cats, others_label)
  cat_levels <- intersect(cat_levels, unique(per_group[[category_col]]))
  per_group[[category_col]] <- factor(per_group[[category_col]],
                                       levels = cat_levels)

  # Group order: alphabetic by default; callers can pre-factor if they
  # want a specific order.
  if (!is.factor(per_group[[group_col]])) {
    per_group[[group_col]] <- factor(per_group[[group_col]])
  }
  per_group$.x <- as.integer(per_group[[group_col]])

  use_stream <- requireNamespace("ggstream", quietly = TRUE)

  if (use_stream) {
    p <- ggplot2::ggplot(
        per_group,
        ggplot2::aes(x    = .data$.x,
                     y    = .data$pct,
                     fill = .data[[category_col]])
      ) +
      ggstream::geom_stream(type = stream_type, bw = 0.65,
                             extra_span = 0.05)
  } else {
    p <- ggplot2::ggplot(
        per_group,
        ggplot2::aes(x    = .data$.x,
                     y    = .data$pct,
                     fill = .data[[category_col]])
      ) +
      ggplot2::geom_area(position = "stack", alpha = 0.95)
  }

  group_levels <- levels(per_group[[group_col]])
  p <- p +
    ggplot2::scale_x_continuous(breaks = seq_along(group_levels),
                                 labels = group_levels,
                                 expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(
      x     = group_col,
      y     = "Relative abundance (%)",
      fill  = category_col,
      title = title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      plot.title      = ggplot2::element_text(hjust = 0.5),
      text            = ggplot2::element_text(size = 12)
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette,
                                          na.value = "grey60")
  }

  p
}

# Convenience wrapper: build + ggsave in one call. Returns the file path
# on success, NULL when the data is empty.
save_stream_composition <- function(df, png_path, ..., width = 9, height = 5.5,
                                     dpi = 300) {
  p <- build_stream_composition(df, ...)
  if (is.null(p)) return(NULL)
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi,
                  bg = "white")
  png_path
}
