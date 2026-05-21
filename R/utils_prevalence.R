# utils_prevalence.R — Shared paired (total | prevalence) plot helper.
#
# Builds a side-by-side patchwork ggplot for any (sample × category) long-form
# data frame:
#   left  panel — horizontal Total <value_col> per category, log10 axis,
#                 K-formatted bar labels, dashed median line.
#   right panel — horizontal n_distinct(sample) per category, dashed median
#                 line, y-axis labels/ticks hidden so the left panel's labels
#                 read across the combined view (matching the GT design).
#
# Both panels share the same y-axis ordering (categories sorted by descending
# Total). Used by R/05 (taxa), R/09 (drug class), R/10 (VF function), R/11
# (MGE gene + family). Returns a patchwork object — caller ggsaves it.
#
# Args:
#   df            long-form data frame with at least sample, category, value cols
#   category_col  string — column name for the category axis (e.g. "name",
#                 "DRUG", "Functions", "GENE", "Replicon_Family")
#   value_col     string — column name to sum for totals (e.g. "count", "TPM")
#   sample_col    string — column name for the sample id (default "sample")
#   palette       optional named vector mapping category levels -> hex
#   value_label   x-axis label for the total panel
#   category_label y-axis label for the total panel (right panel has no label)
#   log_x_total   if TRUE (default) use log10 x-scale on the total panel
#   widths        patchwork width ratio (default c(2, 1) — total wider)

`%||%` <- function(a, b) if (is.null(a)) b else a

build_count_prevalence <- function(df,
                                    category_col,
                                    value_col,
                                    sample_col      = "sample",
                                    palette         = NULL,
                                    value_label     = "Total count",
                                    category_label  = NULL,
                                    log_x_total     = TRUE,
                                    widths          = c(2, 1)) {
  for (pkg in c("patchwork", "ggplot2", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("build_count_prevalence: required package '", pkg, "' not available")
    }
  }
  stopifnot(is.data.frame(df),
            category_col %in% colnames(df),
            value_col    %in% colnames(df),
            sample_col   %in% colnames(df))

  df <- df[!is.na(df[[category_col]]) & nzchar(as.character(df[[category_col]])), ,
           drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  totals <- df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(Total = sum(.data[[value_col]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Total))

  prev <- df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(SampleCount = dplyr::n_distinct(.data[[sample_col]]),
                     .groups = "drop")

  order_vec <- as.character(totals[[category_col]])
  totals[[category_col]] <- factor(totals[[category_col]],
                                   levels = rev(order_vec))
  prev[[category_col]]   <- factor(prev[[category_col]],
                                   levels = rev(order_vec))
  prev <- prev[!is.na(prev[[category_col]]), , drop = FALSE]

  fmt_K <- function(x) ifelse(x >= 1000,
                              paste0(round(x / 1000, 1), "K"),
                              as.character(round(x, 1)))

  med_total <- stats::median(totals$Total, na.rm = TRUE)
  med_prev  <- stats::median(prev$SampleCount, na.rm = TRUE)

  p_total <- ggplot2::ggplot(
      totals,
      ggplot2::aes(y    = .data[[category_col]],
                   x    = .data$Total,
                   fill = .data[[category_col]])
    ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = med_total, linetype = "dashed",
                        colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_K(.data$Total)),
                       hjust = -0.05, size = 3.4) +
    ggplot2::labs(x = value_label, y = category_label %||% category_col) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 12)) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18)),
      trans  = if (isTRUE(log_x_total)) "log10" else "identity"
    )

  p_prev <- ggplot2::ggplot(
      prev,
      ggplot2::aes(y    = .data[[category_col]],
                   x    = .data$SampleCount,
                   fill = .data[[category_col]])
    ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = med_prev, linetype = "dashed",
                        colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = .data$SampleCount),
                       hjust = -0.3, size = 3.4) +
    ggplot2::labs(x = "Sample count", y = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y     = ggplot2::element_blank(),
      axis.ticks.y    = ggplot2::element_blank(),
      text            = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18))
    )

  if (!is.null(palette)) {
    p_total <- p_total +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey60")
    p_prev  <- p_prev +
      ggplot2::scale_fill_manual(values = palette, na.value = "grey60")
  }

  patchwork::wrap_plots(p_total, p_prev, widths = widths)
}

# Convenience wrapper: build + ggsave in one call. Returns the file path on
# success, NULL when the data is empty.
save_count_prevalence <- function(df, png_path, ..., width = 12, height = 6,
                                   dpi = 300) {
  p <- build_count_prevalence(df, ...)
  if (is.null(p)) return(NULL)
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi,
                  bg = "white")
  png_path
}
