# utils_prevalence.R — Abundance × Prevalence dot plot helper.
#
# Redesigned 2026-06-12 (was a paired horizontal-bar plot; see
# templates/panels_ref/TPM_Prevalence.png for the new visual). Builds a
# single scatter ggplot for any (sample × category) long-form data frame:
#
#   X = prevalence (%) — 100 * n_distinct(sample) / total_samples per category
#   Y = log10(Average TPM) per category — Total / SampleCount on positive samples
#   colour = category
#   dashed quadrant lines at configurable thresholds (default 80% prevalence
#     and TPM = 1, i.e. log10(TPM) = 0)
#   text labels for the top-N categories by prevalence × Average TPM
#     (ggrepel if available; falls back to geom_text)
#
# API is unchanged: existing callers in R/05 / R/09 / R/10 / R/11 keep
# working with the same arguments. The `value_label`, `widths`, and
# `log_x_total` args are accepted but ignored (kept for backward compat).
#
# Args:
#   df              long-form data frame with at least sample, category, value cols
#   category_col    string — column name for category (e.g. "name", "DRUG",
#                   "Functions", "GENE", "Replicon_Family")
#   value_col       string — abundance column to summarise (e.g. "TPM", "count")
#   sample_col      string — sample id column (default "sample")
#   palette         optional named vector mapping category levels -> hex
#   category_label  legend / colour-guide title
#   prev_cutoff_pct vertical dashed line position (% prevalence; default 80)
#   tpm_cutoff_log  horizontal dashed line position (log10 TPM; default 0)
#   top_n_labels    number of categories to label by prev × abundance product
#                   (default 10)
#   ...             extra args swallowed for backward compatibility
#
# Returns a single ggplot object — caller ggsaves it.

`%||%` <- function(a, b) if (is.null(a)) b else a

build_count_prevalence <- function(df,
                                    category_col,
                                    value_col,
                                    sample_col      = "sample",
                                    palette         = NULL,
                                    category_label  = NULL,
                                    prev_cutoff_pct = 80,
                                    tpm_cutoff_log  = 0,
                                    top_n_labels    = 10,
                                    ...) {
  for (pkg in c("ggplot2", "dplyr")) {
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

  total_samples <- dplyr::n_distinct(df[[sample_col]])
  if (total_samples == 0) return(NULL)

  per_cat <- df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(
      Total       = sum(.data[[value_col]], na.rm = TRUE),
      SampleCount = dplyr::n_distinct(.data[[sample_col]]),
      .groups     = "drop"
    ) |>
    dplyr::mutate(
      Prevalence  = 100 * .data$SampleCount / total_samples,
      AvgTPM      = ifelse(.data$SampleCount > 0,
                            .data$Total / .data$SampleCount, NA_real_),
      LogAvgTPM   = log10(pmax(.data$AvgTPM, 1e-3, na.rm = FALSE)),
      Score       = .data$Prevalence * .data$AvgTPM
    ) |>
    dplyr::arrange(dplyr::desc(.data$Score))

  if (nrow(per_cat) == 0) return(NULL)

  to_label <- utils::head(per_cat, max(0L, as.integer(top_n_labels)))

  p <- ggplot2::ggplot(
      per_cat,
      ggplot2::aes(x      = .data$Prevalence,
                   y      = .data$LogAvgTPM,
                   colour = .data[[category_col]])
    ) +
    ggplot2::geom_vline(xintercept = prev_cutoff_pct,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = tpm_cutoff_log,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::scale_x_continuous(limits = c(0, 100),
                                 breaks = seq(0, 100, by = 25)) +
    ggplot2::labs(
      x      = "Prevalence (%)",
      y      = expression(log[10] ~ "Average TPM"),
      colour = category_label %||% category_col
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      text = ggplot2::element_text(size = 12)
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_colour_manual(values = palette,
                                            na.value = "grey60")
  }

  # Plain geom_text rather than ggrepel — ggrepel's viewport-bound
  # placement crashes ("Viewport has zero dimension(s)") on log-scaled
  # y-axes with extreme dynamic range. hjust/vjust offsets keep labels
  # off the dots; small size keeps them legible even when categories
  # cluster in the upper-right quadrant.
  if (nrow(to_label) > 0) {
    p <- p + ggplot2::geom_text(
      data = to_label,
      ggplot2::aes(label = .data[[category_col]]),
      hjust = -0.12, vjust = -0.4, size = 3.0,
      check_overlap = TRUE, show.legend = FALSE
    )
  }

  p
}

# Convenience wrapper: build + ggsave in one call. Returns the file path on
# success, NULL when the data is empty.
save_count_prevalence <- function(df, png_path, ..., width = 9, height = 6,
                                   dpi = 300) {
  p <- build_count_prevalence(df, ...)
  if (is.null(p)) return(NULL)
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi,
                  bg = "white")
  png_path
}
