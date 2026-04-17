# 05_relative_abundance.R — stacked bar plots of relative abundance per sample,
# plus an interactive plotly HTML. Mirrors relativeAbundance.R but reads sample
# IDs dynamically from the metadata and pulls the filter threshold from config.

run_relative_abundance <- function(cleaned, cfg) {
  pipeline_log(cfg, "Relative abundance")
  df <- cleaned$noncontaminants
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "relative_abundance")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  min_count <- cfg$filters$min_count_per_sample %||% 17

  df_pct <- df |>
    dplyr::filter(!is.na(name), count > min_count) |>
    dplyr::group_by(sample) |>
    dplyr::mutate(percentage = count / sum(count) * 100) |>
    dplyr::ungroup()

  t <- length(unique(df_pct$name))
  palette <- paletteer::paletteer_d("palettesForR::Named", n = min(t, 255))

  p <- ggplot2::ggplot(df_pct, ggplot2::aes(factor(sample), percentage, fill = name)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)", fill = "Taxa") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(file.path(fig_dir, "relative_abundance.png"), p,
                  width = 10, height = 6, dpi = 300)

  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    htmlwidgets::saveWidget(
      plotly::ggplotly(p),
      file = file.path(fig_dir, "relative_abundance.html"),
      selfcontained = TRUE
    )
  }

  invisible(df_pct)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
