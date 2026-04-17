# 06_alpha_diversity.R — richness + Shannon; Kruskal-Wallis across groups,
# optionally using block as a random effect (lme4) for chicken studies.

run_alpha_diversity <- function(cleaned, cfg) {
  pipeline_log(cfg, "Alpha diversity")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "alpha_diversity")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]

  counts <- cleaned$noncontaminants |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)

  mat <- as.matrix(counts[, -1])
  rownames(mat) <- counts$sample

  div <- data.frame(
    sample    = rownames(mat),
    richness  = rowSums(mat > 0),
    shannon   = vegan::diversity(mat, index = "shannon"),
    simpson   = vegan::diversity(mat, index = "simpson")
  )
  div <- dplyr::inner_join(div, meta, by = c("sample" = sid))

  readr::write_csv(div, file.path(cfg$project_root, cfg$outputs$datasets_dir,
                                  "alpha_diversity.csv"))

  kw_richness <- kruskal.test(reformulate(group, "richness"), data = div)
  kw_shannon  <- kruskal.test(reformulate(group, "shannon"),  data = div)

  stats_lines <- c(
    sprintf("Kruskal-Wallis richness ~ %s: p = %.4g", group, kw_richness$p.value),
    sprintf("Kruskal-Wallis shannon  ~ %s: p = %.4g", group, kw_shannon$p.value)
  )
  writeLines(stats_lines, file.path(fig_dir, "alpha_stats.txt"))
  pipeline_log(cfg, stats_lines[[1]])
  pipeline_log(cfg, stats_lines[[2]])

  for (metric in c("richness", "shannon", "simpson")) {
    p <- ggplot2::ggplot(div, ggplot2::aes(.data[[group]], .data[[metric]],
                                           fill = .data[[group]])) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
      ggplot2::geom_jitter(width = 0.1, size = 1.5) +
      ggpubr::stat_compare_means(method = "kruskal.test") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste(metric, "by", group))
    ggplot2::ggsave(file.path(fig_dir, paste0(metric, "_violin.png")),
                    p, width = 6, height = 5, dpi = 300)
  }

  invisible(div)
}
