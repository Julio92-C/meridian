# utils_ge_profile.R — Shared plumbing for gene-element profile stages
# (R/09 resistome, R/10 virulome, R/11 mobilome). All three stages read the
# TPM-normalised gene table from R/03, filter by a single DATABASE value,
# extract a per-gene category column, and produce a parallel suite of alpha
# diversity / Venn / heatmap / relative-abundance / total-count / beta-PCoA
# outputs. The helpers below encapsulate the structural pieces; each module
# supplies the database name, per-gene category extraction, output filenames,
# and labels.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Read the TPM-normalised gene table produced by R/03. Returns NULL with a
# log line when the file isn't on disk (e.g. cfg$stages$normalisation = false).
ge_load_norm_table <- function(cfg, log_label) {
  path <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                    "genetable_normdata.csv")
  if (!file.exists(path)) {
    pipeline_log(cfg, sprintf(
      "%s: %s not found — enable cfg$stages$normalisation first. Skipping.",
      log_label, path
    ))
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

# R/03 only attaches fixed_effects + random_effect to the normalised table.
# When the requested group column isn't already present, re-join metadata.
ge_ensure_group_column <- function(df, cfg, group, log_label) {
  if (group %in% colnames(df)) return(df)
  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE)
  if (!group %in% colnames(meta)) {
    stop(sprintf("%s: group column '%s' not found in metadata",
                 log_label, group))
  }
  df$sample    <- as.character(df$sample)
  meta[[sid]]  <- as.character(meta[[sid]])
  dplyr::inner_join(df, meta[, c(sid, group), drop = FALSE],
                    by = c("sample" = sid))
}

# Discrete palette of length `n` from a paletteer name; falls back to
# grDevices::hcl.colors when the name is unknown or paletteer is missing.
ge_discrete_palette <- function(name, n) {
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

# Build a named palette for `levels`. If `user_map` is a non-null named list
# (cfg-supplied), use it; otherwise pull `length(levels)` colours from
# `palette_name` via ge_discrete_palette().
ge_named_palette <- function(levels, user_map = NULL,
                              palette_name = "ggsci::default_nejm") {
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- ge_discrete_palette(palette_name, length(levels))
  setNames(pal, levels)
}

# Apply named-list literal-substring renames to a character vector. A hit
# replaces the whole value (mirrors the GT case_when where any substring
# match swaps the GENE / Functions value).
ge_apply_literal_renames <- function(x, renames) {
  if (is.null(renames) || length(renames) == 0) return(as.character(x))
  pats <- names(renames)
  out  <- as.character(x)
  for (i in seq_along(pats)) {
    hits <- !is.na(out) & stringr::str_detect(out, stringr::fixed(pats[[i]]))
    out[hits] <- renames[[i]]
  }
  out
}

# Per-sample richness + Shannon on a log(value + 1) gene matrix.
ge_compute_alpha <- function(df, group, value_col = "TPM") {
  abund <- df |>
    dplyr::mutate(val_log = log(.data[[value_col]] + 1)) |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(val_log = sum(.data$val_log), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "val_log",
                       values_fill = 0)
  mat <- as.matrix(abund[, -1, drop = FALSE])
  rownames(mat) <- abund$sample
  group_per_sample <- df |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group))
  out <- data.frame(
    sample   = rownames(mat),
    richness = rowSums(mat > 0),
    shannon  = vegan::diversity(mat, "shannon"),
    stringsAsFactors = FALSE
  )
  dplyr::left_join(out, group_per_sample, by = "sample")
}

# Kruskal-Wallis test of `metric` ~ `group` on the alpha table. Returns NULL
# (and skips logging) if the metric is missing or there's only one group.
ge_alpha_kw <- function(alpha, group, metric, cfg, log_label) {
  if (!metric %in% colnames(alpha)) {
    pipeline_log(cfg, sprintf("%s alpha: metric '%s' not found",
                              log_label, metric))
    return(NULL)
  }
  sub <- alpha[!is.na(alpha[[metric]]), , drop = FALSE]
  if (dplyr::n_distinct(sub[[group]]) < 2) return(NULL)
  res <- tryCatch(
    kruskal.test(reformulate(group, metric), data = sub),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    pipeline_log(cfg, sprintf("%s alpha %s ~ %s KW p = %.4g",
                              log_label, metric, group, res$p.value))
  }
  res
}

# Min-max scale each column of a matrix independently to [0, 1]. Columns
# with zero range collapse to all-zero. Matches the GT pheatmap recipe.
ge_minmax_per_col <- function(mat) {
  apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
}

# Violin + boxplot + jitter of `metric` ~ `group`, fill by group. KW
# annotation drawn in the upper-right corner when `kw` is non-NULL, so the
# label doesn't collide with the violin distributions.
ge_plot_alpha_violin <- function(alpha, group, metric, kw, pal_group, file) {
  if (!metric %in% colnames(alpha)) return(invisible(NULL))
  p <- ggplot2::ggplot(alpha,
        ggplot2::aes(x = .data[[group]], y = .data[[metric]],
                     fill = .data[[group]])) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          position = ggplot2::position_dodge(0.9)) +
    ggplot2::geom_jitter(width = 0.08, size = 1.4, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = stringr::str_to_title(metric)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 13))
  if (!is.null(kw)) {
    p <- p + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = sprintf("Kruskal-Wallis p = %.2g", kw$p.value),
      hjust = 1.05, vjust = 1.4, size = 4.2, colour = "black"
    )
  }
  ggplot2::ggsave(file, p, width = 7, height = 5, dpi = 300)
}

ge_plot_alpha_bar <- function(alpha, group, metric, kw, pal_group, file) {
  if (!metric %in% colnames(alpha)) return(invisible(NULL))
  metric_label <- stringr::str_to_title(metric)
  title <- if (!is.null(kw))
    sprintf("%s by %s (Kruskal-Wallis p = %.3g)",
            metric_label, group, kw$p.value)
  else
    sprintf("%s by %s", metric_label, group)
  p <- ggplot2::ggplot(alpha,
        ggplot2::aes(x = sample, y = .data[[metric]],
                     fill = .data[[group]])) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = round(.data[[metric]], 1)),
                       vjust = -0.6, size = 3) +
    ggplot2::geom_hline(yintercept = mean(alpha[[metric]], na.rm = TRUE),
                        linetype = "dashed", colour = "red") +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(title = title, x = "Sample", y = metric_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      text            = ggplot2::element_text(size = 13)
    )
  bw <- max(8, 0.5 * nrow(alpha) + 3)
  ggplot2::ggsave(file, p, width = bw, height = 5, dpi = 300)
}

# Per-(sample, gene) log(value + 1) violin + box + jitter coloured by group,
# with Kruskal-Wallis annotation. Complements the richness-based alpha plot
# by showing per-feature abundance distribution differences across groups.
ge_plot_abundance_violin <- function(df, group, pal_group, file,
                                      title, log_label, cfg,
                                      value_col = "TPM",
                                      y_label = "log(TPM)") {
  if (!value_col %in% colnames(df) || nrow(df) == 0) return(invisible(NULL))
  v <- data.frame(
    sample  = df$sample,
    grp_val = df[[group]],
    log_val = log(df[[value_col]] + 1)
  )
  v <- v[is.finite(v$log_val), , drop = FALSE]
  if (dplyr::n_distinct(v$grp_val) < 2 || nrow(v) < 3) {
    pipeline_log(cfg, sprintf(
      "%s abundance: insufficient data for KW — skipping", log_label
    ))
    return(invisible(NULL))
  }
  names(v)[names(v) == "grp_val"] <- group
  kw <- tryCatch(
    kruskal.test(reformulate(group, "log_val"), data = v),
    error = function(e) NULL
  )
  if (!is.null(kw)) {
    pipeline_log(cfg, sprintf("%s log(%s+1) ~ %s KW p = %.4g",
                              log_label, value_col, group, kw$p.value))
  }
  p <- ggplot2::ggplot(v,
        ggplot2::aes(x = .data[[group]], y = .data$log_val,
                     fill = .data[[group]])) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          position = ggplot2::position_dodge(0.9)) +
    ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.35) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = y_label, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold",
                                                       size = 13),
                   text = ggplot2::element_text(size = 13))
  # KW annotation in upper-right corner so it doesn't overlap the
  # log(TPM+1) distributions (which often peak near the panel top).
  if (!is.null(kw)) {
    p <- p + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = sprintf("Kruskal-Wallis p = %.2g", kw$p.value),
      hjust = 1.05, vjust = 1.4, size = 4.2, colour = "black"
    )
  }
  ggplot2::ggsave(file, p, width = 7, height = 5.5, dpi = 300)
}

# Venn diagram of unique category values per group. Supports 2-5 groups
# (VennDiagram::venn.diagram). Category counts are aggregated by `value_col`
# so a row contributes to a set only when total > 0.
ge_plot_category_venn <- function(df, category_col, group, pal_group,
                                   file, main_title, log_label, cfg,
                                   value_col = "TPM") {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    pipeline_log(cfg, sprintf("%s: VennDiagram not available — Venn skipped",
                              log_label))
    return(invisible(NULL))
  }
  per_group <- df |>
    dplyr::group_by(.data[[group]], .data[[category_col]]) |>
    dplyr::summarise(total = sum(.data[[value_col]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(.data$total > 0)
  sets  <- split(per_group[[category_col]], per_group[[group]])
  sets  <- lapply(sets, unique)
  n_grp <- length(sets)
  if (n_grp < 2) {
    pipeline_log(cfg, sprintf("%s: <2 groups — Venn skipped", log_label))
    return(invisible(NULL))
  }
  if (n_grp > 5) {
    pipeline_log(cfg, sprintf(
      "%s: %d groups — Venn supports 2-5, skipped", log_label, n_grp
    ))
    return(invisible(NULL))
  }
  fill_pal <- unname(pal_group[names(sets)])
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    suppressMessages(futile.logger::flog.threshold(
      futile.logger::FATAL, name = "VennDiagramLogger"
    ))
  }
  cat_args <- switch(as.character(n_grp),
    "2" = list(cat.pos = c(-20, 20),       cat.dist = c(0.05, 0.05)),
    "3" = list(cat.pos = c(-25, 25, 180),  cat.dist = c(0.08, 0.08, 0.04)),
    "4" = list(cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22, 0.22, 0.12, 0.12)),
    "5" = list(),
    list()
  )
  do.call(VennDiagram::venn.diagram, c(
    list(
      x         = sets,
      filename  = file,
      imagetype = "png",
      height = 2800, width = 2800, resolution = 300,
      fill = fill_pal, alpha = 0.6,
      cat.cex = 1.25, cat.fontface = "bold",
      cex = 1.5,
      main = main_title,
      main.cex = 1.3, margin = 0.18,
      disable.logging = TRUE
    ),
    cat_args
  ))
}

# Gene-level pheatmap (rows = GENE, cols = sample). Rows ordered by total
# descending then top_n filtered; min-max scaled per column; rows annotated
# by `category_col` and columns annotated by `group`, with gaps between
# group runs. Caller is responsible for any GENE renames in `df`.
ge_plot_gene_heatmap <- function(df, group, pal_group, pal_category,
                                  category_col, file, palette,
                                  top_n             = NULL,
                                  fontsize_row      = 8,
                                  row_height_factor = 0.18,
                                  min_height        = 7) {
  gene_wide <- df |>
    dplyr::group_by(.data$GENE, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("GENE") |>
    as.matrix()
  gene_wide <- gene_wide[order(rowSums(gene_wide), decreasing = TRUE), ,
                          drop = FALSE]
  if (!is.null(top_n) && nrow(gene_wide) > top_n) {
    gene_wide <- gene_wide[seq_len(top_n), , drop = FALSE]
  }
  gene_scaled <- ge_minmax_per_col(gene_wide)

  cat_per_gene <- df |>
    dplyr::distinct(.data$GENE, .keep_all = TRUE) |>
    dplyr::select("GENE", dplyr::all_of(category_col)) |>
    tibble::column_to_rownames("GENE")
  cat_per_gene <- cat_per_gene[rownames(gene_scaled), , drop = FALSE]

  sample_group <- df |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group)) |>
    tibble::column_to_rownames("sample")
  sample_group <- sample_group[order(sample_group[[group]]), , drop = FALSE]
  gene_scaled  <- gene_scaled[, rownames(sample_group), drop = FALSE]
  grp_run      <- as.character(sample_group[[group]])
  grp_lvls     <- unique(grp_run)
  grp_sizes    <- as.integer(table(grp_run)[grp_lvls])
  gaps_col     <- if (length(grp_sizes) > 1) {
    utils::head(cumsum(grp_sizes), -1L)
  } else NULL

  ann_colors <- list()
  ann_colors[[group]] <- pal_group[grp_lvls]
  cat_levels <- unique(stats::na.omit(cat_per_gene[[category_col]]))
  ann_colors[[category_col]] <- pal_category[intersect(names(pal_category),
                                                       cat_levels)]

  hm_height <- max(min_height, row_height_factor * nrow(gene_scaled) + 3)
  hm_width  <- max(8, 0.4 * ncol(gene_scaled) + 4)
  pheatmap::pheatmap(
    gene_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(palette)(100),
    border_color      = NA,
    annotation_row    = cat_per_gene,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = fontsize_row, fontsize_col = 9, fontsize = 10,
    filename = file,
    width = hm_width, height = hm_height
  )
}

# Category-level pheatmap: rows aggregated by `category_col`, columns =
# sample. Rows sorted by total TPM descending; columns sorted by group.
ge_plot_category_heatmap <- function(df, category_col, group, pal_group,
                                      file, palette,
                                      fontsize_row      = 11,
                                      row_height_factor = 0.32,
                                      min_height        = 4) {
  cat_wide <- df |>
    dplyr::group_by(.data[[category_col]], .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames(category_col) |>
    as.matrix()
  cat_wide   <- cat_wide[order(rowSums(cat_wide), decreasing = TRUE), ,
                          drop = FALSE]
  cat_scaled <- ge_minmax_per_col(cat_wide)

  sample_group <- df |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group)) |>
    tibble::column_to_rownames("sample")
  sample_group <- sample_group[order(sample_group[[group]]), , drop = FALSE]
  cat_scaled   <- cat_scaled[, rownames(sample_group), drop = FALSE]
  grp_run      <- as.character(sample_group[[group]])
  grp_lvls     <- unique(grp_run)
  grp_sizes    <- as.integer(table(grp_run)[grp_lvls])
  gaps_col     <- if (length(grp_sizes) > 1) {
    utils::head(cumsum(grp_sizes), -1L)
  } else NULL

  ann_colors <- list()
  ann_colors[[group]] <- pal_group[grp_lvls]

  hm_height <- max(min_height, row_height_factor * nrow(cat_scaled) + 3)
  hm_width  <- max(8, 0.4 * ncol(cat_scaled) + 4)
  pheatmap::pheatmap(
    cat_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(palette)(100),
    border_color      = NA,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = fontsize_row, fontsize_col = 10, fontsize = 11,
    filename = file,
    width = hm_width, height = hm_height
  )
}

# Stacked-bar relative abundance of `category_col` per sample, faceted by
# group. Also writes a category totals CSV and returns the totals tibble
# (caller passes the tibble to ge_plot_total_bar).
ge_plot_category_relative_abundance <- function(df, category_col, group,
                                                  palette, fig_file, ds_file,
                                                  fill_label,
                                                  value_col = "TPM",
                                                  totals_value_name = "Total_TPM") {
  share <- df |>
    dplyr::group_by(.data$sample, .data[[category_col]], .data[[group]]) |>
    dplyr::summarise(val = sum(.data[[value_col]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(Percentage = 100 * .data$val / sum(.data$val)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(share,
        ggplot2::aes(x = sample, y = Percentage,
                     fill = .data[[category_col]])) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                  fill = fill_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  bw <- max(8, 0.55 * dplyr::n_distinct(share$sample) + 3)
  ggplot2::ggsave(fig_file, p, width = bw, height = 6, dpi = 300)

  totals <- df |>
    dplyr::group_by(.data[[category_col]]) |>
    dplyr::summarise(!!totals_value_name :=
                       sum(.data[[value_col]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data[[totals_value_name]]))
  readr::write_csv(totals, ds_file)
  totals
}

# Horizontal log10 total-count bar with K-formatted labels and a median line.
# `category_col` and `value_col` name columns of `totals`.
ge_plot_total_bar <- function(totals, category_col, value_col, palette,
                               file, x_label,
                               y_label = "Total TPM (log10)",
                               width = 7, height = 6) {
  fmt_label <- function(x) {
    ifelse(x >= 1000,
           paste0(round(x / 1000, 1), "K"),
           as.character(round(x, 1)))
  }
  p <- ggplot2::ggplot(totals,
        ggplot2::aes(x = stats::reorder(.data[[category_col]],
                                         .data[[value_col]]),
                     y = .data[[value_col]],
                     fill = .data[[category_col]])) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = stats::median(totals[[value_col]]),
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_label(.data[[value_col]])),
                       hjust = -0.05, size = 3.4) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 11),
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  ggplot2::ggsave(file, p, width = width, height = height, dpi = 300)
}

# Gene-level PCoA + PERMANOVA + PERMDISP on a TPM profile, matching the
# R/07 recipe (Hellinger -> Bray-Curtis -> cmdscale + adonis2 + betadisper).
ge_plot_beta <- function(df, group, pal_group, opts, fig_file,
                          ds_perm_file, ds_disp_file, ds_pcoa_file,
                          log_label, title_prefix, skip_data_label, cfg) {
  wide <- df |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM",
                       values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) < 3) {
    pipeline_log(cfg, sprintf(
      "%s beta: only %d samples with %s — need >= 3, skipping",
      log_label, nrow(mat), skip_data_label
    ))
    return(invisible(NULL))
  }

  transform <- opts$beta_transform %||% "hellinger"
  distance  <- opts$beta_distance  %||% "bray"
  mat_t <- switch(transform,
    none      = mat,
    log       = log1p(mat),
    hellinger = vegan::decostand(mat, method = "hellinger"),
    stop(sprintf("%s beta: unknown transform '%s'", log_label, transform))
  )
  d <- vegan::vegdist(mat_t, method = distance)

  meta <- df |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group))
  meta <- meta[match(rownames(mat_t), meta$sample), , drop = FALSE]

  perms <- cfg$stats$permanova_permutations %||% 9999
  permanova <- vegan::adonis2(
    stats::reformulate(group, "d"),
    data = meta, permutations = perms
  )
  utils::capture.output(permanova, file = ds_perm_file)
  r2 <- permanova$R2[1]
  pv <- permanova$`Pr(>F)`[1]
  pipeline_log(cfg, sprintf("%s PERMANOVA %s: R2 = %.3f, p = %.4g",
                            log_label, group, r2, pv))

  bd_test <- tryCatch({
    bd <- vegan::betadisper(d, factor(meta[[group]]))
    vegan::permutest(bd, permutations = perms)
  }, error = function(e) NULL)
  permdisp_p <- NA_real_
  if (!is.null(bd_test)) {
    utils::capture.output(bd_test, file = ds_disp_file)
    permdisp_p <- bd_test$tab$`Pr(>F)`[1]
    pipeline_log(cfg, sprintf("%s PERMDISP %s: p = %.4g",
                              log_label, group, permdisp_p))
  }

  pcoa <- stats::cmdscale(d, eig = TRUE, k = 2)
  var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
  scores <- data.frame(
    sample = rownames(mat_t),
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2]
  )
  scores <- dplyr::left_join(scores, meta, by = "sample")
  readr::write_csv(scores, ds_pcoa_file)

  annot <- sprintf("PERMANOVA R² = %.3f, p = %.4g\nPERMDISP p = %.4g",
                   r2, pv, permdisp_p)
  ellipse_type <- opts$ellipse_type     %||% "norm"
  ellipse_line <- opts$ellipse_linetype %||% "dashed"

  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal_group) +
    ggplot2::labs(
      title = sprintf("%s PCoA (%s, %s-transformed) — by %s",
                      title_prefix, distance, transform, group),
      colour = group,
      x = sprintf("PC1 (%.1f%%)", var_expl[1]),
      y = sprintf("PC2 (%.1f%%)", var_expl[2])
    ) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = annot,
                      hjust = 1.05, vjust = 1.5,
                      size = 4, colour = "black") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  ggplot2::ggsave(fig_file, p, width = 7, height = 5.5, dpi = 300)

  invisible(list(scores = scores, permanova = permanova,
                 permdisp = bd_test, var_explained = var_expl[1:2]))
}
