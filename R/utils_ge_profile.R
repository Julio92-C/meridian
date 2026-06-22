# utils_ge_profile.R — Shared plumbing for gene-element profile stages
# (R/09 resistome, R/10 virulome, R/11 mobilome). All three stages read the
# TPM-normalised gene table from R/03, filter by a single DATABASE value,
# extract a per-gene category column, and produce a parallel suite of alpha
# diversity / Venn / heatmap / relative-abundance / total-count / beta-PCoA
# outputs. The helpers below encapsulate the structural pieces; each module
# supplies the database name, per-gene category extraction, output filenames,
# and labels.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build a single horizontal legend grob combining the heatmap colour scale
# (continuous gradient) and any number of categorical annotation legends.
# Used by .pheatmap_save_legend_top to replace pheatmap's vertical legends.
.make_gradient_legend_grob <- function(palette, limits = c(0, 1),
                                        title    = "value",
                                        fontsize = 8) {
  df <- data.frame(x = c(0, 1), v = limits)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$x,
                                          fill = .data$v)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = palette,
      limits  = limits,
      name    = title,
      guide   = ggplot2::guide_colorbar(
        direction      = "horizontal",
        title.position = "left",
        title.vjust    = 0.85,
        barwidth       = grid::unit(2.8, "cm"),
        barheight      = grid::unit(0.3, "cm")
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "top",
      legend.text     = ggplot2::element_text(size = fontsize),
      legend.title    = ggplot2::element_text(size = fontsize + 1)
    )
  cowplot::get_legend(p)
}

.make_categorical_legend_grob <- function(pal, title, fontsize = 8) {
  df <- data.frame(x = seq_along(pal), y = 1,
                   cat = factor(names(pal), levels = names(pal)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          fill = .data$cat)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = pal, name = title,
      guide  = ggplot2::guide_legend(
        direction      = "horizontal",
        title.position = "left",
        nrow           = 1,
        keywidth       = grid::unit(0.45, "cm"),
        keyheight      = grid::unit(0.35, "cm")
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "top",
      legend.text     = ggplot2::element_text(size = fontsize),
      legend.title    = ggplot2::element_text(size = fontsize + 1)
    )
  cowplot::get_legend(p)
}

# Render a pheatmap with the colour-scale + annotation-key legends laid out
# horizontally on top, instead of pheatmap's default vertical stack on the
# right. Frees horizontal canvas for the matrix body — needed for the
# figS_diet_effects panel D supps (2026-06-16).
#
# Strategy: caller passes ph (from pheatmap with silent=TRUE) plus the
# colour palette + ann_colors so we can rebuild the legends as proper
# horizontal grobs via ggplot guides. We then gtable_filter pheatmap's
# original (vertical) legends out of the body, stack our horizontal band
# above the body, and grid::grid.draw onto a fresh PNG device.
.pheatmap_save_legend_top <- function(ph, file, width_in, height_in,
                                       palette, ann_colors,
                                       dpi                 = 300,
                                       legend_title        = "value",
                                       fontsize_legend     = 8,
                                       fontsize_annotation = 8) {
  has_deps <- requireNamespace("gtable", quietly = TRUE) &&
              requireNamespace("gridExtra", quietly = TRUE) &&
              requireNamespace("cowplot", quietly = TRUE)
  if (!has_deps) return(invisible(FALSE))

  body <- gtable::gtable_filter(ph$gtable, "legend|annotation_legend",
                                 invert = TRUE)

  gradient_grob <- .make_gradient_legend_grob(palette, limits = c(0, 1),
                                                title    = legend_title,
                                                fontsize = fontsize_legend)
  cat_grobs <- lapply(names(ann_colors), function(nm) {
    .make_categorical_legend_grob(ann_colors[[nm]], nm,
                                  fontsize = fontsize_annotation)
  })

  # Layout policy:
  #   0 categorical: gradient alone (1 row, 1.0cm).
  #   1 categorical: gradient + cat side-by-side (1 row, 1.2cm).
  #   ≥2 categorical: each cat on its own row below the gradient. Single-
  #     row packing breaks once the combined key count exceeds ~10 — the
  #     ARG gene heatmap has 4 (Treatment_Bird) + 13 (DRUG) keys and the
  #     titles collide when arrayed horizontally on one row.
  if (length(cat_grobs) == 0) {
    top_band    <- gradient_grob
    top_band_cm <- 1.0
  } else if (length(cat_grobs) == 1) {
    top_band    <- gridExtra::arrangeGrob(grobs = c(list(gradient_grob),
                                                     cat_grobs),
                                            nrow  = 1)
    top_band_cm <- 1.2
  } else {
    top_band    <- gridExtra::arrangeGrob(
      grobs = c(list(gradient_grob), cat_grobs),
      ncol  = 1
    )
    top_band_cm <- 1.0 + 0.9 * length(cat_grobs)
  }

  combined <- gridExtra::arrangeGrob(
    top_band, body,
    ncol    = 1,
    heights = grid::unit.c(grid::unit(top_band_cm, "cm"),
                            grid::unit(1, "null"))
  )

  grDevices::png(file, width = width_in, height = height_in,
                 units = "in", res = dpi, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.draw(combined)
  invisible(TRUE)
}

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
# Also appends the raw p-value to the cross-module accumulator (surface =
# "ge_alpha_kw") with `family = metric`, so resistome / virulome / mobilome
# per-metric KW tests can be family-adjusted across domains by
# `padj_summary_finalise(cfg, "ge_alpha_kw")` once all GE-side modules
# have run.
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
    tryCatch(
      padj_summary_record(cfg, surface = "ge_alpha_kw",
                          family = metric, key = log_label,
                          raw_p  = res$p.value),
      error = function(e) NULL
    )
  }
  res
}

# Min-max scale each column of a matrix independently to [0, 1]. Columns
# with zero range collapse to all-zero. Matches the GT pheatmap recipe.
#
# Explicitly re-attaches dimnames after apply() — when any column has zero
# range and FUN returns an unnamed rep(0, n), apply silently strips ALL
# rownames from the assembled matrix (because consistency of names across
# columns can't be guaranteed). That break manifested in the VF gene
# heatmap as missing row labels + a dropped annotation_row sidebar (rownames
# = NULL meant cat_per_gene[rownames(.), ] aligned to nothing).
ge_minmax_per_col <- function(mat) {
  result <- apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
  if (is.matrix(result) && all(dim(result) == dim(mat))) {
    dimnames(result) <- dimnames(mat)
  }
  result
}

# Boxplot + jittered dots of `metric` ~ `group`, fill by group. Violin
# was removed 2026-06-12 — the per-sample dots already communicate the
# distribution shape, and the box + dot combo is the reviewer-preferred
# style for small-N alpha-diversity comparisons. KW annotation drawn in
# the upper-right corner when `kw` is non-NULL.
ge_plot_alpha_violin <- function(alpha, group, metric, kw, pal_group, file,
                                  p_adj = NULL, padj_method = "BH") {
  if (!metric %in% colnames(alpha)) return(invisible(NULL))
  p <- ggplot2::ggplot(alpha,
        ggplot2::aes(x = .data[[group]], y = .data[[metric]],
                     fill = .data[[group]])) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    ggplot2::geom_jitter(width = 0.12, size = 1.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = stringr::str_to_title(metric)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 13))
  if (!is.null(kw)) {
    label_str <- if (!is.null(p_adj) && !is.na(p_adj)) {
      sprintf("KW p = %.2g   |   p_adj (%s) = %.2g",
              kw$p.value, padj_method, p_adj)
    } else {
      sprintf("Kruskal-Wallis p = %.2g", kw$p.value)
    }
    p <- p + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = label_str,
      hjust = 1.05, vjust = 1.4, size = 4.2, colour = "black"
    )
  }
  save_panel_ggplot(file, p, width = 7, height = 5, dpi = 300)
}

ge_plot_alpha_bar <- function(alpha, group, metric, kw, pal_group, file,
                               p_adj = NULL, padj_method = "BH") {
  if (!metric %in% colnames(alpha)) return(invisible(NULL))
  metric_label <- stringr::str_to_title(metric)
  title <- if (!is.null(kw)) {
    if (!is.null(p_adj) && !is.na(p_adj)) {
      sprintf("%s by %s (KW p = %.3g | p_adj (%s) = %.3g)",
              metric_label, group, kw$p.value, padj_method, p_adj)
    } else {
      sprintf("%s by %s (Kruskal-Wallis p = %.3g)",
              metric_label, group, kw$p.value)
    }
  } else {
    sprintf("%s by %s", metric_label, group)
  }
  # Per-group means for the in-facet trend line (PIPELINE_V2_GAPS A1).
  # facet_wrap dispatches the geom by matching the group column, so a
  # tibble with one row per group draws one mean line per panel without
  # cross-facet bleed.
  group_means <- alpha |>
    dplyr::group_by(.data[[group]]) |>
    dplyr::summarise(grp_mean = mean(.data[[metric]], na.rm = TRUE),
                     .groups = "drop")
  p <- ggplot2::ggplot(alpha,
        ggplot2::aes(x = sample, y = .data[[metric]],
                     fill = .data[[group]])) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = round(.data[[metric]], 1)),
                       vjust = -0.6, size = 3) +
    ggplot2::geom_hline(yintercept = mean(alpha[[metric]], na.rm = TRUE),
                        linetype = "dashed", colour = "red") +
    ggplot2::geom_hline(data = group_means,
                        ggplot2::aes(yintercept = .data$grp_mean),
                        linetype = "dashed", colour = "black",
                        linewidth = 0.7) +
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
  save_panel_ggplot(file, p, width = bw, height = 5, dpi = 300)
}

# Per-(sample, gene) log(value + 1) abundance plot coloured by group,
# with Kruskal-Wallis global p annotated in the upper-right corner AND
# pairwise Wilcoxon comparisons bracketed above the boxplots
# (asterisks: * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001).
#
# Style: box + jittered dots (no violin envelope) — same as alpha-side
# polish 2026-06-12. Pairwise comparisons run via ggpubr::compare_means
# with cfg$stats$padjust_method applied across all pairs in the family.
# Skipped silently when ggpubr isn't available or < 2 groups present.
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
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
    ggplot2::geom_jitter(width = 0.12, size = 0.6, alpha = 0.3) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = y_label, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold",
                                                       size = 13),
                   text = ggplot2::element_text(size = 13))
  # KW annotation in upper-LEFT corner. Right/centre collide with the
  # pairwise significance brackets added below, which stack vertically
  # above the box-plot from the centre x of each pair.
  if (!is.null(kw)) {
    p <- p + ggplot2::annotate(
      "text", x = -Inf, y = Inf,
      label = sprintf("KW p = %.2g", kw$p.value),
      hjust = -0.05, vjust = 1.4, size = 3.4, colour = "black"
    )
  }
  # Pairwise Wilcoxon comparisons — significant pairs only. Pre-compute
  # the family-adjusted p-values via ggpubr::compare_means (so the BH /
  # Holm / whatever from cfg$stats$padjust_method runs across ALL pairs),
  # filter to p.adj < cfg$stats$alpha, then draw only those via
  # stat_pvalue_manual. Pre-filtering is necessary because ggpubr's
  # stat_compare_means(hide.ns = TRUE) hides the "ns" label but still
  # draws the empty bracket, which clutters panel C with floating lines
  # for every non-significant pair.
  group_levels <- sort(unique(as.character(v[[group]])))
  if (requireNamespace("ggpubr", quietly = TRUE) &&
      length(group_levels) >= 2) {
    pad_method <- padjust_method(cfg)
    alpha_thr  <- cfg$stats$alpha %||% 0.05
    cm <- tryCatch(
      ggpubr::compare_means(
        stats::reformulate(group, "log_val"),
        data            = v,
        method          = "wilcox.test",
        p.adjust.method = pad_method
      ),
      error = function(e) NULL
    )
    cm_sig <- if (!is.null(cm)) cm[cm$p.adj < alpha_thr, , drop = FALSE]
              else cm
    if (!is.null(cm_sig) && nrow(cm_sig) > 0) {
      # Anchor brackets to the per-group UPPER WHISKER (Q3 + 1.5*IQR),
      # not max(v$log_val). Using the raw max means a single high
      # outlier pushes brackets way above the box cloud.
      y_top <- max(tapply(v$log_val, v[[group]], function(x) {
        q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
        as.numeric(q[2] + 1.5 * (q[2] - q[1]))
      }), na.rm = TRUE)
      y_rng <- diff(range(v$log_val, na.rm = TRUE))
      steps <- y_rng * 0.02   # bracket stack step
      cm_sig$y.position <- y_top + steps * seq_len(nrow(cm_sig))
      p <- p + ggpubr::stat_pvalue_manual(
        data       = cm_sig,
        label      = "p.signif",
        tip.length = 0.005,
        vjust      = 0.4,
        size       = 5
      )
      # Headroom so the highest bracket clears the panel top.
      p <- p + ggplot2::expand_limits(
        y = y_top + steps * (nrow(cm_sig) + 1)
      )
    }
  }
  save_panel_ggplot(file, p, width = 7, height = 5.5, dpi = 300)
}

# UpSet plot of unique category values per group. Scales better than a
# Venn when the number of sets or intersections is large. Skipped (with a
# log line) if UpSetR isn't installed. `nintersects` caps the bars shown
# (NA = unlimited). Mirrors ge_plot_category_venn's input shape so the
# same long-form table works for both renders. When `pal_group` is
# supplied (named character vector keyed by group level), the per-set
# size bars (lower-left of the UpSet figure) are coloured by treatment.
ge_plot_category_upset <- function(df, category_col, group, file,
                                    log_label, cfg, value_col = "TPM",
                                    nintersects = NA, pal_group = NULL,
                                    value_label = "Intersection Size") {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    pipeline_log(cfg, sprintf("%s: UpSetR not available — UpSet skipped",
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
  if (length(sets) < 2) {
    pipeline_log(cfg, sprintf("%s: <2 groups — UpSet skipped", log_label))
    return(invisible(NULL))
  }
  df_upset <- UpSetR::fromList(sets)
  set_names_ord <- names(sets)
  # Align palette to set order; fall back to UpSetR's default grey when
  # the palette is incomplete or absent.
  sets_bar_color <- "gray23"
  if (!is.null(pal_group)) {
    aligned <- unname(pal_group[set_names_ord])
    if (length(aligned) == length(set_names_ord) && !any(is.na(aligned))) {
      sets_bar_color <- aligned
    }
  }
  grDevices::png(file, width = 2200, height = 1500, res = 300, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(UpSetR::upset(
    df_upset,
    sets             = set_names_ord,
    keep.order       = TRUE,
    order.by         = "freq",
    nintersects      = if (is.na(nintersects)) NA else as.integer(nintersects),
    sets.bar.color   = sets_bar_color,
    mainbar.y.label  = value_label,
    text.scale       = c(1.4, 1.2, 1.2, 1.2, 1.3, 1.1)
  ))
  invisible(NULL)
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
  # Display-only label cleanup: strip underscores from the group-level
  # category labels (Control_W4 -> Control W4) and from the heading. The
  # underlying `sets` names and the manifest entry are untouched — this
  # only affects what's drawn on the canvas. Reads cleaner in print and
  # matches the figure-caption convention of using spaces in column refs.
  display_names <- gsub("_", " ", names(sets), fixed = TRUE)
  display_title <- gsub("_", " ", main_title,  fixed = TRUE)
  do.call(VennDiagram::venn.diagram, c(
    list(
      x         = sets,
      category.names = display_names,
      filename  = file,
      imagetype = "png",
      height = 2800, width = 2800, resolution = 300,
      fill = fill_pal, alpha = 0.6,
      # cat.cex 1.0 -> 1.5 and cex 1.0 -> 1.4 so the per-set labels and
      # intersection counts match the visual weight of axis text in
      # neighbouring panels (boxplots, bar charts) when rendered into a
      # 2x2 composite.
      cat.cex = 1.5, cat.fontface = "bold",
      cex = 1.4,
      main = display_title,
      # main.cex bumped 1.0 -> 1.6 so the heading carries weight at
      # composite scale; main.pos[2] lowered 1.05 -> 0.98 to close the
      # white-space gap between heading and diagram; margin tightened
      # 0.18 -> 0.08 so the diagram fills more of the canvas.
      main.cex = 1.6, main.pos = c(0.5, 0.98),
      margin = 0.08,
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
                                  top_n               = NULL,
                                  rank_by             = "abundance",
                                  fontsize_row        = 8,
                                  row_height_factor   = 0.15,
                                  min_height          = 6,
                                  width               = NULL,
                                  height              = NULL,
                                  legend_top          = FALSE,
                                  legend_title        = "value",
                                  fontsize_legend     = 8,
                                  fontsize_annotation = 8) {
  gene_wide <- df |>
    dplyr::group_by(.data$GENE, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("GENE") |>
    as.matrix()
  # Row ordering / top_n ranking strategy:
  #   abundance               (default) sum(TPM) across samples
  #   prevalence              count of samples where TPM > 0
  #   prevalence_x_abundance  composite — handles dominant-but-rare AND
  #                           widespread-but-low cases evenly
  rank_score <- switch(rank_by,
    abundance              = rowSums(gene_wide),
    prevalence             = rowSums(gene_wide > 0),
    prevalence_x_abundance = rowSums(gene_wide > 0) * rowSums(gene_wide),
    rowSums(gene_wide)
  )
  gene_wide <- gene_wide[order(rank_score, decreasing = TRUE), , drop = FALSE]
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
  grp_palette <- pal_group[intersect(names(pal_group), grp_lvls)]
  if (length(grp_palette) > 0) {
    ann_colors[[group]] <- grp_palette
  }
  cat_levels <- unique(stats::na.omit(cat_per_gene[[category_col]]))
  cat_in_pal <- intersect(names(pal_category), cat_levels)
  if (length(cat_in_pal) > 0) {
    ann_colors[[category_col]] <- pal_category[cat_in_pal]
  } else {
    # Top-N subset had no genes with a recognised category (e.g. all
    # NA Functions on the VF side after prevalence ranking). Drop the
    # row annotation rather than crashing pheatmap with a 0-length gpar.
    cat_per_gene <- NA
  }

  # Display-only: strip underscores from the annotation header labels
  # (Treatment_Bird -> Treatment Bird, Replicon_Family -> Replicon
  # Family). pheatmap uses the data-frame column name as the legend /
  # annotation-strip title, so renaming the column at the last moment
  # is the simplest path. The data fetching above keeps using actual
  # column names; only the labels passed to pheatmap get cleaned.
  display_group <- gsub("_", " ", group,        fixed = TRUE)
  display_cat   <- gsub("_", " ", category_col, fixed = TRUE)
  names(sample_group)[names(sample_group) == group] <- display_group
  if (is.data.frame(cat_per_gene)) {
    names(cat_per_gene)[names(cat_per_gene) == category_col] <- display_cat
  }
  names(ann_colors)[names(ann_colors) == group]        <- display_group
  names(ann_colors)[names(ann_colors) == category_col] <- display_cat

  hm_width  <- width  %||% max(6, 0.32 * ncol(gene_scaled) + 3)
  hm_height <- height %||% max(min_height, row_height_factor * nrow(gene_scaled) + 2.5)
  ph_args <- list(
    mat               = gene_scaled,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    color             = grDevices::colorRampPalette(palette)(100),
    border_color      = NA,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    annotation_row    = cat_per_gene,
    annotation_col    = sample_group,
    annotation_colors = if (length(ann_colors) > 0) ann_colors else NA,
    annotation_legend = TRUE,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    gaps_col          = gaps_col,
    fontsize_row = fontsize_row, fontsize_col = 8, fontsize = 8
  )
  if (isTRUE(legend_top)) {
    ph <- do.call(pheatmap::pheatmap, c(ph_args, list(silent = TRUE)))
    .pheatmap_save_legend_top(ph, file, hm_width, hm_height,
                               palette    = grDevices::colorRampPalette(palette)(100),
                               ann_colors = if (length(ann_colors) > 0) ann_colors
                                            else list(),
                               legend_title        = legend_title,
                               fontsize_legend     = fontsize_legend,
                               fontsize_annotation = fontsize_annotation)
  } else {
    do.call(pheatmap::pheatmap, c(ph_args,
                                    list(filename = file,
                                         width    = hm_width,
                                         height   = hm_height)))
  }
}

# Category-level pheatmap: rows aggregated by `category_col`, columns =
# sample. Rows sorted by total TPM descending; columns sorted by group.
#
# Defaults retuned 2026-06-16 (option A — source-side shrink) so the raster
# survives the ~3-4x downscale into figS_diet_effects panel D. Width/height
# now locked to fixed defaults so the three panel D heatmaps render at
# uniform composite size regardless of nrow. Pass non-NULL width/height
# to override per caller.
ge_plot_category_heatmap <- function(df, category_col, group, pal_group,
                                      file, palette,
                                      fontsize_row        = 12,
                                      fontsize_col        = NULL,
                                      row_height_factor   = 0.28,
                                      min_height          = 4,
                                      width               = 8.5,
                                      height              = 6,
                                      legend_top          = TRUE,
                                      legend_title        = "value",
                                      fontsize_legend     = 8,
                                      fontsize_annotation = 8) {
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

  # Display-only: strip underscores from the annotation header label
  # (Treatment_Bird -> Treatment Bird). Same rationale as the gene-
  # heatmap counterpart — data fetching keeps the actual column name;
  # only the label passed to pheatmap is cleaned.
  display_group <- gsub("_", " ", group, fixed = TRUE)
  names(sample_group)[names(sample_group) == group] <- display_group
  names(ann_colors)[names(ann_colors) == group]    <- display_group

  hm_width  <- width  %||% max(6, 0.32 * ncol(cat_scaled) + 3)
  hm_height <- height %||% max(min_height, row_height_factor * nrow(cat_scaled) + 2.5)
  # Auto-scale column fontsize when the caller doesn't override it: shrink
  # below 11 once the sample count crosses ~18 so labels at the canvas
  # bottom don't overlap. 32 samples (chicken_batch2) lands at ~6 pt.
  fs_col <- fontsize_col %||% max(4, min(11, round(220 / max(ncol(cat_scaled), 1L))))
  ph_args <- list(
    mat               = cat_scaled,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    color             = grDevices::colorRampPalette(palette)(100),
    border_color      = NA,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = fontsize_row, fontsize_col = fs_col, fontsize = 11
  )
  if (isTRUE(legend_top)) {
    ph <- do.call(pheatmap::pheatmap, c(ph_args, list(silent = TRUE)))
    .pheatmap_save_legend_top(ph, file, hm_width, hm_height,
                               palette             = grDevices::colorRampPalette(palette)(100),
                               ann_colors          = ann_colors,
                               legend_title        = legend_title,
                               fontsize_legend     = fontsize_legend,
                               fontsize_annotation = fontsize_annotation)
  } else {
    do.call(pheatmap::pheatmap, c(ph_args,
                                    list(filename = file,
                                         width    = hm_width,
                                         height   = hm_height)))
  }
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
  save_panel_ggplot(fig_file, p, width = bw, height = 6, dpi = 300)

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
                       hjust = -0.05, size = 2.8) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  save_panel_ggplot(file, p, width = width, height = height, dpi = 300)
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
  pipeline_log(cfg, sprintf("%s PERMANOVA %s: R2 = %.3f, p_raw = %.4g",
                            log_label, group, r2, pv))

  bd_test <- tryCatch({
    bd <- vegan::betadisper(d, factor(meta[[group]]))
    vegan::permutest(bd, permutations = perms)
  }, error = function(e) NULL)
  permdisp_p <- NA_real_
  if (!is.null(bd_test)) {
    utils::capture.output(bd_test, file = ds_disp_file)
    permdisp_p <- bd_test$tab$`Pr(>F)`[1]
    pipeline_log(cfg, sprintf("%s PERMDISP %s: p_raw = %.4g",
                              log_label, group, permdisp_p))
  }

  # Family-adjust the (PERMANOVA, PERMDISP) pair within this domain. Both
  # test beta-diversity differences (location vs dispersion); reporting
  # adjusted alongside raw gives the reader the conservative reading
  # without hiding the underlying numbers.
  pad_method <- padjust_method(cfg)
  raw_betas  <- c(permanova = pv, permdisp = permdisp_p)
  adj_betas  <- padjust_p(raw_betas, cfg)
  pv_adj         <- adj_betas[["permanova"]]
  permdisp_p_adj <- adj_betas[["permdisp"]]
  pipeline_log(cfg, sprintf(
    "%s beta p_adj (%s): PERMANOVA = %.4g, PERMDISP = %.4g",
    log_label, pad_method, pv_adj, permdisp_p_adj
  ))

  pcoa <- stats::cmdscale(d, eig = TRUE, k = 2)
  var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
  scores <- data.frame(
    sample = rownames(mat_t),
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2]
  )
  scores <- dplyr::left_join(scores, meta, by = "sample")
  readr::write_csv(scores, ds_pcoa_file)

  annot <- sprintf(
    "PERMANOVA R² = %.3f, p_adj (%s) = %.4g\nPERMDISP p_adj = %.4g",
    r2, pad_method, pv_adj, permdisp_p_adj
  )
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
  save_panel_ggplot(fig_file, p, width = 7, height = 5.5, dpi = 300)

  invisible(list(scores = scores, permanova = permanova,
                 permdisp = bd_test, var_explained = var_expl[1:2]))
}
