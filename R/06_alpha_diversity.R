# 06_alpha_diversity.R — alpha diversity per sample, tested across groups with
# Kruskal-Wallis and visualised as violin (+ inner boxplot) and an optional
# per-sample bar plot faceted by group. Mirrors alpha_diversity_violinPlot.R
# from the reference chicken_batch1 scripts but parameterised by config — no
# hardcoded sample IDs, treatments, palettes, or metric names.
#
# Data source priority:
#   1. cfg$inputs$diversity_csv — pre-computed indices (e.g. Oxford Nanopore
#      wf-metagenomics output). Each row is a metric, each non-Indices column
#      is a sample, optionally a trailing summary column ("total") that's
#      dropped before joining metadata.
#   2. Fallback: compute richness/shannon/simpson on the fly from
#      cleaned$noncontaminants using vegan.
#
# Selection is forced via cfg$alpha_diversity$source ("auto"|"precomputed"|"computed").

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build a named palette for the group levels. Mirrors group_palette() in R/04
# but lives here to avoid coupling stages. Falls back to ggsci::default_nejm.
alpha_palette <- function(levels, cfg) {
  # 1. Top-level cfg$colors[[group_var]] (centralised resolver)
  group_var <- cfg$metadata$group_cols[[1]]
  top <- resolve_top_level_colors(group_var, levels, cfg)
  if (!is.null(top)) return(top)
  # 2. Legacy per-stage / taxonomy fallback
  user_map <- cfg$alpha_diversity$group_colors %||%
              cfg$taxonomy$treatment_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  name <- cfg$alpha_diversity$palette %||% "ggsci::default_nejm"
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(length(levels), palette = "Dark 3")
  }
  pal <- rep_len(pal, length(levels))
  names(pal) <- levels
  pal
}

# Reshape the precomputed diversity table into a wide per-sample tibble.
# Input shape: first column "Indices" (metric names), remaining columns one
# per sample (+ optional "total" trailing column). Returns NULL if the file
# is missing or the metric column isn't recognisable.
load_precomputed_diversity <- function(cfg) {
  src <- cfg$inputs$diversity_csv
  if (is.null(src)) return(NULL)
  path <- file.path(cfg$project_root, src)
  if (!file.exists(path)) return(NULL)

  raw <- readr::read_csv(path, show_col_types = FALSE)
  idx_col <- intersect(c("Indices", "Index", "Metric"), colnames(raw))
  if (length(idx_col) == 0) return(NULL)
  idx_col <- idx_col[[1]]

  # Drop trailing summary columns (e.g. "total") — not real samples.
  drop_cols <- intersect(c("total", "Total", "TOTAL"), colnames(raw))
  # Also drop negative-control columns. R/02 strips controls from the
  # cleaned counts table during baseline subtraction, but this loader
  # reads the precomputed diversity CSV directly from disk and bypasses
  # that. Leaving controls in produces two problems: (1) wf-metagenomics
  # emits string sentinels like "None" for some indices on controls,
  # which makes pivot_longer fail with a "Can't combine <character> and
  # <double>" vctrs error; (2) even if it didn't, controls don't belong
  # in alpha-diversity comparisons.
  ctrl_cols <- intersect(cfg$metadata$controls %||% character(0),
                         colnames(raw))
  drop_cols <- unique(c(drop_cols, ctrl_cols))
  if (length(drop_cols) > 0) raw <- raw[, setdiff(colnames(raw), drop_cols)]

  long <- raw |>
    tidyr::pivot_longer(cols = -dplyr::all_of(idx_col),
                        names_to = "sample", values_to = "value") |>
    dplyr::mutate(value = suppressWarnings(as.numeric(value)))
  wide <- tidyr::pivot_wider(long, names_from = dplyr::all_of(idx_col),
                             values_from = value)
  wide
}

# Compute richness/shannon/simpson from a long-form count table.
compute_diversity_from_counts <- function(df) {
  mat <- counts_to_matrix(df)
  data.frame(
    sample   = rownames(mat),
    richness = rowSums(mat > 0),
    shannon  = vegan::diversity(mat, index = "shannon"),
    simpson  = vegan::diversity(mat, index = "simpson")
  )
}

# (sample × species) integer matrix from the long-form count table. Shared
# between compute_diversity_from_counts and the Chao1 augmentation step
# (PIPELINE_V2_GAPS C2) so both see the same per-sample integer counts.
counts_to_matrix <- function(df) {
  wide <- df |>
    dplyr::filter(!is.na(name)) |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$sample
  # vegan::estimateR requires non-negative integers.
  storage.mode(mat) <- "integer"
  mat
}

# Per-sample Chao1 from a long-form count table via vegan::estimateR.
# Returns a (sample, chao1) tibble. Used to augment div_wide regardless of
# whether the rest of the alpha metrics came from a precomputed CSV or the
# vegan fallback path — wf-metagenomics doesn't ship Chao1, so this is the
# only way to satisfy the v1.2 fig2_alpha_diversity slot's Chao1 panel.
compute_chao1_from_counts <- function(df) {
  mat <- counts_to_matrix(df)
  est <- tryCatch(
    suppressWarnings(vegan::estimateR(mat)),
    error = function(e) NULL
  )
  if (is.null(est)) return(NULL)
  tibble::tibble(
    sample = colnames(est),
    chao1  = as.numeric(est["S.chao1", ])
  )
}

# Convert "Shannon diversity index" -> "shannon_diversity_index" so metric
# names survive as ggplot2 column refs without quoting.
slugify_metric <- function(x) {
  s <- tolower(trimws(x))
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("^_|_$", "", s)
  s
}

run_alpha_diversity <- function(cleaned, cfg) {
  pipeline_log(cfg, "Alpha diversity")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "alpha_diversity")
  out_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ad_cfg <- cfg$alpha_diversity %||% list()
  source <- ad_cfg$source %||% "auto"
  sid    <- cfg$metadata$sample_id_col
  group  <- cfg$metadata$group_cols[[1]]

  # ---- Pick the data source ----------------------------------------------
  div_wide <- NULL
  if (source %in% c("auto", "precomputed")) {
    div_wide <- load_precomputed_diversity(cfg)
    if (!is.null(div_wide)) {
      pipeline_log(cfg, sprintf(
        "Alpha diversity: using precomputed table (%d samples x %d metrics)",
        nrow(div_wide), ncol(div_wide) - 1
      ))
    } else if (source == "precomputed") {
      pipeline_log(cfg, "Alpha diversity: precomputed forced but unavailable — skipping")
      return(invisible(NULL))
    }
  }
  if (is.null(div_wide)) {
    pipeline_log(cfg, "Alpha diversity: computing richness/shannon/simpson via vegan")
    div_wide <- compute_diversity_from_counts(cleaned$noncontaminants)
  }

  # ---- Augment with Chao1 (PIPELINE_V2_GAPS C2) --------------------------
  # The precomputed wf-metagenomics CSV doesn't ship Chao1; compute it from
  # cleaned counts and merge in so the v1.2 fig2_alpha_diversity slot gets
  # the metric regardless of which source supplied the others.
  chao_existing <- intersect(c("chao1", "Chao1", "chao_1"), colnames(div_wide))
  if (length(chao_existing) == 0 && !is.null(cleaned$noncontaminants)) {
    chao <- compute_chao1_from_counts(cleaned$noncontaminants)
    if (!is.null(chao)) {
      div_wide <- dplyr::left_join(div_wide, chao, by = "sample")
      pipeline_log(cfg, sprintf(
        "Alpha diversity: appended Chao1 estimator (%d/%d samples covered)",
        sum(!is.na(div_wide$chao1)), nrow(div_wide)
      ))
    }
  }

  # ---- Slugify metric column names for safe ggplot2 referencing ----------
  raw_metric_cols <- setdiff(colnames(div_wide), "sample")
  metric_lookup   <- setNames(slugify_metric(raw_metric_cols), raw_metric_cols)
  colnames(div_wide)[match(raw_metric_cols, colnames(div_wide))] <-
    unname(metric_lookup)
  pretty_labels <- setNames(raw_metric_cols, unname(metric_lookup))

  # ---- Choose which metrics to plot --------------------------------------
  available <- setdiff(colnames(div_wide), "sample")
  requested <- ad_cfg$metrics
  if (!is.null(requested)) {
    # Match against either the slug or the pretty label.
    keep <- unique(c(
      intersect(requested, available),
      unname(metric_lookup[intersect(requested, names(metric_lookup))])
    ))
    if (length(keep) == 0) {
      pipeline_log(cfg, sprintf(
        "Alpha diversity: none of cfg$alpha_diversity$metrics matched (%s) — using all available",
        paste(requested, collapse = ", ")
      ))
      keep <- available
    }
    metrics <- keep
  } else {
    metrics <- available
  }

  # ---- Join metadata -----------------------------------------------------
  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  div_wide$sample <- as.character(div_wide$sample)
  meta[[sid]]     <- as.character(meta[[sid]])
  div <- dplyr::inner_join(div_wide, meta, by = c("sample" = sid))
  if (!group %in% colnames(div)) {
    stop(sprintf("Alpha diversity: group column '%s' not in metadata", group))
  }
  if (nrow(div) == 0) {
    pipeline_log(cfg, "Alpha diversity: no samples after metadata join — skipping")
    return(invisible(NULL))
  }

  readr::write_csv(div, file.path(out_dir, "alpha_diversity.csv"))

  # ---- Kruskal-Wallis + plots per metric ---------------------------------
  group_levels <- sort(unique(as.character(div[[group]])))
  pal          <- alpha_palette(group_levels, cfg)
  show_bar     <- ad_cfg$bar_plot %||% TRUE
  stats_lines  <- character()

  for (m in metrics) {
    label <- pretty_labels[[m]] %||% m
    sub   <- div[!is.na(div[[m]]), , drop = FALSE]
    if (nrow(sub) < 2 || dplyr::n_distinct(sub[[group]]) < 2) {
      pipeline_log(cfg, sprintf("Alpha diversity: %s — not enough data for KW", label))
      next
    }
    kw <- tryCatch(
      kruskal.test(reformulate(group, m), data = sub),
      error = function(e) NULL
    )
    p_val <- if (is.null(kw)) NA_real_ else kw$p.value
    label_kw <- sprintf("Kruskal-Wallis p = %s",
                        if (is.na(p_val)) "NA" else format(p_val, digits = 3))
    line <- sprintf("%s ~ %s: p = %s", label, group,
                    if (is.na(p_val)) "NA" else format(p_val, digits = 4))
    stats_lines <- c(stats_lines, line)
    pipeline_log(cfg, paste("Alpha diversity:", line))

    # Violin + inner boxplot. KW annotation parked in the upper-right
    # corner with a small inset so it never overlaps the violin caps.
    p_violin <- ggplot2::ggplot(sub,
        ggplot2::aes(x = .data[[group]], y = .data[[m]], fill = .data[[group]])) +
      ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
      ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                            position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_jitter(width = 0.08, size = 1.4, alpha = 0.8) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = label_kw,
                        hjust = 1.05, vjust = 1.4,
                        size = 4.2, colour = "black") +
      ggplot2::labs(x = group, y = label, title = paste(label, "by", group)) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "top",
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     text       = ggplot2::element_text(size = 13))

    ggplot2::ggsave(file.path(fig_dir, paste0(m, "_violin.png")),
                    p_violin, width = 7, height = 5, dpi = 300)

    # Per-sample bar plot, faceted by group, with overall mean line +
    # per-group mean line (PIPELINE_V2_GAPS A1). facet_wrap matches the
    # group column on the geom_hline tibble so each panel gets its own
    # group mean.
    if (isTRUE(show_bar)) {
      mean_val <- mean(sub[[m]], na.rm = TRUE)
      group_means <- sub |>
        dplyr::group_by(.data[[group]]) |>
        dplyr::summarise(grp_mean = mean(.data[[m]], na.rm = TRUE),
                         .groups = "drop")
      p_bar <- ggplot2::ggplot(sub,
          ggplot2::aes(x = sample, y = .data[[m]], fill = .data[[group]])) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(yintercept = mean_val,
                            linetype = "dashed", colour = "red") +
        ggplot2::geom_hline(data = group_means,
                            ggplot2::aes(yintercept = .data$grp_mean),
                            linetype = "dashed", colour = "black",
                            linewidth = 0.7) +
        ggplot2::geom_text(
          ggplot2::aes(label = round(.data[[m]], 2)),
          vjust = -0.5, size = 3
        ) +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::labs(
          title = sprintf("%s by %s (%s)", label, group, label_kw),
          x = "Sample", y = label
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = "none",
          plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
          text            = ggplot2::element_text(size = 13)
        ) +
        ggplot2::facet_wrap(stats::as.formula(paste("~", group)),
                            scales = "free_x", nrow = 1)
      bw <- max(8, 0.5 * nrow(sub) + 3)
      ggplot2::ggsave(file.path(fig_dir, paste0(m, "_bar.png")),
                      p_bar, width = bw, height = 5, dpi = 300)
    }
  }

  writeLines(stats_lines, file.path(fig_dir, "alpha_stats.txt"))

  # ---- Rarefaction curves (PIPELINE_V2_GAPS C1) ---------------------------
  # vegan::rarecurve subsamples each sample down to step-spaced sizes and
  # computes expected species counts. Curves that have plateaued at actual
  # sequencing depth indicate adequate coverage; rising curves mean the
  # sample is under-sequenced. Coloured by primary grouping; endpoints
  # marked so the real per-sample depth is visible.
  if (!isFALSE(ad_cfg$rarefaction %||% TRUE) &&
      requireNamespace("vegan", quietly = TRUE) &&
      !is.null(cleaned$noncontaminants)) {
    mat_rc <- counts_to_matrix(cleaned$noncontaminants)
    mat_rc <- mat_rc[intersect(rownames(mat_rc), as.character(div$sample)), ,
                      drop = FALSE]
    if (nrow(mat_rc) >= 2) {
      depths <- rowSums(mat_rc)
      step   <- max(1L, as.integer(round(min(depths) / 100)))
      rc_long <- tryCatch(
        suppressWarnings(vegan::rarecurve(mat_rc, step = step, tidy = TRUE)),
        error = function(e) {
          pipeline_log(cfg, sprintf(
            "Rarefaction curves: vegan::rarecurve failed (%s) — skipping",
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (!is.null(rc_long)) {
        rc_long$Site <- as.character(rc_long$Site)
        meta_sub <- meta[, c(sid, group)] |> dplyr::distinct()
        meta_sub[[sid]] <- as.character(meta_sub[[sid]])
        rc_long <- dplyr::inner_join(rc_long, meta_sub,
                                      by = c("Site" = sid))
        endpoints <- rc_long |>
          dplyr::group_by(.data$Site) |>
          dplyr::slice_max(.data$Sample, n = 1) |>
          dplyr::ungroup()
        p_rc <- ggplot2::ggplot(rc_long,
                  ggplot2::aes(x = .data$Sample, y = .data$Species,
                                colour = .data[[group]],
                                group  = .data$Site)) +
          ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
          ggplot2::geom_point(data = endpoints, size = 2) +
          ggplot2::scale_colour_manual(values = pal) +
          ggplot2::labs(x = "Reads subsampled", y = "Observed species",
                        colour = group) +
          ggplot2::theme_classic() +
          ggplot2::theme(legend.position = "top",
                         text = ggplot2::element_text(size = 12))
        ggplot2::ggsave(file.path(fig_dir, "rarefaction_curves.png"),
                        p_rc, width = 9, height = 6, dpi = 300)
        pipeline_log(cfg, sprintf(
          "Rarefaction curves: %d samples, step=%d reads",
          nrow(mat_rc), step
        ))
      }
    } else {
      pipeline_log(cfg,
        "Rarefaction curves: <2 samples in count matrix — skipped")
    }
  }

  invisible(div)
}
