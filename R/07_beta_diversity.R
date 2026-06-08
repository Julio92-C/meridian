# 07_beta_diversity.R — sample-level beta diversity: Hellinger-transformed
# Bray-Curtis -> PCoA, PERMANOVA (adonis2), PERMDISP (betadisper) homogeneity
# check, and PCoA scatter with ellipses annotated by R^2 / p. Mirrors
# taxaNormProfile_PCoA.R from the reference chicken_batch1 scripts but is
# parameterised by config — no hardcoded sample IDs, treatments, or palettes.
#
# Data source priority (configurable via cfg$beta_diversity$source):
#   1. "bracken"        — cleaned$merged (Bracken table with metadata)
#   2. "noncontaminants" — cleaned$noncontaminants (Re-centrifuge counts)
#   "auto" picks bracken when available, otherwise noncontaminants.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Named palette for the group levels. Falls back to ggsci::default_igv (the
# palette GT used in PCoA.R / taxaProfile_PCoA.R).
beta_palette <- function(levels, cfg) {
  # 1. Top-level cfg$colors[[group_var]] (centralised resolver)
  group_var <- cfg$metadata$group_cols[[1]]
  top <- resolve_top_level_colors(group_var, levels, cfg)
  if (!is.null(top)) return(top)
  # 2. Legacy per-stage / taxonomy fallback
  user_map <- cfg$beta_diversity$group_colors %||%
              cfg$taxonomy$treatment_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  name <- cfg$beta_diversity$palette %||% "ggsci::default_igv"
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

# Pick the long-form (sample, name, count) table the user wants to ordinate.
# Bracken's per-sample read counts live under `sampleCount` in cleaned$merged.
beta_source_table <- function(cleaned, cfg) {
  src <- cfg$beta_diversity$source %||% "auto"
  if (src %in% c("auto", "bracken") && !is.null(cleaned$merged)) {
    df <- cleaned$merged
    if (!"sampleCount" %in% colnames(df)) {
      stop("Beta diversity: cleaned$merged has no `sampleCount` column")
    }
    df <- dplyr::transmute(df, sample, name, count = as.numeric(sampleCount))
    return(list(df = df, source = "bracken"))
  }
  list(df = cleaned$noncontaminants, source = "noncontaminants")
}

run_beta_diversity <- function(cleaned, cfg) {
  pipeline_log(cfg, "Beta diversity (PCoA + PERMANOVA + PERMDISP)")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "beta_diversity")
  out_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bd_cfg <- cfg$beta_diversity %||% list()
  sid    <- cfg$metadata$sample_id_col
  group  <- cfg$metadata$group_cols[[1]]

  # ---- 1. Pick source + clean taxa names ---------------------------------
  picked <- beta_source_table(cleaned, cfg)
  df <- clean_taxa_names(picked$df, cfg)
  pipeline_log(cfg, sprintf("Beta diversity: source=%s, %d rows after clean_taxa_names",
                            picked$source, nrow(df)))
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Beta diversity: nothing left after clean_taxa_names — skipping")
    return(invisible(NULL))
  }

  # ---- 2. Per-(taxon, sample) count floor --------------------------------
  min_count <- bd_cfg$min_count_per_sample %||%
               cfg$filters$min_count_per_sample %||% 5
  df <- dplyr::filter(df, !is.na(name), count > min_count)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Beta diversity: nothing left after count filter — skipping")
    return(invisible(NULL))
  }

  # ---- 3. Sample x taxon wide matrix -------------------------------------
  wide <- df |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = name, values_from = count, values_fill = 0)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) < 3) {
    pipeline_log(cfg, sprintf(
      "Beta diversity: only %d samples after filtering — need >= 3 for PCoA", nrow(mat)
    ))
    return(invisible(NULL))
  }

  # ---- 4. Transform + distance ------------------------------------------
  transform <- bd_cfg$transform %||% "hellinger"
  distance  <- bd_cfg$distance  %||% "bray"
  mat_t <- switch(transform,
    none      = mat,
    log       = log1p(mat),
    hellinger = vegan::decostand(mat, method = "hellinger"),
    stop(sprintf("Beta diversity: unknown transform '%s'", transform))
  )
  pipeline_log(cfg, sprintf("Beta diversity: transform=%s, distance=%s",
                            transform, distance))

  d <- vegan::vegdist(mat_t, method = distance)

  # ---- 5. PCoA (cmdscale matches GT) ------------------------------------
  pcoa <- stats::cmdscale(d, eig = TRUE, k = 2)
  var_explained <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
  scores <- data.frame(
    sample = rownames(mat_t),
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2]
  )

  # ---- 6. Join metadata --------------------------------------------------
  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  meta[[sid]]   <- as.character(meta[[sid]])
  scores$sample <- as.character(scores$sample)
  scores <- dplyr::inner_join(scores, meta, by = c("sample" = sid))
  if (!group %in% colnames(scores)) {
    stop(sprintf("Beta diversity: group column '%s' not in metadata", group))
  }
  scores[[group]] <- as.character(scores[[group]])

  readr::write_csv(scores, file.path(out_dir, "beta_diversity_pcoa.csv"))

  # ---- 7. PERMANOVA ------------------------------------------------------
  meta_aligned <- meta[match(rownames(mat_t), meta[[sid]]), , drop = FALSE]
  permanova <- vegan::adonis2(
    stats::reformulate(group, "d"),
    data         = meta_aligned,
    permutations = cfg$stats$permanova_permutations %||% 9999
  )
  utils::capture.output(permanova,
                        file = file.path(fig_dir, "permanova.txt"))
  r2_val <- permanova$R2[1]
  p_val  <- permanova$`Pr(>F)`[1]
  pipeline_log(cfg, sprintf("PERMANOVA %s: R2 = %.3f, p = %.4g", group, r2_val, p_val))

  # ---- 8. PERMDISP (homogeneity of dispersion) --------------------------
  betadisp <- vegan::betadisper(d, factor(meta_aligned[[group]]))
  betadisp_test <- vegan::permutest(
    betadisp,
    permutations = cfg$stats$permanova_permutations %||% 9999
  )
  utils::capture.output(betadisp_test,
                        file = file.path(fig_dir, "permdisp.txt"))
  permdisp_p <- betadisp_test$tab$`Pr(>F)`[1]
  pipeline_log(cfg, sprintf("PERMDISP %s: p = %.4g", group, permdisp_p))

  # ---- 9. PCoA scatter + ellipses --------------------------------------
  group_levels  <- sort(unique(scores[[group]]))
  pal           <- beta_palette(group_levels, cfg)
  ellipse_type  <- bd_cfg$ellipse_type     %||% "norm"
  ellipse_line  <- bd_cfg$ellipse_linetype %||% "dashed"

  annot_label <- sprintf(
    "PERMANOVA R² = %.3f, p = %.4g\nPERMDISP p = %.4g",
    r2_val, p_val, permdisp_p
  )

  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(
      title  = sprintf("%s PCoA (%s, %s-transformed) — by %s",
                       toupper(distance), distance, transform, group),
      colour = group,
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    ggplot2::annotate(
      "text",
      x = Inf, y = Inf,
      label = annot_label,
      hjust = 1.05, vjust = 1.5,
      size  = 4, colour = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   plot.title = ggplot2::element_text(size = 12, hjust = 0.5))

  ggplot2::ggsave(file.path(fig_dir, "pcoa.png"),
                  p, width = 7, height = 5.5, dpi = 300)

  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    # See note in R/05 — suppress benign re-run dir.create() warnings.
    suppressWarnings(
      htmlwidgets::saveWidget(
        plotly::ggplotly(p),
        file = file.path(fig_dir, "pcoa.html"),
        selfcontained = TRUE
      )
    )
  }

  invisible(list(
    scores      = scores,
    permanova   = permanova,
    permdisp    = betadisp_test,
    var_explained = var_explained[1:2]
  ))
}
