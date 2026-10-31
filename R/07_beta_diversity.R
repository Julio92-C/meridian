# 07_beta_diversity.R — sample-level beta diversity: Hellinger-transformed
# Bray-Curtis -> PCoA, PERMANOVA (adonis2), PERMDISP (betadisper) homogeneity
# check, and PCoA scatter with ellipses annotated by R^2 / p. Mirrors
# taxaNormProfile_PCoA.R from the reference chicken_batch1 scripts but is
# parameterised by config — no hardcoded sample IDs, treatments, or palettes.
#
# Data source (cfg$beta_diversity$source):
#   "noncontaminants" (default) — cleaned$noncontaminants, the community taxa
#                                 table (every taxon, every sample).
#   "bracken" / "merged"        — cleaned$merged: Bracken sampleCount for taxa
#                                 that carry an ARG/VF/plasmid hit only, one row
#                                 per gene hit. A resistome-linked view, opt-in.
#   "auto" is kept for old configs and means "bracken" when available.

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

# Pick the long-form (sample, name, count) table to ordinate. Beta diversity of
# taxa is a COMMUNITY-level question, so the default is the community taxa table
# (as for taxa DAA in R/08). cleaned$merged is an inner join of the gene-hit
# table with Bracken: it keeps only gene-carrying taxa and repeats each taxon's
# count once per gene hit, which distorts the ordination (on the lung cohort it
# used 172 of 2,881 taxa and hid a batch difference). Opt in with
# `beta_diversity.source: bracken` (or `merged`).
beta_source_table <- function(cleaned, cfg) {
  src <- cfg$beta_diversity$source %||% "noncontaminants"
  if (src %in% c("auto", "bracken", "merged") && !is.null(cleaned$merged)) {
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
  if (picked$source == "noncontaminants") {
    picked$df <- community_filter(picked$df, cfg, "Beta diversity")
  }
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

  # Ordinate only metadata samples with a group value. The community table can
  # carry columns that are not samples (e.g. Re-centrifuge EXCLUSIVE/CTRL/SHARED
  # summary columns on chicken_batch2) or samples absent from the metadata;
  # adonis2 fails on their NA group.
  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  meta[[sid]] <- as.character(meta[[sid]])
  if (!group %in% colnames(meta)) {
    stop(sprintf("Beta diversity: group column '%s' not in metadata", group))
  }
  keep_ids <- meta[[sid]][!is.na(meta[[group]])]
  dropped  <- setdiff(rownames(mat), keep_ids)
  if (length(dropped) > 0) {
    pipeline_log(cfg, sprintf(
      "Beta diversity: dropping %d sample(s) with NA %s (table column not in metadata or blank group value)",
      length(dropped), group))
  }
  mat <- mat[rownames(mat) %in% keep_ids, , drop = FALSE]
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  # Optional minimum depth per sample (reads summed over the taxa left after the
  # count floor). Off by default. Very shallow samples give sparse, extreme
  # profiles that can drive a PERMANOVA on their own (lung: <100-read sputum
  # samples, mostly from young patients, produced an age-band effect).
  min_reads <- bd_cfg$min_sample_reads %||% 0
  if (min_reads > 0) {
    shallow <- rownames(mat)[rowSums(mat) < min_reads]
    if (length(shallow) > 0) {
      pipeline_log(cfg, sprintf(
        "Beta diversity: dropped %d sample(s) with < %s reads (min_sample_reads): %s",
        length(shallow), format(min_reads), paste(shallow, collapse = ", ")))
    }
    mat <- mat[rowSums(mat) >= min_reads, , drop = FALSE]
  }
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
  if (ncol(pcoa$points) < 2) {
    pipeline_log(cfg, sprintf(
      "Beta diversity: %s distances give < 2 PCoA axes (degenerate/identical profiles) — skipping beta diversity",
      distance
    ))
    return(invisible(NULL))
  }
  var_explained <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
  scores <- data.frame(
    sample = rownames(mat_t),
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2]
  )

  # ---- 6. Join metadata --------------------------------------------------
  scores$sample <- as.character(scores$sample)
  scores <- dplyr::inner_join(scores, meta, by = c("sample" = sid))
  if (!group %in% colnames(scores)) {
    stop(sprintf("Beta diversity: group column '%s' not in metadata", group))
  }
  scores[[group]] <- as.character(scores[[group]])

  readr::write_csv(scores, file.path(out_dir, "beta_diversity_pcoa.csv"))

  # ---- 7. PERMANOVA ------------------------------------------------------
  # Group-contrast stats (PERMANOVA / PERMDISP, and the Jaccard pair below)
  # need >= 2 populated group levels among the ordinated samples. Some
  # cohorts collapse to a single level here even when the full metadata has
  # two — e.g. a resistome-linked source table where one group carries ~no
  # ARG hits (hospital_wastewater_eg: tap water). Degrade gracefully like
  # taxonomy/relative_abundance instead of letting adonis2 -> model.matrix
  # throw "contrasts can be applied only to factors with 2 or more levels";
  # the PCoA ordination + CSV are still emitted below.
  meta_aligned <- meta[match(rownames(mat_t), meta[[sid]]), , drop = FALSE]
  n_grp <- length(unique(stats::na.omit(meta_aligned[[group]])))
  can_contrast <- n_grp >= 2
  r2_val <- NA_real_; p_val <- NA_real_; permdisp_p <- NA_real_
  permanova <- NULL; betadisp_test <- NULL   # stay defined for the return list
  # Patient-aware mode: same pseudo-F / R2 / dispersion F, p from a
  # cluster-respecting permutation null (R/utils_cluster_perm.R).
  use_cl   <- cluster_tests_enabled(cfg)
  clusters <- if (use_cl) cluster_vector(cfg, rownames(mat_t)) else NULL
  if (can_contrast && use_cl) {
    permanova <- cluster_permanova(d, meta_aligned[[group]], clusters,
                                   nperm = cluster_nperm(cfg), seed = cluster_seed(cfg))
    cluster_write_beta_txt(permanova, file.path(fig_dir, "permanova.txt"),
                           group, cfg$metadata$random_effect, "PERMANOVA (Bray-Curtis)")
    r2_val <- permanova$R2
    p_val  <- permanova$p.value
    pipeline_log(cfg, sprintf("PERMANOVA %s: R2 = %.3f, p_raw = %.4g (%s)",
                              group, r2_val, p_val, cluster_tag(permanova)))

    betadisp_test <- cluster_permdisp(d, meta_aligned[[group]], clusters,
                                      nperm = cluster_nperm(cfg), seed = cluster_seed(cfg))
    cluster_write_beta_txt(betadisp_test, file.path(fig_dir, "permdisp.txt"),
                           group, cfg$metadata$random_effect, "PERMDISP (Bray-Curtis)")
    permdisp_p <- betadisp_test$p.value
    pipeline_log(cfg, sprintf("PERMDISP %s: p_raw = %.4g (%s)",
                              group, permdisp_p, cluster_tag(betadisp_test)))
  } else if (can_contrast) {
    permanova <- vegan::adonis2(
      stats::reformulate(group, "d"),
      data         = meta_aligned,
      permutations = cfg$stats$permanova_permutations %||% 9999
    )
    utils::capture.output(permanova,
                          file = file.path(fig_dir, "permanova.txt"))
    r2_val <- permanova$R2[1]
    p_val  <- permanova$`Pr(>F)`[1]
    pipeline_log(cfg, sprintf("PERMANOVA %s: R2 = %.3f, p_raw = %.4g",
                              group, r2_val, p_val))

    # ---- 8. PERMDISP (homogeneity of dispersion) ------------------------
    betadisp <- vegan::betadisper(d, factor(meta_aligned[[group]]))
    betadisp_test <- vegan::permutest(
      betadisp,
      permutations = cfg$stats$permanova_permutations %||% 9999
    )
    utils::capture.output(betadisp_test,
                          file = file.path(fig_dir, "permdisp.txt"))
    permdisp_p <- betadisp_test$tab$`Pr(>F)`[1]
    pipeline_log(cfg, sprintf("PERMDISP %s: p_raw = %.4g", group, permdisp_p))
  } else {
    pipeline_log(cfg, sprintf(
      "Beta diversity: only %d level of '%s' among the %d ordinated samples — PERMANOVA/PERMDISP skipped (PCoA still written)",
      n_grp, group, nrow(mat_t)))
  }

  # ---- 9. PCoA scatter + ellipses --------------------------------------
  # Bray rendering is deferred until after the Jaccard branch runs so we
  # can family-adjust Bray + Jaccard PERMANOVA p-values (and likewise
  # PERMDISP) together before annotating either panel. Bray-only
  # adjustment is degenerate when Jaccard is disabled; we still surface
  # the raw p so the title stays readable.
  group_levels  <- sort(unique(scores[[group]]))
  pal           <- beta_palette(group_levels, cfg)
  ellipse_type  <- bd_cfg$ellipse_type     %||% "norm"
  ellipse_line  <- bd_cfg$ellipse_linetype %||% "dashed"

  # ---- 9.5. Jaccard PERMANOVA + PERMDISP (computed BEFORE rendering) -----
  # Compute Jaccard stats up front so the {Bray, Jaccard} PERMANOVA pair
  # (and the {Bray, Jaccard} PERMDISP pair) can be family-adjusted with
  # cfg$stats$padjust_method before either PCoA gets annotated.
  has_jac <- !isFALSE(bd_cfg$jaccard %||% TRUE)
  scores_j <- NULL; var_j <- NULL
  p_j <- NA_real_; r2_j <- NA_real_; permdisp_j_p <- NA_real_
  if (has_jac) {
    mat_pa  <- vegan::decostand(mat, method = "pa")
    # (Jaccard adonis2/betadisper are gated on can_contrast below — a
    # single-level group would otherwise throw the same contrasts error.)
    d_jac   <- vegan::vegdist(mat_pa, method = "jaccard", binary = TRUE)
    pcoa_j  <- stats::cmdscale(d_jac, eig = TRUE, k = 2)
    if (ncol(pcoa_j$points) < 2) {
      pipeline_log(cfg, "Beta diversity: Jaccard distances give < 2 PCoA axes (identical presence/absence profiles) — skipping Jaccard PCoA")
      has_jac <- FALSE
    }
  }
  if (has_jac) {
    var_j   <- pcoa_j$eig / sum(pcoa_j$eig[pcoa_j$eig > 0]) * 100
    scores_j <- data.frame(
      sample = rownames(mat_pa),
      PC1    = pcoa_j$points[, 1],
      PC2    = pcoa_j$points[, 2]
    )
    scores_j$sample <- as.character(scores_j$sample)
    scores_j <- dplyr::inner_join(scores_j, meta, by = c("sample" = sid))
    scores_j[[group]] <- as.character(scores_j[[group]])
    readr::write_csv(scores_j,
                     file.path(out_dir, "beta_diversity_jaccard_pcoa.csv"))

    if (can_contrast && use_cl) {
      perm_j <- cluster_permanova(d_jac, meta_aligned[[group]], clusters,
                                  nperm = cluster_nperm(cfg), seed = cluster_seed(cfg))
      cluster_write_beta_txt(perm_j, file.path(fig_dir, "jaccard_permanova.txt"),
                             group, cfg$metadata$random_effect, "PERMANOVA (Jaccard)")
      r2_j <- perm_j$R2
      p_j  <- perm_j$p.value
      pipeline_log(cfg, sprintf("PERMANOVA Jaccard %s: R2 = %.3f, p_raw = %.4g (%s)",
                                group, r2_j, p_j, cluster_tag(perm_j)))

      bd_j_test <- cluster_permdisp(d_jac, meta_aligned[[group]], clusters,
                                    nperm = cluster_nperm(cfg), seed = cluster_seed(cfg))
      cluster_write_beta_txt(bd_j_test, file.path(fig_dir, "jaccard_permdisp.txt"),
                             group, cfg$metadata$random_effect, "PERMDISP (Jaccard)")
      permdisp_j_p <- bd_j_test$p.value
      pipeline_log(cfg, sprintf("PERMDISP Jaccard %s: p_raw = %.4g (%s)",
                                group, permdisp_j_p, cluster_tag(bd_j_test)))
    } else if (can_contrast) {
      perm_j <- vegan::adonis2(
        stats::reformulate(group, "d_jac"),
        data         = meta_aligned,
        permutations = cfg$stats$permanova_permutations %||% 9999
      )
      utils::capture.output(perm_j,
                            file = file.path(fig_dir, "jaccard_permanova.txt"))
      r2_j <- perm_j$R2[1]
      p_j  <- perm_j$`Pr(>F)`[1]
      pipeline_log(cfg, sprintf("PERMANOVA Jaccard %s: R2 = %.3f, p_raw = %.4g",
                                group, r2_j, p_j))

      bd_j      <- vegan::betadisper(d_jac, factor(meta_aligned[[group]]))
      bd_j_test <- vegan::permutest(
        bd_j, permutations = cfg$stats$permanova_permutations %||% 9999
      )
      utils::capture.output(bd_j_test,
                            file = file.path(fig_dir, "jaccard_permdisp.txt"))
      permdisp_j_p <- bd_j_test$tab$`Pr(>F)`[1]
      pipeline_log(cfg, sprintf("PERMDISP Jaccard %s: p_raw = %.4g",
                                group, permdisp_j_p))
    }
  }

  # ---- 9.6. Family-adjust PERMANOVA + PERMDISP pairs ---------------------
  pad_method <- padjust_method(cfg)
  if (has_jac) {
    permanova_padj <- padjust_p(c(bray = p_val, jaccard = p_j), cfg)
    permdisp_padj  <- padjust_p(c(bray = permdisp_p,
                                   jaccard = permdisp_j_p), cfg)
  } else {
    permanova_padj <- c(bray = p_val,    jaccard = NA_real_)
    permdisp_padj  <- c(bray = permdisp_p, jaccard = NA_real_)
  }
  p_val_adj        <- permanova_padj[["bray"]]
  permdisp_p_adj   <- permdisp_padj[["bray"]]
  p_j_adj          <- permanova_padj[["jaccard"]]
  permdisp_j_p_adj <- permdisp_padj[["jaccard"]]

  beta_padj_summary <- data.frame(
    distance = c("bray-curtis", "jaccard"),
    test     = rep(c("PERMANOVA", "PERMDISP"), each = 1),
    stringsAsFactors = FALSE
  )
  beta_padj_summary <- rbind(
    data.frame(distance = "bray-curtis", test = "PERMANOVA",
               r2 = r2_val, p_raw = p_val, p_adj = p_val_adj,
               method = pad_method, stringsAsFactors = FALSE),
    data.frame(distance = "bray-curtis", test = "PERMDISP",
               r2 = NA_real_, p_raw = permdisp_p, p_adj = permdisp_p_adj,
               method = pad_method, stringsAsFactors = FALSE),
    data.frame(distance = "jaccard",     test = "PERMANOVA",
               r2 = r2_j, p_raw = p_j, p_adj = p_j_adj,
               method = pad_method, stringsAsFactors = FALSE),
    data.frame(distance = "jaccard",     test = "PERMDISP",
               r2 = NA_real_, p_raw = permdisp_j_p, p_adj = permdisp_j_p_adj,
               method = pad_method, stringsAsFactors = FALSE)
  )
  readr::write_csv(beta_padj_summary,
                    file.path(out_dir, "beta_diversity_padj_summary.csv"))

  annot_label <- if (can_contrast) {
    sprintf(
      "PERMANOVA R² = %.3f, p_adj (%s) = %.4g\nPERMDISP p_adj = %.4g",
      r2_val, pad_method, p_val_adj, permdisp_p_adj
    )
  } else {
    sprintf("Single '%s' level — PERMANOVA/PERMDISP not applicable", group)
  }
  cl_note <- if (use_cl && can_contrast)
    sprintf("\n(p: %s-aware permutation)", cfg$metadata$random_effect) else ""
  annot_label <- paste0(annot_label, cl_note)

  # ---- 9. Bray-Curtis PCoA scatter + ellipses --------------------------
  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(
      title  = sprintf("%s PCoA (%s, %s-transformed) — by %s",
                       toupper(distance), distance, transform, tolower(pretty_var(group))),
      colour = pretty_var(group),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    ggplot2::scale_y_continuous(expand = stats_headroom(n_text_lines(annot_label))) +
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

  save_panel_ggplot(file.path(fig_dir, "pcoa.png"),
                  p, width = 7, height = 5.5, dpi = 300)

  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    # See note in R/05 — suppress benign re-run dir.create() warnings.
    save_widget_safe(
      plotly::ggplotly(p),
      file = file.path(fig_dir, "pcoa.html"),
      selfcontained = TRUE, cfg = cfg
    )
  }

  # ---- 10. Jaccard PCoA (PIPELINE_V2_GAPS C3) -----------------------------
  # Render using family-adjusted PERMANOVA / PERMDISP p-values computed
  # above. Presence/absence distance captures rare-biosphere structure that
  # Bray-Curtis can miss when a handful of abundant taxa dominate.
  if (has_jac) {
    annot_j <- if (can_contrast) {
      sprintf(
        "PERMANOVA R² = %.3f, p_adj (%s) = %.4g\nPERMDISP p_adj = %.4g",
        r2_j, pad_method, p_j_adj, permdisp_j_p_adj
      )
    } else {
      sprintf("Single '%s' level — PERMANOVA/PERMDISP not applicable", group)
    }
    annot_j <- paste0(annot_j, cl_note)
    p_jac <- ggplot2::ggplot(scores_j,
              ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
      ggplot2::geom_point(size = 3) +
      ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                            linetype = ellipse_line) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::labs(
        title  = sprintf("Jaccard PCoA (presence/absence) — by %s", tolower(pretty_var(group))),
        colour = pretty_var(group),
        x = sprintf("PC1 (%.1f%%)", var_j[1]),
        y = sprintf("PC2 (%.1f%%)", var_j[2])
      ) +
      ggplot2::scale_y_continuous(expand = stats_headroom(n_text_lines(annot_j))) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = annot_j,
                        hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = 13),
                     plot.title = ggplot2::element_text(size = 12, hjust = 0.5))
    save_panel_ggplot(file.path(fig_dir, "jaccard_pcoa.png"),
                    p_jac, width = 7, height = 5.5, dpi = 300)
  }

  invisible(list(
    scores      = scores,
    permanova   = permanova,
    permdisp    = betadisp_test,
    var_explained = var_explained[1:2]
  ))
}
