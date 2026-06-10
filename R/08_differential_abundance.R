# 08_differential_abundance.R — ALDEx2 differential abundance.
# Mirrors aldex_overall_and_pairwise.R (taxa, Bracken_SpeciesCounts) and
# aldexGENE_overall_and_pairwise.R (gene, geneTable_rawCounts) from the
# reference chicken_batch1 scripts, but parameterised by config — no
# hardcoded sample IDs, treatments, or palettes.
#
# Two passes, both gated by cfg$differential_abundance$levels (default
# c("taxa", "gene")):
#   - "taxa": counts come from cleaned$merged$sampleCount (Bracken) when
#     available, falling back to cleaned$noncontaminants.
#   - "gene": counts derived from cleaned$abri_kraken2 by counting unique
#     contig hits per (sample, GENE) — matches the geneTable_rawCounts.csv
#     shape used by the GT script.
#
# For each pass we:
#   1. Build a feature x sample count matrix and filter low-prevalence rows.
#   2. Run an overall Kruskal-Wallis ALDEx2 (omnibus) across all groups.
#   3. Run pairwise Welch's t ALDEx2 for every group pair, save per-pair TSV.
#   4. Merge pairwise (prefixed columns) + overall + per-group raw count
#      mean ± SE into a single merged_all table.
#   5. Per-comparison publication CSV of taxa/genes passing effect/overlap
#      thresholds (falls back to top-by-|effect| when none pass).
#   6. Optional plots (cfg$differential_abundance$make_plots): per-comparison
#      ALDEx2 built-in MA/MW/volcano, custom effect-vs-q volcano, and
#      per-feature bar (raw counts faceted by group) + violin (log10 with KW
#      annotation and pairwise Wilcoxon brackets).
#
# Outputs go under <datasets_dir>/aldex2/<level>/ and
# <figures_dir>/differential_abundance/<level>/.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build a named palette for the group levels (mirrors helpers in R/05–R/07).
da_palette <- function(levels, cfg) {
  # 1. Top-level cfg$colors[[group_var]] (centralised resolver)
  group_var <- cfg$metadata$group_cols[[1]]
  top <- resolve_top_level_colors(group_var, levels, cfg)
  if (!is.null(top)) return(top)
  # 2. Legacy per-stage / taxonomy fallback
  user_map <- cfg$differential_abundance$group_colors %||%
              cfg$taxonomy$treatment_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  name <- cfg$differential_abundance$palette %||% "ggsci::default_nejm"
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

# Taxa counts matrix: prefer Bracken sampleCount via cleaned$merged, fall back
# to Re-centrifuge cleaned$noncontaminants. Returns features x samples matrix.
da_taxa_matrix <- function(cleaned, cfg) {
  src <- cfg$differential_abundance$source %||% "auto"
  if (src %in% c("auto", "bracken") && !is.null(cleaned$merged) &&
      "sampleCount" %in% colnames(cleaned$merged)) {
    df <- cleaned$merged
    df <- dplyr::transmute(df, sample, name, count = as.numeric(sampleCount))
    src_used <- "bracken"
  } else {
    df <- cleaned$noncontaminants
    src_used <- "noncontaminants"
  }
  df <- clean_taxa_names(df, cfg)
  list(df = df, source = src_used)
}

# Gene counts matrix: count unique contig hits per (sample, GENE) from
# cleaned$abri_kraken2. Each row in abri_kraken2 is already a unique
# (sequence, START, END) hit from R/02, so n() per group is the raw count.
da_gene_matrix <- function(cleaned, cfg) {
  df <- cleaned$abri_kraken2
  if (is.null(df) || !"GENE" %in% colnames(df)) return(NULL)
  df <- df |>
    dplyr::filter(!is.na(GENE), GENE != "") |>
    dplyr::group_by(sample, GENE) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::rename(name = GENE)
  list(df = df, source = "abri_kraken2")
}

# Long-form (sample, name, count) -> features x samples count matrix, filtered
# to rowSums(>0) >= min_prevalence.
da_pivot_matrix <- function(df, samples, min_prevalence) {
  wide <- df |>
    dplyr::filter(sample %in% samples) |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0)
  feat <- wide$name
  mat  <- as.matrix(wide[, -1, drop = FALSE])
  mat[is.na(mat)] <- 0
  rownames(mat) <- feat
  # Add any missing sample columns as zeros so column order matches `samples`.
  missing_cols <- setdiff(samples, colnames(mat))
  if (length(missing_cols) > 0) {
    pad <- matrix(0, nrow = nrow(mat), ncol = length(missing_cols),
                  dimnames = list(rownames(mat), missing_cols))
    mat <- cbind(mat, pad)
  }
  mat <- mat[, samples, drop = FALSE]
  keep <- rowSums(mat > 0) >= min_prevalence
  mat[keep, , drop = FALSE]
}

# Build a model.matrix from a formula string + metadata. Character columns
# are coerced to factors so model.matrix doesn't drop them. Used by the
# optional aldex.glm path (cfg$differential_abundance$design).
da_build_design <- function(design_str, meta_use, sid, level) {
  form <- tryCatch(
    stats::as.formula(design_str),
    error = function(e) stop(sprintf(
      "DA[%s] design: cannot parse formula '%s' — %s",
      level, design_str, conditionMessage(e)
    ))
  )
  vars <- all.vars(form)
  missing_vars <- setdiff(vars, colnames(meta_use))
  if (length(missing_vars) > 0) {
    stop(sprintf("DA[%s] design: metadata missing columns: %s",
                 level, paste(missing_vars, collapse = ", ")))
  }
  meta_sub <- meta_use[, c(sid, vars), drop = FALSE]
  for (v in vars) {
    if (is.character(meta_sub[[v]])) meta_sub[[v]] <- factor(meta_sub[[v]])
  }
  list(form = form, vars = vars,
       mm = stats::model.matrix(form, data = meta_sub))
}

# aldex.glm path — fits the supplied design across all samples in one shot.
# Writes glm_results.tsv (full coefficient table) and one
# candidates_<coef>.csv per non-intercept coefficient with BH-adjusted
# p < alpha.
da_run_glm <- function(level, mat, meta_use, sid, design_str,
                       mc, denom, alpha, ds_dir, feature_col, cfg) {
  ds <- da_build_design(design_str, meta_use, sid, level)
  vars <- ds$vars
  # Drop samples with NA in any design variable so model.matrix and the
  # counts matrix line up; warn if anything was dropped.
  ok <- stats::complete.cases(meta_use[, vars, drop = FALSE])
  if (!all(ok)) {
    dropped <- as.character(meta_use[[sid]])[!ok]
    pipeline_log(cfg, sprintf(
      "DA[%s] glm: dropping %d sample(s) with NA in design vars: %s",
      level, length(dropped), paste(dropped, collapse = ", ")
    ))
    meta_use <- meta_use[ok, , drop = FALSE]
    keep_cols <- intersect(colnames(mat), as.character(meta_use[[sid]]))
    mat <- mat[, keep_cols, drop = FALSE]
    ds <- da_build_design(design_str, meta_use, sid, level)
  }
  mm <- ds$mm
  if (nrow(mm) != ncol(mat)) {
    stop(sprintf("DA[%s] design: model.matrix has %d rows but counts have %d samples",
                 level, nrow(mm), ncol(mat)))
  }
  if (ncol(mm) < 2) {
    stop(sprintf("DA[%s] design: model.matrix has no non-intercept terms (only '%s')",
                 level, paste(colnames(mm), collapse = ", ")))
  }
  pipeline_log(cfg, sprintf(
    "DA[%s] glm: design=%s (coefficients: %s)",
    level, design_str, paste(colnames(mm), collapse = ", ")
  ))

  clr <- ALDEx2::aldex.clr(
    reads = mat, conds = mm, mc.samples = mc,
    denom = denom, verbose = FALSE
  )
  res <- ALDEx2::aldex.glm(clr, mm, fdr.method = "BH")
  res_df <- tibble::rownames_to_column(as.data.frame(res), feature_col)
  readr::write_tsv(res_df, file.path(ds_dir, "glm_results.tsv"))
  pipeline_log(cfg, sprintf("DA[%s] glm → glm_results.tsv (%d features)",
                            level, nrow(res_df)))

  # ALDEx2 names columns "<coef>:Est", "<coef>:pval", "<coef>:pval.padj"
  # (BH-adjusted by default when fdr.method = "BH"). The intercept is
  # written as "Intercept::Est" but "(Intercept):pval" — handle both.
  coef_names <- setdiff(colnames(mm), "(Intercept)")
  for (cf in coef_names) {
    safe <- gsub("[^A-Za-z0-9]+", "_", cf)
    est_col <- paste0(cf, ":Est")
    p_col   <- paste0(cf, ":pval")
    q_col   <- paste0(cf, ":pval.padj")
    if (!est_col %in% colnames(res_df) || !q_col %in% colnames(res_df)) {
      pipeline_log(cfg, sprintf(
        "DA[%s] glm: no Est/pval.padj columns for coefficient '%s' — skipping CSV",
        level, cf
      ))
      next
    }
    cand <- res_df[!is.na(res_df[[q_col]]) & res_df[[q_col]] < alpha,
                   , drop = FALSE]
    cand <- cand[order(cand[[q_col]]), , drop = FALSE]
    pub_cols <- intersect(c(feature_col, est_col, p_col, q_col),
                          colnames(cand))
    readr::write_csv(cand[, pub_cols, drop = FALSE],
                     file.path(ds_dir, paste0("candidates_", safe, ".csv")))
    pipeline_log(cfg, sprintf(
      "DA[%s] glm %s: %d features with BH-adj p < %.3g",
      level, cf, nrow(cand), alpha
    ))
  }
  invisible(NULL)
}

# One ALDEx2 pair call with the args common to overall+pairwise.
da_aldex_call <- function(mat, conds, test, mc, denom) {
  ALDEx2::aldex(
    reads = mat, conditions = conds,
    mc.samples = mc, test = test,
    effect = (test == "t"), CI = FALSE,
    include.sample.summary = FALSE, verbose = FALSE,
    paired.test = FALSE, denom = denom,
    iterate = FALSE, gamma = NULL
  )
}

# Per-comparison raw-count stats — mean ± SE per group as both a combined
# string column and numeric mean_/se_ columns (mirrors GT merged_all2).
da_raw_count_stats <- function(df, meta, sid, group, feature_col) {
  long <- df |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::inner_join(meta[, c(sid, group)], by = c("sample" = sid))
  names(long)[names(long) == group] <- ".group"

  stats <- long |>
    dplyr::group_by(name, .group) |>
    dplyr::summarise(
      mean_count = mean(count, na.rm = TRUE),
      se_count   = ifelse(dplyr::n() > 1,
                          sd(count, na.rm = TRUE) / sqrt(dplyr::n()),
                          NA_real_),
      .groups = "drop"
    )

  combined <- stats |>
    dplyr::mutate(
      mean_r  = round(mean_count, 2),
      se_r    = ifelse(is.na(se_count), NA_real_, round(se_count, 2)),
      mean_se = ifelse(is.na(se_r), as.character(mean_r),
                       paste0(mean_r, " ± ", se_r))
    ) |>
    dplyr::select(name, .group, mean_se) |>
    tidyr::pivot_wider(names_from = .group, values_from = mean_se,
                       values_fill = "")

  numeric <- stats |>
    dplyr::select(name, .group, mean_count, se_count) |>
    tidyr::pivot_wider(names_from = .group,
                       values_from = c(mean_count, se_count),
                       names_sep   = "_", values_fill = NA_real_)

  out <- dplyr::full_join(combined, numeric, by = "name")
  names(out)[names(out) == "name"] <- feature_col
  out
}

# Run a single level pass (taxa or gene). Returns invisible(NULL).
da_run_level <- function(level, picked, cleaned, cfg, meta, sid, group,
                         mc, denom, min_prev, effect_thresh, overlap_thresh,
                         top_n_plots, make_plots, da_cfg,
                         design_str = NULL, alpha = 0.05) {
  feature_col <- if (level == "taxa") "taxon" else "GENE"
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "aldex2", level)
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir,
                       "differential_abundance", level)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(make_plots)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  df <- picked$df
  if (is.null(df) || nrow(df) == 0) {
    pipeline_log(cfg, sprintf("DA[%s]: no data — skipping", level))
    return(invisible(NULL))
  }

  # Align samples with metadata (only keep samples that have a group label).
  meta_use <- meta[!is.na(meta[[group]]), , drop = FALSE]
  samples  <- intersect(as.character(meta_use[[sid]]),
                        as.character(df$sample))
  if (length(samples) < 4) {
    pipeline_log(cfg, sprintf(
      "DA[%s]: only %d aligned samples — need >= 4 — skipping",
      level, length(samples)
    ))
    return(invisible(NULL))
  }
  meta_use <- meta_use[match(samples, as.character(meta_use[[sid]])), , drop = FALSE]
  conds    <- as.character(meta_use[[group]])

  mat <- da_pivot_matrix(df, samples, min_prev)
  if (nrow(mat) < 2) {
    pipeline_log(cfg, sprintf(
      "DA[%s]: only %d features after min_prevalence=%d — skipping",
      level, nrow(mat), min_prev
    ))
    return(invisible(NULL))
  }
  group_levels <- sort(unique(conds))
  group_sizes  <- table(conds)
  pipeline_log(cfg, sprintf(
    "DA[%s]: source=%s, %d features x %d samples, %d group(s) (%s), mc=%d, denom=%s",
    level, picked$source, nrow(mat), ncol(mat),
    length(group_levels),
    paste(sprintf("%s=%d", names(group_sizes), as.integer(group_sizes)),
          collapse = ", "),
    mc, denom
  ))

  # ---- Optional aldex.glm path -------------------------------------------
  # If cfg$differential_abundance$design is set (formula string), run
  # ALDEx2's model-based test once and skip the KW + pairwise path. Useful
  # for adjusting on covariates / random effects (e.g. "~ Treatment + Block")
  # or interaction designs ("~ Treatment * Tissue").
  if (!is.null(design_str) && nzchar(design_str)) {
    return(da_run_glm(level, mat, meta_use, sid, design_str,
                      mc, denom, alpha, ds_dir, feature_col, cfg))
  }

  if (length(group_levels) < 2) {
    pipeline_log(cfg, sprintf("DA[%s]: < 2 groups present — nothing to test", level))
    return(invisible(NULL))
  }

  # ---- 1. Overall (omnibus) ----------------------------------------------
  # With only 2 groups the omnibus is just the single pairwise t — running
  # ALDEx2's KW on top would duplicate the pairwise output (and waste
  # ~mc.samples seconds). Skip it; pairwise handles the 2-group case.
  if (length(group_levels) >= 3) {
    overall <- tryCatch(
      da_aldex_call(mat, conds, test = "kw", mc = mc, denom = denom),
      error = function(e) { pipeline_log(cfg, sprintf("DA[%s] overall failed: %s",
                                                      level, conditionMessage(e))); NULL }
    )
    if (!is.null(overall)) {
      overall_df <- tibble::rownames_to_column(as.data.frame(overall), feature_col)
      readr::write_tsv(overall_df, file.path(ds_dir, "overall_kw.tsv"))
    } else {
      overall_df <- NULL
    }
  } else {
    pipeline_log(cfg, sprintf(
      "DA[%s]: 2 groups → skipping omnibus KW (pairwise t-test is the omnibus)",
      level
    ))
    overall_df <- NULL
  }

  # ---- 2. Pairwise t-test for every group pair ---------------------------
  pairs <- utils::combn(group_levels, 2, simplify = FALSE)
  comp_name <- function(p) gsub("[^A-Za-z0-9]+", "_", paste(p, collapse = "_vs_"))

  raw_results   <- list()
  prefixed_list <- list()
  for (pr in pairs) {
    cn   <- comp_name(pr)
    # ALDEx2's t-test needs >= 2 samples per group; otherwise it errors
    # out inside aldex.effect(). Skip the pair early instead.
    pair_sizes <- as.integer(group_sizes[pr])
    if (any(pair_sizes < 2)) {
      pipeline_log(cfg, sprintf(
        "DA[%s] %s: group sizes %s — need >= 2 per group; skipping",
        level, cn,
        paste(sprintf("%s=%d", pr, pair_sizes), collapse = ", ")
      ))
      next
    }
    keep <- conds %in% pr
    sub  <- mat[, keep, drop = FALSE]
    sub  <- sub[rowSums(sub > 0) >= min_prev, , drop = FALSE]
    if (nrow(sub) < 2) {
      pipeline_log(cfg, sprintf("DA[%s] %s: < 2 features after filter — skipping",
                                level, cn))
      next
    }
    res <- tryCatch(
      da_aldex_call(sub, conds[keep], test = "t", mc = mc, denom = denom),
      error = function(e) { pipeline_log(cfg, sprintf("DA[%s] %s failed: %s",
                                                      level, cn, conditionMessage(e))); NULL }
    )
    if (is.null(res)) next
    raw_results[[cn]] <- res
    df_pr <- tibble::rownames_to_column(as.data.frame(res), feature_col)
    readr::write_tsv(df_pr, file.path(ds_dir, paste0(cn, ".tsv")))
    pipeline_log(cfg, sprintf("DA[%s] %s → %s.tsv (%d features)",
                              level, cn, cn, nrow(df_pr)))

    # Prefix non-feature columns with comparison name for merged_all.
    prefixed <- df_pr
    other    <- setdiff(colnames(prefixed), feature_col)
    names(prefixed)[match(other, names(prefixed))] <- paste0(cn, "_", other)
    prefixed_list[[cn]] <- prefixed
  }

  if (length(prefixed_list) == 0) {
    pipeline_log(cfg, sprintf("DA[%s]: no successful pairwise comparisons", level))
    return(invisible(NULL))
  }

  # ---- 3. Merge pairwise + overall + raw count stats ---------------------
  merged <- Reduce(function(a, b) dplyr::full_join(a, b, by = feature_col),
                   prefixed_list)
  if (!is.null(overall_df)) {
    keep_overall <- intersect(c("kw.ep", "kw.eBH"), colnames(overall_df))
    if (length(keep_overall) == 0) {
      keep_overall <- setdiff(colnames(overall_df), feature_col)
    }
    merged <- dplyr::left_join(
      merged,
      overall_df[, c(feature_col, keep_overall), drop = FALSE],
      by = feature_col
    )
  }

  raw_stats <- da_raw_count_stats(df, meta_use, sid, group, feature_col)
  merged    <- dplyr::left_join(merged, raw_stats, by = feature_col)
  readr::write_tsv(merged, file.path(ds_dir, "merged_pairwise_with_overall.tsv"))

  # ---- 4. Per-comparison publication CSVs (candidates) -------------------
  per_comp_candidates <- list()
  for (cn in names(raw_results)) {
    res_df <- tibble::rownames_to_column(as.data.frame(raw_results[[cn]]),
                                         feature_col)
    qcol   <- intersect(c("wi.eBH", "we.eBH"), colnames(res_df))[1]
    pcol   <- intersect(c("wi.ep",  "we.ep"),  colnames(res_df))[1]
    res_df$abs_effect <- abs(res_df$effect)
    cand <- res_df |>
      dplyr::filter(!is.na(effect),
                    abs_effect >= effect_thresh,
                    overlap    <= overlap_thresh) |>
      dplyr::arrange(dplyr::desc(abs_effect))
    if (nrow(cand) == 0) {
      cand <- res_df |>
        dplyr::filter(!is.na(effect)) |>
        dplyr::arrange(dplyr::desc(abs_effect)) |>
        dplyr::slice_head(n = top_n_plots)
      pipeline_log(cfg, sprintf(
        "DA[%s] %s: no features pass |effect|>=%.2f & overlap<=%.2f — using top %d by |effect|",
        level, cn, effect_thresh, overlap_thresh, nrow(cand)
      ))
    }
    rab_cols <- grep("^rab\\.win\\.", colnames(cand), value = TRUE)
    pub_cols <- unique(c(feature_col, "effect", "overlap",
                         pcol, qcol, rab_cols))
    pub_cols <- intersect(pub_cols, colnames(cand))
    pub <- cand[, pub_cols, drop = FALSE]
    if ("effect"  %in% pub_cols) pub$effect  <- round(pub$effect, 3)
    if ("overlap" %in% pub_cols) pub$overlap <- round(pub$overlap, 4)
    if (!is.na(qcol) && qcol %in% pub_cols) pub[[qcol]] <- round(pub[[qcol]], 4)
    readr::write_csv(pub, file.path(ds_dir, paste0("candidates_", cn, ".csv")))
    per_comp_candidates[[cn]] <- head(cand[[feature_col]], top_n_plots)
  }

  # ---- 5. Optional plots --------------------------------------------------
  if (!isTRUE(make_plots)) {
    return(invisible(NULL))
  }

  pal <- da_palette(group_levels, cfg)
  pdf_path <- file.path(fig_dir, "aldex_plots.pdf")
  grDevices::pdf(pdf_path, width = 14, height = 8)
  pdf_open <- TRUE
  on.exit(if (pdf_open) grDevices::dev.off(), add = TRUE)

  for (cn in names(raw_results)) {
    res_obj <- raw_results[[cn]]
    res_df  <- tibble::rownames_to_column(as.data.frame(res_obj), feature_col)

    # ALDEx2 builtins (MA/MW/volcano) — wrap each in tryCatch in case one fails.
    graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    for (pt in c("MA", "MW", "volcano")) {
      tryCatch(
        ALDEx2::aldex.plot(res_obj, type = pt, test = "welch",
                           main = sprintf("%s — %s", cn, pt)),
        error = function(e) {
          graphics::plot.new()
          graphics::title(main = sprintf("%s — %s failed", cn, pt))
        }
      )
    }

    # Custom volcano: effect vs -log10(q).
    qcol <- intersect(c("wi.eBH", "we.eBH"), colnames(res_df))[1]
    if (!is.na(qcol)) {
      graphics::par(mfrow = c(1, 1), mar = c(5, 4, 3, 1))
      graphics::plot(
        res_df$effect, -log10(res_df[[qcol]]), pch = 20,
        xlab = "ALDEx2 effect (CLR diff)",
        ylab = sprintf("-log10(%s)", qcol),
        main = sprintf("%s — effect vs -log10(q)", cn)
      )
      graphics::abline(v = c(-effect_thresh, effect_thresh), col = "grey", lty = 2)
      graphics::abline(h = -log10(0.05), col = "red", lty = 2)
    }

    # Per-feature bar + violin plots intentionally omitted from this PDF —
    # the cross-comparison top_candidates.pdf below shows each pooled
    # candidate plotted against all groups, which is more useful than
    # repeating the same features per pair.
  }

  grDevices::dev.off(); pdf_open <- FALSE

  # ---- 5b. Per-pair volcano PNG (effect vs -log10(q)) --------------------
  # Quarto dashboards can't inline a multi-page PDF, so emit a flat PNG per
  # pair for the report. Same data as the PDF's custom volcano above.
  alpha <- cfg$stats$alpha %||% 0.05
  for (cn in names(raw_results)) {
    res_obj <- raw_results[[cn]]
    res_df  <- tibble::rownames_to_column(as.data.frame(res_obj), feature_col)
    qcol <- intersect(c("wi.eBH", "we.eBH"), colnames(res_df))[1]
    if (is.na(qcol)) next
    res_df$neglogq <- -log10(pmax(res_df[[qcol]], 1e-300))
    res_df$sig     <- res_df[[qcol]] < alpha
    p <- ggplot2::ggplot(res_df,
          ggplot2::aes(x = .data$effect, y = .data$neglogq,
                       color = .data$sig)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::geom_vline(xintercept = c(-effect_thresh, effect_thresh),
                          linetype = "dashed", colour = "grey60") +
      ggplot2::geom_hline(yintercept = -log10(alpha),
                          linetype = "dashed", colour = "red") +
      ggplot2::scale_color_manual(values = c(`TRUE` = "#d62728",
                                              `FALSE` = "#7f7f7f"),
                                   labels = c("ns", sprintf("q < %.2f", alpha))) +
      ggplot2::labs(
        title = sprintf("%s — effect vs -log10(q)", cn),
        x = "ALDEx2 effect (CLR diff)",
        y = sprintf("-log10(%s)", qcol),
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     text = ggplot2::element_text(size = 12))
    ggplot2::ggsave(file.path(fig_dir, paste0("volcano_", cn, ".png")),
                    p, width = 7, height = 5, dpi = 300)
  }

  # ---- 5c. Per-pair ALDEx2 MA-plot PNG (PIPELINE_V2_GAPS C5) ------------
  # Same MA data as the multi-page aldex_plots.pdf but written one PNG per
  # pair so the dashboard / manuscript can embed individual panels. Taxa
  # level only per spec (the gene-level path skips this to avoid bloat).
  if (level == "taxa") {
    for (cn in names(raw_results)) {
      res_obj <- raw_results[[cn]]
      grDevices::png(file.path(fig_dir, paste0("aldex2_maplot_", cn, ".png")),
                     width = 1800, height = 1400, res = 220, bg = "white")
      tryCatch({
        graphics::par(mar = c(4, 4, 3, 1))
        ALDEx2::aldex.plot(res_obj, type = "MA", test = "welch",
                           main = sprintf("ALDEx2 MA — %s", cn))
      }, error = function(e) {
        pipeline_log(cfg, sprintf("DA[%s] %s MA-plot failed: %s",
                                  level, cn, conditionMessage(e)))
      })
      grDevices::dev.off()
    }
  }

  # ---- 6. Cross-comparison summary ---------------------------------------
  # Pool candidates across pairs, dedupe, rank by max |effect|, then plot
  # each surviving feature against ALL groups (not just one pair). Mirrors
  # GT's aldex_top_taxa_bar_violin.pdf + aldex_top_candidates_summary.csv.
  all_cands <- unique(unlist(per_comp_candidates))
  if (length(all_cands) > 0) {
    # Rank by the max absolute pairwise effect across comparisons.
    effect_cols <- grep("_effect$", colnames(merged), value = TRUE)
    if (length(effect_cols) > 0) {
      eff_mat <- as.matrix(merged[, effect_cols, drop = FALSE])
      rownames(eff_mat) <- merged[[feature_col]]
      eff_mat[is.na(eff_mat)] <- 0
      max_abs <- apply(abs(eff_mat[all_cands, , drop = FALSE]), 1, max)
      all_cands <- names(sort(max_abs, decreasing = TRUE))
    }
    cap <- top_n_plots * max(1L, length(raw_results))
    if (length(all_cands) > cap) all_cands <- all_cands[seq_len(cap)]
    pipeline_log(cfg, sprintf(
      "DA[%s] cross-comp summary: %d unique candidates",
      level, length(all_cands)
    ))

    sum_pdf <- file.path(fig_dir, "top_candidates.pdf")
    grDevices::pdf(sum_pdf, width = 12, height = 8)
    pdf_open <- TRUE
    on.exit(if (pdf_open) grDevices::dev.off(), add = TRUE)
    for (ft in all_cands) {
      panel <- da_feature_panel(ft, mat, meta_use, sid, group,
                                feature_col, restrict_groups = NULL,
                                pal = pal,
                                title_suffix = "across all groups")
      if (!is.null(panel)) { print(panel$bar); print(panel$viol) }
    }
    grDevices::dev.off(); pdf_open <- FALSE

    # Top-30 sample-level heatmap (Frontiers Fig S3, taxa level only).
    # Picks the 30 features ranked by max |effect| across pairs (same
    # ranking as all_cands but without the top_n_plots cap), log10(count+1)
    # transformed and per-row z-scored so a single colour ramp covers all
    # rows. Columns are samples ordered by treatment, with gaps + a
    # treatment annotation strip on top.
    if (level == "taxa" &&
        requireNamespace("pheatmap", quietly = TRUE) &&
        length(effect_cols) > 0) {
      max_abs_all <- apply(abs(eff_mat), 1, max)
      top30 <- names(sort(max_abs_all, decreasing = TRUE))
      top30 <- intersect(top30, rownames(mat))
      top30 <- top30[seq_len(min(30L, length(top30)))]
      if (length(top30) >= 2) {
        samp_order <- order(conds)
        sub_mat    <- mat[top30, samp_order, drop = FALSE]
        grp_run    <- conds[samp_order]
        grp_lvls   <- unique(grp_run)
        grp_sizes  <- as.integer(table(grp_run)[grp_lvls])
        gaps_col   <- if (length(grp_sizes) > 1) {
          utils::head(cumsum(grp_sizes), -1L)
        } else NULL
        z <- t(scale(t(log10(sub_mat + 1))))
        z[is.na(z)] <- 0
        ann_col    <- data.frame(Treatment = grp_run,
                                  row.names = colnames(sub_mat))
        ann_colors <- list(Treatment = pal[grp_lvls])
        pheatmap::pheatmap(
          z,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          color = grDevices::colorRampPalette(
            c("#0612bd", "#bbbbbd", "#bd0606")
          )(100),
          border_color      = NA,
          annotation_col    = ann_col,
          annotation_colors = ann_colors,
          gaps_col          = gaps_col,
          fontsize_row = 8, fontsize_col = 9, fontsize = 10,
          filename = file.path(fig_dir, "top30_daa_heatmap.png"),
          width  = max(8, 0.4  * ncol(sub_mat) + 4),
          height = max(6, 0.22 * nrow(sub_mat) + 2.5)
        )
        pipeline_log(cfg, sprintf(
          "DA[%s] top30_daa_heatmap.png (%d features x %d samples)",
          level, nrow(sub_mat), ncol(sub_mat)
        ))
      }
    }

    # Summary CSV: union of top candidates with overall + per-pair stats.
    keep_cols <- c(feature_col,
                   intersect(c("kw.ep", "kw.eBH"), colnames(merged)),
                   grep("_(effect|overlap|wi\\.eBH|we\\.eBH)$",
                        colnames(merged), value = TRUE))
    keep_cols <- intersect(unique(keep_cols), colnames(merged))
    sum_df <- merged[merged[[feature_col]] %in% all_cands, keep_cols,
                     drop = FALSE]
    sum_df <- sum_df[match(all_cands, sum_df[[feature_col]]), , drop = FALSE]
    num <- vapply(sum_df, is.numeric, logical(1))
    sum_df[num] <- lapply(sum_df[num], function(x) round(x, 4))
    readr::write_csv(sum_df, file.path(ds_dir, "top_candidates_summary.csv"))

    # ---- 6b. Cross-comparison dot plot (PIPELINE_V2_GAPS C6) -----------
    # Forest-style summary: rows = features significant in >=1 pair (capped
    # at top-30 by max |effect|), columns = pairwise comparisons (facets),
    # x = ALDEx2 effect, point shape = significance, colour = direction,
    # BH-adjusted padj as star annotation. Single PNG so the dashboard can
    # embed it as fig5 panel D.
    if (level == "taxa" && length(effect_cols) > 0) {
      alpha_dt  <- cfg$stats$alpha %||% 0.05
      padj_cols <- grep("_(we|wi)\\.eBH$", colnames(merged), value = TRUE)
      if (length(padj_cols) > 0) {
        sig_mat <- as.matrix(merged[, padj_cols, drop = FALSE]) < alpha_dt
        sig_mat[is.na(sig_mat)] <- FALSE
        feature_sig <- rowSums(sig_mat) > 0
      } else {
        feature_sig <- rep(FALSE, nrow(merged))
      }
      if (any(feature_sig)) {
        max_abs <- apply(abs(eff_mat), 1, max)
        keep_idx <- which(feature_sig)
        keep_idx <- keep_idx[order(-max_abs[keep_idx])]
        keep_idx <- utils::head(keep_idx, 30L)
        sub_m    <- merged[keep_idx, , drop = FALSE]

        long_parts <- lapply(effect_cols, function(ec) {
          cn       <- sub("_effect$", "", ec)
          padj_col <- intersect(c(paste0(cn, "_we.eBH"),
                                  paste0(cn, "_wi.eBH")),
                                colnames(sub_m))[1]
          data.frame(
            feature = sub_m[[feature_col]],
            pair    = cn,
            effect  = sub_m[[ec]],
            padj    = if (!is.na(padj_col)) sub_m[[padj_col]] else NA_real_,
            stringsAsFactors = FALSE
          )
        })
        long_dt <- do.call(rbind, long_parts)
        long_dt$sig   <- !is.na(long_dt$padj) & long_dt$padj < alpha_dt
        long_dt$stars <- dplyr::case_when(
          is.na(long_dt$padj)      ~ "",
          long_dt$padj < 0.001     ~ "***",
          long_dt$padj < 0.01      ~ "**",
          long_dt$padj < alpha_dt  ~ "*",
          TRUE                      ~ ""
        )
        # Preserve max-|effect| ordering (most extreme at top of y-axis).
        long_dt$feature <- factor(long_dt$feature,
                                  levels = rev(sub_m[[feature_col]]))
        long_dt$direction <- ifelse(long_dt$effect > 0, "Up", "Down")

        p_dot <- ggplot2::ggplot(long_dt,
                  ggplot2::aes(x = .data$effect, y = .data$feature,
                                colour = .data$direction,
                                shape  = .data$sig)) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                              colour = "grey60") +
          ggplot2::geom_point(size = 3, alpha = 0.85) +
          ggplot2::geom_text(ggplot2::aes(label = .data$stars),
                             hjust = -0.6, vjust = 0.5,
                             size = 3, colour = "black",
                             show.legend = FALSE) +
          ggplot2::scale_colour_manual(values = c(Up   = "#d62728",
                                                   Down = "#1f77b4"),
                                        na.value = "grey60") +
          ggplot2::scale_shape_manual(values = c(`TRUE`  = 16,
                                                  `FALSE` = 1),
                                       labels = c("ns",
                                                  sprintf("q < %.2g",
                                                          alpha_dt))) +
          ggplot2::facet_wrap(~ pair, nrow = 1) +
          ggplot2::labs(x = "ALDEx2 effect (CLR diff)", y = NULL,
                        colour = "Direction", shape = "Significance") +
          ggplot2::theme_classic() +
          ggplot2::theme(
            legend.position = "top",
            strip.text      = ggplot2::element_text(face = "bold"),
            text            = ggplot2::element_text(size = 11)
          )
        nf <- nrow(sub_m)
        dh <- max(6, 0.28 * nf + 2.5)
        dw <- max(9, 3.5 * length(effect_cols) + 2)
        ggplot2::ggsave(file.path(fig_dir, "aldex2_dotplot_summary.png"),
                        p_dot, width = dw, height = dh, dpi = 300)
        pipeline_log(cfg, sprintf(
          "DA[%s] aldex2_dotplot_summary.png (%d features x %d pairs)",
          level, nf, length(effect_cols)
        ))
      } else {
        pipeline_log(cfg, sprintf(
          "DA[%s] aldex2_dotplot: no features significant in >=1 pair — skipping",
          level
        ))
      }
    }
  }

  invisible(NULL)
}

# Build a bar + violin panel for one feature. `restrict_groups` (or NULL
# for all groups) controls which samples are included. Uses the long-form
# `mat` already filtered to aligned samples.
da_feature_panel <- function(feat, mat, meta_use, sid, group, feature_col,
                             restrict_groups = NULL, pal,
                             title_suffix = NULL) {
  if (!feat %in% rownames(mat)) return(NULL)

  meta_df <- as.data.frame(meta_use[, c(sid, group), drop = FALSE])
  names(meta_df) <- c("sample", ".group")
  meta_df$sample <- as.character(meta_df$sample)

  use_samples <- meta_df$sample
  if (!is.null(restrict_groups)) {
    use_samples <- meta_df$sample[meta_df$.group %in% restrict_groups]
  }
  use_samples <- intersect(use_samples, colnames(mat))
  if (length(use_samples) < 2) return(NULL)

  cls <- data.frame(
    sample = use_samples,
    count  = as.numeric(mat[feat, use_samples])
  )
  cls <- merge(cls, meta_df, by = "sample", all.x = TRUE)
  cls$log_count <- log10(cls$count + 1)

  kw_p <- tryCatch(
    kruskal.test(count ~ .group, data = cls)$p.value,
    error = function(e) NA_real_
  )
  kw_lab <- sprintf("Kruskal-Wallis p = %s",
                    if (is.na(kw_p)) "NA" else format(kw_p, digits = 3))

  sample_stats <- cls |>
    dplyr::group_by(sample, .group) |>
    dplyr::summarise(
      mean_count = mean(count, na.rm = TRUE),
      se_count   = ifelse(dplyr::n() > 1,
                          sd(count, na.rm = TRUE) / sqrt(dplyr::n()),
                          NA_real_),
      .groups = "drop"
    )

  bar_title <- if (is.null(title_suffix)) {
    sprintf("%s — counts per sample (%s)", feat, kw_lab)
  } else {
    sprintf("%s — counts per sample (%s, %s)", feat, kw_lab, title_suffix)
  }

  barp <- ggplot2::ggplot(cls,
        ggplot2::aes(x = sample, y = count, fill = .group)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = mean(cls$count, na.rm = TRUE),
                        linetype = "dashed", colour = "red") +
    ggplot2::geom_errorbar(
      data = sample_stats,
      ggplot2::aes(x = sample,
                   ymin = pmax(0, mean_count - se_count),
                   ymax = mean_count + se_count),
      width = 0.2, colour = "black", inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = sample_stats,
      ggplot2::aes(x = sample, y = mean_count),
      colour = "red", size = 2, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::facet_wrap(~ .group, scales = "free_x", nrow = 1) +
    ggplot2::labs(title = bar_title,
                  x = "Sample", y = "Raw count", fill = group) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title  = ggplot2::element_text(hjust = 0.5)
    )

  present_groups <- unique(cls$.group)
  pair_combos <- if (length(present_groups) >= 2) {
    utils::combn(sort(present_groups), 2, simplify = FALSE)
  } else list()

  viol_title <- if (is.null(title_suffix)) feat else sprintf("%s — %s", feat, title_suffix)
  violp <- ggplot2::ggplot(cls,
        ggplot2::aes(x = .group, y = log_count, fill = .group)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width") +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          position = ggplot2::position_dodge(0.9)) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = group, y = "log10(count + 1)", title = viol_title) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5))
  if (length(pair_combos) > 0 && requireNamespace("ggpubr", quietly = TRUE)) {
    violp <- violp + ggpubr::stat_compare_means(
      comparisons = pair_combos, method = "wilcox.test", label = "p.signif"
    )
  }
  violp <- violp + ggplot2::annotate(
    "text",
    x = (length(present_groups) + 1) / 2,
    y = max(cls$log_count, na.rm = TRUE) +
        if (length(pair_combos) > 0) 0.6 * length(pair_combos) + 0.5 else 0.5,
    label = kw_lab, size = 4
  )

  list(bar = barp, viol = violp)
}

run_aldex <- function(cleaned, cfg) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    pipeline_log(cfg, "ALDEx2 not installed — skipping differential abundance")
    return(invisible(NULL))
  }
  pipeline_log(cfg, "Differential abundance (ALDEx2)")

  da_cfg         <- cfg$differential_abundance %||% list()
  levels         <- da_cfg$levels         %||% c("taxa", "gene")
  mc             <- da_cfg$mc_samples     %||% cfg$stats$aldex_mc_samples   %||% 128
  denom          <- da_cfg$denom          %||% "all"
  min_prev       <- da_cfg$min_prevalence %||% cfg$filters$min_prevalence_aldex %||% 2
  effect_thresh  <- da_cfg$effect_threshold  %||% 0.8
  overlap_thresh <- da_cfg$overlap_threshold %||% 0.20
  top_n_plots    <- da_cfg$top_n_plots       %||% 5
  make_plots     <- da_cfg$make_plots        %||% TRUE
  seed           <- da_cfg$seed              %||% 12345
  design_str     <- da_cfg$design                                # NULL = KW+pairwise mode
  alpha          <- da_cfg$alpha %||% cfg$stats$alpha %||% 0.05
  set.seed(seed)

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid   <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]
  if (is.null(group) || !group %in% colnames(meta)) {
    stop(sprintf("DA: group column '%s' missing from metadata", group))
  }
  meta[[sid]] <- as.character(meta[[sid]])

  if ("taxa" %in% levels) {
    picked <- da_taxa_matrix(cleaned, cfg)
    da_run_level("taxa", picked, cleaned, cfg, meta, sid, group,
                 mc, denom, min_prev, effect_thresh, overlap_thresh,
                 top_n_plots, make_plots, da_cfg,
                 design_str = design_str, alpha = alpha)
  }
  if ("gene" %in% levels) {
    picked <- da_gene_matrix(cleaned, cfg)
    if (is.null(picked)) {
      pipeline_log(cfg, "DA[gene]: cleaned$abri_kraken2 unavailable — skipping")
    } else {
      da_run_level("gene", picked, cleaned, cfg, meta, sid, group,
                   mc, denom, min_prev, effect_thresh, overlap_thresh,
                   top_n_plots, make_plots, da_cfg,
                   design_str = design_str, alpha = alpha)
    }
  }
  invisible(NULL)
}
