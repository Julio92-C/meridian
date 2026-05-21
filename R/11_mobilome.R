# 11_mobilome.R — MGE / plasmid profile analysis (PlasmidFinder hits).
# Mirrors the GT scripts MGEsNorm_relativeAbundance.R + MGEsNorm_pHeatmap.R.
# Reads the TPM-normalised gene table written by R/03 and classifies each
# replicon GENE into a Replicon_Family via configurable regex patterns
# (GT defaults: ^Col -> "Col-like", ^IncF -> "IncF", ^IncX -> "IncX",
# ^Inc -> "Other Inc", else -> "Unknown/Other"). Driven entirely by cfg.
#
# Outputs (under <project>/<figures_dir>/mobilome/):
#   venn_replicons.png             Plasmid replicons (GENEs) shared across groups
#   alpha_<metric>_violin.png      Violin + boxplot + KW annotation
#   alpha_<metric>_bar.png         Per-sample bar plot faceted by group
#   mge_abundance_violin.png       Per-(sample,gene) log(TPM) by group + KW
#   pcoa.png                       MGE-profile PCoA, PERMANOVA / PERMDISP
#   pheatmap_genes.png             Replicon × sample, Family row sidebar
#   pheatmap_replicon_families.png Family × sample (GT yellow-cyan-magenta)
#   mge_relative_abundance.png     Stacked replicon % per sample (GT-style, by GENE)
#   mge_family_abundance.png       Stacked family % per sample
#   mge_total_count.png            Total TPM per replicon (log10, hbar, GT-style)
#   mge_family_total_count.png     Total TPM per family (log10, hbar)
#
# Datasets (under <project>/<datasets_dir>/mobilome/):
#   alpha_diversity.csv, mge_total_TPM.csv, replicon_family_total_TPM.csv,
#   kw_<metric>.txt, beta_permanova.txt, beta_permdisp.txt,
#   beta_pcoa_scores.csv

`%||%` <- function(a, b) if (is.null(a)) b else a

run_mobilome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Mobilome (PlasmidFinder)")
  mcfg     <- cfg$mobilome %||% list()
  database <- tolower(mcfg$database %||% "plasmidfinder")

  norm <- .mobilome_load_norm_table(cfg)
  if (is.null(norm)) return(invisible(NULL))

  pfdb <- norm[tolower(norm$DATABASE) %in% database, , drop = FALSE]
  if (nrow(pfdb) == 0) {
    pipeline_log(cfg, sprintf(
      "Mobilome: no rows in DATABASE == '%s' — skipping",
      paste(database, collapse = "/")
    ))
    return(invisible(NULL))
  }

  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "mobilome")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir,  "mobilome")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  group        <- cfg$metadata$group_cols[[1]]
  pfdb         <- .mobilome_ensure_group_column(pfdb, cfg, group)
  group_levels <- sort(unique(as.character(pfdb[[group]])))
  pal_group    <- .mobilome_group_palette(group_levels, mcfg)

  # ---- Replicon family classification ----------------------------------
  patterns <- mcfg$family_patterns %||% .default_replicon_patterns()
  pfdb$Replicon_Family <- classify_replicon_family(pfdb$GENE, patterns)

  pfdb_fam     <- dplyr::filter(pfdb, !is.na(.data$Replicon_Family))
  family_levels <- sort(unique(pfdb_fam$Replicon_Family))
  pal_family    <- .mobilome_family_palette(family_levels, mcfg)
  gene_levels   <- sort(unique(pfdb$GENE))
  pal_gene      <- .mobilome_gene_palette(gene_levels, mcfg)

  pipeline_log(cfg, sprintf(
    "Mobilome: %d PlasmidFinder rows, %d replicons, %d families, %d samples",
    nrow(pfdb), length(gene_levels),
    length(family_levels), dplyr::n_distinct(pfdb$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- mcfg$alpha_metric %||% "richness"
  alpha <- .mobilome_compute_alpha(pfdb, group)
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  kw <- .mobilome_alpha_kw(alpha, group, alpha_metric, cfg)
  if (!is.null(kw)) {
    capture.output(
      kw,
      file = file.path(ds_dir, sprintf("kw_%s.txt", alpha_metric))
    )
  }
  .mobilome_plot_alpha_violin(alpha, group, alpha_metric, kw,
                              pal_group, fig_dir)
  .mobilome_plot_alpha_bar   (alpha, group, alpha_metric, kw,
                              pal_group, fig_dir)

  # ---- (1b) MGE abundance per Treatment (log TPM violin) ---------------
  .mobilome_plot_abundance_violin(pfdb, group, pal_group, fig_dir, cfg)

  # ---- (2) Venn of plasmid replicons per group -------------------------
  .mobilome_plot_replicon_venn(pfdb, group, pal_group, fig_dir, cfg)

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(mcfg$beta %||% TRUE)) {
    .mobilome_plot_beta(pfdb, group, pal_group, mcfg, fig_dir, ds_dir, cfg)
  }

  # ---- (3) Gene-level pheatmap with Family row annotation --------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    .mobilome_plot_gene_heatmap(pfdb, group, pal_group, pal_family,
                                mcfg, fig_dir)
  } else {
    pipeline_log(cfg, "Mobilome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Family-level pheatmap ---------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(pfdb_fam) > 0) {
    .mobilome_plot_family_heatmap(pfdb_fam, group, pal_group, mcfg, fig_dir)
  }

  # ---- (5) Gene & family relative abundance + total bars ---------------
  gene_totals <- .mobilome_plot_gene_relative_abundance(
    pfdb, group, pal_gene, fig_dir, ds_dir
  )
  .mobilome_plot_gene_total_bar(gene_totals, pal_gene, fig_dir)
  save_count_prevalence(
    pfdb,
    file.path(fig_dir, "mge_count_prevalence.png"),
    category_col   = "GENE",
    value_col      = "TPM",
    palette        = pal_gene,
    value_label    = "Total TPM",
    category_label = "Plasmid replicon",
    width = 13, height = 7
  )

  family_totals <- NULL
  if (nrow(pfdb_fam) > 0) {
    family_totals <- .mobilome_plot_family_relative_abundance(
      pfdb_fam, group, pal_family, fig_dir, ds_dir
    )
    .mobilome_plot_family_total_bar(family_totals, pal_family, fig_dir)
    save_count_prevalence(
      pfdb_fam,
      file.path(fig_dir, "mge_family_count_prevalence.png"),
      category_col   = "Replicon_Family",
      value_col      = "TPM",
      palette        = pal_family,
      value_label    = "Total TPM",
      category_label = "Replicon family",
      width = 11, height = 4
    )
  } else {
    pipeline_log(cfg, "Mobilome: no family-classified rows — family RA skipped")
  }

  invisible(list(alpha = alpha,
                 gene_totals = gene_totals,
                 family_totals = family_totals))
}

# ============================================================================
# Helpers
# ============================================================================

.mobilome_load_norm_table <- function(cfg) {
  path <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                    "genetable_normdata.csv")
  if (!file.exists(path)) {
    pipeline_log(cfg, sprintf(
      "Mobilome: %s not found — enable cfg$stages$normalisation first. Skipping.",
      path
    ))
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

.mobilome_ensure_group_column <- function(df, cfg, group) {
  if (group %in% colnames(df)) return(df)
  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE)
  if (!group %in% colnames(meta)) {
    stop(sprintf("Mobilome: group column '%s' not found in metadata", group))
  }
  df$sample    <- as.character(df$sample)
  meta[[sid]]  <- as.character(meta[[sid]])
  dplyr::inner_join(df, meta[, c(sid, group), drop = FALSE],
                    by = c("sample" = sid))
}

# Ordered list of (family -> regex) rules, matched first-hit-wins. Mirrors
# the GT classify_plasmids case_when chain. Override via cfg$mobilome$family_
# patterns (a named list of family -> regex; unmatched genes fall back to
# "Unknown/Other").
.default_replicon_patterns <- function() {
  list(
    "Col-like"  = "^Col",
    "IncF"      = "^IncF",
    "IncX"      = "^IncX",
    "Other Inc" = "^Inc"
  )
}

# Classify each GENE name into a replicon family using ordered regex rules.
# First matching rule wins; anything unmatched lands in "Unknown/Other".
classify_replicon_family <- function(genes, patterns) {
  out <- rep(NA_character_, length(genes))
  remaining <- !is.na(genes)
  for (fam in names(patterns)) {
    hits <- remaining & stringr::str_detect(genes, patterns[[fam]])
    out[hits] <- fam
    remaining[hits] <- FALSE
  }
  out[is.na(out) & !is.na(genes)] <- "Unknown/Other"
  out
}

.mobilome_discrete_palette <- function(name, n) {
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

.mobilome_group_palette <- function(levels, mcfg) {
  user_map <- mcfg$group_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .mobilome_discrete_palette(
    mcfg$group_palette %||% "ggsci::default_nejm",
    length(levels)
  )
  setNames(pal, levels)
}

.mobilome_family_palette <- function(levels, mcfg) {
  user_map <- mcfg$family_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .mobilome_discrete_palette(
    mcfg$family_palette %||% "ggsci::default_nejm",
    length(levels)
  )
  setNames(pal, levels)
}

.mobilome_gene_palette <- function(levels, mcfg) {
  user_map <- mcfg$gene_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .mobilome_discrete_palette(
    mcfg$gene_palette %||% "ggsci::default_ucscgb",
    length(levels)
  )
  setNames(pal, levels)
}

.mobilome_compute_alpha <- function(pfdb, group) {
  abund <- pfdb |>
    dplyr::mutate(TPM_log = log(.data$TPM + 1)) |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM_log = sum(.data$TPM_log), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM_log",
                       values_fill = 0)
  mat <- as.matrix(abund[, -1, drop = FALSE])
  rownames(mat) <- abund$sample
  group_per_sample <- pfdb |>
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

.mobilome_alpha_kw <- function(alpha, group, metric, cfg) {
  if (!metric %in% colnames(alpha)) {
    pipeline_log(cfg, sprintf("Mobilome alpha: metric '%s' not found", metric))
    return(NULL)
  }
  sub <- alpha[!is.na(alpha[[metric]]), , drop = FALSE]
  if (dplyr::n_distinct(sub[[group]]) < 2) return(NULL)
  res <- tryCatch(
    kruskal.test(reformulate(group, metric), data = sub),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    pipeline_log(cfg, sprintf("Mobilome alpha %s ~ %s KW p = %.4g",
                              metric, group, res$p.value))
  }
  res
}

.mobilome_plot_alpha_violin <- function(alpha, group, metric, kw,
                                         pal_group, fig_dir) {
  if (!metric %in% colnames(alpha)) return(invisible(NULL))
  y_max <- max(alpha[[metric]], na.rm = TRUE)
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
      "text",
      x = (dplyr::n_distinct(alpha[[group]]) + 1) / 2,
      y = y_max * 1.08,
      label = sprintf("Kruskal-Wallis p = %.2g", kw$p.value),
      size = 4.5
    )
  }
  ggplot2::ggsave(file.path(fig_dir, sprintf("alpha_%s_violin.png", metric)),
                  p, width = 7, height = 5, dpi = 300)
}

.mobilome_plot_alpha_bar <- function(alpha, group, metric, kw,
                                      pal_group, fig_dir) {
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
  ggplot2::ggsave(file.path(fig_dir, sprintf("alpha_%s_bar.png", metric)),
                  p, width = bw, height = 5, dpi = 300)
}

.mobilome_plot_abundance_violin <- function(pfdb, group, pal_group,
                                             fig_dir, cfg) {
  if (!"TPM" %in% colnames(pfdb) || nrow(pfdb) == 0) return(invisible(NULL))
  df <- data.frame(
    sample  = pfdb$sample,
    group   = pfdb[[group]],
    log_TPM = log(pfdb$TPM + 1)
  )
  df <- df[is.finite(df$log_TPM), , drop = FALSE]
  if (dplyr::n_distinct(df$group) < 2 || nrow(df) < 3) {
    pipeline_log(cfg,
                 "Mobilome abundance: insufficient data for KW — skipping")
    return(invisible(NULL))
  }
  names(df)[names(df) == "group"] <- group
  kw <- tryCatch(
    kruskal.test(reformulate(group, "log_TPM"), data = df),
    error = function(e) NULL
  )
  if (!is.null(kw)) {
    pipeline_log(cfg,
                 sprintf("Mobilome log(TPM+1) ~ %s KW p = %.4g",
                         group, kw$p.value))
  }
  y_max <- max(df$log_TPM, na.rm = TRUE)
  p <- ggplot2::ggplot(df,
        ggplot2::aes(x = .data[[group]], y = .data$log_TPM,
                     fill = .data[[group]])) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          position = ggplot2::position_dodge(0.9)) +
    ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.45) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = "log(TPM)",
                  title = "MGE abundance per Treatment group") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5,
                                                       face = "bold",
                                                       size = 13),
                   text = ggplot2::element_text(size = 13))
  if (!is.null(kw)) {
    p <- p + ggplot2::annotate(
      "text",
      x = (dplyr::n_distinct(df[[group]]) + 1) / 2,
      y = y_max * 1.05,
      label = sprintf("Kruskal-Wallis p = %.2g", kw$p.value),
      size = 4.5
    )
  }
  ggplot2::ggsave(file.path(fig_dir, "mge_abundance_violin.png"),
                  p, width = 7, height = 5.5, dpi = 300)
}

.mobilome_plot_replicon_venn <- function(pfdb, group, pal_group,
                                          fig_dir, cfg) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    pipeline_log(cfg, "Mobilome: VennDiagram not available — Venn skipped")
    return(invisible(NULL))
  }
  rep_by_group <- pfdb |>
    dplyr::group_by(.data[[group]], .data$GENE) |>
    dplyr::summarise(total = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(.data$total > 0)
  sets  <- split(rep_by_group$GENE, rep_by_group[[group]])
  sets  <- lapply(sets, unique)
  n_grp <- length(sets)
  if (n_grp < 2) {
    pipeline_log(cfg, "Mobilome: <2 groups — Venn skipped")
    return(invisible(NULL))
  }
  if (n_grp > 5) {
    pipeline_log(cfg, sprintf(
      "Mobilome: %d groups — Venn supports 2-5, skipped", n_grp
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
      filename  = file.path(fig_dir, "venn_replicons.png"),
      imagetype = "png",
      height = 2800, width = 2800, resolution = 300,
      fill = fill_pal, alpha = 0.6,
      cat.cex = 1.25, cat.fontface = "bold",
      cex = 1.5,
      main = sprintf("Plasmid replicons shared across %s groups", group),
      main.cex = 1.3, margin = 0.18,
      disable.logging = TRUE
    ),
    cat_args
  ))
}

.mobilome_minmax_per_col <- function(mat) {
  apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
}

.mobilome_plot_gene_heatmap <- function(pfdb, group, pal_group, pal_family,
                                         mcfg, fig_dir) {
  gene_wide <- pfdb |>
    dplyr::group_by(.data$GENE, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("GENE") |>
    as.matrix()
  # GT sorts rows by total abundance descending before min-max.
  gene_wide   <- gene_wide[order(rowSums(gene_wide), decreasing = TRUE), ,
                            drop = FALSE]
  top_n <- mcfg$top_n_genes
  if (!is.null(top_n) && nrow(gene_wide) > top_n) {
    gene_wide <- gene_wide[seq_len(top_n), , drop = FALSE]
  }
  gene_scaled <- .mobilome_minmax_per_col(gene_wide)

  fam_per_gene <- pfdb |>
    dplyr::distinct(.data$GENE, .keep_all = TRUE) |>
    dplyr::select("GENE", "Replicon_Family") |>
    tibble::column_to_rownames("GENE")
  fam_per_gene <- fam_per_gene[rownames(gene_scaled), , drop = FALSE]

  sample_group <- pfdb |>
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
  fam_levels <- unique(stats::na.omit(fam_per_gene$Replicon_Family))
  ann_colors$Replicon_Family <- pal_family[intersect(names(pal_family),
                                                       fam_levels)]

  hm_height <- max(5, 0.32 * nrow(gene_scaled) + 3)
  hm_width  <- max(8, 0.4  * ncol(gene_scaled) + 4)
  pheatmap::pheatmap(
    gene_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      mcfg$gene_heatmap_palette %||% c("#ebc3c0", "#bd0606", "#d98c07")
    )(100),
    border_color      = NA,
    annotation_row    = fam_per_gene,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 10, fontsize_col = 10, fontsize = 11,
    filename = file.path(fig_dir, "pheatmap_genes.png"),
    width = hm_width, height = hm_height
  )
}

.mobilome_plot_family_heatmap <- function(pfdb_fam, group, pal_group,
                                           mcfg, fig_dir) {
  fam_wide <- pfdb_fam |>
    dplyr::group_by(.data$Replicon_Family, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("Replicon_Family") |>
    as.matrix()
  fam_wide   <- fam_wide[order(rowSums(fam_wide), decreasing = TRUE), ,
                          drop = FALSE]
  fam_scaled <- .mobilome_minmax_per_col(fam_wide)

  sample_group <- pfdb_fam |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group)) |>
    tibble::column_to_rownames("sample")
  sample_group <- sample_group[order(sample_group[[group]]), , drop = FALSE]
  fam_scaled   <- fam_scaled[, rownames(sample_group), drop = FALSE]
  grp_run      <- as.character(sample_group[[group]])
  grp_lvls     <- unique(grp_run)
  grp_sizes    <- as.integer(table(grp_run)[grp_lvls])
  gaps_col     <- if (length(grp_sizes) > 1) {
    utils::head(cumsum(grp_sizes), -1L)
  } else NULL

  ann_colors <- list()
  ann_colors[[group]] <- pal_group[grp_lvls]

  hm_height <- max(4, 0.4 * nrow(fam_scaled) + 3)
  hm_width  <- max(8, 0.4 * ncol(fam_scaled) + 4)
  pheatmap::pheatmap(
    fam_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      mcfg$family_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )(100),
    border_color      = NA,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 11, fontsize_col = 10, fontsize = 11,
    filename = file.path(fig_dir, "pheatmap_replicon_families.png"),
    width = hm_width, height = hm_height
  )
}

.mobilome_plot_gene_relative_abundance <- function(pfdb, group, pal_gene,
                                                    fig_dir, ds_dir) {
  share <- pfdb |>
    dplyr::group_by(.data$sample, .data$GENE, .data[[group]]) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(Percentage = 100 * .data$TPM / sum(.data$TPM)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(share,
        ggplot2::aes(x = sample, y = Percentage, fill = GENE)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = pal_gene) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                  fill = "Plasmid replicon") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  bw <- max(8, 0.55 * dplyr::n_distinct(share$sample) + 3)
  ggplot2::ggsave(file.path(fig_dir, "mge_relative_abundance.png"), p,
                  width = bw, height = 6, dpi = 300)

  gene_totals <- pfdb |>
    dplyr::group_by(.data$GENE) |>
    dplyr::summarise(Total_TPM = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Total_TPM))
  readr::write_csv(gene_totals, file.path(ds_dir, "mge_total_TPM.csv"))
  gene_totals
}

.mobilome_plot_family_relative_abundance <- function(pfdb_fam, group,
                                                      pal_family,
                                                      fig_dir, ds_dir) {
  share <- pfdb_fam |>
    dplyr::group_by(.data$sample, .data$Replicon_Family, .data[[group]]) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(Percentage = 100 * .data$TPM / sum(.data$TPM)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(share,
        ggplot2::aes(x = sample, y = Percentage, fill = Replicon_Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = pal_family) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                  fill = "Replicon family") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  bw <- max(8, 0.55 * dplyr::n_distinct(share$sample) + 3)
  ggplot2::ggsave(file.path(fig_dir, "mge_family_abundance.png"), p,
                  width = bw, height = 6, dpi = 300)

  family_totals <- pfdb_fam |>
    dplyr::group_by(.data$Replicon_Family) |>
    dplyr::summarise(Total_TPM = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Total_TPM))
  readr::write_csv(family_totals,
                   file.path(ds_dir, "replicon_family_total_TPM.csv"))
  family_totals
}

.mobilome_plot_gene_total_bar <- function(gene_totals, pal_gene, fig_dir) {
  fmt_label <- function(x) {
    ifelse(x >= 1000,
           paste0(round(x / 1000, 1), "K"),
           as.character(round(x, 1)))
  }
  p <- ggplot2::ggplot(gene_totals,
        ggplot2::aes(x = stats::reorder(GENE, Total_TPM),
                     y = Total_TPM, fill = GENE)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = stats::median(gene_totals$Total_TPM),
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_label(Total_TPM)),
                       hjust = -0.05, size = 3.4) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = pal_gene) +
    ggplot2::labs(x = "Plasmid replicon", y = "Total TPM (log10)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 11),
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  ggplot2::ggsave(file.path(fig_dir, "mge_total_count.png"), p,
                  width = 7, height = 6, dpi = 300)
}

.mobilome_plot_family_total_bar <- function(family_totals, pal_family, fig_dir) {
  fmt_label <- function(x) {
    ifelse(x >= 1000,
           paste0(round(x / 1000, 1), "K"),
           as.character(round(x, 1)))
  }
  p <- ggplot2::ggplot(family_totals,
        ggplot2::aes(x = stats::reorder(Replicon_Family, Total_TPM),
                     y = Total_TPM, fill = Replicon_Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = stats::median(family_totals$Total_TPM),
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_label(Total_TPM)),
                       hjust = -0.05, size = 3.4) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = pal_family) +
    ggplot2::labs(x = "Replicon family", y = "Total TPM (log10)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 11),
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  ggplot2::ggsave(file.path(fig_dir, "mge_family_total_count.png"), p,
                  width = 7, height = 5, dpi = 300)
}

.mobilome_plot_beta <- function(pfdb, group, pal_group, mcfg,
                                 fig_dir, ds_dir, cfg) {
  wide <- pfdb |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM",
                       values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) < 3) {
    pipeline_log(cfg, sprintf(
      "Mobilome beta: only %d samples with MGE hits — need >= 3, skipping",
      nrow(mat)
    ))
    return(invisible(NULL))
  }

  transform <- mcfg$beta_transform %||% "hellinger"
  distance  <- mcfg$beta_distance  %||% "bray"
  mat_t <- switch(transform,
    none      = mat,
    log       = log1p(mat),
    hellinger = vegan::decostand(mat, method = "hellinger"),
    stop(sprintf("Mobilome beta: unknown transform '%s'", transform))
  )
  d <- vegan::vegdist(mat_t, method = distance)

  meta <- pfdb |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group))
  meta <- meta[match(rownames(mat_t), meta$sample), , drop = FALSE]

  perms <- cfg$stats$permanova_permutations %||% 9999
  permanova <- vegan::adonis2(
    stats::reformulate(group, "d"),
    data = meta, permutations = perms
  )
  utils::capture.output(permanova,
                        file = file.path(ds_dir, "beta_permanova.txt"))
  r2 <- permanova$R2[1]
  pv <- permanova$`Pr(>F)`[1]
  pipeline_log(cfg, sprintf("Mobilome PERMANOVA %s: R2 = %.3f, p = %.4g",
                            group, r2, pv))

  bd_test <- tryCatch({
    bd <- vegan::betadisper(d, factor(meta[[group]]))
    vegan::permutest(bd, permutations = perms)
  }, error = function(e) NULL)
  permdisp_p <- NA_real_
  if (!is.null(bd_test)) {
    utils::capture.output(bd_test,
                          file = file.path(ds_dir, "beta_permdisp.txt"))
    permdisp_p <- bd_test$tab$`Pr(>F)`[1]
    pipeline_log(cfg, sprintf("Mobilome PERMDISP %s: p = %.4g",
                              group, permdisp_p))
  }

  pcoa <- stats::cmdscale(d, eig = TRUE, k = 2)
  var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
  scores <- data.frame(
    sample = rownames(mat_t),
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2]
  )
  scores <- dplyr::left_join(scores, meta, by = "sample")
  readr::write_csv(scores, file.path(ds_dir, "beta_pcoa_scores.csv"))

  annot <- sprintf("PERMANOVA R² = %.3f, p = %.4g\nPERMDISP p = %.4g",
                   r2, pv, permdisp_p)
  ellipse_type <- mcfg$ellipse_type     %||% "norm"
  ellipse_line <- mcfg$ellipse_linetype %||% "dashed"

  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal_group) +
    ggplot2::labs(
      title = sprintf("MGE profile PCoA (%s, %s-transformed) — by %s",
                      distance, transform, group),
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
  ggplot2::ggsave(file.path(fig_dir, "pcoa.png"), p,
                  width = 7, height = 5.5, dpi = 300)

  invisible(list(scores = scores, permanova = permanova,
                 permdisp = bd_test, var_explained = var_expl[1:2]))
}
