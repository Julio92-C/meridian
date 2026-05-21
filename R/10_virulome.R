# 10_virulome.R — VF profile analysis (VFDB hits). Mirrors the GT scripts
# VFsNorm_relativeAbundance.R and VFsNorm_pHeatmap.R. Reads the TPM-normalised
# gene table written by R/03 and extracts the VFDB functional category from
# the PRODUCT column with a regex (the GT extract_productFunction helper).
# Driven entirely by cfg — no hardcoded sample IDs, treatments, palettes, or
# function-name fixes (cfg$virulome overrides defaults).
#
# Outputs (under <project>/<figures_dir>/virulome/):
#   venn_vf_functions.png         VF functions shared across groups
#   alpha_<metric>_violin.png     Violin + boxplot + KW annotation
#   alpha_<metric>_bar.png        Per-sample bar plot faceted by group
#   vf_abundance_violin.png       Per-(sample,gene) log(TPM) by group + KW
#   pcoa.png                      VF-profile PCoA, ellipses, PERMANOVA/PERMDISP
#   pheatmap_genes.png            Gene × sample, min-max scaled, Function sidebar
#   pheatmap_functions.png        Function × sample, min-max scaled
#   vf_relative_abundance.png     Stacked function % per sample
#   vf_total_count.png            Total TPM per function (log10, hbar)
#
# Datasets (under <project>/<datasets_dir>/virulome/):
#   alpha_diversity.csv, vf_total_TPM.csv, kw_<metric>.txt,
#   beta_permanova.txt, beta_permdisp.txt, beta_pcoa_scores.csv

`%||%` <- function(a, b) if (is.null(a)) b else a

run_virulome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Virulome (VFDB)")
  vcfg     <- cfg$virulome %||% list()
  database <- tolower(vcfg$database %||% "vfdb")

  norm <- .virulome_load_norm_table(cfg)
  if (is.null(norm)) return(invisible(NULL))

  vfdb <- norm[tolower(norm$DATABASE) %in% database, , drop = FALSE]
  if (nrow(vfdb) == 0) {
    pipeline_log(cfg, sprintf(
      "Virulome: no rows in DATABASE == '%s' — skipping",
      paste(database, collapse = "/")
    ))
    return(invisible(NULL))
  }
  if (!"PRODUCT" %in% colnames(vfdb)) {
    pipeline_log(cfg, "Virulome: PRODUCT column missing — skipping")
    return(invisible(NULL))
  }

  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "virulome")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir,  "virulome")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  group        <- cfg$metadata$group_cols[[1]]
  vfdb         <- .virulome_ensure_group_column(vfdb, cfg, group)
  group_levels <- sort(unique(as.character(vfdb[[group]])))
  pal_group    <- .virulome_group_palette(group_levels, vcfg)

  # ---- Function extraction from PRODUCT --------------------------------
  vfdb$Functions <- vapply(vfdb$PRODUCT, extract_vf_function,
                           FUN.VALUE = character(1))
  fn_renames <- vcfg$function_renames %||% .default_virulome_function_renames()
  vfdb$Functions <- .apply_literal_renames(vfdb$Functions, fn_renames)

  vfdb_fun     <- dplyr::filter(vfdb, !is.na(.data$Functions),
                                 .data$Functions != "")
  function_levels <- sort(unique(vfdb_fun$Functions))
  pal_function    <- .virulome_function_palette(function_levels, vcfg)

  pipeline_log(cfg, sprintf(
    "Virulome: %d VFDB rows, %d genes, %d VF functions, %d samples",
    nrow(vfdb), dplyr::n_distinct(vfdb$GENE),
    length(function_levels), dplyr::n_distinct(vfdb$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- vcfg$alpha_metric %||% "richness"
  alpha <- .virulome_compute_alpha(vfdb, group)
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  kw <- .virulome_alpha_kw(alpha, group, alpha_metric, cfg)
  if (!is.null(kw)) {
    capture.output(
      kw,
      file = file.path(ds_dir, sprintf("kw_%s.txt", alpha_metric))
    )
  }
  .virulome_plot_alpha_violin(alpha, group, alpha_metric, kw,
                              pal_group, fig_dir)
  .virulome_plot_alpha_bar   (alpha, group, alpha_metric, kw,
                              pal_group, fig_dir)

  # ---- (1b) VF abundance per Treatment (log TPM violin) ----------------
  .virulome_plot_abundance_violin(vfdb, group, pal_group, fig_dir, cfg)

  # ---- (2) Venn of VF functions per group ------------------------------
  .virulome_plot_function_venn(vfdb_fun, group, pal_group, fig_dir, cfg)

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(vcfg$beta %||% TRUE)) {
    .virulome_plot_beta(vfdb, group, pal_group, vcfg, fig_dir, ds_dir, cfg)
  }

  # ---- (3) Gene-level pheatmap with Function row annotation ------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    .virulome_plot_gene_heatmap(vfdb, group, pal_group, pal_function,
                                vcfg, fig_dir)
  } else {
    pipeline_log(cfg, "Virulome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Function-level pheatmap -------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(vfdb_fun) > 0) {
    .virulome_plot_function_heatmap(vfdb_fun, group, pal_group, vcfg, fig_dir)
  }

  # ---- (5) Function relative abundance + total bar ---------------------
  vf_totals <- NULL
  if (nrow(vfdb_fun) > 0) {
    vf_totals <- .virulome_plot_function_relative_abundance(
      vfdb_fun, group, pal_function, fig_dir, ds_dir
    )
    .virulome_plot_function_total_bar(vf_totals, pal_function, fig_dir)
  } else {
    pipeline_log(cfg, "Virulome: no rows with a function — RA/total skipped")
  }

  invisible(list(alpha = alpha, vf_totals = vf_totals))
}

# ============================================================================
# Helpers
# ============================================================================

.virulome_load_norm_table <- function(cfg) {
  path <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                    "genetable_normdata.csv")
  if (!file.exists(path)) {
    pipeline_log(cfg, sprintf(
      "Virulome: %s not found — enable cfg$stages$normalisation first. Skipping.",
      path
    ))
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

.virulome_ensure_group_column <- function(df, cfg, group) {
  if (group %in% colnames(df)) return(df)
  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE)
  if (!group %in% colnames(meta)) {
    stop(sprintf("Virulome: group column '%s' not found in metadata", group))
  }
  df$sample    <- as.character(df$sample)
  meta[[sid]]  <- as.character(meta[[sid]])
  dplyr::inner_join(df, meta[, c(sid, group), drop = FALSE],
                    by = c("sample" = sid))
}

# Extract the VFDB functional category from a PRODUCT string. VFDB PRODUCT
# values look like:
#   "(chuA) Outer membrane heme/hemoglobin receptor ChuA
#    [Chu (VF0227) - Nutritional/Metabolic factor (VFC0272)] [...]"
# The function category is the substring between " - " and " (" inside the
# bracketed VF block. Mirrors the GT extract_productFunction helper.
extract_vf_function <- function(s) {
  if (is.null(s) || is.na(s) || s == "") return(NA_character_)
  m <- regmatches(s, regexpr("- ([^\\(]+) \\(", s))
  if (length(m) == 0 || !nzchar(m)) return(NA_character_)
  trimws(gsub("- | \\(", "", m))
}

# GT abbreviates one long function label so the pheatmap row sidebar stays
# legible. Override via cfg$virulome$function_renames (named list of literal
# substring -> replacement).
.default_virulome_function_renames <- function() {
  list("Antimicrobial activity/Competitive advantage" =
         "Antimicrobial activity/CA*")
}

# Apply literal-substring renames to a character vector. A hit replaces the
# whole value (matches the GT case_when pattern used for ARG renames).
.apply_literal_renames <- function(x, renames) {
  if (is.null(renames) || length(renames) == 0) return(as.character(x))
  pats <- names(renames)
  out  <- as.character(x)
  for (i in seq_along(pats)) {
    hits      <- !is.na(out) & stringr::str_detect(out, stringr::fixed(pats[[i]]))
    out[hits] <- renames[[i]]
  }
  out
}

.virulome_discrete_palette <- function(name, n) {
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

.virulome_group_palette <- function(levels, vcfg) {
  user_map <- vcfg$group_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .virulome_discrete_palette(
    vcfg$group_palette %||% "ggsci::default_nejm",
    length(levels)
  )
  setNames(pal, levels)
}

.virulome_function_palette <- function(levels, vcfg) {
  user_map <- vcfg$function_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .virulome_discrete_palette(
    vcfg$function_palette %||% "ggsci::default_nejm",
    length(levels)
  )
  setNames(pal, levels)
}

.virulome_compute_alpha <- function(vfdb, group) {
  abund <- vfdb |>
    dplyr::mutate(TPM_log = log(.data$TPM + 1)) |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM_log = sum(.data$TPM_log), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM_log",
                       values_fill = 0)
  mat <- as.matrix(abund[, -1, drop = FALSE])
  rownames(mat) <- abund$sample
  group_per_sample <- vfdb |>
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

.virulome_alpha_kw <- function(alpha, group, metric, cfg) {
  if (!metric %in% colnames(alpha)) {
    pipeline_log(cfg, sprintf("Virulome alpha: metric '%s' not found", metric))
    return(NULL)
  }
  sub <- alpha[!is.na(alpha[[metric]]), , drop = FALSE]
  if (dplyr::n_distinct(sub[[group]]) < 2) return(NULL)
  res <- tryCatch(
    kruskal.test(reformulate(group, metric), data = sub),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    pipeline_log(cfg, sprintf("Virulome alpha %s ~ %s KW p = %.4g",
                              metric, group, res$p.value))
  }
  res
}

.virulome_plot_alpha_violin <- function(alpha, group, metric, kw,
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

.virulome_plot_alpha_bar <- function(alpha, group, metric, kw,
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

# Per-(sample, gene) log(TPM+1) violin + boxplot + jitter colored by group,
# with KW annotation. Complements the richness-based alpha plot by showing
# per-VF abundance distribution differences. Y-axis labelled log(TPM).
.virulome_plot_abundance_violin <- function(vfdb, group, pal_group,
                                             fig_dir, cfg) {
  if (!"TPM" %in% colnames(vfdb) || nrow(vfdb) == 0) return(invisible(NULL))
  df <- data.frame(
    sample  = vfdb$sample,
    group   = vfdb[[group]],
    log_TPM = log(vfdb$TPM + 1)
  )
  df <- df[is.finite(df$log_TPM), , drop = FALSE]
  if (dplyr::n_distinct(df$group) < 2 || nrow(df) < 3) {
    pipeline_log(cfg,
                 "Virulome abundance: insufficient data for KW — skipping")
    return(invisible(NULL))
  }
  names(df)[names(df) == "group"] <- group
  kw <- tryCatch(
    kruskal.test(reformulate(group, "log_TPM"), data = df),
    error = function(e) NULL
  )
  if (!is.null(kw)) {
    pipeline_log(cfg,
                 sprintf("Virulome log(TPM+1) ~ %s KW p = %.4g",
                         group, kw$p.value))
  }
  y_max <- max(df$log_TPM, na.rm = TRUE)
  p <- ggplot2::ggplot(df,
        ggplot2::aes(x = .data[[group]], y = .data$log_TPM,
                     fill = .data[[group]])) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          position = ggplot2::position_dodge(0.9)) +
    ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.35) +
    ggplot2::scale_fill_manual(values = pal_group) +
    ggplot2::labs(x = group, y = "log(TPM)",
                  title = "VF abundance per Treatment group") +
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
  ggplot2::ggsave(file.path(fig_dir, "vf_abundance_violin.png"),
                  p, width = 7, height = 5.5, dpi = 300)
}

.virulome_plot_function_venn <- function(vfdb_fun, group, pal_group,
                                          fig_dir, cfg) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    pipeline_log(cfg, "Virulome: VennDiagram not available — Venn skipped")
    return(invisible(NULL))
  }
  fun_by_group <- vfdb_fun |>
    dplyr::group_by(.data[[group]], .data$Functions) |>
    dplyr::summarise(total = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(.data$total > 0)
  sets  <- split(fun_by_group$Functions, fun_by_group[[group]])
  sets  <- lapply(sets, unique)
  n_grp <- length(sets)
  if (n_grp < 2) {
    pipeline_log(cfg, "Virulome: <2 groups — Venn skipped")
    return(invisible(NULL))
  }
  if (n_grp > 5) {
    pipeline_log(cfg, sprintf(
      "Virulome: %d groups — Venn supports 2-5, skipped", n_grp
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
      filename  = file.path(fig_dir, "venn_vf_functions.png"),
      imagetype = "png",
      height = 2800, width = 2800, resolution = 300,
      fill = fill_pal, alpha = 0.6,
      cat.cex = 1.25, cat.fontface = "bold",
      cex = 1.5,
      main = sprintf("VF functions shared across %s groups", group),
      main.cex = 1.3, margin = 0.18,
      disable.logging = TRUE
    ),
    cat_args
  ))
}

.virulome_minmax_per_col <- function(mat) {
  apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
}

.virulome_plot_gene_heatmap <- function(vfdb, group, pal_group, pal_function,
                                         vcfg, fig_dir) {
  gene_renames <- vcfg$gene_renames %||% list()
  vfdb_r <- vfdb |>
    dplyr::mutate(GENE = .apply_literal_renames(.data$GENE, gene_renames))

  gene_wide <- vfdb_r |>
    dplyr::group_by(.data$GENE, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("GENE") |>
    as.matrix()

  top_n <- vcfg$top_n_genes
  if (!is.null(top_n) && nrow(gene_wide) > top_n) {
    keep <- names(sort(rowSums(gene_wide), decreasing = TRUE))[seq_len(top_n)]
    gene_wide <- gene_wide[keep, , drop = FALSE]
  }
  gene_scaled <- .virulome_minmax_per_col(gene_wide)

  fun_per_gene <- vfdb_r |>
    dplyr::distinct(.data$GENE, .keep_all = TRUE) |>
    dplyr::select("GENE", "Functions") |>
    tibble::column_to_rownames("GENE")
  fun_per_gene <- fun_per_gene[rownames(gene_scaled), , drop = FALSE]

  sample_group <- vfdb_r |>
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
  fun_levels <- unique(stats::na.omit(fun_per_gene$Functions))
  ann_colors$Functions <- pal_function[intersect(names(pal_function), fun_levels)]

  hm_height <- max(7, 0.14 * nrow(gene_scaled) + 3)
  hm_width  <- max(8, 0.4  * ncol(gene_scaled) + 4)
  pheatmap::pheatmap(
    gene_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      vcfg$gene_heatmap_palette %||% c("#0612bd", "#bbbbbd", "#bd0606")
    )(100),
    border_color      = NA,
    annotation_row    = fun_per_gene,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 6, fontsize_col = 9, fontsize = 10,
    filename = file.path(fig_dir, "pheatmap_genes.png"),
    width = hm_width, height = hm_height
  )
}

.virulome_plot_function_heatmap <- function(vfdb_fun, group, pal_group,
                                             vcfg, fig_dir) {
  fun_wide <- vfdb_fun |>
    dplyr::group_by(.data$Functions, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("Functions") |>
    as.matrix()
  fun_wide   <- fun_wide[order(rowSums(fun_wide), decreasing = TRUE), ,
                          drop = FALSE]
  fun_scaled <- .virulome_minmax_per_col(fun_wide)

  sample_group <- vfdb_fun |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group)) |>
    tibble::column_to_rownames("sample")
  sample_group <- sample_group[order(sample_group[[group]]), , drop = FALSE]
  fun_scaled   <- fun_scaled[, rownames(sample_group), drop = FALSE]
  grp_run      <- as.character(sample_group[[group]])
  grp_lvls     <- unique(grp_run)
  grp_sizes    <- as.integer(table(grp_run)[grp_lvls])
  gaps_col     <- if (length(grp_sizes) > 1) {
    utils::head(cumsum(grp_sizes), -1L)
  } else NULL

  ann_colors <- list()
  ann_colors[[group]] <- pal_group[grp_lvls]

  hm_height <- max(4, 0.32 * nrow(fun_scaled) + 3)
  hm_width  <- max(8, 0.4  * ncol(fun_scaled) + 4)
  pheatmap::pheatmap(
    fun_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      vcfg$function_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )(100),
    border_color      = NA,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 11, fontsize_col = 10, fontsize = 11,
    filename = file.path(fig_dir, "pheatmap_functions.png"),
    width = hm_width, height = hm_height
  )
}

.virulome_plot_function_relative_abundance <- function(vfdb_fun, group,
                                                        pal_function,
                                                        fig_dir, ds_dir) {
  fun_share <- vfdb_fun |>
    dplyr::group_by(.data$sample, .data$Functions, .data[[group]]) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(Percentage = 100 * .data$TPM / sum(.data$TPM)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(fun_share,
        ggplot2::aes(x = sample, y = Percentage, fill = Functions)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = pal_function) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                  fill = "Virulent function") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  bw <- max(8, 0.55 * dplyr::n_distinct(fun_share$sample) + 3)
  ggplot2::ggsave(file.path(fig_dir, "vf_relative_abundance.png"), p,
                  width = bw, height = 6, dpi = 300)

  vf_totals <- vfdb_fun |>
    dplyr::group_by(.data$Functions) |>
    dplyr::summarise(Total_TPM = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Total_TPM))
  readr::write_csv(vf_totals, file.path(ds_dir, "vf_total_TPM.csv"))
  vf_totals
}

.virulome_plot_beta <- function(vfdb, group, pal_group, vcfg,
                                 fig_dir, ds_dir, cfg) {
  wide <- vfdb |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM",
                       values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) < 3) {
    pipeline_log(cfg, sprintf(
      "Virulome beta: only %d samples with VF hits — need >= 3, skipping",
      nrow(mat)
    ))
    return(invisible(NULL))
  }

  transform <- vcfg$beta_transform %||% "hellinger"
  distance  <- vcfg$beta_distance  %||% "bray"
  mat_t <- switch(transform,
    none      = mat,
    log       = log1p(mat),
    hellinger = vegan::decostand(mat, method = "hellinger"),
    stop(sprintf("Virulome beta: unknown transform '%s'", transform))
  )
  d <- vegan::vegdist(mat_t, method = distance)

  meta <- vfdb |>
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
  pipeline_log(cfg, sprintf("Virulome PERMANOVA %s: R2 = %.3f, p = %.4g",
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
    pipeline_log(cfg, sprintf("Virulome PERMDISP %s: p = %.4g",
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
  ellipse_type <- vcfg$ellipse_type     %||% "norm"
  ellipse_line <- vcfg$ellipse_linetype %||% "dashed"

  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal_group) +
    ggplot2::labs(
      title = sprintf("VF profile PCoA (%s, %s-transformed) — by %s",
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

.virulome_plot_function_total_bar <- function(vf_totals, pal_function, fig_dir) {
  fmt_label <- function(x) {
    ifelse(x >= 1000,
           paste0(round(x / 1000, 1), "K"),
           as.character(round(x, 1)))
  }
  p <- ggplot2::ggplot(vf_totals,
        ggplot2::aes(x = stats::reorder(Functions, Total_TPM),
                     y = Total_TPM, fill = Functions)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = stats::median(vf_totals$Total_TPM),
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_label(Total_TPM)),
                       hjust = -0.05, size = 3.4) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = pal_function) +
    ggplot2::labs(x = "Virulent function", y = "Total TPM (log10)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 11),
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  ggplot2::ggsave(file.path(fig_dir, "vf_total_count.png"), p,
                  width = 7, height = 6, dpi = 300)
}
