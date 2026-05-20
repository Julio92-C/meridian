# 09_resistome.R — ARG profile analysis (CARD hits). Mirrors the GT scripts
# ARGNorm_VennDiagram.R, ARGNorm_alphaDiversity.R, ARGNorm_pHeatmap.R,
# ARGNorm_relativeAbundance.R, and DrugNorm_pHeatmap.R (drug-class heatmap).
# Reads the TPM-normalised gene table written by R/03 and re-classifies the
# RESISTANCE strings into per-gene drug classes (single class / Multi-drug /
# optional MLS rollup). Driven entirely by cfg — no hardcoded sample IDs,
# treatments, palettes, or gene-name fixes (cfg$resistome overrides defaults).
#
# Outputs (under <project>/<figures_dir>/resistome/):
#   venn_drug_classes.png         Drug classes shared across groups
#   alpha_<metric>_violin.png     Violin + boxplot + KW annotation
#   alpha_<metric>_bar.png        Per-sample bar plot faceted by group
#   pcoa.png                      ARG-profile PCoA, ellipses, PERMANOVA/PERMDISP
#   pheatmap_genes.png            Gene × sample, min-max scaled, DRUG sidebar
#   pheatmap_drug_classes.png     Drug-class × sample, min-max scaled
#   drug_relative_abundance.png   Stacked drug-class % per sample
#   drug_total_count.png          Total TPM per drug class (log10, hbar)
#
# Datasets (under <project>/<datasets_dir>/resistome/):
#   alpha_diversity.csv, drug_total_TPM.csv, kw_<metric>.txt,
#   beta_permanova.txt, beta_permdisp.txt, beta_pcoa_scores.csv

`%||%` <- function(a, b) if (is.null(a)) b else a

run_resistome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Resistome (CARD)")
  rcfg     <- cfg$resistome %||% list()
  database <- tolower(rcfg$database %||% "card")

  norm <- .resistome_load_norm_table(cfg)
  if (is.null(norm)) return(invisible(NULL))

  card <- norm[tolower(norm$DATABASE) %in% database, , drop = FALSE]
  if (nrow(card) == 0) {
    pipeline_log(cfg, sprintf(
      "Resistome: no rows in DATABASE == '%s' — skipping",
      paste(database, collapse = "/")
    ))
    return(invisible(NULL))
  }
  if (!"RESISTANCE" %in% colnames(card)) {
    pipeline_log(cfg, "Resistome: RESISTANCE column missing — skipping")
    return(invisible(NULL))
  }

  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "resistome")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir,  "resistome")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  group        <- cfg$metadata$group_cols[[1]]
  card         <- .resistome_ensure_group_column(card, cfg, group)
  group_levels <- sort(unique(as.character(card[[group]])))
  pal_group    <- .resistome_group_palette(group_levels, rcfg)

  # ---- Drug-class classification ---------------------------------------
  card$RESISTANCE_raw <- card$RESISTANCE
  if (isTRUE(rcfg$mls_rollup %||% TRUE)) {
    card$RESISTANCE <- stringr::str_replace_all(
      card$RESISTANCE,
      stringr::fixed(
        "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B"
      ),
      "MLS"
    )
  }
  card$DRUG <- vapply(card$RESISTANCE, classify_resistance,
                      FUN.VALUE = character(1))
  card$DRUG <- ifelse(card$DRUG %in% c("Mls", "mls"), "MLS", card$DRUG)

  card_drug   <- dplyr::filter(card, !is.na(.data$DRUG))
  drug_levels <- sort(unique(card_drug$DRUG))
  pal_drug    <- .resistome_drug_palette(drug_levels, rcfg)

  pipeline_log(cfg, sprintf(
    "Resistome: %d CARD rows, %d genes, %d drug classes, %d samples",
    nrow(card), dplyr::n_distinct(card$GENE),
    length(drug_levels), dplyr::n_distinct(card$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- rcfg$alpha_metric %||% "richness"
  alpha <- .resistome_compute_alpha(card, group)
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  kw <- .resistome_alpha_kw(alpha, group, alpha_metric, cfg)
  if (!is.null(kw)) {
    capture.output(
      kw,
      file = file.path(ds_dir, sprintf("kw_%s.txt", alpha_metric))
    )
  }
  .resistome_plot_alpha_violin(alpha, group, alpha_metric, kw,
                               pal_group, fig_dir)
  .resistome_plot_alpha_bar   (alpha, group, alpha_metric, kw,
                               pal_group, fig_dir)

  # ---- (2) Venn of drug classes per group ------------------------------
  .resistome_plot_drug_venn(card_drug, group, pal_group, fig_dir, cfg)

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(rcfg$beta %||% TRUE)) {
    .resistome_plot_beta(card, group, pal_group, rcfg, fig_dir, ds_dir, cfg)
  }

  # ---- (3) Gene-level pheatmap with DRUG row annotation ----------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    .resistome_plot_gene_heatmap(card, group, pal_group, pal_drug,
                                 rcfg, fig_dir)
  } else {
    pipeline_log(cfg, "Resistome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Drug-class pheatmap -----------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(card_drug) > 0) {
    .resistome_plot_drug_heatmap(card_drug, group, pal_group, rcfg, fig_dir)
  }

  # ---- (5) Drug-class relative abundance + total bar -------------------
  drug_totals <- NULL
  if (nrow(card_drug) > 0) {
    drug_totals <- .resistome_plot_drug_relative_abundance(
      card_drug, group, pal_drug, fig_dir, ds_dir
    )
    .resistome_plot_drug_total_bar(drug_totals, pal_drug, fig_dir)
  } else {
    pipeline_log(cfg, "Resistome: no rows with a drug class — RA/total skipped")
  }

  invisible(list(alpha = alpha, drug_totals = drug_totals))
}

# ============================================================================
# Helpers
# ============================================================================

# Read the TPM-normalised gene table produced by R/03. Returns NULL with a
# log line if the file isn't on disk (e.g. cfg$stages$normalisation = false).
.resistome_load_norm_table <- function(cfg) {
  path <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                    "genetable_normdata.csv")
  if (!file.exists(path)) {
    pipeline_log(cfg, sprintf(
      "Resistome: %s not found — enable cfg$stages$normalisation first. Skipping.",
      path
    ))
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

# R/03 only attaches fixed_effects + random_effect to the normalised table.
# When the requested group column isn't already present, re-join metadata.
.resistome_ensure_group_column <- function(df, cfg, group) {
  if (group %in% colnames(df)) return(df)
  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE)
  if (!group %in% colnames(meta)) {
    stop(sprintf("Resistome: group column '%s' not found in metadata", group))
  }
  df$sample    <- as.character(df$sample)
  meta[[sid]]  <- as.character(meta[[sid]])
  dplyr::inner_join(df, meta[, c(sid, group), drop = FALSE],
                    by = c("sample" = sid))
}

# Turn a semicolon-separated RESISTANCE string into a single drug-class label.
# Single class -> that class, multiple -> "Multi-drug", empty/NA -> NA.
# Mirrors the GT helper of the same name in ARGNorm_*.R.
classify_resistance <- function(res_string) {
  if (is.null(res_string) || is.na(res_string) || res_string == "") {
    return(NA_character_)
  }
  pretty  <- stringr::str_to_title(stringr::str_replace_all(res_string,
                                                            "_", " "))
  classes <- trimws(unlist(strsplit(pretty, ";")))
  classes <- unique(classes[classes != ""])
  if (length(classes) == 0) return(NA_character_)
  if (length(classes) == 1) return(classes)
  "Multi-drug"
}

# Built-in literal-substring -> short-label rewrites for ARG gene names.
# Mirrors the case_when in the GT pheatmap / RA scripts. Override with
# cfg$resistome$gene_renames (also matched as literal substrings).
.default_resistome_gene_renames <- function() {
  list(
    "vanW_gene_in_vanB_cluster"                            = "vanWB",
    "vanR_gene_in_vanB_cluster"                            = "vanRB",
    "vanX_gene_in_vanB_cluster"                            = "vanXB",
    "vanY_gene_in_vanB_cluster"                            = "vanYB",
    "vanS_gene_in_vanB_cluster"                            = "vanSB",
    "vanH_gene_in_vanB_cluster"                            = "vanHB",
    "Escherichia_coli_ampC_beta-lactamase"                 = "ampC",
    "Escherichia_coli_mdfA"                                = "mdfA",
    "Escherichia_coli_emrE"                                = "emrE",
    "Escherichia_coli_acrA"                                = "acrA",
    "Vibrio_anguillarum_chloramphenicol_acetyltransferase" = "vacA*",
    "Bifidobacterium_adolescentis_rpoB_mutants_conferring_resistance_to_rifampicin" = "rpoB_mutants*",
    "Bifidobacterium_bifidum_ileS_conferring_resistance_to_mupirocin"               = "ileS*",
    "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein"          = "aac(6')-Ie/aph(2'')-Ia"
  )
}

# Apply literal-substring renames to a gene vector. A hit replaces the whole
# value (mirrors the GT case_when where any substring match swaps the GENE).
.apply_gene_renames <- function(genes, renames) {
  if (is.null(renames) || length(renames) == 0) return(as.character(genes))
  pats <- names(renames)
  out  <- as.character(genes)
  for (i in seq_along(pats)) {
    hits      <- stringr::str_detect(out, stringr::fixed(pats[[i]]))
    out[hits] <- renames[[i]]
  }
  out
}

# Discrete palette from `name` of length `n`, recycled if needed.
.resistome_discrete_palette <- function(name, n) {
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

.resistome_group_palette <- function(levels, rcfg) {
  user_map <- rcfg$group_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .resistome_discrete_palette(
    rcfg$group_palette %||% "ggsci::default_nejm",
    length(levels)
  )
  setNames(pal, levels)
}

.resistome_drug_palette <- function(levels, rcfg) {
  user_map <- rcfg$drug_colors
  if (!is.null(user_map)) {
    pal <- unlist(user_map[levels])
    if (length(pal) == length(levels) && all(!is.na(pal))) {
      names(pal) <- levels
      return(pal)
    }
  }
  pal <- .resistome_discrete_palette(
    rcfg$drug_palette %||% "ggsci::default_ucscgb",
    length(levels)
  )
  setNames(pal, levels)
}

# Per-sample richness + Shannon on a log(TPM+1) gene matrix. Matches GT.
.resistome_compute_alpha <- function(card, group) {
  abund <- card |>
    dplyr::mutate(TPM_log = log(.data$TPM + 1)) |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM_log = sum(.data$TPM_log), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM_log",
                       values_fill = 0)
  mat <- as.matrix(abund[, -1, drop = FALSE])
  rownames(mat) <- abund$sample
  group_per_sample <- card |>
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

.resistome_alpha_kw <- function(alpha, group, metric, cfg) {
  if (!metric %in% colnames(alpha)) {
    pipeline_log(cfg, sprintf("Resistome alpha: metric '%s' not found", metric))
    return(NULL)
  }
  sub <- alpha[!is.na(alpha[[metric]]), , drop = FALSE]
  if (dplyr::n_distinct(sub[[group]]) < 2) return(NULL)
  res <- tryCatch(
    kruskal.test(reformulate(group, metric), data = sub),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    pipeline_log(cfg, sprintf("Resistome alpha %s ~ %s KW p = %.4g",
                              metric, group, res$p.value))
  }
  res
}

.resistome_plot_alpha_violin <- function(alpha, group, metric, kw,
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

.resistome_plot_alpha_bar <- function(alpha, group, metric, kw,
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

.resistome_plot_drug_venn <- function(card_drug, group, pal_group,
                                       fig_dir, cfg) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    pipeline_log(cfg, "Resistome: VennDiagram not available — Venn skipped")
    return(invisible(NULL))
  }
  drug_by_group <- card_drug |>
    dplyr::group_by(.data[[group]], .data$DRUG) |>
    dplyr::summarise(total = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(.data$total > 0)
  sets  <- split(drug_by_group$DRUG, drug_by_group[[group]])
  sets  <- lapply(sets, unique)
  n_grp <- length(sets)
  if (n_grp < 2) {
    pipeline_log(cfg, "Resistome: <2 groups — Venn skipped")
    return(invisible(NULL))
  }
  if (n_grp > 5) {
    pipeline_log(cfg, sprintf(
      "Resistome: %d groups — Venn supports 2-5, skipped", n_grp
    ))
    return(invisible(NULL))
  }
  fill_pal <- unname(pal_group[names(sets)])
  # disable.logging silences the file log but futile.logger still chatters the
  # call arguments to the console; raise its threshold for the duration.
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    suppressMessages(futile.logger::flog.threshold(
      futile.logger::FATAL, name = "VennDiagramLogger"
    ))
  }
  # Push category labels off the circles so long strings ("Reference diet",
  # "Soyabean meal") don't collide with the diagram or get clipped at the
  # canvas edge; also widen the canvas + bump the margin a bit.
  cat_args <- switch(as.character(n_grp),
    "2" = list(cat.pos = c(-20, 20),      cat.dist = c(0.05, 0.05)),
    "3" = list(cat.pos = c(-25, 25, 180), cat.dist = c(0.08, 0.08, 0.04)),
    "4" = list(cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22, 0.22, 0.12, 0.12)),
    "5" = list(),
    list()
  )
  do.call(VennDiagram::venn.diagram, c(
    list(
      x         = sets,
      filename  = file.path(fig_dir, "venn_drug_classes.png"),
      imagetype = "png",
      height = 2800, width = 2800, resolution = 300,
      fill = fill_pal, alpha = 0.6,
      cat.cex = 1.25, cat.fontface = "bold",
      cex = 1.5,
      main = sprintf("Drug classes shared across %s groups", group),
      main.cex = 1.3, margin = 0.18,
      disable.logging = TRUE
    ),
    cat_args
  ))
}

# Min-max scale each column independently. Mirrors the GT pheatmap recipe so
# each sample is shown on its own 0..1 colour scale.
.resistome_minmax_per_col <- function(mat) {
  apply(mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
}

.resistome_plot_gene_heatmap <- function(card, group, pal_group, pal_drug,
                                          rcfg, fig_dir) {
  renames <- rcfg$gene_renames %||% .default_resistome_gene_renames()
  card_r <- card |>
    dplyr::mutate(GENE = .apply_gene_renames(.data$GENE, renames))

  gene_wide <- card_r |>
    dplyr::group_by(.data$GENE, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("GENE") |>
    as.matrix()

  top_n <- rcfg$top_n_genes
  if (!is.null(top_n) && nrow(gene_wide) > top_n) {
    keep <- names(sort(rowSums(gene_wide), decreasing = TRUE))[seq_len(top_n)]
    gene_wide <- gene_wide[keep, , drop = FALSE]
  }
  gene_scaled <- .resistome_minmax_per_col(gene_wide)

  drug_per_gene <- card_r |>
    dplyr::distinct(.data$GENE, .keep_all = TRUE) |>
    dplyr::select("GENE", "DRUG") |>
    tibble::column_to_rownames("GENE")
  drug_per_gene <- drug_per_gene[rownames(gene_scaled), , drop = FALSE]

  sample_group <- card_r |>
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
  drug_levels <- unique(stats::na.omit(drug_per_gene$DRUG))
  ann_colors$DRUG <- pal_drug[intersect(names(pal_drug), drug_levels)]

  hm_height <- max(7, 0.18 * nrow(gene_scaled) + 3)
  hm_width  <- max(8, 0.4  * ncol(gene_scaled) + 4)
  pheatmap::pheatmap(
    gene_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      rcfg$gene_heatmap_palette %||% c("#0612bd", "#bbbbbd", "#bd0606")
    )(100),
    border_color      = NA,
    annotation_row    = drug_per_gene,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 8, fontsize_col = 9, fontsize = 10,
    filename = file.path(fig_dir, "pheatmap_genes.png"),
    width = hm_width, height = hm_height
  )
}

.resistome_plot_drug_heatmap <- function(card_drug, group, pal_group,
                                          rcfg, fig_dir) {
  drug_wide <- card_drug |>
    dplyr::group_by(.data$DRUG, .data$sample) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "TPM",
                       values_fill = 0) |>
    tibble::column_to_rownames("DRUG") |>
    as.matrix()
  drug_wide   <- drug_wide[order(rowSums(drug_wide), decreasing = TRUE), ,
                            drop = FALSE]
  drug_scaled <- .resistome_minmax_per_col(drug_wide)

  sample_group <- card_drug |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::select("sample", dplyr::all_of(group)) |>
    tibble::column_to_rownames("sample")
  sample_group <- sample_group[order(sample_group[[group]]), , drop = FALSE]
  drug_scaled  <- drug_scaled[, rownames(sample_group), drop = FALSE]
  grp_run      <- as.character(sample_group[[group]])
  grp_lvls     <- unique(grp_run)
  grp_sizes    <- as.integer(table(grp_run)[grp_lvls])
  gaps_col     <- if (length(grp_sizes) > 1) {
    utils::head(cumsum(grp_sizes), -1L)
  } else NULL

  ann_colors <- list()
  ann_colors[[group]] <- pal_group[grp_lvls]

  hm_height <- max(4, 0.32 * nrow(drug_scaled) + 3)
  hm_width  <- max(8, 0.4  * ncol(drug_scaled) + 4)
  pheatmap::pheatmap(
    drug_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      rcfg$drug_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )(100),
    border_color      = NA,
    annotation_col    = sample_group,
    annotation_colors = ann_colors,
    gaps_col          = gaps_col,
    fontsize_row = 11, fontsize_col = 10, fontsize = 11,
    filename = file.path(fig_dir, "pheatmap_drug_classes.png"),
    width = hm_width, height = hm_height
  )
}

.resistome_plot_drug_relative_abundance <- function(card_drug, group, pal_drug,
                                                     fig_dir, ds_dir) {
  drug_share <- card_drug |>
    dplyr::group_by(.data$sample, .data$DRUG, .data[[group]]) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(Percentage = 100 * .data$TPM / sum(.data$TPM)) |>
    dplyr::ungroup()
  p <- ggplot2::ggplot(drug_share,
        ggplot2::aes(x = sample, y = Percentage, fill = DRUG)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = pal_drug) +
    ggplot2::facet_wrap(stats::reformulate(group),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)",
                  fill = "Drug class") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 13),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  bw <- max(8, 0.55 * dplyr::n_distinct(drug_share$sample) + 3)
  ggplot2::ggsave(file.path(fig_dir, "drug_relative_abundance.png"), p,
                  width = bw, height = 6, dpi = 300)

  drug_totals <- card_drug |>
    dplyr::group_by(.data$DRUG) |>
    dplyr::summarise(Total_TPM = sum(.data$TPM, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$Total_TPM))
  readr::write_csv(drug_totals, file.path(ds_dir, "drug_total_TPM.csv"))
  drug_totals
}

# Gene-level PCoA + PERMANOVA + PERMDISP on the ARG TPM profile. Mirrors
# R/07's recipe (Hellinger transform -> Bray-Curtis -> cmdscale + adonis2 +
# betadisper) but on resistome features instead of taxa, so the question is
# "do treatment groups have distinct ARG community structure?".
.resistome_plot_beta <- function(card, group, pal_group, rcfg,
                                  fig_dir, ds_dir, cfg) {
  wide <- card |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "GENE", values_from = "TPM",
                       values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$sample
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) < 3) {
    pipeline_log(cfg, sprintf(
      "Resistome beta: only %d samples with ARG hits — need >= 3, skipping",
      nrow(mat)
    ))
    return(invisible(NULL))
  }

  transform <- rcfg$beta_transform %||% "hellinger"
  distance  <- rcfg$beta_distance  %||% "bray"
  mat_t <- switch(transform,
    none      = mat,
    log       = log1p(mat),
    hellinger = vegan::decostand(mat, method = "hellinger"),
    stop(sprintf("Resistome beta: unknown transform '%s'", transform))
  )
  d <- vegan::vegdist(mat_t, method = distance)

  meta <- card |>
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
  pipeline_log(cfg, sprintf("Resistome PERMANOVA %s: R2 = %.3f, p = %.4g",
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
    pipeline_log(cfg, sprintf("Resistome PERMDISP %s: p = %.4g",
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
  ellipse_type <- rcfg$ellipse_type     %||% "norm"
  ellipse_line <- rcfg$ellipse_linetype %||% "dashed"

  p <- ggplot2::ggplot(scores,
        ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(type = ellipse_type, linewidth = 0.8,
                          linetype = ellipse_line) +
    ggplot2::scale_color_manual(values = pal_group) +
    ggplot2::labs(
      title = sprintf("ARG profile PCoA (%s, %s-transformed) — by %s",
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

.resistome_plot_drug_total_bar <- function(drug_totals, pal_drug, fig_dir) {
  fmt_label <- function(x) {
    ifelse(x >= 1000,
           paste0(round(x / 1000, 1), "K"),
           as.character(round(x, 1)))
  }
  p <- ggplot2::ggplot(drug_totals,
        ggplot2::aes(x = stats::reorder(DRUG, Total_TPM),
                     y = Total_TPM, fill = DRUG)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = stats::median(drug_totals$Total_TPM),
                        linetype = "dashed", colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = fmt_label(Total_TPM)),
                       hjust = -0.05, size = 3.4) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = pal_drug) +
    ggplot2::labs(x = "Drug class", y = "Total TPM (log10)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_text(size = 11),
                   text = ggplot2::element_text(size = 13)) +
    ggplot2::coord_flip()
  ggplot2::ggsave(file.path(fig_dir, "drug_total_count.png"), p,
                  width = 7, height = 6, dpi = 300)
}
