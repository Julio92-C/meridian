# 10_virulome.R — VF profile analysis (VFDB hits). Mirrors the GT scripts
# VFsNorm_relativeAbundance.R and VFsNorm_pHeatmap.R. Reads the TPM-normalised
# gene table written by R/03 and extracts the VFDB functional category from
# the PRODUCT column with a regex (the GT extract_productFunction helper).
# Driven entirely by cfg — no hardcoded sample IDs, treatments, palettes, or
# function-name fixes (cfg$virulome overrides defaults). Structural plumbing
# lives in R/utils_ge_profile.R; this module supplies the VFDB-specific
# function extractor, default function renames, and orchestration.
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
#   vf_count_prevalence.png       Paired total / sample-count panels (function)
#   gene_count_prevalence.png     Paired total / sample-count panels (gene)
#
# Datasets (under <project>/<datasets_dir>/virulome/):
#   alpha_diversity.csv, vf_total_TPM.csv, kw_<metric>.txt,
#   beta_permanova.txt, beta_permdisp.txt, beta_pcoa_scores.csv

`%||%` <- function(a, b) if (is.null(a)) b else a

run_virulome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Virulome (VFDB)")
  vcfg     <- cfg$virulome %||% list()
  database <- tolower(vcfg$database %||% "vfdb")
  label    <- "Virulome"

  norm <- ge_load_norm_table(cfg, label)
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
  vfdb         <- ge_ensure_group_column(vfdb, cfg, group, label)
  group_levels <- sort(unique(as.character(vfdb[[group]])))
  pal_group    <- resolve_top_level_colors(group, group_levels, cfg) %||%
                  ge_named_palette(group_levels, vcfg$group_colors,
                                    vcfg$group_palette %||% "ggsci::default_nejm")

  # ---- Function extraction from PRODUCT --------------------------------
  vfdb$Functions <- vapply(vfdb$PRODUCT, extract_vf_function,
                           FUN.VALUE = character(1))
  fn_renames <- vcfg$function_renames %||% .default_virulome_function_renames()
  vfdb$Functions <- ge_apply_literal_renames(vfdb$Functions, fn_renames)

  vfdb_fun        <- dplyr::filter(vfdb, !is.na(.data$Functions),
                                    .data$Functions != "")
  function_levels <- sort(unique(vfdb_fun$Functions))
  pal_function    <- ge_named_palette(function_levels, vcfg$function_colors,
                                       vcfg$function_palette %||% "ggsci::default_nejm")

  pipeline_log(cfg, sprintf(
    "Virulome: %d VFDB rows, %d genes, %d VF functions, %d samples",
    nrow(vfdb), dplyr::n_distinct(vfdb$GENE),
    length(function_levels), dplyr::n_distinct(vfdb$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- vcfg$alpha_metric %||% "richness"
  alpha <- ge_compute_alpha(vfdb, group)
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  # Pre-pass: compute both per-metric KW p-values before rendering so we
  # can family-adjust across {richness, shannon} within the domain, then
  # annotate raw + adjusted p in the plot titles. Cross-module accumulator
  # still records the raw values.
  kw_metrics <- unique(c(alpha_metric, "shannon"))
  kw_objects <- setNames(vector("list", length(kw_metrics)), kw_metrics)
  for (m in kw_metrics) {
    kw_objects[[m]] <- ge_alpha_kw(alpha, group, m, cfg, label)
  }
  raw_ps <- vapply(kw_metrics, function(m) {
    k <- kw_objects[[m]]
    if (is.null(k)) NA_real_ else k$p.value
  }, numeric(1))
  pad_method <- padjust_method(cfg)
  adj_ps     <- padjust_p(raw_ps, cfg)

  for (m in kw_metrics) {
    kw <- kw_objects[[m]]
    if (!is.null(kw)) {
      capture.output(kw, file = file.path(ds_dir, sprintf("kw_%s.txt", m)))
    }
    ge_plot_alpha_violin(
      alpha, group, m, kw, pal_group,
      file.path(fig_dir, sprintf("alpha_%s_violin.png", m)),
      p_adj = adj_ps[[m]], padj_method = pad_method
    )
    ge_plot_alpha_bar(
      alpha, group, m, kw, pal_group,
      file.path(fig_dir, sprintf("alpha_%s_bar.png", m)),
      p_adj = adj_ps[[m]], padj_method = pad_method
    )
  }

  # ---- (1b) VF abundance per Treatment (log TPM violin) ----------------
  ge_plot_abundance_violin(
    vfdb, group, pal_group,
    file = file.path(fig_dir, "vf_abundance_violin.png"),
    title = "VF abundance per Treatment group",
    log_label = label, cfg = cfg
  )

  # ---- (2) Venn of VF functions per group ------------------------------
  ge_plot_category_venn(
    vfdb_fun, category_col = "Functions", group = group, pal_group = pal_group,
    file = file.path(fig_dir, "venn_vf_functions.png"),
    main_title = sprintf("VF functions shared across %s groups", group),
    log_label = label, cfg = cfg
  )
  # Gene-level VF Venn — added 2026-06-12 for the diet-effects
  # supplementary which wants treatment-overlap of individual VF genes
  # (matching the ARG and MGE rows which already operate at GENE level).
  # Falls through R/13's auto-classifier as kind `ge_venn_vfs`.
  ge_plot_category_venn(
    vfdb, category_col = "GENE", group = group, pal_group = pal_group,
    file = file.path(fig_dir, "venn_vfs.png"),
    main_title = sprintf("VF genes shared across %s groups", group),
    log_label = label, cfg = cfg
  )

  # VF UpSet (gene-level) — UpSet counterpart to the function-level Venn.
  ge_plot_category_upset(
    vfdb, category_col = "GENE", group = group,
    file = file.path(fig_dir, "vf_upset_treatments.png"),
    log_label = label, cfg = cfg, pal_group = pal_group
  )

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(vcfg$beta %||% TRUE)) {
    ge_plot_beta(
      vfdb, group, pal_group, opts = vcfg,
      fig_file     = file.path(fig_dir, "pcoa.png"),
      ds_perm_file = file.path(ds_dir,  "beta_permanova.txt"),
      ds_disp_file = file.path(ds_dir,  "beta_permdisp.txt"),
      ds_pcoa_file = file.path(ds_dir,  "beta_pcoa_scores.csv"),
      log_label = label, title_prefix = "VF profile",
      skip_data_label = "VF hits", cfg = cfg
    )
  }

  # ---- (3) Gene-level pheatmap with Function row annotation ------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    gene_renames <- vcfg$gene_renames %||% list()
    # Restrict to genes with a recognised Functions assignment so the
    # row annotation sidebar always renders (genes whose PRODUCT didn't
    # match any of the function patterns get filtered out — they'd
    # produce an unlabelled row and an empty palette entry otherwise).
    vfdb_r <- vfdb |>
      dplyr::filter(!is.na(.data$Functions), nzchar(.data$Functions)) |>
      dplyr::mutate(GENE = ge_apply_literal_renames(.data$GENE, gene_renames))
    ge_plot_gene_heatmap(
      vfdb_r, group, pal_group, pal_function,
      category_col = "Functions",
      file    = file.path(fig_dir, "pheatmap_genes.png"),
      palette = vcfg$gene_heatmap_palette %||% c("#0612bd", "#bbbbbd", "#bd0606"),
      top_n   = vcfg$top_n_genes %||% 30L,
      rank_by = vcfg$gene_rank_by %||% "prevalence",
      fontsize_row = 8, row_height_factor = 0.18, min_height = 7
    )
  } else {
    pipeline_log(cfg, "Virulome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Function-level pheatmap -------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(vfdb_fun) > 0) {
    ge_plot_category_heatmap(
      vfdb_fun, category_col = "Functions", group = group, pal_group = pal_group,
      file    = file.path(fig_dir, "pheatmap_functions.png"),
      palette = vcfg$function_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )
  }

  # ---- (5) Function relative abundance + total bar ---------------------
  vf_totals <- NULL
  if (nrow(vfdb_fun) > 0) {
    # Per-sample VF function TPM sums — feeds the v1.2 Mantel triangle
    # (R/12 reads this when cfg$network$mantel$level == "category") and
    # mirrors R/09's drug_class_per_sample_TPM.csv pattern.
    vfdb_fun |>
      dplyr::group_by(.data$sample, .data$Functions) |>
      dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE),
                       .groups = "drop") |>
      readr::write_csv(file.path(ds_dir, "vf_function_per_sample_TPM.csv"))

    vf_totals <- ge_plot_category_relative_abundance(
      vfdb_fun, category_col = "Functions", group = group,
      palette = pal_function,
      fig_file = file.path(fig_dir, "vf_relative_abundance.png"),
      ds_file  = file.path(ds_dir,  "vf_total_TPM.csv"),
      fill_label = "VF function"
    )
    # Streamgraph companion (2026-06-12 polish) — per-treatment view for
    # composite panel D.
    save_stream_composition(
      vfdb_fun,
      file.path(fig_dir, "vf_stream_abundance.png"),
      category_col = "Functions",
      group_col    = group,
      value_col    = "TPM",
      palette      = pal_function,
      fill_label   = "VF function",
      title        = "VF function composition by treatment"
    )
    ge_plot_total_bar(
      vf_totals, category_col = "Functions", value_col = "Total_TPM",
      palette = pal_function,
      file    = file.path(fig_dir, "vf_total_count.png"),
      x_label = "VF function"
    )
    save_count_prevalence(
      vfdb_fun,
      file.path(fig_dir, "vf_count_prevalence.png"),
      category_col   = "Functions",
      value_col      = "TPM",
      palette        = pal_function,
      value_label    = "Total TPM",
      category_label = "VF function",
      width = 12, height = 6
    )
  } else {
    pipeline_log(cfg, "Virulome: no rows with a function — RA/total skipped")
  }

  # ---- (6) Gene-level paired count + prevalence ------------------------
  if ("TPM" %in% colnames(vfdb) && nrow(vfdb) > 0) {
    gh <- max(6, 0.18 * dplyr::n_distinct(vfdb$GENE) + 2)
    save_count_prevalence(
      vfdb,
      file.path(fig_dir, "gene_count_prevalence.png"),
      category_col   = "GENE",
      value_col      = "TPM",
      value_label    = "Total TPM",
      category_label = "VF gene",
      width = 13, height = gh
    )
  }

  # ---- (7) VF x ARG cross-domain correlation (PIPELINE_V2_GAPS C9) -----
  # Spearman correlation between per-sample VF-function TPM sums and
  # ARG-drug-class TPM sums. R/09 saves drug_class_per_sample_TPM.csv;
  # we read it here (virulome runs after resistome) and join on sample.
  # Significance: per-cell Spearman p-value, BH-adjusted across the
  # whole matrix; stars annotated on the heatmap cells.
  arg_csv <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                       "resistome/drug_class_per_sample_TPM.csv")
  if (file.exists(arg_csv) && nrow(vfdb_fun) > 0 &&
      requireNamespace("pheatmap", quietly = TRUE)) {
    arg_long <- readr::read_csv(arg_csv, show_col_types = FALSE)
    vf_long  <- vfdb_fun |>
      dplyr::group_by(.data$sample, .data$Functions) |>
      dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE), .groups = "drop")

    vf_mat  <- vf_long  |>
      tidyr::pivot_wider(names_from = "Functions", values_from = "TPM",
                         values_fill = 0) |>
      tibble::column_to_rownames("sample") |>
      as.matrix()
    arg_mat <- arg_long |>
      tidyr::pivot_wider(names_from = "DRUG", values_from = "TPM",
                         values_fill = 0) |>
      tibble::column_to_rownames("sample") |>
      as.matrix()

    # Backfill missing samples with zero abundance on both sides so the
    # correlation is computed across the full metadata sample set
    # (rather than the smaller intersect of "samples with detections in
    # both domains"). Zero = "no detection" which is itself the relevant
    # signal for omics correlation, and consistent N across cross-domain
    # analyses is the methodological norm.
    sid_col   <- cfg$metadata$sample_id_col
    meta_full <- tryCatch(
      readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                       show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    controls  <- cfg$metadata$controls %||% character(0)
    all_samples <- if (!is.null(meta_full)) {
      s <- as.character(meta_full[[sid_col]])
      sort(setdiff(s, controls))
    } else {
      sort(union(rownames(vf_mat), rownames(arg_mat)))
    }
    .reindex_zero <- function(mat, samples) {
      out <- matrix(0, nrow = length(samples), ncol = ncol(mat),
                     dimnames = list(samples, colnames(mat)))
      hit <- intersect(samples, rownames(mat))
      if (length(hit) > 0) out[hit, ] <- mat[hit, , drop = FALSE]
      out
    }
    vf_mat  <- .reindex_zero(vf_mat,  all_samples)
    arg_mat <- .reindex_zero(arg_mat, all_samples)
    common <- all_samples
    if (length(common) >= 3 && ncol(vf_mat) >= 2 && ncol(arg_mat) >= 2) {

      r_mat <- matrix(NA_real_, nrow = ncol(vf_mat), ncol = ncol(arg_mat),
                      dimnames = list(colnames(vf_mat), colnames(arg_mat)))
      p_mat <- r_mat
      for (i in seq_len(ncol(vf_mat))) {
        for (j in seq_len(ncol(arg_mat))) {
          ct <- suppressWarnings(stats::cor.test(
            vf_mat[, i], arg_mat[, j], method = "spearman", exact = FALSE
          ))
          r_mat[i, j] <- as.numeric(ct$estimate)
          p_mat[i, j] <- ct$p.value
        }
      }
      # Flatten + adjust across ALL (VF category × drug class) cells as
      # one family. Family-wide method comes from cfg$stats$padjust_method
      # (default BH; flip to BY for arbitrary-dependence FDR).
      padj <- matrix(padjust_p(as.vector(p_mat), cfg),
                     nrow = nrow(p_mat), ncol = ncol(p_mat),
                     dimnames = dimnames(p_mat))
      stars <- ifelse(is.na(padj), "",
        ifelse(padj < 0.001, "***",
        ifelse(padj < 0.01,  "**",
        ifelse(padj < 0.05,  "*", ""))))

      n_vf  <- nrow(r_mat); n_arg <- ncol(r_mat)
      pheatmap::pheatmap(
        r_mat,
        cluster_rows     = TRUE, cluster_cols = TRUE,
        clustering_method = "complete",
        color = grDevices::colorRampPalette(
          c("#2166AC", "white", "#B2182B"))(100),
        breaks            = seq(-1, 1, length.out = 101),
        display_numbers   = stars,
        number_color      = "black",
        fontsize_number   = 10,
        border_color      = "grey80",
        fontsize_row = 10, fontsize_col = 10, fontsize = 10,
        # Title omitted — panel composites strip titles, and the
        # standalone caption belongs in the manuscript figure legend,
        # not baked into the raster.
        main = NA,
        filename = file.path(fig_dir, "vf_arg_correlation_heatmap.png"),
        width  = max(8,  0.55 * n_arg + 4),
        height = max(6,  0.50 * n_vf  + 3)
      )
      pipeline_log(cfg, sprintf(
        "Virulome: vf_arg_correlation_heatmap.png (%d VF functions x %d drug classes, %d samples)",
        n_vf, n_arg, length(common)
      ))
    } else {
      pipeline_log(cfg, sprintf(
        "Virulome: insufficient overlap for VFxARG correlation (samples=%d, vf=%d, arg=%d)",
        length(common), ncol(vf_mat), ncol(arg_mat)
      ))
    }
  } else {
    pipeline_log(cfg,
      "Virulome: ARG drug-class CSV not found or pheatmap missing — VFxARG corr skipped")
  }

  invisible(list(alpha = alpha, vf_totals = vf_totals))
}

# ============================================================================
# Virulome-specific helpers (used here and re-used by R/12 network)
# ============================================================================

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
