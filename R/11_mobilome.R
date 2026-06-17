# 11_mobilome.R — MGE / plasmid profile analysis (PlasmidFinder hits).
# Mirrors the GT scripts MGEsNorm_relativeAbundance.R + MGEsNorm_pHeatmap.R.
# Reads the TPM-normalised gene table written by R/03 and classifies each
# replicon GENE into a Replicon_Family via configurable regex patterns
# (GT defaults: ^Col -> "Col-like", ^IncF -> "IncF", ^IncX -> "IncX",
# ^Inc -> "Other Inc", else -> "Unknown/Other"). Driven entirely by cfg.
# Structural plumbing lives in R/utils_ge_profile.R; this module supplies
# the replicon-family classifier, default family patterns, and orchestration.
#
# Output suite extends R/09 / R/10 with both gene-level AND family-level
# views because GT only stratifies by GENE but family-level summaries are
# useful for noisier studies.
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
#   mge_count_prevalence.png       Paired total / sample-count (gene)
#   mge_family_count_prevalence.png Paired total / sample-count (family)
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
  label    <- "Mobilome"

  norm <- ge_load_norm_table(cfg, label)
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
  pfdb         <- ge_ensure_group_column(pfdb, cfg, group, label)
  group_levels <- sort(unique(as.character(pfdb[[group]])))
  pal_group    <- resolve_top_level_colors(group, group_levels, cfg) %||%
                  ge_named_palette(group_levels, mcfg$group_colors,
                                    mcfg$group_palette %||% "ggsci::default_nejm")

  # ---- Replicon family classification ----------------------------------
  patterns <- mcfg$family_patterns %||% .default_replicon_patterns()
  pfdb$Replicon_Family <- classify_replicon_family(pfdb$GENE, patterns)

  pfdb_fam      <- dplyr::filter(pfdb, !is.na(.data$Replicon_Family))
  family_levels <- sort(unique(pfdb_fam$Replicon_Family))
  pal_family    <- ge_named_palette(family_levels, mcfg$family_colors,
                                     mcfg$family_palette %||% "ggsci::default_nejm")
  gene_levels   <- sort(unique(pfdb$GENE))
  pal_gene      <- ge_named_palette(gene_levels, mcfg$gene_colors,
                                     mcfg$gene_palette %||% "ggsci::default_ucscgb")

  pipeline_log(cfg, sprintf(
    "Mobilome: %d PlasmidFinder rows, %d replicons, %d families, %d samples",
    nrow(pfdb), length(gene_levels),
    length(family_levels), dplyr::n_distinct(pfdb$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- mcfg$alpha_metric %||% "richness"
  alpha <- ge_compute_alpha(pfdb, group)
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

  # ---- (1b) MGE abundance per Treatment (log TPM violin) ---------------
  ge_plot_abundance_violin(
    pfdb, group, pal_group,
    file = file.path(fig_dir, "mge_abundance_violin.png"),
    title = "MGE abundance per Treatment group",
    log_label = label, cfg = cfg
  )

  # ---- (2) Venn of plasmid replicons (GENE) per group ------------------
  ge_plot_category_venn(
    pfdb, category_col = "GENE", group = group, pal_group = pal_group,
    file = file.path(fig_dir, "venn_replicons.png"),
    main_title = sprintf("Plasmid replicons shared across %s groups", group),
    log_label = label, cfg = cfg
  )

  # MGE UpSet (replicon-level) — UpSet counterpart to venn_replicons.
  ge_plot_category_upset(
    pfdb, category_col = "GENE", group = group,
    file = file.path(fig_dir, "mge_upset_treatments.png"),
    log_label = label, cfg = cfg, pal_group = pal_group
  )

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(mcfg$beta %||% TRUE)) {
    ge_plot_beta(
      pfdb, group, pal_group, opts = mcfg,
      fig_file     = file.path(fig_dir, "pcoa.png"),
      ds_perm_file = file.path(ds_dir,  "beta_permanova.txt"),
      ds_disp_file = file.path(ds_dir,  "beta_permdisp.txt"),
      ds_pcoa_file = file.path(ds_dir,  "beta_pcoa_scores.csv"),
      log_label = label, title_prefix = "MGE profile",
      skip_data_label = "MGE hits", cfg = cfg
    )
  }

  # ---- (3) Gene-level pheatmap with Replicon_Family row annotation -----
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    gene_renames <- mcfg$gene_renames %||% list()
    pfdb_r <- pfdb |>
      dplyr::mutate(GENE = ge_apply_literal_renames(.data$GENE, gene_renames))
    ge_plot_gene_heatmap(
      pfdb_r, group, pal_group, pal_family,
      category_col = "Replicon_Family",
      file    = file.path(fig_dir, "pheatmap_genes.png"),
      palette = mcfg$gene_heatmap_palette %||% c("#ebc3c0", "#bd0606", "#d98c07"),
      top_n   = mcfg$top_n_genes %||% 30L,
      rank_by = mcfg$gene_rank_by %||% "prevalence",
      fontsize_row = 8, row_height_factor = 0.32, min_height = 5,
      # Lock to uniform dims so MGE panel D in figS_diet_effects_mge
      # matches the ARG/VF category-heatmap canvas (8.5 x 6). Legend on
      # top to free horizontal room (same as category heatmaps).
      width = 8.5, height = 6, legend_top = TRUE
    )
  } else {
    pipeline_log(cfg, "Mobilome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Family-level pheatmap ---------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(pfdb_fam) > 0) {
    ge_plot_category_heatmap(
      pfdb_fam, category_col = "Replicon_Family", group = group,
      pal_group = pal_group,
      file    = file.path(fig_dir, "pheatmap_replicon_families.png"),
      palette = mcfg$family_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )
  }

  # ---- (5) Gene & family relative abundance + total bars ---------------
  gene_totals <- ge_plot_category_relative_abundance(
    pfdb, category_col = "GENE", group = group, palette = pal_gene,
    fig_file = file.path(fig_dir, "mge_relative_abundance.png"),
    ds_file  = file.path(ds_dir,  "mge_total_TPM.csv"),
    fill_label = "Plasmid replicon"
  )
  # Streamgraph companions (2026-06-12 polish) — per-treatment view for
  # the MGE composite. The user asked 2026-06-12 to use plasmid replicon
  # GENES (IncF, IncX, ColE1 ...) rather than the family rollup for the
  # main MGE panels, so we emit both: the gene-level stream is what the
  # YAML targets in fig04 / fig08, the family stream stays for
  # supplementary or alternative layouts.
  save_stream_composition(
    pfdb,
    file.path(fig_dir, "mge_stream_abundance.png"),
    category_col = "GENE",
    group_col    = group,
    value_col    = "TPM",
    palette      = pal_gene,
    fill_label   = "Plasmid replicon",
    title        = "MGE plasmid replicon (gene) composition by treatment"
  )
  save_stream_composition(
    pfdb_fam,
    file.path(fig_dir, "mge_family_stream_abundance.png"),
    category_col = "Replicon_Family",
    group_col    = group,
    value_col    = "TPM",
    palette      = pal_family,
    fill_label   = "Replicon family",
    title        = "MGE replicon family composition by treatment"
  )
  ge_plot_total_bar(
    gene_totals, category_col = "GENE", value_col = "Total_TPM",
    palette = pal_gene,
    file    = file.path(fig_dir, "mge_total_count.png"),
    x_label = "Plasmid replicon"
  )
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
    # Per-sample replicon-family TPM sums — feeds the v1.2 Mantel triangle
    # (R/12 reads this when cfg$network$mantel$level == "category") and
    # mirrors R/09's drug_class_per_sample_TPM.csv pattern.
    pfdb_fam |>
      dplyr::group_by(.data$sample, .data$Replicon_Family) |>
      dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE),
                       .groups = "drop") |>
      readr::write_csv(file.path(ds_dir,
                                  "replicon_family_per_sample_TPM.csv"))

    family_totals <- ge_plot_category_relative_abundance(
      pfdb_fam, category_col = "Replicon_Family", group = group,
      palette = pal_family,
      fig_file = file.path(fig_dir, "mge_family_abundance.png"),
      ds_file  = file.path(ds_dir,  "replicon_family_total_TPM.csv"),
      fill_label = "Replicon family"
    )
    ge_plot_total_bar(
      family_totals, category_col = "Replicon_Family", value_col = "Total_TPM",
      palette = pal_family,
      file    = file.path(fig_dir, "mge_family_total_count.png"),
      x_label = "Replicon family",
      height  = 5
    )
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
# Mobilome-specific helpers
# ============================================================================

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
