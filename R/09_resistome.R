# 09_resistome.R — ARG profile analysis (CARD hits). Mirrors the GT scripts
# ARGNorm_VennDiagram.R, ARGNorm_alphaDiversity.R, ARGNorm_pHeatmap.R,
# ARGNorm_relativeAbundance.R, and DrugNorm_pHeatmap.R (drug-class heatmap).
# Reads the TPM-normalised gene table written by R/03 and re-classifies the
# RESISTANCE strings into per-gene drug classes (single class / Multi-drug /
# optional MLS rollup). Driven entirely by cfg — no hardcoded sample IDs,
# treatments, palettes, or gene-name fixes (cfg$resistome overrides defaults).
# Structural plumbing (load_norm_table, alpha/Venn/heatmap/PCoA helpers) is
# delegated to R/utils_ge_profile.R; this module supplies only the resistome-
# specific drug-class classification, default gene renames, and orchestration.
#
# Outputs (under <project>/<figures_dir>/resistome/):
#   venn_drug_classes.png         Drug classes shared across groups
#   alpha_<metric>_violin.png     Violin + boxplot + KW annotation
#   alpha_<metric>_bar.png        Per-sample bar plot faceted by group
#   arg_abundance_violin.png      Per-(sample,gene) log(TPM+1) by group + KW
#   pcoa.png                      ARG-profile PCoA, ellipses, PERMANOVA/PERMDISP
#   pheatmap_genes.png            Gene × sample, min-max scaled, DRUG sidebar
#   pheatmap_drug_classes.png     Drug-class × sample, min-max scaled
#   drug_relative_abundance.png   Stacked drug-class % per sample
#   drug_total_count.png          Total TPM per drug class (log10, hbar)
#   drug_count_prevalence.png     Paired total / sample-count panels
#   gene_count_prevalence.png     Paired total / sample-count panels (gene)
#
# Datasets (under <project>/<datasets_dir>/resistome/):
#   alpha_diversity.csv, drug_total_TPM.csv, kw_<metric>.txt,
#   beta_permanova.txt, beta_permdisp.txt, beta_pcoa_scores.csv

`%||%` <- function(a, b) if (is.null(a)) b else a

run_resistome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Resistome (CARD)")
  rcfg     <- cfg$resistome %||% list()
  database <- tolower(rcfg$database %||% "card")
  label    <- "Resistome"

  norm <- ge_load_norm_table(cfg, label)
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
  card         <- ge_ensure_group_column(card, cfg, group, label)
  group_levels <- sort(unique(as.character(card[[group]])))
  pal_group    <- resolve_top_level_colors(group, group_levels, cfg) %||%
                  ge_named_palette(group_levels, rcfg$group_colors,
                                    rcfg$group_palette %||% "ggsci::default_nejm")

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
  pal_drug    <- ge_named_palette(drug_levels, rcfg$drug_colors,
                                   rcfg$drug_palette %||% "ggsci::default_ucscgb")

  pipeline_log(cfg, sprintf(
    "Resistome: %d CARD rows, %d genes, %d drug classes, %d samples",
    nrow(card), dplyr::n_distinct(card$GENE),
    length(drug_levels), dplyr::n_distinct(card$sample)
  ))

  # ---- (1) Alpha diversity on gene-level TPM ---------------------------
  alpha_metric <- rcfg$alpha_metric %||% "richness"
  alpha <- ge_compute_alpha(card, group)
  readr::write_csv(alpha, file.path(ds_dir, "alpha_diversity.csv"))

  kw <- ge_alpha_kw(alpha, group, alpha_metric, cfg, label)
  if (!is.null(kw)) {
    capture.output(kw,
                   file = file.path(ds_dir, sprintf("kw_%s.txt", alpha_metric)))
  }
  ge_plot_alpha_violin(alpha, group, alpha_metric, kw, pal_group,
                       file.path(fig_dir, sprintf("alpha_%s_violin.png",
                                                  alpha_metric)))
  ge_plot_alpha_bar(alpha, group, alpha_metric, kw, pal_group,
                    file.path(fig_dir, sprintf("alpha_%s_bar.png",
                                               alpha_metric)))

  # ge_compute_alpha returns both richness and shannon; always emit the
  # shannon companion so the v1.1 manifest contract (ge_alpha_shannon_bar
  # for Fig S4/S7/S12) doesn't depend on cfg$resistome$alpha_metric.
  if (alpha_metric != "shannon") {
    kw_sh <- ge_alpha_kw(alpha, group, "shannon", cfg, label)
    if (!is.null(kw_sh)) {
      capture.output(kw_sh, file = file.path(ds_dir, "kw_shannon.txt"))
    }
    ge_plot_alpha_violin(alpha, group, "shannon", kw_sh, pal_group,
                         file.path(fig_dir, "alpha_shannon_violin.png"))
    ge_plot_alpha_bar(alpha, group, "shannon", kw_sh, pal_group,
                      file.path(fig_dir, "alpha_shannon_bar.png"))
  }

  # ---- (1b) ARG abundance per Treatment (log TPM violin) ---------------
  ge_plot_abundance_violin(
    card, group, pal_group,
    file = file.path(fig_dir, "arg_abundance_violin.png"),
    title = "ARG abundance per Treatment group",
    log_label = label, cfg = cfg
  )

  # ---- (2) Venn of drug classes per group ------------------------------
  ge_plot_category_venn(
    card_drug, category_col = "DRUG", group = group, pal_group = pal_group,
    file = file.path(fig_dir, "venn_drug_classes.png"),
    main_title = sprintf("Drug classes shared across %s groups", group),
    log_label = label, cfg = cfg
  )

  # Gene-level counterpart for v1.1 (Fig 4 panel D): ARGs (genes) shared
  # vs unique across treatment groups, complementing the drug-class Venn.
  ge_plot_category_venn(
    card, category_col = "GENE", group = group, pal_group = pal_group,
    file = file.path(fig_dir, "venn_args.png"),
    main_title = sprintf("ARG genes shared across %s groups", group),
    log_label = label, cfg = cfg
  )

  # ARG UpSet (PIPELINE_V2_GAPS C7) — same set membership as venn_args but
  # rendered as an UpSet plot, which scales beyond 4-5 groups where Venns
  # become unreadable.
  ge_plot_category_upset(
    card, category_col = "GENE", group = group,
    file = file.path(fig_dir, "arg_upset_treatments.png"),
    log_label = label, cfg = cfg
  )

  # ---- (2b) Beta diversity on the gene-level TPM matrix ----------------
  if (isTRUE(rcfg$beta %||% TRUE)) {
    ge_plot_beta(
      card, group, pal_group, opts = rcfg,
      fig_file     = file.path(fig_dir, "pcoa.png"),
      ds_perm_file = file.path(ds_dir,  "beta_permanova.txt"),
      ds_disp_file = file.path(ds_dir,  "beta_permdisp.txt"),
      ds_pcoa_file = file.path(ds_dir,  "beta_pcoa_scores.csv"),
      log_label = label, title_prefix = "ARG profile",
      skip_data_label = "ARG hits", cfg = cfg
    )
  }

  # ---- (3) Gene-level pheatmap with DRUG row annotation ----------------
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    renames <- rcfg$gene_renames %||% .default_resistome_gene_renames()
    card_r  <- card |>
      dplyr::mutate(GENE = ge_apply_literal_renames(.data$GENE, renames))
    ge_plot_gene_heatmap(
      card_r, group, pal_group, pal_drug,
      category_col = "DRUG",
      file    = file.path(fig_dir, "pheatmap_genes.png"),
      palette = rcfg$gene_heatmap_palette %||% c("#0612bd", "#bbbbbd", "#bd0606"),
      top_n   = rcfg$top_n_genes,
      fontsize_row = 8, row_height_factor = 0.18, min_height = 7
    )
  } else {
    pipeline_log(cfg, "Resistome: pheatmap not available — heatmaps skipped")
  }

  # ---- (4) Drug-class pheatmap -----------------------------------------
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(card_drug) > 0) {
    ge_plot_category_heatmap(
      card_drug, category_col = "DRUG", group = group, pal_group = pal_group,
      file    = file.path(fig_dir, "pheatmap_drug_classes.png"),
      palette = rcfg$drug_heatmap_palette %||% c("#f5f06e", "#3df5e9", "#f53d8d")
    )
  }

  # ---- (5) Drug-class relative abundance + total bar -------------------
  drug_totals <- NULL
  if (nrow(card_drug) > 0) {
    drug_totals <- ge_plot_category_relative_abundance(
      card_drug, category_col = "DRUG", group = group, palette = pal_drug,
      fig_file = file.path(fig_dir, "drug_relative_abundance.png"),
      ds_file  = file.path(ds_dir,  "drug_total_TPM.csv"),
      fill_label = "Drug class"
    )
    ge_plot_total_bar(
      drug_totals, category_col = "DRUG", value_col = "Total_TPM",
      palette = pal_drug,
      file    = file.path(fig_dir, "drug_total_count.png"),
      x_label = "Drug class"
    )
    save_count_prevalence(
      card_drug,
      file.path(fig_dir, "drug_count_prevalence.png"),
      category_col   = "DRUG",
      value_col      = "TPM",
      palette        = pal_drug,
      value_label    = "Total TPM",
      category_label = "Drug class",
      width = 12, height = 6
    )
  } else {
    pipeline_log(cfg, "Resistome: no rows with a drug class — RA/total skipped")
  }

  # ---- (6) Gene-level paired count + prevalence ------------------------
  if ("TPM" %in% colnames(card) && nrow(card) > 0) {
    gh <- max(6, 0.20 * dplyr::n_distinct(card$GENE) + 2)
    save_count_prevalence(
      card,
      file.path(fig_dir, "gene_count_prevalence.png"),
      category_col   = "GENE",
      value_col      = "TPM",
      value_label    = "Total TPM",
      category_label = "ARG",
      width = 13, height = gh
    )
  }

  invisible(list(alpha = alpha, drug_totals = drug_totals))
}

# ============================================================================
# Resistome-specific helpers (used here and re-used by R/12 network)
# ============================================================================

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
