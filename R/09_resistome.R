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

  # Pre-pass: compute both per-metric KW p-values BEFORE rendering so we
  # can family-adjust across {richness, shannon} within the domain. The
  # plotters then annotate raw + adjusted p in the title / corner. The
  # cross-module accumulator (R/00_setup.R) still records the raw values
  # and finalises a cross-domain summary at manifest time.
  kw_metrics  <- unique(c(alpha_metric, "shannon"))
  kw_objects  <- setNames(vector("list", length(kw_metrics)), kw_metrics)
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
  # become unreadable. Per-set bars coloured by treatment via pal_group.
  ge_plot_category_upset(
    card, category_col = "GENE", group = group,
    file = file.path(fig_dir, "arg_upset_treatments.png"),
    log_label = label, cfg = cfg, pal_group = pal_group,
    value_label = "ARGs count"
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
      top_n   = rcfg$top_n_genes %||% 30L,
      rank_by = rcfg$gene_rank_by %||% "prevalence",
      fontsize_row = 8, row_height_factor = 0.18, min_height = 7,
      legend_top = TRUE, legend_title = "Scaled TPM (0-1)"
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
    # Per-sample drug-class TPM sums — feeds the v1.2 VF×ARG correlation
    # heatmap computed in R/10 (PIPELINE_V2_GAPS C9).
    card_drug |>
      dplyr::group_by(.data$sample, .data$DRUG) |>
      dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE),
                       .groups = "drop") |>
      readr::write_csv(file.path(ds_dir, "drug_class_per_sample_TPM.csv"))

    drug_totals <- ge_plot_category_relative_abundance(
      card_drug, category_col = "DRUG", group = group, palette = pal_drug,
      fig_file = file.path(fig_dir, "drug_relative_abundance.png"),
      ds_file  = file.path(ds_dir,  "drug_total_TPM.csv"),
      fill_label = "Drug class"
    )
    # Streamgraph companion (2026-06-12 polish) — per-treatment view for
    # composite panel D. Stacked bar above stays for the per-sample
    # standalone slot.
    save_stream_composition(
      card_drug,
      file.path(fig_dir, "drug_stream_abundance.png"),
      category_col = "DRUG",
      group_col    = group,
      value_col    = "TPM",
      palette      = pal_drug,
      fill_label   = "Drug class",
      title        = "Drug class composition by treatment"
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

    # ---- (5b) ARG drug-class circos (PIPELINE_V2_GAPS C8) -------------
    # Distinct from R/12's sample->taxon->gene->drug-class chord. This is
    # a 2-tier chord: outer = drug classes, inner = treatment groups,
    # link width = summed TPM. Shows at a glance which drug classes
    # dominate which treatments without scrolling through the heatmap.
    if (requireNamespace("circlize", quietly = TRUE)) {
      # Long-form (DRUG | group | TPM) — passing this to chordDiagram makes the
      # from / to / value semantics explicit. An n×2 matrix tripped circlize
      # 0.4.18 into reading the cell values as sector labels (lung/hospital
      # showed numeric labels like "109321.92" instead of "Aminoglycoside");
      # the data-frame form sidesteps that codepath. complete() preserves
      # every (drug, group) combination as a zero so the matrix-form's
      # values_fill = 0 behaviour is unchanged.
      drug_grp_long <- card_drug |>
        dplyr::group_by(.data$DRUG, .data[[group]]) |>
        dplyr::summarise(TPM = sum(.data$TPM, na.rm = TRUE),
                         .groups = "drop") |>
        tidyr::complete(DRUG, !!rlang::sym(group),
                        fill = list(TPM = 0))
      drug_levels_used  <- sort(unique(drug_grp_long$DRUG))
      group_levels_used <- sort(unique(as.character(drug_grp_long[[group]])))
      if (length(drug_levels_used) >= 2 && length(group_levels_used) >= 2) {
        grid_col <- c(pal_drug[drug_levels_used],
                      pal_group[group_levels_used])
        circos_png <- file.path(fig_dir, "arg_circos_drugclass.png")
        # 3200px canvas + circle.margin only (no explicit canvas.xlim/ylim)
        # mirrors R/12's per-treatment chord renderer so the long drug-class
        # labels (Aminocoumarin, Fluoroquinolone, Streptogramin, ...) don't
        # clip and the composite quadrant matches the chord panels.
        # Retuned 2026-06-16.
        grDevices::png(circos_png, width = 3200, height = 3200,
                       res = 300, bg = "white")
        tryCatch({
          circlize::circos.clear()
          circlize::circos.par(
            start.degree   = 90,
            gap.degree     = 3,
            circle.margin  = c(0.30, 0.30, 0.30, 0.30),
            unit.circle.segments = 500
          )
          circlize::chordDiagram(
            drug_grp_long,
            grid.col        = grid_col,
            transparency    = 0.4,
            annotationTrack = "grid",
            preAllocateTracks = list(track.height = 0.06)
          )
          circlize::circos.trackPlotRegion(
            track.index = 1,
            panel.fun = function(x, y) {
              sector_idx <- circlize::get.cell.meta.data("sector.index")
              circlize::circos.text(
                circlize::get.cell.meta.data("xcenter"),
                circlize::get.cell.meta.data("ylim")[1],
                sector_idx, facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.7
              )
            }, bg.border = NA
          )
          circlize::circos.clear()
        }, error = function(e) {
          pipeline_log(cfg, sprintf("Resistome arg_circos failed: %s",
                                    conditionMessage(e)))
        })
        grDevices::dev.off()
        pipeline_log(cfg, sprintf(
          "Resistome: arg_circos_drugclass.png (%d drug classes x %d groups)",
          length(drug_levels_used), length(group_levels_used)
        ))
      }
    } else {
      pipeline_log(cfg,
        "Resistome: circlize not available — arg_circos skipped")
    }
  } else {
    pipeline_log(cfg, "Resistome: no rows with a drug class — RA/total skipped")
  }

  # ---- (6) Gene-level paired count + prevalence ------------------------
  if ("TPM" %in% colnames(card) && nrow(card) > 0) {
    gh <- min(48, max(6, 0.20 * dplyr::n_distinct(card$GENE) + 2))  # cap: ggsave aborts >50in on dense gene sets
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

  # ---- (7) Per-organism deep-dive (PIPELINE_V2_GAPS B6) -----------------
  # Two panels per organism: per-sample count by treatment (with KW p)
  # and gggenes-style ARG neighbourhood map (Y=sample, X=contig coords).
  # Hardcoded organism list — matches the v1.1 manuscript spec
  # (Frontiers Figs S6 + S11). Override via cfg$resistome$organisms.
  if (isTRUE(rcfg$species_deepdive %||% TRUE)) {
    organisms <- rcfg$organisms %||% .default_resistome_organisms()
    .resistome_organism_deepdive(cleaned, cfg, group, pal_group, organisms,
                                  fig_dir)
  }

  invisible(list(alpha = alpha, drug_totals = drug_totals))
}

# ============================================================================
# Resistome-specific helpers (used here and re-used by R/12 network)
# ============================================================================

# Default per-organism list for the B6 deep-dive. Each entry has:
#   label : human-readable name (used in figure titles + manifest organism field)
#   match : regex applied to cleaned$noncontaminants$name (species column)
#   slug  : filename-safe slug used in the PNG basename
# Override via cfg$resistome$organisms.
.default_resistome_organisms <- function() {
  list(
    list(
      label = "Clostridioides difficile",
      match = "^Clostridioides difficile",
      slug  = "c_difficile"
    ),
    list(
      label = "Enterobacteriaceae",
      match = "^(Escherichia|Klebsiella|Salmonella|Enterobacter|Citrobacter|Shigella|Proteus|Serratia|Yersinia|Pantoea|Cronobacter|Hafnia|Morganella)",
      slug  = "enterobacteriaceae"
    )
  )
}

# Per-organism deep-dive renderer. Two PNGs per organism (count + gene_map);
# both skip cleanly if the organism is absent or has no ARG hits.
.resistome_organism_deepdive <- function(cleaned, cfg, group, pal_group,
                                          organisms, fig_dir) {
  if (is.null(cleaned$noncontaminants)) return(invisible(NULL))

  sid    <- cfg$metadata$sample_id_col
  noncon <- cleaned$noncontaminants
  if (!group %in% colnames(noncon)) {
    meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                             show_col_types = FALSE)
    noncon <- dplyr::inner_join(noncon, meta[, c(sid, group)],
                                 by = c("sample" = sid))
  }
  abri <- cleaned$abri_kraken2

  for (org in organisms) {
    org_rows <- dplyr::filter(noncon, grepl(org$match, .data$name))
    if (nrow(org_rows) == 0) {
      pipeline_log(cfg, sprintf(
        "Resistome %s: no taxa rows match — deep-dive skipped", org$slug
      ))
      next
    }
    per_sample <- org_rows |>
      dplyr::group_by(.data$sample, .data[[group]]) |>
      dplyr::summarise(count = sum(.data$count, na.rm = TRUE),
                       .groups = "drop")

    kw <- if (dplyr::n_distinct(per_sample[[group]]) >= 2) {
      tryCatch(
        kruskal.test(stats::reformulate(group, "count"), data = per_sample),
        error = function(e) NULL
      )
    } else NULL
    kw_label <- if (!is.null(kw))
      sprintf("Kruskal-Wallis p = %.3g", kw$p.value)
    else "Kruskal-Wallis p = NA"

    p_count <- ggplot2::ggplot(per_sample,
                ggplot2::aes(x = sample, y = .data$count,
                              fill = .data[[group]])) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(stats::reformulate(group),
                           scales = "free_x", nrow = 1) +
      ggplot2::scale_fill_manual(values = pal_group) +
      ggplot2::labs(
        title = sprintf("%s — per-sample count (%s)", org$label, kw_label),
        x = "Sample", y = "Read count"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        plot.title  = ggplot2::element_text(hjust = 0.5, size = 12),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        text        = ggplot2::element_text(size = 11)
      )
    bw <- max(7, 0.5 * dplyr::n_distinct(per_sample$sample) + 3)
    save_panel_ggplot(
      file.path(fig_dir, sprintf("%s_count_per_treatment.png", org$slug)),
      p_count, width = bw, height = 5, dpi = 300
    )
    pipeline_log(cfg, sprintf(
      "Resistome %s: count_per_treatment.png (%d samples)",
      org$slug, dplyr::n_distinct(per_sample$sample)
    ))

    if (!requireNamespace("gggenes", quietly = TRUE)) {
      pipeline_log(cfg, sprintf(
        "Resistome %s: gggenes not available — gene_map skipped",
        org$slug
      ))
      next
    }
    if (is.null(abri) || nrow(abri) == 0) next
    org_taxa <- unique(org_rows$name)
    abri_org <- dplyr::filter(abri,
                               .data$name %in% org_taxa,
                               tolower(.data$DATABASE) == "card",
                               !is.na(.data$sequence),
                               !is.na(.data$START), !is.na(.data$END))
    if (nrow(abri_org) == 0) {
      pipeline_log(cfg, sprintf(
        "Resistome %s: no ARG hits on this organism's contigs — gene_map skipped",
        org$slug
      ))
      next
    }
    abri_org$strand_sign <- if ("STRAND" %in% colnames(abri_org)) {
      ifelse(tolower(as.character(abri_org$STRAND)) %in% c("+", "1", "plus"),
             1L, -1L)
    } else {
      1L
    }

    p_genes <- ggplot2::ggplot(abri_org,
                ggplot2::aes(xmin = .data$START, xmax = .data$END,
                              y = .data$sample, fill = .data$GENE,
                              forward = .data$strand_sign > 0)) +
      gggenes::geom_gene_arrow(arrowhead_height = grid::unit(3, "mm"),
                                arrowhead_width  = grid::unit(2, "mm")) +
      ggplot2::scale_fill_viridis_d(option = "turbo") +
      ggplot2::labs(
        title = sprintf("%s — ARG gene neighbourhoods", org$label),
        x = "Contig coordinate (bp)", y = "Sample", fill = "ARG"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.text     = ggplot2::element_text(size = 9),
        text            = ggplot2::element_text(size = 11)
      )
    gh <- max(4, 0.35 * dplyr::n_distinct(abri_org$sample) + 2)
    save_panel_ggplot(
      file.path(fig_dir, sprintf("%s_gene_map.png", org$slug)),
      p_genes, width = 12, height = gh, dpi = 300
    )
    pipeline_log(cfg, sprintf(
      "Resistome %s: gene_map.png (%d hits, %d samples)",
      org$slug, nrow(abri_org), dplyr::n_distinct(abri_org$sample)
    ))
  }
  invisible(NULL)
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
