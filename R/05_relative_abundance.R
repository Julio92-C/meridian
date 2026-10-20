# 05_relative_abundance.R — relative abundance per sample (stacked bar)
# (top-N taxa, remainder collapsed into "Others"), plus a horizontal total-
# count bar for the same taxa. Mirrors relativeAbundance.R and
# bracken_relativeAbundance.R from the reference chicken_batch1 scripts but
# parameterised by config — no hardcoded sample IDs, treatments, or palettes.
#
# Steps
#   1. Apply shared taxa-name cleanup (clean_taxa_names from utils_taxa.R):
#      drop clade-group noise, keep species-level only, strip NCBI brackets,
#      abbreviate genus, truncate long names.
#   2. Apply a per-(taxon, sample) count floor.
#   3. Optionally join metadata so the relative abundance bar can be facetted by group.
#   4. Compute per-sample relative abundance %.
#   5. Identify the top-N taxa by MEAN per-sample relative abundance (not total
#      counts, which favour the deepest samples); collapse everything else into
#      "Others" so the palette only has to cover top_n + 1 levels.
#   6. Relative abundance stacked bar (PNG + plotly HTML) and total-count plot
#      (PNG, log10 scale with median line and count labels). A per-group
#      stacked bar (mean of the per-sample %) is the main figure; the per-sample
#      bar is the supplementary one.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build a palette sized to exactly `n` levels. Caller picks the palette name
# via cfg$relative_abundance$palette; falls back to ggsci::default_igv (matches
# GT bracken_relativeAbundance.R) and finally to hcl.colors if paletteer fails.
relab_palette <- function(n, cfg) {
  name <- cfg$relative_abundance$palette %||% "ggsci::default_igv"
  pal <- tryCatch(
    as.character(paletteer::paletteer_d(name, n = n)),
    error = function(e) NULL
  )
  if (is.null(pal) || length(pal) == 0) {
    # paletteer_d with n > length(palette) errors on some palettes; cycle a
    # broader default rather than failing the stage.
    pal <- tryCatch(
      as.character(paletteer::paletteer_d(name)),
      error = function(e) NULL
    )
  }
  if (is.null(pal) || length(pal) == 0) {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  rep_len(pal, n)
}

# "Others" is always this light grey, in every R/05 figure; taxon colours are
# drawn so they never use it.
RA_OTHERS_COLOUR <- "#BDBDBD"

# Named colour vector over `taxa_levels` (in that order): taxa get relab_palette
# colours minus the Others grey, `others_label` gets RA_OTHERS_COLOUR.
ra_taxa_palette <- function(taxa_levels, others_label, cfg) {
  taxa <- setdiff(taxa_levels, others_label)
  pal  <- relab_palette(length(taxa) + 3L, cfg)
  pal  <- pal[!toupper(pal) %in% toupper(RA_OTHERS_COLOUR)]
  out  <- setNames(rep_len(pal, length(taxa)), taxa)
  if (others_label %in% taxa_levels) out[others_label] <- RA_OTHERS_COLOUR
  out[taxa_levels]
}

# Shared top-N + Others recipe for every R/05 composition figure. Takes
# long-form (sample, name, count, [facet_by]) and returns
#   df_pct      per-(sample, name[, facet]) count + % with non-top taxa
#               collapsed into others_label, `name` a factor in stack order
#   top_taxa    the top_n taxa, ranked by MEAN per-sample % (a sample where the
#               taxon is absent counts as 0, so the ranking is not driven by
#               the deepest samples)
#   taxa_levels legend / stack order (largest total % first, Others last)
ra_prepare_composition <- function(df_in, facet_by, ra_cfg) {
  top_n        <- ra_cfg$top_n        %||% 15
  others_label <- ra_cfg$others_label %||% "Others"

  df_pct <- df_in |>
    dplyr::group_by(sample) |>
    dplyr::mutate(percentage = count / sum(count, na.rm = TRUE) * 100) |>
    dplyr::ungroup()
  n_samples <- dplyr::n_distinct(df_pct$sample)

  top_taxa <- df_pct |>
    dplyr::group_by(name) |>
    dplyr::summarise(mean_pct = sum(percentage, na.rm = TRUE) / n_samples,
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mean_pct), name) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(name)

  df_pct <- df_pct |>
    dplyr::mutate(name = ifelse(name %in% top_taxa, name, others_label)) |>
    dplyr::group_by(sample, name, dplyr::across(dplyr::any_of(facet_by))) |>
    dplyr::summarise(
      count      = sum(count, na.rm = TRUE),
      percentage = sum(percentage, na.rm = TRUE),
      .groups    = "drop"
    )

  # Order taxa globally by total % so the legend matches the largest stacks.
  taxa_levels <- df_pct |>
    dplyr::group_by(name) |>
    dplyr::summarise(total = sum(percentage, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::pull(name)
  # "Others" pinned to the bottom of the stack regardless of size.
  taxa_levels <- c(setdiff(taxa_levels, others_label),
                   intersect(others_label, taxa_levels))
  df_pct$name <- factor(df_pct$name, levels = taxa_levels)

  list(df_pct = df_pct, top_taxa = top_taxa, taxa_levels = taxa_levels,
       others_label = others_label, top_n = top_n)
}

# Per-group mean of the per-sample % (top-N + Others frame from
# ra_prepare_composition). A taxon absent from a sample counts as 0, so each
# group's bar sums to 100. `label` carries "<group> (n=<samples>)".
ra_group_means <- function(df_pct, facet_by) {
  d <- df_pct |>
    dplyr::mutate(group = as.character(.data[[facet_by]]))
  n_per_group <- d |>
    dplyr::distinct(sample, group) |>
    dplyr::count(group, name = "n_samples")
  d |>
    dplyr::group_by(group, name) |>
    dplyr::summarise(total = sum(percentage, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(n_per_group, by = "group") |>
    dplyr::mutate(mean_percentage = total / n_samples,
                  label = sprintf("%s (n=%d)", group, n_samples)) |>
    dplyr::select(group, name, n_samples, mean_percentage, label)
}

# Compact stacked-bar render for filtered species partitions (v1.1
# unique/shared composition plots, Frontiers Fig S1/S2). Mirrors the
# main 6a render in run_relative_abundance() but skips the plotly HTML
# and dataset CSV — callers only need the PNG. Takes long-form
# (sample, name, count, [facet_by]) and applies the same top-N + Others
# + taxa ordering recipe (ra_prepare_composition) so the legend style matches.
ra_render_stacked_bar <- function(df_in, facet_by, ra_cfg, cfg,
                                   file_path, log_label) {
  if (nrow(df_in) == 0) {
    pipeline_log(cfg, sprintf("%s: no rows — skipping %s", log_label,
                              basename(file_path)))
    return(invisible(NULL))
  }
  prep         <- ra_prepare_composition(df_in, facet_by, ra_cfg)
  df_pct       <- prep$df_pct
  taxa_levels  <- prep$taxa_levels

  n_levels    <- length(taxa_levels)
  palette     <- ra_taxa_palette(taxa_levels, prep$others_label, cfg)
  show_legend <- ra_cfg$show_legend %||% TRUE
  legend_pos  <- if (isTRUE(show_legend)) "top" else "none"

  p <- ggplot2::ggplot(df_pct,
                       ggplot2::aes(x = factor(sample), y = percentage, fill = name)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)", fill = "Taxa") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = legend_pos,
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      text            = ggplot2::element_text(size = 12)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 10))
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                                 scales = "free_x", nrow = 1)
  }
  n_samples <- length(unique(df_pct$sample))
  pw <- max(10, 0.5 * n_samples + 3)
  ph <- if (isTRUE(show_legend)) max(7, 6 + ceiling(n_levels / 10) * 0.4) else 6
  save_panel_ggplot(file_path, p, width = pw, height = ph, dpi = 300)
  pipeline_log(cfg, sprintf("%s: %d species, %d samples → %s",
                            log_label, dplyr::n_distinct(df_in$name),
                            n_samples, basename(file_path)))
  invisible(df_pct)
}

run_relative_abundance <- function(cleaned, cfg) {
  pipeline_log(cfg, "Relative abundance")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "relative_abundance")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir,  "relative_abundance")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  # ---- 1. Species rank only, then clean taxa names -----------------------
  # R/02 carries each row's kraken2 rank code in `tax_rank`. Keep "S" only so
  # strain / subspecies rows (S1, S2, ...) don't double-count their species.
  raw_nc <- cleaned$noncontaminants
  if ("tax_rank" %in% colnames(raw_nc)) {
    n_rows <- nrow(raw_nc)
    raw_nc <- dplyr::filter(raw_nc, .data$tax_rank %in% "S")
    pipeline_log(cfg, sprintf(
      "Relative abundance: %d/%d rows at species rank (kraken2 rank S)",
      nrow(raw_nc), n_rows
    ))
  } else {
    pipeline_log(cfg,
      "Relative abundance: no tax_rank column — strain rows not filtered")
  }
  df <- clean_taxa_names(raw_nc, cfg)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Relative abundance: no rows after clean_taxa_names — skipping")
    return(invisible(NULL))
  }

  # ---- 2. Per-(taxon, sample) count floor ---------------------------------
  ra_cfg    <- cfg$relative_abundance %||% list()
  min_count <- ra_cfg$min_count_per_sample %||%
               cfg$filters$min_count_per_sample %||% 5
  n_before  <- nrow(df)
  df <- dplyr::filter(df, !is.na(name), count > min_count)
  pipeline_log(cfg, sprintf(
    "Relative abundance: %d/%d (taxon, sample) rows kept with count > %d",
    nrow(df), n_before, min_count
  ))
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Relative abundance: nothing left after count filter — skipping")
    return(invisible(NULL))
  }

  # ---- 3. Optional metadata join for facetting ----------------------------
  facet_by <- ra_cfg$facet_by %||% cfg$metadata$group_cols[[1]]
  if (!is.null(facet_by) && nzchar(facet_by)) {
    sid  <- cfg$metadata$sample_id_col
    meta <- readr::read_csv(
      file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
    )
    if (facet_by %in% colnames(meta)) {
      df <- dplyr::inner_join(df, meta[, c(sid, facet_by)],
                              by = c("sample" = sid))
    } else {
      pipeline_log(cfg, sprintf(
        "Relative abundance: facet column '%s' not in metadata — single panel",
        facet_by
      ))
      facet_by <- NULL
    }
  }

  # ---- 4-5. Per-sample %, top-N by mean % + collapse rest into "Others" --
  prep         <- ra_prepare_composition(df, facet_by, ra_cfg)
  df_pct       <- prep$df_pct
  top_taxa     <- prep$top_taxa
  taxa_levels  <- prep$taxa_levels
  top_n        <- prep$top_n
  others_label <- prep$others_label

  # Long-form tidy table for the dashboard composition card. Kept separate
  # from the "Others"-collapsed plotting frame so the report can re-rank by
  # cfg$report$rank_taxa_by without being locked to top_n from this stage.
  readr::write_csv(df_pct, file.path(ds_dir, "composition_long.csv"))

  n_levels <- length(taxa_levels)
  palette  <- ra_taxa_palette(taxa_levels, others_label, cfg)
  pipeline_log(cfg, sprintf(
    "Relative abundance: %d top taxa + %s (palette sized to %d)",
    min(top_n, length(top_taxa)),
    if (others_label %in% taxa_levels) sprintf("'%s'", others_label) else "no overflow",
    n_levels
  ))

  # ---- 6a. Relative abundance stacked bar --------------------------------
  show_legend <- ra_cfg$show_species_legend %||% TRUE
  legend_pos  <- if (isTRUE(show_legend)) "top" else "none"

  # Per-project fontsize knobs. R/13 scales this standalone source PNG down
  # into the fig05 composite cell, so cohorts with high sample counts (lung
  # 42, hospital 64) need smaller at-source axis text to remain legible
  # after the scale-down. Defaults preserve chicken_batch1's tuned look.
  bar_axis_size   <- ra_cfg$bar_axis_text_size   %||% 12
  bar_legend_size <- ra_cfg$bar_legend_text_size %||% 12

  p <- ggplot2::ggplot(df_pct,
                       ggplot2::aes(x = factor(sample), y = percentage, fill = name)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Sample", y = "Relative abundance (%)", fill = "Taxa") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = legend_pos,
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                              size = bar_axis_size),
      legend.text     = ggplot2::element_text(size = bar_legend_size),
      legend.title    = ggplot2::element_text(size = bar_legend_size + 1),
      text            = ggplot2::element_text(size = 12)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 10))

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                                 scales = "free_x", nrow = 1)
  }

  # Size canvas to sample count so labels stay readable across studies.
  n_samples <- length(unique(df_pct$sample))
  pw <- max(10, 0.5 * n_samples + 3)
  ph <- if (isTRUE(show_legend)) max(7, 6 + ceiling(n_levels / 10) * 0.4) else 6

  save_panel_ggplot(file.path(fig_dir, "relative_abundance.png"), p,
                  width = pw, height = ph, dpi = 300)

  # ---- 6a-group. Per-group mean composition (main figure) ----------------
  # One stacked bar per level of facet_by: mean of the per-sample % over that
  # group's samples, same top-N taxa and colours as the per-sample bar above.
  if (!is.null(facet_by) && facet_by %in% colnames(df_pct)) {
    grp_df <- ra_group_means(df_pct, facet_by)
    grp_df$name  <- factor(grp_df$name, levels = taxa_levels)
    grp_lvls     <- unique(grp_df$label[order(grp_df$group)])
    grp_df$label <- factor(grp_df$label, levels = grp_lvls)
    readr::write_csv(
      dplyr::transmute(grp_df, !!facet_by := group, taxon = as.character(name),
                       n_samples, mean_percentage),
      file.path(ds_dir, "relative_abundance_group.csv")
    )
    pg <- ggplot2::ggplot(grp_df,
                          ggplot2::aes(x = label, y = mean_percentage, fill = name)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::labs(x = facet_by, y = "Mean relative abundance (%)",
                    fill = "Taxa") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = if (isTRUE(show_legend)) "right" else "none",
        # Rotate the "<group> (n=..)" labels only when there are enough groups
        # for them to collide.
        axis.text.x     = ggplot2::element_text(
          angle = if (length(grp_lvls) > 4) 45 else 0,
          hjust = if (length(grp_lvls) > 4) 1 else 0.5,
          size  = bar_axis_size),
        legend.text     = ggplot2::element_text(size = bar_legend_size),
        legend.title    = ggplot2::element_text(size = bar_legend_size + 1),
        text            = ggplot2::element_text(size = 12)
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
    n_grp <- length(grp_lvls)
    save_panel_ggplot(file.path(fig_dir, "relative_abundance_group.png"), pg,
                      width = max(7, 0.9 * n_grp + 5),
                      height = max(6, 0.3 * n_levels + 2), dpi = 300)
    pipeline_log(cfg, sprintf(
      "Relative abundance: relative_abundance_group.png (%d groups)", n_grp))
  } else {
    pipeline_log(cfg,
      "Relative abundance: no facet column — group composition skipped")
  }

  # ---- 6a-stream. Per-treatment streamgraph (2026-06-12 polish) ---------
  # Mirrors templates/panels_ref/relativate_abundace.png — per-(treatment,
  # category) mean % as a smoothed stacked area. Lives in composite panel
  # F where the small panel size makes the per-sample stacked bar above
  # hard to read. Species-level here; per-rank streams emitted in the
  # 6a-bis loop below.
  #
  # Palette: project-wide relab_palette keyed by taxon name so the stream
  # colors match the standalone stacked bar exactly. Categories outside
  # the top_n (and the Others bin) still get a colour via the named
  # vector — scale_fill_manual ignores extras.
  if (!is.null(facet_by)) {
    stream_palette <- setNames(palette, taxa_levels)
    save_stream_composition(
      df,
      file.path(fig_dir, "relative_abundance_stream.png"),
      category_col = "name",
      group_col    = facet_by,
      value_col    = "count",
      palette      = stream_palette,
      categories   = top_taxa,
      top_n        = top_n,
      others_label = others_label,
      fill_label   = "Species",
      title        = paste("Species composition by", facet_by)
    )
  }

  # ---- 6a-bis. Per-rank composition stacked bars (2026-06-12 polish) ----
  # Roll up the species-level df to class / order / family / genus via the
  # kraken2 ancestry columns now carried on cleaned$noncontaminants
  # (build_taxid_ancestry, R/02). Each rank gets its own
  # relative_abundance_<rank>.png plus an entry in the manifest tagged
  # with `rank = <rank>` so the panels stage can resolve a specific rank
  # for the composite-D streamgraph or for figS4-style breakdowns.
  for (rank_col in c("class", "order", "family", "genus")) {
    if (!rank_col %in% colnames(df)) next
    df_rank <- df |>
      dplyr::filter(!is.na(.data[[rank_col]]),
                    nzchar(as.character(.data[[rank_col]]))) |>
      dplyr::group_by(sample, .data[[rank_col]],
                       dplyr::across(dplyr::any_of(facet_by))) |>
      dplyr::summarise(count = sum(.data$count, na.rm = TRUE),
                       .groups = "drop") |>
      dplyr::rename(name = !!rank_col)
    if (nrow(df_rank) == 0 ||
        dplyr::n_distinct(df_rank$name) < 2) {
      pipeline_log(cfg, sprintf(
        "Relative abundance (%s): not enough categories — skipping", rank_col
      ))
      next
    }
    ra_render_stacked_bar(
      df_rank,
      facet_by = facet_by,
      ra_cfg   = ra_cfg,
      cfg      = cfg,
      file_path = file.path(fig_dir,
                             sprintf("relative_abundance_%s.png", rank_col)),
      log_label = sprintf("Relative abundance (%s)", rank_col)
    )
    # Per-rank streamgraph companion. Uses the per-rank df_rank (already
    # rolled up to (sample, rank_value) and with facet_by attached) so
    # the rank stream matches the rank stacked bar exactly. Palette is
    # built per rank via relab_palette (the project core palette) so
    # colours are consistent with the rest of the pipeline.
    if (!is.null(facet_by) && facet_by %in% colnames(df_rank)) {
      rank_prep    <- ra_prepare_composition(df_rank, facet_by, ra_cfg)
      rank_palette <- ra_taxa_palette(rank_prep$taxa_levels,
                                      rank_prep$others_label, cfg)
      save_stream_composition(
        df_rank,
        file.path(fig_dir,
                   sprintf("relative_abundance_stream_%s.png", rank_col)),
        category_col = "name",
        group_col    = facet_by,
        value_col    = "count",
        palette      = rank_palette,
        categories   = rank_prep$top_taxa,
        top_n        = rank_prep$top_n,
        others_label = rank_prep$others_label,
        # Capitalised display label for fig01 F (taxa class stream).
        fill_label   = paste0(toupper(substr(rank_col, 1, 1)),
                                substr(rank_col, 2, nchar(rank_col))),
        title        = sprintf("Composition by %s (%s level)", facet_by, rank_col)
      )
    }
  }

  if (requireNamespace("plotly", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    # selfcontained=TRUE still extracts plotly assets into <name>_files/ during
    # render; on re-run dir.create() warns that the libdir already exists.
    # Warnings are benign housekeeping — suppress to keep the log clean.
    suppressWarnings(
      htmlwidgets::saveWidget(
        plotly::ggplotly(p),
        file = file.path(fig_dir, "relative_abundance.html"),
        selfcontained = TRUE
      )
    )
  }

  # ---- 6b. Horizontal total-count plot ------------------------------------
  min_for_species <- ra_cfg$min_count_for_species %||%
                     cfg$filters$min_count_for_species %||% 2000

  df_species <- df_pct |>
    dplyr::filter(name != others_label) |>
    dplyr::group_by(name) |>
    dplyr::summarise(totCount = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(totCount >= min_for_species) |>
    dplyr::arrange(dplyr::desc(totCount))

  if (nrow(df_species) > 0) {
    name_colors <- setNames(
      palette[match(as.character(df_species$name), taxa_levels)],
      df_species$name
    )

    sp <- ggplot2::ggplot(
        df_species,
        ggplot2::aes(x = stats::reorder(name, totCount), y = totCount, fill = name)
      ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_hline(yintercept = stats::median(df_species$totCount),
                          linetype = "dashed", colour = "black") +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(round(totCount / 1000, 1), "K")),
        hjust = -0.05, size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = name_colors) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Taxon", y = "Total count", fill = "Taxon") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 12)) +
      ggplot2::coord_flip()

    sph <- max(5, 0.25 * nrow(df_species) + 2)
    save_panel_ggplot(file.path(fig_dir, "species_count.png"), sp,
                    width = 9, height = sph, dpi = 300)

    # ---- Paired count + prevalence (total | sample prevalence) ---------
    df_long_species <- df_pct |>
      dplyr::filter(name %in% as.character(df_species$name))
    save_count_prevalence(
      df_long_species,
      file.path(fig_dir, "species_count_prevalence.png"),
      category_col   = "name",
      value_col      = "count",
      palette        = name_colors,
      value_label    = "Total count",
      category_label = "Taxon",
      width = 13, height = sph
    )
  } else {
    pipeline_log(cfg, sprintf(
      "Relative abundance: no taxa with total count >= %d — species_count.png skipped",
      min_for_species
    ))
  }

  # ---- 7. Group-membership partitions (Frontiers Fig S1 / S2) -------------
  # Filtered composition plots: species observed in exactly one treatment
  # group ("unique") and species observed in every group ("shared"). Both
  # are derived from `df` (post-clean, post-count-filter, post-metadata
  # join) before the Others collapse. Skipped when no facet column.
  if (!is.null(facet_by) && facet_by %in% colnames(df)) {
    group_species <- split(as.character(df$name), as.character(df[[facet_by]]))
    group_species <- lapply(group_species, unique)
    all_groups    <- names(group_species)
    if (length(all_groups) >= 2) {
      freq           <- table(unlist(group_species))
      unique_species <- names(freq[freq == 1])
      shared_species <- Reduce(intersect, group_species)

      ra_render_stacked_bar(
        dplyr::filter(df, name %in% unique_species),
        facet_by, ra_cfg, cfg,
        file_path = file.path(fig_dir, "unique_species_relative_abundance.png"),
        log_label = "Relative abundance (unique species)"
      )
      ra_render_stacked_bar(
        dplyr::filter(df, name %in% shared_species),
        facet_by, ra_cfg, cfg,
        file_path = file.path(fig_dir, "shared_species_relative_abundance.png"),
        log_label = "Relative abundance (shared species)"
      )
    } else {
      pipeline_log(cfg, sprintf(
        "Relative abundance: only %d group(s) — unique/shared composition skipped",
        length(all_groups)
      ))
    }
  } else {
    pipeline_log(cfg,
      "Relative abundance: no facet column — unique/shared composition skipped")
  }

  # ---- 8. Top-30 genus heatmap (PIPELINE_V2_GAPS C4) ----------------------
  # Genus rollup uses the kraken2-derived `genus` column attached in R/02
  # via `build_taxid_ancestry` — every row's genus is the bare-G ancestor
  # kraken2 placed it under. Higher-rank rows (phylum / class / order /
  # family) and unclassified rows have genus = NA and are dropped here, so
  # the heatmap shows only true kraken2 genera (no first-word string
  # heuristics leaking phyla / families through). Heatmap values are
  # log10(% + 0.01) for colour scale; clustering uses raw percentages with
  # Bray-Curtis distance + complete linkage on both axes per the v1.2 spec.
  if (!isFALSE(ra_cfg$genus_heatmap %||% TRUE) &&
      requireNamespace("pheatmap", quietly = TRUE) &&
      !is.null(cleaned$noncontaminants) &&
      "genus" %in% colnames(cleaned$noncontaminants)) {
    raw <- cleaned$noncontaminants
    raw <- dplyr::filter(raw, !is.na(.data$name), .data$count > min_count)
    raw <- dplyr::filter(raw, !is.na(.data$genus), nzchar(.data$genus))

    meta_g <- tryCatch(
      readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                       show_col_types = FALSE),
      error = function(e) NULL
    )
    has_group <- !is.null(facet_by) && !is.null(meta_g) &&
                 facet_by %in% colnames(meta_g)
    if (has_group) {
      raw <- dplyr::inner_join(
        raw,
        meta_g[, c(cfg$metadata$sample_id_col, facet_by)],
        by = c("sample" = cfg$metadata$sample_id_col)
      )
    }

    if (nrow(raw) > 0 && dplyr::n_distinct(raw$genus) >= 2) {
      samp_total <- raw |>
        dplyr::group_by(.data$sample) |>
        dplyr::summarise(samp_total = sum(.data$count, na.rm = TRUE),
                         .groups = "drop")
      g_long <- raw |>
        dplyr::group_by(.data$sample, .data$genus) |>
        dplyr::summarise(count = sum(.data$count, na.rm = TRUE),
                         .groups = "drop") |>
        dplyr::inner_join(samp_total, by = "sample") |>
        dplyr::mutate(pct = 100 * .data$count / .data$samp_total)

      top_g <- g_long |>
        dplyr::group_by(.data$genus) |>
        dplyr::summarise(mean_pct = mean(.data$pct, na.rm = TRUE),
                         .groups = "drop") |>
        dplyr::arrange(dplyr::desc(.data$mean_pct)) |>
        dplyr::slice_head(n = 30L) |>
        dplyr::pull(.data$genus)

      g_wide <- g_long |>
        dplyr::filter(.data$genus %in% top_g) |>
        dplyr::select("sample", "genus", "pct") |>
        tidyr::pivot_wider(names_from = "sample", values_from = "pct",
                           values_fill = 0) |>
        tibble::column_to_rownames("genus") |>
        as.matrix()

      d_rows <- vegan::vegdist(g_wide,     method = "bray")
      d_cols <- vegan::vegdist(t(g_wide),  method = "bray")
      z      <- log10(g_wide + 0.01)

      samp_ann <- NA
      ann_colors <- NA
      if (has_group) {
        samp_ann <- raw |>
          dplyr::distinct(.data$sample, .keep_all = TRUE) |>
          dplyr::select("sample", dplyr::all_of(facet_by)) |>
          tibble::column_to_rownames("sample")
        samp_ann <- samp_ann[colnames(z), , drop = FALSE]
        ann_lvls <- sort(unique(as.character(samp_ann[[facet_by]])))
        ann_pal  <- resolve_top_level_colors(facet_by, ann_lvls, cfg)
        if (is.null(ann_pal)) {
          ann_pal <- relab_palette(length(ann_lvls), cfg)
          names(ann_pal) <- ann_lvls
        }
        ann_colors <- setNames(list(ann_pal), facet_by)
      }

      pheatmap::pheatmap(
        z,
        clustering_distance_rows = d_rows,
        clustering_distance_cols = d_cols,
        clustering_method        = "complete",
        color = grDevices::colorRampPalette(
          c("white", "#fee08b", "#d73027"))(100),
        annotation_col    = samp_ann,
        annotation_colors = ann_colors,
        border_color      = NA,
        fontsize_row = 8, fontsize_col = 8, fontsize = 8,
        filename = file.path(fig_dir, "genus_heatmap_top30.png"),
        width    = max(8, 0.4  * ncol(z) + 4),
        height   = max(6, 0.25 * nrow(z) + 2.5)
      )
      pipeline_log(cfg, sprintf(
        "Relative abundance: genus_heatmap_top30.png (%d genera x %d samples)",
        nrow(z), ncol(z)
      ))
    } else {
      pipeline_log(cfg,
        "Relative abundance: not enough genus data — heatmap skipped")
    }
  }

  invisible(df_pct)
}
