# 14_manifest.R — emit a single manifest.json describing what the pipeline
# produced this run. The downstream `metaomics-scribe` manuscript agent reads
# this file (and nothing else from the pipeline) so the two repos couple
# through a stable JSON contract rather than file-naming conventions.
#
# Schema: ../metaomics-scribe/docs/MANIFEST_SCHEMA.md
# Manifest is written next to pipeline_report.html in the report directory
# (`dirname(cfg$outputs$report_html)`). Paths inside stages stay relative to
# project_root; `outputs.project_root` gives the agent the relative walk
# from the manifest's own directory back to project_root.

# Pretty alpha-metric labels keyed by the slugified column names that
# R/06_alpha_diversity.R writes to Datasets/alpha_diversity.csv. The pipeline
# round-trips between pretty and slug forms in several places; this map is
# the canonical inverse so the manifest can quote the original metric names.
# Stored as a list so missing keys return NULL (not "subscript out of
# bounds") under `[[slug]]` lookup.
.ALPHA_METRIC_LABELS <- list(
  berger_parker_index         = "Berger Parker index",
  chao1                       = "Chao1 estimator",
  effective_number_of_species = "Effective number of species",
  fisher_s_alpha              = "Fisher's alpha",
  inverse_simpson_s_index     = "Inverse Simpson's index",
  pielou_s_evenness           = "Pielou's evenness",
  richness                    = "Richness",
  shannon_diversity_index     = "Shannon diversity index",
  simpson_s_index             = "Simpson's index",
  total_counts                = "Total counts",
  # Short-form aliases for the vegan fallback path (R/06 compute_diversity_from_counts
  # emits richness/shannon/simpson directly). Mapping them to the same canonical
  # labels the precomputed wf-metagenomics path uses keeps the manifest "metric"
  # key uniform, so slot YAML filters (e.g. fig01_taxa_overview panel B) match
  # regardless of which alpha source the project used.
  shannon                     = "Shannon diversity index",
  simpson                     = "Simpson's index"
)

# Stage primary-artifact paths (relative to project_root). The manifest
# writer treats a stage as `failed` if enabled in cfg$stages but missing its
# primary artifact AND missing a stage_times entry; presence of either is
# enough to mark `complete`.
.STAGE_PRIMARY <- list(
  taxonomy               = c("figures_dir", "taxonomy/taxa_heatmap.png"),
  relative_abundance     = c("datasets_dir", "relative_abundance/composition_long.csv"),
  alpha_diversity        = c("datasets_dir", "alpha_diversity.csv"),
  beta_diversity         = c("datasets_dir", "beta_diversity_pcoa.csv"),
  differential_abundance = c("datasets_dir", "aldex2/taxa/merged_pairwise_with_overall.tsv"),
  resistome              = c("datasets_dir", "resistome/alpha_diversity.csv"),
  virulome               = c("datasets_dir", "virulome/alpha_diversity.csv"),
  mobilome               = c("datasets_dir", "mobilome/alpha_diversity.csv"),
  network                = c("datasets_dir", "network/gephi_nodes.csv"),
  panels                 = c("datasets_dir", "panels/supplementary_tables.xlsx")
)

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

# ISO-8601 with colon-separated timezone offset ("+01:00" not "+0100").
iso8601 <- function(t) {
  if (is.null(t) || is.na(t)) return(NA_character_)
  s <- format(as.POSIXct(t), "%Y-%m-%dT%H:%M:%S%z")
  sub("([+-]\\d{2})(\\d{2})$", "\\1:\\2", s)
}

# Resolve a path inside `project_root`. Returns NULL if the file is missing
# so the caller can drop the entry from the manifest.
project_path <- function(cfg, ...) {
  p <- file.path(cfg$project_root, ...)
  if (file.exists(p)) p else NULL
}

# Path-under-the-configured-roots helpers. `cfg$outputs$datasets_dir` and
# `cfg$outputs$figures_dir` are "Datasets"/"Figures" in the default template
# and "test_run/Datasets"/"test_run/Figures" in sandbox runs — every stage
# builder must go through these so the manifest follows the config, not a
# hardcoded layout.
ds_path  <- function(cfg, ...) file.path(cfg$outputs$datasets_dir, ...)
fig_path <- function(cfg, ...) file.path(cfg$outputs$figures_dir, ...)

# Path relative to project_root. Used for every `path` field inside stages
# so the manifest stays portable when project_root moves. Always emits
# forward slashes for cross-platform JSON.
rel_to_project_root <- function(cfg, abs_path) {
  if (is.null(abs_path)) return(NULL)
  root <- normalizePath(cfg$project_root, winslash = "/", mustWork = FALSE)
  ap   <- normalizePath(abs_path,         winslash = "/", mustWork = FALSE)
  out  <- sub(paste0("^", root, "/?"), "", ap)
  gsub("\\\\", "/", out)
}

# Read the first few rows of a CSV/TSV and return a flat {col: type} map.
# Falls back to "string" for anything readr can't auto-type.
infer_table_schema <- function(path, sep = ",") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  reader <- if (sep == "\t") readr::read_tsv else readr::read_csv
  raw <- tryCatch(
    suppressWarnings(reader(path, n_max = 50, show_col_types = FALSE,
                            progress = FALSE)),
    error = function(e) NULL
  )
  if (is.null(raw)) return(NULL)
  vapply(raw, function(col) {
    if (is.character(col))           "string"
    else if (is.integer(col))        "integer"
    else if (is.numeric(col))        "double"
    else if (is.logical(col))        "boolean"
    else                              class(col)[[1]]
  }, character(1))
}

# Count rows of a CSV/TSV without loading everything. Used for `row_count`.
table_row_count <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_integer_)
  n <- tryCatch(
    length(readr::read_lines(path, progress = FALSE)) - 1L, # minus header
    error = function(e) NA_integer_
  )
  if (is.na(n) || n < 0) NA_integer_ else as.integer(n)
}

# Build a tables[] entry. Returns NULL if the file is missing so the writer
# can keep the rest of the manifest valid.
table_entry <- function(cfg, rel_path, kind, description, sep = ",") {
  abs <- project_path(cfg, rel_path)
  if (is.null(abs)) return(NULL)
  fmt <- if (sep == "\t") "tsv" else tools::file_ext(rel_path)
  list(
    path        = rel_to_project_root(cfg, abs),
    kind        = kind,
    format      = fmt,
    schema      = as.list(infer_table_schema(abs, sep)),
    row_count   = table_row_count(abs),
    description = description
  )
}

# Build a figures[] entry. Extra annotations (metric/pair/groups/caption_seed)
# are passed as `...` and merged in.
figure_entry <- function(cfg, rel_path, kind, ...) {
  abs <- project_path(cfg, rel_path)
  if (is.null(abs)) return(NULL)
  base <- list(
    path   = rel_to_project_root(cfg, abs),
    kind   = kind,
    format = tools::file_ext(rel_path)
  )
  c(base, list(...))
}

# Glob inside `project_root` and return the matched relative paths (sorted).
project_glob <- function(cfg, pattern) {
  hits <- Sys.glob(file.path(cfg$project_root, pattern))
  if (length(hits) == 0) return(character(0))
  vapply(hits, function(h) rel_to_project_root(cfg, h), character(1),
         USE.NAMES = FALSE) |> sort()
}

# Get a short git SHA for the pipeline checkout. Returns NA if git is not
# available or the directory isn't a repo.
git_short_sha <- function(repo_root) {
  if (is.null(repo_root) || !dir.exists(repo_root)) return(NA_character_)
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  ))
  if (length(out) == 0 || any(is.na(out))) NA_character_ else out[[1]]
}

# Stage status: "complete" if stage_times has an entry OR the stage's primary
# artifact exists on disk, "skipped" if disabled in cfg$stages, "failed"
# otherwise. The filesystem fallback lets the manifest be re-emitted after a
# previous run without rerunning the full pipeline.
#
# Per-stage defaults must mirror run_pipeline.R: `panels` is enabled when the
# cfg key is absent (the runner uses `%||% TRUE` for it). Without this, a
# config that omits `stages: panels` runs the panels stage to disk but the
# manifest silently records zero panel_composite entries.
.STAGE_DEFAULT_ENABLED <- list(panels = TRUE)

stage_status <- function(cfg, name, stage_times) {
  default <- isTRUE(.STAGE_DEFAULT_ENABLED[[name]])
  enabled <- isTRUE(cfg$stages[[name]] %||% default)
  if (!enabled) return("skipped")
  if (!is.null(stage_times[[name]])) return("complete")
  spec <- .STAGE_PRIMARY[[name]]
  if (!is.null(spec)) {
    base <- cfg$outputs[[ spec[[1]] ]]
    if (!is.null(base)) {
      if (file.exists(file.path(cfg$project_root, base, spec[[2]]))) {
        return("complete")
      }
    }
  }
  "failed"
}

# ---------------------------------------------------------------------------
# Per-stage builders. Each returns a stage entry list (or NULL).
# ---------------------------------------------------------------------------

build_stage_taxonomy <- function(cfg, stage_times) {
  status <- stage_status(cfg, "taxonomy", stage_times)
  if (status == "skipped") return(NULL)
  groups <- cfg_group_levels(cfg)
  list(
    status     = status,
    duration_s = stage_times[["taxonomy"]],
    tables     = list(),
    figures    = Filter(Negate(is.null), list(
      figure_entry(cfg, fig_path(cfg, "taxonomy/taxa_venn.png"),
                   kind = "venn", groups = as.list(groups),
                   caption_seed = "Venn diagram of taxa shared and unique across the primary grouping."),
      figure_entry(cfg, fig_path(cfg, "taxonomy/taxa_heatmap.png"),
                   kind = "taxa_heatmap", groups = as.list(groups),
                   caption_seed = "Heatmap of top taxa by mean relative abundance across the primary grouping.")
    ))
  )
}

build_stage_relative_abundance <- function(cfg, stage_times) {
  status <- stage_status(cfg, "relative_abundance", stage_times)
  if (status == "skipped") return(NULL)
  groups <- cfg_group_levels(cfg)
  list(
    status     = status,
    duration_s = stage_times[["relative_abundance"]],
    tables     = Filter(Negate(is.null), list(
      table_entry(cfg, ds_path(cfg, "relative_abundance/composition_long.csv"),
                  kind = "composition_long",
                  description = "Long-form per-sample taxon relative abundance, filtered to count > min_count_per_sample.")
    )),
    figures    = Filter(Negate(is.null), list(
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance.png"),
                   kind = "composition_stacked_bar", groups = as.list(groups),
                   rank = "species",
                   caption_seed = "Stacked-bar composition of the top taxa per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_class.png"),
                   kind = "composition_stacked_bar", groups = as.list(groups),
                   rank = "class",
                   caption_seed = "Stacked-bar composition at class level per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_order.png"),
                   kind = "composition_stacked_bar", groups = as.list(groups),
                   rank = "order",
                   caption_seed = "Stacked-bar composition at order level per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_family.png"),
                   kind = "composition_stacked_bar", groups = as.list(groups),
                   rank = "family",
                   caption_seed = "Stacked-bar composition at family level per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_genus.png"),
                   kind = "composition_stacked_bar", groups = as.list(groups),
                   rank = "genus",
                   caption_seed = "Stacked-bar composition at genus level per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_stream.png"),
                   kind = "composition_stream", groups = as.list(groups),
                   rank = "species",
                   caption_seed = "Smoothed stacked-area composition by treatment group (per-group mean %), species level."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_stream_class.png"),
                   kind = "composition_stream", groups = as.list(groups),
                   rank = "class",
                   caption_seed = "Smoothed stacked-area composition by treatment group, class level."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_stream_order.png"),
                   kind = "composition_stream", groups = as.list(groups),
                   rank = "order",
                   caption_seed = "Smoothed stacked-area composition by treatment group, order level."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_stream_family.png"),
                   kind = "composition_stream", groups = as.list(groups),
                   rank = "family",
                   caption_seed = "Smoothed stacked-area composition by treatment group, family level."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/relative_abundance_stream_genus.png"),
                   kind = "composition_stream", groups = as.list(groups),
                   rank = "genus",
                   caption_seed = "Smoothed stacked-area composition by treatment group, genus level."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/species_count.png"),
                   kind = "species_count",
                   caption_seed = "Distinct species observed per sample, grouped by primary grouping."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/species_count_prevalence.png"),
                   kind = "species_prevalence",
                   caption_seed = "Per-sample distinct-species counts paired with feature prevalence."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/unique_species_relative_abundance.png"),
                   kind = "composition_unique_species", groups = as.list(groups),
                   caption_seed = "Relative abundance of species unique to each treatment group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/shared_species_relative_abundance.png"),
                   kind = "composition_shared_species", groups = as.list(groups),
                   caption_seed = "Relative abundance of species shared across all treatment groups."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/genus_heatmap_top30.png"),
                   kind = "genus_heatmap", groups = as.list(groups),
                   caption_seed = "Heatmap of the top 30 genera by mean relative abundance, hierarchically clustered on both axes (Bray-Curtis, complete linkage).")
    ))
  )
}

build_stage_alpha_diversity <- function(cfg, stage_times) {
  status <- stage_status(cfg, "alpha_diversity", stage_times)
  if (status == "skipped") return(NULL)
  groups <- cfg_group_levels(cfg)

  # Discover per-metric figures via glob; the slug is everything before
  # the final _violin / _bar.
  metric_figs <- list()
  for (kind_suffix in c("violin", "bar")) {
    hits <- project_glob(cfg, fig_path(cfg, paste0("alpha_diversity/*_", kind_suffix, ".png")))
    for (rel in hits) {
      slug <- sub(paste0("_", kind_suffix, "$"), "",
                  tools::file_path_sans_ext(basename(rel)))
      label <- .ALPHA_METRIC_LABELS[[slug]] %||% slug
      metric_figs[[length(metric_figs) + 1]] <- figure_entry(
        cfg, rel,
        kind = paste0("alpha_", kind_suffix),
        metric = label,
        groups = as.list(groups),
        caption_seed = sprintf("%s across the primary grouping (%s).", label, kind_suffix)
      )
    }
  }

  list(
    status     = status,
    duration_s = stage_times[["alpha_diversity"]],
    tables     = Filter(Negate(is.null), list(
      table_entry(cfg, ds_path(cfg, "alpha_diversity.csv"),
                  kind = "alpha_diversity_per_sample",
                  description = "Per-sample alpha-diversity indices (negative controls excluded).")
    )),
    figures    = c(
      metric_figs,
      Filter(Negate(is.null), list(
        figure_entry(cfg, fig_path(cfg, "alpha_diversity/rarefaction_curves.png"),
                     kind = "rarefaction_curves",
                     groups = as.list(groups),
                     caption_seed = "Rarefaction curves for all samples coloured by primary grouping; endpoints mark actual sequencing depth.")
      ))
    ),
    stats_text = Filter(Negate(is.null), list(
      stats_entry(cfg, fig_path(cfg, "alpha_diversity/alpha_stats.txt"),
                  kind = "kruskal_wallis",
                  description = "Per-metric Kruskal-Wallis tests across the primary grouping.")
    ))
  )
}

build_stage_beta_diversity <- function(cfg, stage_times) {
  status <- stage_status(cfg, "beta_diversity", stage_times)
  if (status == "skipped") return(NULL)
  groups <- cfg_group_levels(cfg)
  list(
    status     = status,
    duration_s = stage_times[["beta_diversity"]],
    tables     = Filter(Negate(is.null), list(
      table_entry(cfg, ds_path(cfg, "beta_diversity_pcoa.csv"),
                  kind = "beta_pcoa_scores",
                  description = "PCoA scores from Bray-Curtis on the filtered relative-abundance table."),
      table_entry(cfg, ds_path(cfg, "beta_diversity_jaccard_pcoa.csv"),
                  kind = "beta_pcoa_scores",
                  description = "PCoA scores from Jaccard (presence/absence) on the filtered relative-abundance table."),
      table_entry(cfg, ds_path(cfg, "beta_diversity_padj_summary.csv"),
                  kind = "beta_padj_summary",
                  description = "Family-adjusted PERMANOVA + PERMDISP p-values for {Bray, Jaccard} via cfg$stats$padjust_method."),
      table_entry(cfg, ds_path(cfg, "stats/ge_alpha_kw_padj_summary.csv"),
                  kind = "ge_alpha_kw_padj_summary",
                  description = "Cross-domain GE-side alpha-diversity KW p-values, family-adjusted per metric across resistome / virulome / mobilome via cfg$stats$padjust_method.")
    )),
    figures    = Filter(Negate(is.null), list(
      figure_entry(cfg, fig_path(cfg, "beta_diversity/pcoa.png"),
                   kind = "pcoa_scatter", groups = as.list(groups),
                   caption_seed = "PCoA of Bray-Curtis distances coloured by primary grouping; 95% confidence ellipses."),
      figure_entry(cfg, fig_path(cfg, "beta_diversity/jaccard_pcoa.png"),
                   kind = "pcoa_jaccard", groups = as.list(groups),
                   caption_seed = "PCoA of Jaccard (presence/absence) dissimilarity coloured by primary grouping; PERMANOVA R² and p-value annotated.")
    )),
    stats_text = Filter(Negate(is.null), list(
      stats_entry(cfg, fig_path(cfg, "beta_diversity/permanova.txt"),
                  kind = "permanova",
                  description = "Adonis PERMANOVA on Bray-Curtis distances."),
      stats_entry(cfg, fig_path(cfg, "beta_diversity/permdisp.txt"),
                  kind = "permdisp",
                  description = "PERMDISP homogeneity-of-dispersion test."),
      stats_entry(cfg, fig_path(cfg, "beta_diversity/jaccard_permanova.txt"),
                  kind = "permanova",
                  description = "Adonis PERMANOVA on Jaccard (presence/absence) distances."),
      stats_entry(cfg, fig_path(cfg, "beta_diversity/jaccard_permdisp.txt"),
                  kind = "permdisp",
                  description = "PERMDISP on Jaccard distances.")
    ))
  )
}

build_stage_differential_abundance <- function(cfg, stage_times) {
  status <- stage_status(cfg, "differential_abundance", stage_times)
  if (status == "skipped") return(NULL)

  tables <- list()
  figures <- list()
  for (domain in c("taxa", "gene")) {
    tables <- c(tables, list(
      table_entry(cfg, ds_path(cfg, sprintf("aldex2/%s/merged_pairwise_with_overall.tsv", domain)),
                  kind = "daa_pairwise", sep = "\t",
                  description = sprintf("ALDEx2 pairwise results (%s) merged with overall Kruskal-Wallis.", domain)),
      table_entry(cfg, ds_path(cfg, sprintf("aldex2/%s/overall_kw.tsv", domain)),
                  kind = "daa_overall_kw", sep = "\t",
                  description = sprintf("Overall Kruskal-Wallis BH-adjusted p-values (%s).", domain)),
      table_entry(cfg, ds_path(cfg, sprintf("aldex2/%s/top_candidates_summary.csv", domain)),
                  kind = "daa_top_candidates",
                  description = sprintf("Top DAA candidates ranked across all pairwise contrasts (%s).", domain))
    ))

    # Volcano plots: filenames carry the pair as volcano_<lvlA>_vs_<lvlB>.png
    for (rel in project_glob(cfg, fig_path(cfg, sprintf("differential_abundance/%s/volcano_*.png", domain)))) {
      m <- regmatches(basename(rel),
                      regexec("^volcano_(.+)_vs_(.+)\\.png$", basename(rel)))[[1]]
      pair <- if (length(m) == 3) list(m[[2]], m[[3]]) else NULL
      figures[[length(figures) + 1]] <- figure_entry(
        cfg, rel,
        kind = "volcano",
        domain = domain,
        pair = pair,
        caption_seed = sprintf("Volcano plot for %s contrast (%s vs %s).",
                               domain, pair[[1]] %||% "?", pair[[2]] %||% "?")
      )
    }
    # ALDEx2 MA-plot PNGs (PIPELINE_V2_GAPS C5; taxa only — gene path
    # currently skips MA-plot PNG emission). Same pair-from-filename
    # parser as volcano.
    for (rel in project_glob(cfg, fig_path(cfg, sprintf("differential_abundance/%s/aldex2_maplot_*.png", domain)))) {
      m <- regmatches(basename(rel),
                      regexec("^aldex2_maplot_(.+)_vs_(.+)\\.png$", basename(rel)))[[1]]
      pair <- if (length(m) == 3) list(m[[2]], m[[3]]) else NULL
      figures[[length(figures) + 1]] <- figure_entry(
        cfg, rel,
        kind = "aldex2_maplot",
        domain = domain,
        pair = pair,
        caption_seed = sprintf("ALDEx2 MA plot for %s contrast (%s vs %s). Significant features (BH p < 0.05) coloured by direction.",
                               domain, pair[[1]] %||% "?", pair[[2]] %||% "?")
      )
    }
    figures[[length(figures) + 1]] <- figure_entry(
      cfg, fig_path(cfg, sprintf("differential_abundance/%s/top_candidates.pdf", domain)),
      kind = "top_candidates", domain = domain,
      caption_seed = sprintf("Per-sample abundance for top DAA candidates (%s).", domain))
    figures[[length(figures) + 1]] <- figure_entry(
      cfg, fig_path(cfg, sprintf("differential_abundance/%s/aldex_plots.pdf", domain)),
      kind = "aldex_plot", domain = domain,
      caption_seed = sprintf("ALDEx2 MA / effect plots (%s).", domain))
    # Top-30 DAA sample heatmap (Frontiers Fig S3; taxa only — gene-level
    # PNG isn't currently emitted so figure_entry returns NULL there).
    if (domain == "taxa") {
      figures[[length(figures) + 1]] <- figure_entry(
        cfg, fig_path(cfg, "differential_abundance/taxa/top30_daa_heatmap.png"),
        kind = "taxa_daa_heatmap", domain = "taxa",
        caption_seed = "Sample-level heatmap of the top 30 differentially-abundant taxa across treatment groups.")
      # ALDEx2 cross-comparison dotplot summary (PIPELINE_V2_GAPS C6;
      # taxa only). Figure_entry returns NULL if the PNG isn't on disk
      # (e.g. no features significant in any pair).
      figures[[length(figures) + 1]] <- figure_entry(
        cfg, fig_path(cfg, "differential_abundance/taxa/aldex2_dotplot_summary.png"),
        kind = "aldex2_dotplot", domain = "taxa",
        caption_seed = "Effect-size summary across all pairwise ALDEx2 comparisons; rows are features significant in at least one comparison, columns are pairs.")
    }
  }

  list(
    status     = status,
    duration_s = stage_times[["differential_abundance"]],
    tables     = Filter(Negate(is.null), tables),
    figures    = Filter(Negate(is.null), figures)
  )
}

build_stage_ge_domain <- function(cfg, domain, stage_times) {
  status <- stage_status(cfg, domain, stage_times)
  if (status == "skipped") return(NULL)
  tpm_glob <- project_glob(cfg, ds_path(cfg, sprintf("%s/*_total_TPM.csv", domain)))
  list(
    status     = status,
    duration_s = stage_times[[domain]],
    tables     = c(
      Filter(Negate(is.null), list(
        table_entry(cfg, ds_path(cfg, sprintf("%s/alpha_diversity.csv", domain)),
                    kind = "ge_alpha_diversity",
                    description = sprintf("Per-sample alpha diversity of the %s gene set.", domain)),
        table_entry(cfg, ds_path(cfg, sprintf("%s/beta_pcoa_scores.csv", domain)),
                    kind = "beta_pcoa_scores",
                    description = sprintf("PCoA scores for the %s gene set.", domain))
      )),
      lapply(tpm_glob, function(rel) table_entry(
        cfg, rel, kind = "ge_tpm_totals",
        description = sprintf("Per-sample TPM totals (%s).", domain)
      )) |> Filter(Negate(is.null), x = _)
    ),
    figures    = lapply(
      project_glob(cfg, fig_path(cfg, sprintf("%s/*.png", domain))),
      function(rel) {
        base <- tools::file_path_sans_ext(basename(rel))
        # Filename-dispatch for kinds that don't follow the ge_<filename>
        # convention. Keep the override list short; everything else falls
        # through to the auto-classified ge_* kind.
        if (domain == "resistome" && base == "arg_upset_treatments") {
          return(figure_entry(
            cfg, rel,
            kind   = "arg_upset",
            domain = "resistome",
            groups = as.list(cfg_group_levels(cfg)),
            caption_seed = "UpSet plot of ARG gene families shared and unique across dietary treatment groups."
          ))
        }
        if (domain == "virulome" && base == "vf_upset_treatments") {
          return(figure_entry(
            cfg, rel,
            kind   = "vf_upset",
            domain = "virulome",
            groups = as.list(cfg_group_levels(cfg)),
            caption_seed = "UpSet plot of virulence-factor genes shared and unique across dietary treatment groups."
          ))
        }
        if (domain == "mobilome" && base == "mge_upset_treatments") {
          return(figure_entry(
            cfg, rel,
            kind   = "mge_upset",
            domain = "mobilome",
            groups = as.list(cfg_group_levels(cfg)),
            caption_seed = "UpSet plot of plasmid-replicon genes shared and unique across dietary treatment groups."
          ))
        }
        if (domain == "resistome" && base == "arg_circos_drugclass") {
          return(figure_entry(
            cfg, rel,
            kind   = "arg_circos",
            domain = "resistome",
            caption_seed = "Circos plot showing proportional contribution of antibiotic drug classes to total ARG TPM across dietary treatment groups."
          ))
        }
        # Per-organism deep-dive (B6, Frontiers Figs S6/S11). Hardcoded
        # filenames + organism labels — kept in sync with R/09's default
        # organism list. The `organism` annotation lets the agent target
        # the right panel without parsing filenames.
        if (domain == "resistome" &&
            base %in% c("c_difficile_count_per_treatment",
                        "c_difficile_gene_map")) {
          return(figure_entry(
            cfg, rel,
            kind     = "species_count_genmap",
            domain   = "resistome",
            organism = "Clostridioides difficile",
            caption_seed = if (grepl("count_per_treatment", base))
              "Per-sample Clostridioides difficile count grouped by treatment with Kruskal-Wallis p-value."
            else
              "Per-sample ARG neighbourhood map for Clostridioides difficile."
          ))
        }
        if (domain == "resistome" &&
            base %in% c("enterobacteriaceae_count_per_treatment",
                        "enterobacteriaceae_gene_map")) {
          return(figure_entry(
            cfg, rel,
            kind     = "species_count_genmap",
            domain   = "resistome",
            organism = "Enterobacteriaceae",
            caption_seed = if (grepl("count_per_treatment", base))
              "Per-sample Enterobacteriaceae count grouped by treatment with Kruskal-Wallis p-value."
            else
              "Per-sample ARG neighbourhood map for Enterobacteriaceae."
          ))
        }
        # VFxARG correlation heatmap is cross-domain; rendered from R/10
        # so its file sits under virulome/ but the kind is domain-neutral.
        if (domain == "virulome" && base == "vf_arg_correlation_heatmap") {
          return(figure_entry(
            cfg, rel,
            kind = "vf_arg_corr_heatmap",
            caption_seed = "Spearman correlation heatmap between virulence-factor functional categories and antibiotic drug classes; significance stars from BH-adjusted p-values."
          ))
        }
        figure_entry(
          cfg, rel,
          kind = paste0("ge_", base),
          domain = domain,
          caption_seed = sprintf("%s — %s.", domain, base)
        )
      }
    ),
    stats_text = Filter(Negate(is.null), list(
      stats_entry(cfg, ds_path(cfg, sprintf("%s/beta_permanova.txt", domain)),
                  kind = "permanova",
                  description = sprintf("PERMANOVA on %s beta diversity.", domain)),
      stats_entry(cfg, ds_path(cfg, sprintf("%s/beta_permdisp.txt", domain)),
                  kind = "permdisp",
                  description = sprintf("PERMDISP on %s beta diversity.", domain)),
      stats_entry(cfg, ds_path(cfg, sprintf("%s/kw_richness.txt", domain)),
                  kind = "kruskal_wallis",
                  description = sprintf("Kruskal-Wallis on %s richness across the primary grouping.", domain))
    ))
  )
}

build_stage_network <- function(cfg, stage_times) {
  status <- stage_status(cfg, "network", stage_times)
  if (status == "skipped") return(NULL)
  list(
    status     = status,
    duration_s = stage_times[["network"]],
    tables     = Filter(Negate(is.null), list(
      table_entry(cfg, ds_path(cfg, "network/gephi_nodes.csv"),
                  kind = "network_nodes",
                  description = "Network nodes with degree and module assignment."),
      table_entry(cfg, ds_path(cfg, "network/gephi_edges.csv"),
                  kind = "network_edges",
                  description = "Network edges with correlation weights."),
      table_entry(cfg, ds_path(cfg, "network/topology.csv"),
                  kind = "network_topology",
                  description = "Whole-network topology metrics."),
      table_entry(cfg, ds_path(cfg, "network/chord_long.csv"),
                  kind = "network_chord_long",
                  description = "Long-form chord data: source/target/weight triples."),
      table_entry(cfg, ds_path(cfg, "network/sankey_long.csv"),
                  kind = "network_sankey_long",
                  description = "Long-form sankey data."),
      table_entry(cfg, ds_path(cfg, "network/sample_clusters.csv"),
                  kind = "network_sample_clusters",
                  description = "Per-sample cluster / module assignment."),
      table_entry(cfg, ds_path(cfg, "network/mantel_correlation_triangle.csv"),
                  kind = "mantel_triangle_matrix",
                  description = "Pairwise Mantel correlations (Spearman) between Bray-Curtis distance matrices of the 4 omics layers."),
      table_entry(cfg, ds_path(cfg, "network/mantel_pairwise_padj.csv"),
                  kind = "mantel_pairwise_padj",
                  description = "Per-pair Mantel r, raw permutation p, family-adjusted p (via cfg$stats$padjust_method), and adjustment method name."),
      table_entry(cfg, ds_path(cfg, "network/mobile_arg_fraction_per_sample.csv"),
                  kind = "mobile_arg_fraction_per_sample",
                  description = "Per-sample mobile vs non-mobile ARG TPM totals; an ARG is mobile when its contig co-harbours a PlasmidFinder hit."),
      table_entry(cfg, ds_path(cfg, "network/sankey_taxon_arg_mge_long.csv"),
                  kind = "network_sankey_taxon_arg_mge_long",
                  description = "Long-form (phylum, ARG drug class, MGE replicon family, weight) flow data for the taxon→ARG→MGE Sankey.")
    )),
    figures    = c(
      Filter(Negate(is.null), list(
        figure_entry(cfg, fig_path(cfg, "network/network.png"),
                     kind = "network_graph",
                     caption_seed = "Co-occurrence network coloured by module; edge weight by correlation."),
        figure_entry(cfg, fig_path(cfg, "network/chord_overall.png"),
                     kind = "chord",
                     caption_seed = "Overall chord diagram of taxon co-occurrence."),
        figure_entry(cfg, fig_path(cfg, "network/sankey_overall.html"),
                     kind = "sankey",
                     caption_seed = "Overall sankey diagram of taxon flows."),
        figure_entry(cfg, fig_path(cfg, "network/sankey_overall.png"),
                     kind = "sankey_png",
                     caption_seed = "VF taxa-to-function sankey (PNG render of the interactive supplementary)."),
        figure_entry(cfg, fig_path(cfg, "network/degree_distribution.png"),
                     kind = "degree_distribution",
                     caption_seed = "Network degree distribution."),
        figure_entry(cfg, fig_path(cfg, "network/connectivity_venn_taxa.png"),
                     kind = "connectivity_venn_taxa",
                     groups = as.list(cfg_group_levels(cfg)),
                     caption_seed = "Network-connected taxa shared and unique across treatment groups."),
        figure_entry(cfg, fig_path(cfg, "network/connectivity_venn_genesets.png"),
                     kind = "connectivity_venn_genesets",
                     groups = list("ARGs", "VFs", "MGEs"),
                     caption_seed = "Network-connected gene elements (ARGs / VFs / MGEs)."),
        figure_entry(cfg, fig_path(cfg, "network/mantel_correlation_triangle.png"),
                     kind = "mantel_triangle",
                     caption_seed = "Mantel test correlation triangle showing pairwise community-level Spearman correlations between taxonomy, resistome, virulome, and mobilome distance matrices."),
        figure_entry(cfg, fig_path(cfg, "network/mobile_arg_fraction_bar.png"),
                     kind = "mobile_fraction_bar",
                     groups = as.list(cfg_group_levels(cfg)),
                     caption_seed = "Proportion of ARG TPM carried on predicted mobile contigs (co-harbouring a PlasmidFinder hit) by treatment group; Kruskal-Wallis annotated."),
        figure_entry(cfg, fig_path(cfg, "network/sankey_taxon_arg_mge.png"),
                     kind = "sankey_taxon_arg_mge",
                     caption_seed = "Sankey of mobile-ARG flow from bacterial phylum to ARG drug class to MGE replicon family; ribbons weighted by CARD TPM summed over mobile contigs (CARD + PlasmidFinder co-occurrence).")
      )),
      lapply(project_glob(cfg, fig_path(cfg, "network/chord_*.png")), function(rel) {
        if (basename(rel) == "chord_overall.png") return(NULL)
        grp <- sub("^chord_", "", tools::file_path_sans_ext(basename(rel)))
        figure_entry(cfg, rel, kind = "chord", group = grp,
                     caption_seed = sprintf("Chord diagram for group %s.", grp))
      }) |> Filter(Negate(is.null), x = _),
      lapply(project_glob(cfg, fig_path(cfg, "network/sankey_*.html")), function(rel) {
        if (basename(rel) == "sankey_overall.html") return(NULL)
        grp <- sub("^sankey_", "", tools::file_path_sans_ext(basename(rel)))
        figure_entry(cfg, rel, kind = "sankey", group = grp,
                     caption_seed = sprintf("Sankey diagram for group %s.", grp))
      }) |> Filter(Negate(is.null), x = _)
    )
  )
}

# ---------------------------------------------------------------------------
# Panels stage builder. Scans <figures_dir>/panels/{main,supplementary}/ for
# every emitted composite (TIFF / PNG / PDF) and emits a `figures` array; the
# supplementary tables XLSX lands in `tables`. R/13 runs BEFORE this writer,
# so this builder catalogues what was actually rendered, not predictions.
# ---------------------------------------------------------------------------

build_stage_panels <- function(cfg, stage_times) {
  status <- stage_status(cfg, "panels", stage_times)
  if (status == "skipped") return(NULL)

  # Build a slot -> subsection lookup from the slots YAML so each emitted
  # panel_composite figure carries the manuscript subsection id it belongs
  # to (community_overview, resistome, virulome, mobilome, ...). The
  # downstream agent (metaomics-scribe) groups composites by this field when
  # filling the results-section figure references — see
  # metaomics-scribe/docs/PIPELINE_SUBSECTION_GAP.md. Slots without an
  # explicit `subsection:` simply get NA and the agent surfaces them under
  # an "awaiting subsection assignment" block.
  .manifest_repo_root <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1) {
      return(normalizePath(dirname(dirname(sub("^--file=", "", file_arg[1])))))
    }
    normalizePath(getwd())
  }
  slots_yaml_path <- (cfg$panels %||% list())$slots_yaml %||%
    file.path(.manifest_repo_root(), "templates", "frontiers_v2_slots.yaml")

  subsection_lookup <- character(0)
  if (file.exists(slots_yaml_path)) {
    slots_doc <- yaml::read_yaml(slots_yaml_path)
    combined <- c(slots_doc$main %||% list(),
                  slots_doc$supplementary %||% list())
    subsection_lookup <- vapply(combined, function(spec) {
      val <- spec$subsection
      if (is.null(val)) NA_character_ else as.character(val)
    }, character(1))
    names(subsection_lookup) <- names(combined)
  } else {
    message(sprintf(
      "Panels manifest: slots YAML not found at %s; emitting composites without `subsection`",
      slots_yaml_path
    ))
  }

  panel_globs <- function(section) {
    # One glob per format so the relative paths stay sorted deterministically
    # by format (png, tiff, pdf) within each section.
    fmts <- c("png", "tiff", "pdf")
    out <- character(0)
    for (fmt in fmts) {
      out <- c(out, project_glob(
        cfg, fig_path(cfg, sprintf("panels/%s/*.%s", section, fmt))
      ))
    }
    out
  }

  panel_figures <- function(section) {
    rels <- panel_globs(section)
    lapply(rels, function(rel) {
      slot <- tools::file_path_sans_ext(basename(rel))
      sub  <- subsection_lookup[slot]
      extras <- list(
        kind    = "panel_composite",
        section = section,
        slot    = slot
      )
      # Only include `subsection` when the slots YAML actually carries one —
      # an absent field is more honest than a null in the JSON.
      if (length(sub) == 1L && !is.na(sub) && nzchar(sub)) {
        extras$subsection <- unname(sub)
      }
      extras$caption_seed <- sprintf(
        "Publication composite — slot %s (%s figure).",
        slot, section
      )
      do.call(figure_entry, c(list(cfg, rel), extras))
    })
  }

  list(
    status     = status,
    duration_s = stage_times[["panels"]],
    tables     = Filter(Negate(is.null), list(
      table_entry(cfg, ds_path(cfg, "panels/supplementary_tables.xlsx"),
                  kind = "supplementary_tables",
                  description = paste(
                    "Multi-tab supplementary tables — one sheet per",
                    "table_entry kind across the run, plus an Index sheet."
                  ))
    )),
    figures    = c(panel_figures("main"), panel_figures("supplementary"))
  )
}

# ---------------------------------------------------------------------------
# Stages-index builder. Walks every per-stage builder and returns the
# `stages = { taxonomy = {...}, ... }` dict. Used by `write_manifest` to fill
# manifest.json AND by R/13_panels.R to build its figures_index in-process
# (panels runs before manifest in the current ordering, so it can't read
# manifest.json — it asks for the same data structure directly).
# ---------------------------------------------------------------------------

build_stages_index <- function(cfg, stage_times = list()) {
  Filter(Negate(is.null), list(
    taxonomy               = build_stage_taxonomy(cfg, stage_times),
    relative_abundance     = build_stage_relative_abundance(cfg, stage_times),
    alpha_diversity        = build_stage_alpha_diversity(cfg, stage_times),
    beta_diversity         = build_stage_beta_diversity(cfg, stage_times),
    differential_abundance = build_stage_differential_abundance(cfg, stage_times),
    resistome              = build_stage_ge_domain(cfg, "resistome", stage_times),
    virulome               = build_stage_ge_domain(cfg, "virulome",  stage_times),
    mobilome               = build_stage_ge_domain(cfg, "mobilome",  stage_times),
    network                = build_stage_network(cfg, stage_times),
    panels                 = build_stage_panels(cfg, stage_times)
  ))
}

# ---------------------------------------------------------------------------
# Top-level entry: gather everything and write `<project_root>/manifest.json`.
# ---------------------------------------------------------------------------

# stats_text entries follow the same shape as table_entry minus the schema.
stats_entry <- function(cfg, rel_path, kind, description) {
  abs <- project_path(cfg, rel_path)
  if (is.null(abs)) return(NULL)
  list(
    path        = rel_to_project_root(cfg, abs),
    kind        = kind,
    primary_var = cfg$metadata$group_cols[[1]],
    description = description
  )
}

# Resolve the ordered factor levels for the primary grouping, excluding
# negative controls. Used by figures/groups annotations.
cfg_group_levels <- function(cfg) {
  meta_path <- file.path(cfg$project_root, cfg$metadata$file)
  if (!file.exists(meta_path)) return(character(0))
  meta <- tryCatch(
    readr::read_csv(meta_path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(meta)) return(character(0))
  group_col <- cfg$metadata$group_cols[[1]]
  sid_col   <- cfg$metadata$sample_id_col
  controls  <- cfg$metadata$controls %||% character(0)
  real <- meta[!meta[[sid_col]] %in% controls, , drop = FALSE]
  sort(unique(as.character(real[[group_col]])))
}

write_manifest <- function(cfg,
                           stage_times        = list(),
                           pipeline_start     = NULL,
                           pipeline_end       = NULL,
                           pipeline_repo_root = NULL) {
  if (!isTRUE(cfg$stages$manifest %||% TRUE)) {
    pipeline_log(cfg, "Manifest stage disabled — skipping manifest.json")
    return(invisible(NULL))
  }
  pipeline_log(cfg, "Writing manifest.json")

  # Cross-module padj accumulators get finalised once, here, after every
  # GE-side module has had a chance to append. Currently surfaces the
  # ge_alpha_kw family — per-metric KW p-values across resistome /
  # virulome / mobilome adjusted within each metric.
  padj_summary_finalise(cfg, "ge_alpha_kw")

  meta_path <- file.path(cfg$project_root, cfg$metadata$file)
  meta <- tryCatch(
    readr::read_csv(meta_path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  sid_col  <- cfg$metadata$sample_id_col
  controls <- cfg$metadata$controls %||% character(0)
  n_samples <- if (!is.null(meta)) sum(!meta[[sid_col]] %in% controls) else NA_integer_

  manifest <- list(
    manifest_version = "1.2",
    study = list(
      id                = cfg$study$id,
      name              = cfg$study$name,
      description       = cfg$study$description,
      primary_group_col = cfg$metadata$group_cols[[1]],
      random_effect     = cfg$metadata$random_effect,
      fixed_effects     = as.list(cfg$metadata$fixed_effects %||% list()),
      group_levels      = as.list(cfg_group_levels(cfg)),
      n_samples         = n_samples,
      n_controls        = length(controls),
      control_ids       = as.list(controls)
    ),
    config = list(
      filters = cfg$filters,
      stats   = cfg$stats
    ),
    outputs = list(
      # Relative walk from the manifest's directory back to project_root.
      # The agent resolves stage paths as `manifest_dir / project_root / path`.
      project_root = manifest_project_root_rel(cfg),
      datasets_dir = cfg$outputs$datasets_dir,
      figures_dir  = cfg$outputs$figures_dir,
      log_file     = cfg$outputs$log_file
    ),
    stages = build_stages_index(cfg, stage_times),
    pipeline = list(
      name             = "Metagenomics_pipeline_automation",
      repo             = "https://github.com/Julio92-C/Metagenomics_pipeline_automation",
      version          = git_short_sha(pipeline_repo_root),
      run_started_at   = iso8601(pipeline_start),
      run_finished_at  = iso8601(pipeline_end %||% Sys.time())
    )
  )

  out_dir  <- file.path(cfg$project_root, dirname(cfg$outputs$report_html))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "manifest.json")
  jsonlite::write_json(
    manifest, out_path,
    pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"
  )
  pipeline_log(cfg, sprintf("Manifest written to %s", out_path))
  invisible(out_path)
}

# Relative path from the manifest's directory back to project_root. The
# manifest is placed inside `dirname(cfg$outputs$report_html)`, so the walk
# is one ".." per path component in that subdir.
manifest_project_root_rel <- function(cfg) {
  sub_path <- dirname(cfg$outputs$report_html %||% "")
  if (sub_path %in% c("", ".")) return(".")
  parts <- strsplit(sub_path, "[/\\\\]+")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) "." else paste(rep("..", length(parts)), collapse = "/")
}
