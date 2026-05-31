# 13_manifest.R — emit a single manifest.json describing what the pipeline
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
  effective_number_of_species = "Effective number of species",
  fisher_s_alpha              = "Fisher's alpha",
  inverse_simpson_s_index     = "Inverse Simpson's index",
  pielou_s_evenness           = "Pielou's evenness",
  richness                    = "Richness",
  shannon_diversity_index     = "Shannon diversity index",
  simpson_s_index             = "Simpson's index",
  total_counts                = "Total counts"
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
  network                = c("datasets_dir", "network/gephi_nodes.csv")
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
stage_status <- function(cfg, name, stage_times) {
  enabled <- isTRUE(cfg$stages[[name]])
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
                   caption_seed = "Stacked-bar composition of the top taxa per sample, faceted by group."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/species_count.png"),
                   kind = "species_count",
                   caption_seed = "Distinct species observed per sample, grouped by primary grouping."),
      figure_entry(cfg, fig_path(cfg, "relative_abundance/species_count_prevalence.png"),
                   kind = "species_prevalence",
                   caption_seed = "Per-sample distinct-species counts paired with feature prevalence.")
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
    figures    = metric_figs,
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
                  description = "PCoA scores from Bray-Curtis on the filtered relative-abundance table.")
    )),
    figures    = Filter(Negate(is.null), list(
      figure_entry(cfg, fig_path(cfg, "beta_diversity/pcoa.png"),
                   kind = "pcoa_scatter", groups = as.list(groups),
                   caption_seed = "PCoA of Bray-Curtis distances coloured by primary grouping; 95% confidence ellipses.")
    )),
    stats_text = Filter(Negate(is.null), list(
      stats_entry(cfg, fig_path(cfg, "beta_diversity/permanova.txt"),
                  kind = "permanova",
                  description = "Adonis PERMANOVA on Bray-Curtis distances."),
      stats_entry(cfg, fig_path(cfg, "beta_diversity/permdisp.txt"),
                  kind = "permdisp",
                  description = "PERMDISP homogeneity-of-dispersion test.")
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
    figures[[length(figures) + 1]] <- figure_entry(
      cfg, fig_path(cfg, sprintf("differential_abundance/%s/top_candidates.pdf", domain)),
      kind = "top_candidates", domain = domain,
      caption_seed = sprintf("Per-sample abundance for top DAA candidates (%s).", domain))
    figures[[length(figures) + 1]] <- figure_entry(
      cfg, fig_path(cfg, sprintf("differential_abundance/%s/aldex_plots.pdf", domain)),
      kind = "aldex_plot", domain = domain,
      caption_seed = sprintf("ALDEx2 MA / effect plots (%s).", domain))
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
      function(rel) figure_entry(
        cfg, rel,
        kind = paste0("ge_", tools::file_path_sans_ext(basename(rel))),
        domain = domain,
        caption_seed = sprintf("%s — %s.", domain,
                               tools::file_path_sans_ext(basename(rel)))
      )
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
                  description = "Per-sample cluster / module assignment.")
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
        figure_entry(cfg, fig_path(cfg, "network/degree_distribution.png"),
                     kind = "degree_distribution",
                     caption_seed = "Network degree distribution.")
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

  meta_path <- file.path(cfg$project_root, cfg$metadata$file)
  meta <- tryCatch(
    readr::read_csv(meta_path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  sid_col  <- cfg$metadata$sample_id_col
  controls <- cfg$metadata$controls %||% character(0)
  n_samples <- if (!is.null(meta)) sum(!meta[[sid_col]] %in% controls) else NA_integer_

  manifest <- list(
    manifest_version = "1.0",
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
    stages = Filter(Negate(is.null), list(
      taxonomy               = build_stage_taxonomy(cfg, stage_times),
      relative_abundance     = build_stage_relative_abundance(cfg, stage_times),
      alpha_diversity        = build_stage_alpha_diversity(cfg, stage_times),
      beta_diversity         = build_stage_beta_diversity(cfg, stage_times),
      differential_abundance = build_stage_differential_abundance(cfg, stage_times),
      resistome              = build_stage_ge_domain(cfg, "resistome", stage_times),
      virulome               = build_stage_ge_domain(cfg, "virulome",  stage_times),
      mobilome               = build_stage_ge_domain(cfg, "mobilome",  stage_times),
      network                = build_stage_network(cfg, stage_times)
    )),
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
