# utils_figures.R — which figures a run writes (cfg$figures).
#
#   figures:
#     set: core            # core (default) | full
#     include: []          # file stems to write even when the set skips them, e.g. [venn_args]
#     exclude: []          # file stems to skip even when the set writes them
#
# `core` skips one figure of each duplicate pair (same data drawn twice; agreed
# with Julio 2026-10-27): per-sample bars (violins kept), total-count bars
# (count + prevalence kept), stacked composition bars of taxa and genes (stream
# kept) and Venns that have an UpSet twin. `full` writes every
# figure, as before. Skipped figures are not deleted from the code: they come
# back with `set: full` or `include:`.
#
# The set is read once by load_config() (set_figure_set) and checked where the
# files are written (save_panel_ggplot, ge_plot_category_venn), so the stage
# modules need no per-figure switches.

# Regexes on "<stage dir>/<file name>".
FIGURES_FULL_ONLY <- c(
  "^alpha_diversity/[a-z0-9]+_bar\\.png$",
  "^(resistome|virulome|mobilome)/alpha_[a-z0-9]+_bar\\.png$",
  "^(resistome|virulome|mobilome)/(drug|vf|mge|mge_family)_total_count\\.png$",
  "^relative_abundance/species_count\\.png$",
  "^relative_abundance/relative_abundance(_(phylum|class|order|family|genus|species))?\\.png$",
  "^(resistome|virulome|mobilome)/((drug|vf|mge)_relative_abundance|mge_family_abundance)\\.png$",
  "^(resistome|virulome|mobilome)/venn_(args|vfs|replicons)\\.png$",
  # Second cut (Julio 2026-10-29, lung contact sheet): figures nobody reads
  # in the standard set. None of them is the source of a composite slot.
  "^alpha_diversity/(chao1|simpson)_violin\\.png$",
  "^beta_diversity/jaccard_pcoa\\.png$",
  "^gene/volcano_.+\\.png$",
  "^taxa/(aldex2_dotplot_summary|aldex2_maplot_.+|top30_daa_heatmap)\\.png$",
  "^mobilome/mge_family_(count_prevalence|stream_abundance)\\.png$",
  "^network/(connectivity_venn_genesets|mantel_correlation_triangle|mobile_arg_fraction_bar|sankey_taxon_arg_mge)\\.png$",
  "^relative_abundance/(genus_heatmap_top30|relative_abundance_stream_(family|order)|(shared|unique)_species_relative_abundance)\\.png$",
  "^resistome/[a-z0-9_]+_(count_per_treatment|gene_map)\\.png$",
  "^(resistome|virulome)/gene_count_prevalence\\.png$",
  "^resistome/venn_drug_classes\\.png$",
  "^virulome/venn_vf_functions\\.png$"
)

set_figure_set <- function(cfg) {
  fc  <- cfg$figures %||% list()
  set <- fc$set %||% "core"
  if (!set %in% c("core", "full")) {
    stop(sprintf("figures.set must be 'core' or 'full', not '%s'", set))
  }
  options(meridian.figures = list(set = set,
                                  include = unlist(fc$include),
                                  exclude = unlist(fc$exclude)))
  invisible(set)
}

# TRUE when the figure at `path` should be written. With no set configured
# (helpers used outside a pipeline run, e.g. in tests) everything is written.
figure_wanted <- function(path) {
  fs <- getOption("meridian.figures")
  if (is.null(fs)) return(TRUE)
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
  key   <- paste(utils::tail(parts, 2), collapse = "/")
  stem  <- tools::file_path_sans_ext(basename(path))
  if (stem %in% fs$exclude) return(FALSE)
  if (stem %in% fs$include) return(TRUE)
  if (identical(fs$set, "full")) return(TRUE)
  !any(vapply(FIGURES_FULL_ONLY, grepl, logical(1), x = key))
}
