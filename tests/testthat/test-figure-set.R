# Tests for the figure set (cfg$figures$set: core | full). `core` skips the
# duplicate of each pair (bar vs violin, total count vs count + prevalence,
# stacked composition bar vs stream, Venn vs UpSet); `full` writes everything. The skipped figures' code is kept.
`%||%` <- function(a, b) if (is.null(a)) b else a
R_DIR <- file.path("..", "..", "R")
source(file.path(R_DIR, "utils_figures.R"))
source(file.path(R_DIR, "utils_panel_io.R"))

core_cut <- c(
  "alpha_diversity/richness_bar.png", "alpha_diversity/chao1_bar.png",
  "resistome/alpha_shannon_bar.png", "mobilome/alpha_richness_bar.png",
  "resistome/drug_total_count.png", "virulome/vf_total_count.png",
  "mobilome/mge_total_count.png", "mobilome/mge_family_total_count.png",
  "relative_abundance/species_count.png",
  "relative_abundance/relative_abundance.png", "relative_abundance/relative_abundance_genus.png",
  "relative_abundance/relative_abundance_class.png",
  "resistome/drug_relative_abundance.png", "virulome/vf_relative_abundance.png",
  "mobilome/mge_relative_abundance.png", "mobilome/mge_family_abundance.png",
  "resistome/venn_args.png", "virulome/venn_vfs.png", "mobilome/venn_replicons.png",
  # second cut (lung contact sheet)
  "alpha_diversity/chao1_violin.png", "beta_diversity/jaccard_pcoa.png",
  "differential_abundance/gene/volcano_A_vs_B.png",
  "differential_abundance/taxa/aldex2_maplot_A_vs_B.png",
  "differential_abundance/taxa/top30_daa_heatmap.png",
  "mobilome/mge_family_stream_abundance.png", "network/mantel_correlation_triangle.png",
  "network/sankey_taxon_arg_mge.png", "relative_abundance/relative_abundance_stream_order.png",
  "relative_abundance/unique_species_relative_abundance.png",
  "resistome/c_difficile_count_per_treatment.png", "resistome/enterobacteriaceae_gene_map.png",
  "virulome/gene_count_prevalence.png", "resistome/venn_drug_classes.png"
)
core_kept <- c(
  "alpha_diversity/richness_violin.png", "resistome/alpha_shannon_violin.png",
  "resistome/drug_count_prevalence.png", "relative_abundance/species_count_prevalence.png",
  "relative_abundance/relative_abundance_stream.png",
  "relative_abundance/relative_abundance_stream_genus.png",
  "relative_abundance/relative_abundance_group.png",
  "resistome/drug_stream_abundance.png", "mobilome/mge_stream_abundance.png",
  "resistome/arg_upset_treatments.png",
  "taxonomy/taxa_venn.png",
  # composite sources and figures kept after the second cut
  "resistome/pheatmap_genes.png", "network/sankey_overall.png",
  "differential_abundance/taxa/volcano_A_vs_B.png", "network/connectivity_venn_taxa.png",
  "relative_abundance/relative_abundance_stream_class.png", "alpha_diversity/shannon_violin.png",
  "resistome/drug_count_prevalence.png", "virulome/vf_upset_treatments.png"
)

test_that("core skips one figure of each duplicate pair and keeps the other", {
  set_figure_set(list(figures = list(set = "core")))
  for (f in core_cut)  expect_false(figure_wanted(file.path("C:/x/Figures", f)), info = f)
  for (f in core_kept) expect_true(figure_wanted(file.path("C:/x/Figures", f)), info = f)
})

test_that("core is the default and full writes everything", {
  set_figure_set(list())
  expect_false(figure_wanted("Figures/alpha_diversity/richness_bar.png"))
  set_figure_set(list(figures = list(set = "full")))
  for (f in core_cut) expect_true(figure_wanted(file.path("Figures", f)), info = f)
})

test_that("figures.include brings a cut figure back; figures.exclude drops a kept one", {
  set_figure_set(list(figures = list(set = "core", include = list("venn_args"),
                                     exclude = list("pcoa"))))
  expect_true(figure_wanted("Figures/resistome/venn_args.png"))
  expect_false(figure_wanted("Figures/resistome/pcoa.png"))
})

test_that("an unknown figure set is an error", {
  expect_error(set_figure_set(list(figures = list(set = "some"))), "core.*full")
})

test_that("save_panel_ggplot writes nothing for a skipped figure", {
  set_figure_set(list(figures = list(set = "core")))
  d <- file.path(tempfile(), "alpha_diversity"); dir.create(d, recursive = TRUE)
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
  save_panel_ggplot(file.path(d, "richness_bar.png"), p, width = 2, height = 2)
  save_panel_ggplot(file.path(d, "richness_violin.png"), p, width = 2, height = 2)
  expect_equal(sort(list.files(d)), c("richness_violin.png", "richness_violin.rds"))
})

# Leave no figure set behind for the other test files.
options(meridian.figures = NULL)
