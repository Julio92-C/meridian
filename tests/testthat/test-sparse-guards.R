# Regression tests for the sparse / degenerate-input hardening (low-biomass
# data): each test feeds a helper or stage function a tiny degenerate input
# that used to abort the whole pipeline and asserts that it now skips (with a
# log line) or renders instead.
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
R_DIR <- file.path("..", "..", "R")
source(file.path(R_DIR, "00_setup.R"))          # padjust_*, pipeline_log
source(file.path(R_DIR, "utils_palette.R"))
source(file.path(R_DIR, "utils_panel_io.R"))
source(file.path(R_DIR, "utils_taxa.R"))
source(file.path(R_DIR, "utils_ge_profile.R"))
source(file.path(R_DIR, "utils_cluster_perm.R"))   # cluster_tests_enabled() etc.
source(file.path(R_DIR, "04_taxonomy.R"))
source(file.path(R_DIR, "07_beta_diversity.R"))
source(file.path(R_DIR, "08_differential_abundance.R"))

logged <- character(0)
pipeline_log <<- function(cfg, msg) logged <<- c(logged, msg)
reset_log <- function() logged <<- character(0)

# ---- tiny project fixture ------------------------------------------------
make_proj <- function(meta) {
  root <- tempfile("sparse_proj_"); dir.create(root)
  readr::write_csv(meta, file.path(root, "meta.csv"))
  list(
    project_root = root,
    metadata = list(file = "meta.csv", sample_id_col = "sample",
                    group_cols = list("grp")),
    outputs  = list(figures_dir = "figs", datasets_dir = "ds",
                    log_file = "log.txt"),
    stats    = list(permanova_permutations = 99),
    beta_diversity = list(min_count_per_sample = 1)
  )
}

# ---- utils_ge_profile ------------------------------------------------------
test_that("ge_compute_alpha gives each sample its own Shannon (n x 1 matrix)", {
  df <- data.frame(sample = c("s1", "s2", "s3"), GENE = "g1",
                   TPM = c(5, 50, 500), grp = c("A", "A", "B"))
  a <- ge_compute_alpha(df, "grp")
  expect_equal(a$shannon, c(0, 0, 0))     # one gene -> H = 0 for every sample
})

test_that("gene/category heatmaps render a constant (single-gene) matrix", {
  skip_if_not_installed("pheatmap")
  df <- data.frame(sample = c("s1", "s2", "s3"), GENE = "g1", CAT = "c1",
                   TPM = c(5, 5, 5), grp = c("A", "A", "B"))
  pal_g <- c(A = "red", B = "blue"); pal_c <- c(c1 = "green")
  f1 <- tempfile(fileext = ".png"); f2 <- tempfile(fileext = ".png")
  f3 <- tempfile(fileext = ".png")
  expect_no_error(ge_plot_gene_heatmap(df, "grp", pal_g, pal_c, "CAT", f1,
                                       palette = c("white", "red")))
  expect_true(file.exists(f1))
  expect_no_error(ge_plot_category_heatmap(df, "CAT", "grp", pal_g, f2,
                                           palette = c("white", "red")))
  expect_true(file.exists(f2))
  expect_no_error(ge_plot_category_heatmap(df, "CAT", "grp", pal_g, f3,
                                           palette = c("white", "red"),
                                           legend_top = FALSE))
})

test_that("ge_ensure_group_column drops NA-group rows and logs", {
  reset_log()
  df <- data.frame(sample = c("s1", "s2", "s3"), GENE = "g", grp = c("A", NA, "B"))
  out <- ge_ensure_group_column(df, list(), "grp", "Resistome")
  expect_equal(out$sample, c("s1", "s3"))
  expect_true(any(grepl("dropping 1 row.*NA grp", logged)))
})

beta_args <- function(df, tmp = tempfile()) {
  dir.create(tmp)
  list(df = df, group = "grp", pal_group = c(A = "red", B = "blue"),
       opts = list(), fig_file = file.path(tmp, "pcoa.png"),
       ds_perm_file = file.path(tmp, "perm.txt"),
       ds_disp_file = file.path(tmp, "disp.txt"),
       ds_pcoa_file = file.path(tmp, "scores.csv"),
       log_label = "Resistome", title_prefix = "ARG",
       skip_data_label = "ARG hits",
       cfg = list(stats = list(permanova_permutations = 99)))
}

test_that("ge_plot_beta skips (not subscript-out-of-bounds) on < 2 PCoA axes", {
  reset_log()
  df <- data.frame(sample = paste0("s", 1:4), GENE = "g1",
                   TPM = c(5, 5, 5, 5), grp = c("A", "A", "B", "B"))
  a <- beta_args(df)
  r <- NULL
  expect_no_error(r <- do.call(ge_plot_beta, a))
  expect_null(r)
  expect_true(any(grepl("< 2 PCoA axes", logged)))
  expect_false(file.exists(a$fig_file))
})

test_that("ge_plot_beta drops NA-group samples instead of failing in adonis2", {
  reset_log()
  df <- data.frame(
    sample = rep(paste0("s", 1:5), each = 3),
    GENE   = rep(c("g1", "g2", "g3"), 5),
    TPM    = c(9,1,2, 8,2,1, 1,9,3, 2,8,1, 5,5,5),
    grp    = rep(c("A", "A", "B", "B", NA), each = 3)
  )
  a <- beta_args(df)
  r <- NULL
  expect_no_error(r <- do.call(ge_plot_beta, a))
  expect_true(any(grepl("dropping 1 sample.*NA grp", logged)))
  expect_false(any(is.na(r$scores$grp)))
  expect_equal(nrow(r$scores), 4L)
})

# ---- save_widget_safe ----------------------------------------------------
test_that("save_widget_safe falls back to selfcontained = FALSE", {
  skip_if_not_installed("htmlwidgets")
  reset_log()
  calls <- c()
  testthat::local_mocked_bindings(
    saveWidget = function(widget, file, selfcontained = TRUE, ...) {
      calls <<- c(calls, selfcontained)
      if (selfcontained) stop("pandoc is required")
      invisible(NULL)
    }, .package = "htmlwidgets")
  expect_no_error(ok <- save_widget_safe(list(), tempfile(fileext = ".html"),
                                         selfcontained = TRUE, cfg = list()))
  expect_true(ok)
  expect_equal(calls, c(TRUE, FALSE))
  expect_true(any(grepl("retrying with selfcontained = FALSE", logged)))
})

test_that("save_widget_safe logs and returns FALSE when both attempts fail", {
  skip_if_not_installed("htmlwidgets")
  reset_log()
  testthat::local_mocked_bindings(
    saveWidget = function(...) stop("boom"), .package = "htmlwidgets")
  expect_no_error(ok <- save_widget_safe(list(), tempfile(fileext = ".html"),
                                         cfg = list()))
  expect_false(ok)
  expect_true(any(grepl("skipping html", logged)))
})

# ---- R/07 beta diversity -----------------------------------------------------
# R/07 ordinates the community taxa table (cleaned$noncontaminants) by
# default since 12b191c, so the fixtures feed that table (sample, name, count).
beta_community <- function(samples, taxa, counts) {
  expand.grid(sample = samples, name = taxa, stringsAsFactors = FALSE) |>
    dplyr::mutate(count = counts)
}

test_that("run_beta_diversity skips when Bray PCoA has < 2 axes", {
  reset_log()
  s <- paste0("s", 1:4)
  cfg <- make_proj(data.frame(sample = s, grp = c("A", "A", "B", "B")))
  cleaned <- list(noncontaminants = beta_community(s, "Escherichia coli", rep(100, 4)))
  r <- NULL
  expect_no_error(r <- run_beta_diversity(cleaned, cfg))
  expect_null(r)
  expect_true(any(grepl("< 2 PCoA axes", logged)))
})

test_that("run_beta_diversity drops NA-group samples and skips only Jaccard when degenerate", {
  skip_if_not_installed("plotly")
  reset_log()
  s <- paste0("s", 1:6)
  taxa <- c("Escherichia coli", "Bacillus subtilis", "Listeria innocua")
  # same taxa in every sample (identical presence/absence -> Jaccard d = 0)
  # but varied abundances (Bray stays informative)
  cnt <- c(90, 80, 10, 20, 50, 40,   5, 10, 90, 70, 40, 30,   8, 12, 20, 60, 30, 70)
  cfg <- make_proj(data.frame(sample = s,
                              grp = c("A", "A", "A", "B", "B", NA)))
  cleaned <- list(noncontaminants = beta_community(s, taxa, cnt))
  r <- NULL
  expect_no_error(r <- run_beta_diversity(cleaned, cfg))
  expect_true(any(grepl("dropping 1 sample.*NA grp", logged)))
  expect_true(any(grepl("Jaccard distances give < 2 PCoA axes", logged)))
  figs <- file.path(cfg$project_root, "figs", "beta_diversity")
  expect_true(file.exists(file.path(figs, "pcoa.png")))
  expect_false(file.exists(file.path(figs, "jaccard_pcoa.png")))
  expect_equal(nrow(r$scores), 5L)
})

# ---- R/04 taxonomy -------------------------------------------------------------
test_that("run_taxonomy skips a constant heatmap and drops NA-group samples", {
  skip_if_not_installed("pheatmap")
  reset_log()
  s <- paste0("s", 1:4)
  cfg <- make_proj(data.frame(sample = s, grp = c("A", "A", "B", NA)))
  cfg$taxonomy <- list(min_count_per_sample = 5)
  cleaned <- list(noncontaminants = data.frame(
    sample = s, name = "Escherichia coli", count = 100, stringsAsFactors = FALSE))
  expect_no_error(run_taxonomy(cleaned, cfg))
  expect_true(any(grepl("dropping 1 row.*NA grp", logged)))
  expect_true(any(grepl("heatmap matrix is constant", logged)))
  expect_false(file.exists(file.path(cfg$project_root, "figs", "taxonomy",
                                     "taxa_heatmap.png")))
})

# ---- R/08 differential abundance: single-level design ------------------------
test_that("da_run_glm skips a design whose factor has a single level", {
  reset_log()
  meta_use <- data.frame(sample = paste0("s", 1:4), grp = "A")
  mat <- matrix(c(10, 5, 8, 9, 3, 4, 6, 2), nrow = 2,
                dimnames = list(c("f1", "f2"), meta_use$sample))
  r <- NULL
  expect_no_error(r <- da_run_glm("taxa", mat, meta_use, "sample", "~ grp",
                                  mc = 2, denom = "all", alpha = 0.05,
                                  ds_dir = tempdir(), feature_col = "name",
                                  cfg = list()))
  expect_null(r)
  expect_true(any(grepl("have < 2 levels", logged)))
})
