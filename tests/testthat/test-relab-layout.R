# R/05 layout rules (per-sample grid, group_order, main_panel source) and the
# single-row heatmap scaling in utils_ge_profile.R.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("..", "..", "R", "05_relative_abundance.R"))
source(file.path("..", "..", "R", "utils_ge_profile.R"))

cfgp <- function(...) list(relative_abundance = list(per_sample = list(...)))

test_that("auto layout follows the number of groups", {
  expect_equal(ra_auto_layout(1, 10, list()), list(nrow = 1L, ncol = 2L))
  expect_equal(ra_auto_layout(2, 10, list()), list(nrow = 1L, ncol = 2L))
  expect_equal(ra_auto_layout(3, 10, list()), list(nrow = 2L, ncol = 2L))
  expect_equal(ra_auto_layout(4, 10, list()), list(nrow = 2L, ncol = 2L))
  expect_equal(ra_auto_layout(5, 10, list()), list(nrow = 3L, ncol = 2L))
  expect_equal(ra_auto_layout(6, 10, list()), list(nrow = 3L, ncol = 2L))
  expect_null(ra_auto_layout(7, 10, list()))
})

test_that("max_samples boundary (default 40, configurable)", {
  expect_false(is.null(ra_auto_layout(2, 40, list())))
  expect_null(ra_auto_layout(2, 41, list()))
  expect_false(is.null(ra_auto_layout(2, 60, cfgp(max_samples = 60))))
  expect_null(ra_auto_layout(2, 61, cfgp(max_samples = 60)))
  plan <- ra_per_sample_plan(2, 41, list())
  expect_match(plan$reason, "max_samples")
  expect_match(ra_per_sample_plan(7, 5, list())$reason, "7 groups")
})

test_that("explicit layout, enabled true / false", {
  expect_equal(ra_auto_layout(2, 10, cfgp(layout = "2x2")),
               list(nrow = 2L, ncol = 2L))
  expect_equal(ra_auto_layout(5, 10, cfgp(layout = "3X2")),
               list(nrow = 3L, ncol = 2L))
  # a grid too small for the facets is grown, not truncated
  expect_equal(ra_auto_layout(5, 10, cfgp(layout = "1x2"))$nrow, 3L)
  # false never draws; true forces past both limits
  expect_null(ra_auto_layout(2, 10, cfgp(enabled = FALSE)))
  expect_false(is.null(ra_auto_layout(8, 500, cfgp(enabled = TRUE))))
  expect_equal(ra_auto_layout(8, 500, cfgp(enabled = TRUE))$ncol, 2L)
  # YAML strings are accepted too
  expect_null(ra_auto_layout(2, 10, cfgp(enabled = "false")))
  # the block is also accepted without the relative_abundance wrapper
  expect_null(ra_auto_layout(2, 10,
                             list(per_sample = list(enabled = FALSE))))
})

test_that("group_order is honoured, unknown levels ignored, rest alphabetical", {
  g <- c("B", "A", "C", "A")
  expect_equal(ra_group_levels(g), c("A", "B", "C"))
  expect_equal(ra_group_levels(g, c("C", "A", "B")), c("C", "A", "B"))
  expect_equal(ra_group_levels(g, c("C", "Z")), c("C", "A", "B"))
})

test_that("samples are ordered by dominant taxon, then its %", {
  mk <- function(s, a, b, o = 0) {
    data.frame(sample = s, name = c("T1", "T2", "Others"),
               percentage = c(a, b, o))
  }
  df <- dplyr::bind_rows(mk("s1", 10, 60, 30), mk("s2", 70, 20, 10),
                         mk("s3", 50, 40, 10), mk("s4", 5, 5, 90))
  df$name <- factor(df$name, levels = c("T1", "T2", "Others"))
  # T1-dominant first (s2 70 > s3 50), then T2-dominant (s1); s4 is mostly
  # Others, so its dominant *taxon* is T1 (5, tie with T2 -> legend order)
  expect_equal(ra_order_samples(df), c("s2", "s3", "s4", "s1"))
})

test_that("main_panel source: wanted first, then group -> stream -> sample", {
  all3 <- c("sample", "stream", "group")
  cf <- function(x) list(relative_abundance = list(main_panel = x))
  expect_equal(ra_main_panel_source(list(), all3), "group")           # default
  expect_equal(ra_main_panel_source(cf("stream"), all3), "stream")
  expect_equal(ra_main_panel_source(cf("sample"), all3), "sample")
  expect_equal(ra_main_panel_source(cf("sample"), c("group", "stream")), "group")
  expect_equal(ra_main_panel_source(cf("sample"), "stream"), "stream")
  expect_equal(ra_main_panel_source(cf("stream"), c("sample", "group")), "group")
  expect_equal(ra_main_panel_source(cf("group"), "sample"), "sample")
  expect_null(ra_main_panel_source(cf("group"), character(0)))
  expect_equal(ra_main_panel_source(cf("bogus"), all3), "group")      # invalid
})

test_that("split level file names are safe and never clobber built-ins", {
  expect_equal(ra_safe_level("Lung A/B"), "Lung_A_B")
  expect_equal(ra_safe_level("genus"), "genus_level")
  expect_equal(ra_safe_level("stream_x"), "stream_x_level")
  expect_equal(ra_safe_level("Batch1"), "Batch1")
})

test_that("per-sample plot: Others is the last legend level, single column, strips carry n", {
  pal <- c(T1 = "#E41A1C", T2 = "#377EB8", Others = "#BDBDBD")
  mk <- function(s, g, a) {
    data.frame(sample = s, grp = g, name = c("T1", "T2", "Others"),
               percentage = c(a, 100 - a - 10, 10))
  }
  df <- dplyr::bind_rows(mk("a1", "A", 60), mk("a2", "A", 20),
                         mk("b1", "B", 40))
  b <- ra_build_sample_plot(df, "grp", pal, list(nrow = 1L, ncol = 2L), list(),
                            group_order = c("B", "A"))
  built <- ggplot2::ggplot_build(b$plot)
  expect_equal(levels(b$plot$data$name), c("T1", "T2", "Others"))
  lim <- built$plot$scales$get_scales("fill")$get_limits()
  expect_equal(tail(lim, 1), "Others")
  expect_equal(levels(b$plot$data$.facet), c("B (n=1)", "A (n=2)"))
  expect_null(b$plot$labels$x)                       # no "Sample" x title
  g <- ggplot2::ggplotGrob(b$plot)
  expect_true(any(grepl("guide-box", g$layout$name)))
  expect_equal(b$plot$theme$legend.position, "right")
  expect_s3_class(b$plot$facet, "FacetGrid")          # single row
  b2 <- ra_build_sample_plot(df, "grp", pal, list(nrow = 2L, ncol = 2L), list())
  expect_s3_class(b2$plot$facet, "FacetWrap")
  expect_no_error(ggplot2::ggplotGrob(b$plot))
})

test_that("single-gene heatmap matrix is scaled across samples, not all zero", {
  m <- matrix(c(2, 10, 6), nrow = 1, dimnames = list("g1", c("s1", "s2", "s3")))
  r <- ge_scale_for_heatmap(m, "value")
  expect_equal(unname(r$mat[1, ]), c(0, 1, 0.5))
  expect_equal(dimnames(r$mat), dimnames(m))
  expect_equal(r$legend_title, "Scaled TPM (0-1, across samples)")
  # constant single row keeps the old (all zero) behaviour and title
  k <- matrix(c(4, 4, 4), nrow = 1, dimnames = list("g1", c("s1", "s2", "s3")))
  rk <- ge_scale_for_heatmap(k, "value")
  expect_equal(unname(rk$mat[1, ]), c(0, 0, 0))
  expect_equal(rk$legend_title, "value")
  # multi-row matrices are untouched (per-column min-max)
  m2 <- matrix(c(1, 2, 3, 4, 4, 4), ncol = 2,
               dimnames = list(letters[1:3], c("x", "y")))
  r2 <- ge_scale_for_heatmap(m2, "value")
  expect_equal(r2$mat, ge_minmax_per_col(m2))
  expect_equal(r2$legend_title, "value")
})

test_that("category and gene heatmap helpers render a single row without error", {
  skip_if_not_installed("pheatmap")
  df <- data.frame(
    GENE = "tet(M)", Class = "Tetracycline",
    sample = c("s1", "s2", "s3"), TPM = c(1, 5, 3),
    grp = c("A", "A", "B"), stringsAsFactors = FALSE)
  f1 <- tempfile(fileext = ".png"); f2 <- tempfile(fileext = ".png")
  pal_g <- c(A = "#E41A1C", B = "#377EB8")
  expect_no_error(ge_plot_category_heatmap(df, "Class", "grp", pal_g, f1,
                                           palette = c("white", "red")))
  expect_no_error(ge_plot_gene_heatmap(df, "grp", pal_g, c(Tetracycline = "#999999"),
                                       "Class", f2, palette = c("white", "red")))
  expect_true(file.exists(f1)); expect_true(file.exists(f2))
})
