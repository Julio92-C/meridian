# The chord diagram must render when sectors are many (137 samples + taxa +
# genes + classes): a fixed 2-degree gap failed in circlize with
# "gap.degree is too large ... no space to allocate sectors".
suppressPackageStartupMessages({ library(dplyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
pipeline_log <- function(cfg, msg) invisible(NULL)
source(file.path("..", "..", "R", "12_network.R"), local = TRUE)

test_that(".chord_gap_degree keeps 2 degrees for small diagrams and shrinks for large ones", {
  expect_equal(.chord_gap_degree(20), 2)
  expect_lt(.chord_gap_degree(250) * 250, 360)
})

test_that(".network_render_chord draws 250 sectors without error", {
  skip_if_not_installed("circlize")
  set.seed(1)
  df <- tibble::tibble(
    sample = sprintf("S%03d", rep(1:150, each = 2)),
    name   = sprintf("taxon %02d", sample.int(60, 300, replace = TRUE)),
    GENE   = sprintf("gene%02d", sample.int(30, 300, replace = TRUE)),
    RESISTANCE_class = sprintf("class%d", sample.int(10, 300, replace = TRUE))
  )
  png <- tempfile(fileext = ".png")
  expect_no_error(.network_render_chord(df, png, "ggsci::default_ucscgb", 152))
  expect_true(file.exists(png))
})
