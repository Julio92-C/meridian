# R/05 composition helpers: top-N ranks by MEAN per-sample %, group means sum to
# 100, and "Others" is always light grey (never a taxon colour).
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("..", "..", "R", "05_relative_abundance.R"))

# deep: 10000 reads; shallow: 100 reads. "Deep" has a huge raw count of D1 but
# a low share; "S1" is dominant in the shallow samples only.
toy <- function() {
  mk <- function(sample, group, counts) {
    data.frame(sample = sample, group = group, name = names(counts),
               count = as.numeric(counts), stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(
    mk("deep1",    "A", c(D1 = 3000, X = 3500, Y = 3500)),
    mk("shallow1", "A", c(S1 = 60,   X = 20,   Y = 20)),
    mk("shallow2", "B", c(S1 = 60,   X = 20,   Y = 20)),
    mk("shallow3", "B", c(S1 = 60,   X = 20,   Y = 20))
  )
}

test_that("top-N is ranked by mean per-sample %, not total raw counts", {
  d <- toy()
  # By raw counts D1 (3000) beats S1 (180); by mean % S1 (~45%) beats D1 (7.5%).
  by_counts <- d |> group_by(name) |> summarise(t = sum(count)) |>
    arrange(desc(t)) |> head(3) |> pull(name)
  expect_true("D1" %in% by_counts)
  prep <- ra_prepare_composition(d, "group", list(top_n = 3))
  expect_true("S1" %in% prep$top_taxa)
  expect_false("D1" %in% prep$top_taxa)
  expect_equal(length(prep$top_taxa), 3)
  expect_true("Others" %in% prep$taxa_levels)
  expect_equal(tail(prep$taxa_levels, 1), "Others")
})

test_that("default top_n is 15 and ra_cfg$top_n overrides it", {
  d <- data.frame(sample = "s1", group = "A",
                  name = sprintf("T%02d", 1:30), count = 30:1)
  expect_equal(length(ra_prepare_composition(d, "group", list())$top_taxa), 15)
  expect_equal(length(ra_prepare_composition(d, "group", list(top_n = 5))$top_taxa), 5)
})

test_that("group means sum to 100 per group and carry n", {
  prep <- ra_prepare_composition(toy(), "group", list(top_n = 2))
  gm <- ra_group_means(prep$df_pct, "group")
  totals <- gm |> group_by(group) |> summarise(s = sum(mean_percentage))
  expect_equal(totals$s, c(100, 100), tolerance = 1e-8)
  expect_equal(sort(unique(gm$label)), c("A (n=2)", "B (n=2)"))
})

test_that("Others is always #BDBDBD and taxa never use it", {
  cfg <- list(relative_abundance = list(palette = "ggsci::default_igv"))
  lv  <- c(sprintf("T%02d", 1:60), "Others")
  pal <- ra_taxa_palette(lv, "Others", cfg)
  expect_equal(unname(pal["Others"]), "#BDBDBD")
  expect_false(any(toupper(pal[setdiff(lv, "Others")]) == "#BDBDBD"))
  expect_equal(names(pal), lv)
  # custom others_label honoured
  pal2 <- ra_taxa_palette(c("a", "b", "Rest"), "Rest", cfg)
  expect_equal(unname(pal2["Rest"]), "#BDBDBD")
})
