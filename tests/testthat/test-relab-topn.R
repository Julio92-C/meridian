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

# ---- species-rank filter: unknown ranks must never make taxa vanish --------
if (!exists("pipeline_log", envir = globalenv()))
  assign("pipeline_log", function(cfg, msg) invisible(NULL), envir = globalenv())

test_that("all-NA tax_rank (blank lvl_type) keeps every row", {
  d <- data.frame(name = c("a", "b", "c"), count = 1:3, tax_rank = NA_character_)
  expect_equal(nrow(ra_filter_species_rank(d, list())), 3)
  d$tax_rank <- ""
  expect_equal(nrow(ra_filter_species_rank(d, list())), 3)
})

test_that("mixed S / S1 / NA keeps S and NA, drops S1", {
  d <- data.frame(name = c("sp", "strain", "unk", "genus"), count = 1:4,
                  tax_rank = c("S", "S1", NA, "G"))
  out <- ra_filter_species_rank(d, list())
  expect_equal(out$name, c("sp", "unk"))
})

test_that("a taxid missing from kraken2 (NA rank) is kept after the R/02 join", {
  # mimic R/02: left join of rcf rows onto a kraken2 rank lookup
  rcf <- data.frame(taxid = c(1, 2, 3), name = c("a", "b", "c"), count = 1:3)
  k2  <- data.frame(taxid = c(1, 2), tax_rank = c("S", "S1"))
  d <- dplyr::left_join(rcf, k2, by = "taxid")
  expect_true(is.na(d$tax_rank[d$taxid == 3]))
  out <- ra_filter_species_rank(d, list())
  expect_equal(out$name, c("a", "c"))
})

test_that("missing tax_rank column leaves the table untouched", {
  d <- data.frame(name = "a", count = 1)
  expect_equal(ra_filter_species_rank(d, list()), d)
})
