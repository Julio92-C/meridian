# sf_prepare_inputs() must work at species rank even when taxa_filter.rank is
# not configured: strain / subspecies rows (kraken2 S1, S2, ...) would otherwise
# be counted as separate "species" next to their own species, and collapse onto
# the same display label in the widespread figure ("factor level is duplicated").
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
pipeline_log <- function(cfg, msg) invisible(NULL)   # stubs live here: sources below use local = TRUE
source(file.path("..", "..", "R", "utils_taxa.R"), local = TRUE)
source(file.path("..", "..", "R", "utils_shared_fraction.R"), local = TRUE)
source(file.path("..", "..", "R", "05_relative_abundance.R"), local = TRUE)
source(file.path("..", "..", "R", "04_shared_fraction.R"), local = TRUE)

test_that("sf_prepare_inputs keeps species rank only when taxa_filter.rank is unset", {
  root <- tempfile("sfrank"); dir.create(root)
  readr::write_csv(tibble::tibble(sample = c("S1", "S2"), grp = c("A", "B")),
                   file.path(root, "meta.csv"))
  nc <- tibble::tibble(
    sample   = rep(c("S1", "S2"), each = 3),
    name     = rep(c("Bifidobacterium pullorum",
                     "Bifidobacterium pullorum subsp. gallinarum",
                     "Alistipes finegoldii"), 2),
    tax_rank = rep(c("S", "S1", "S"), 2),
    count    = 100L
  )
  cfg <- list(project_root = root,
              metadata = list(file = "meta.csv", sample_id_col = "sample", group_cols = list("grp")))
  inp <- sf_prepare_inputs(list(noncontaminants = nc), cfg)
  expect_setequal(unique(inp$long$name), c("Bifidobacterium pullorum", "Alistipes finegoldii"))
  # display labels of what remains are unique
  expect_false(anyDuplicated(sf_short(unique(inp$long$name))) > 0)
})
