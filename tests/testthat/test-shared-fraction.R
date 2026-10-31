# Unit tests for the shared-fraction helpers (R/utils_shared_fraction.R) on a
# synthetic three-group table with known membership patterns.
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("..", "..", "R", "utils_shared_fraction.R"))

# Groups A (10 samples), B (10), C (5). Species:
#   wide    in 6 A, 5 B, 3 C  -> shared:all, widespread (>=3 samples, >=10%)
#   rec     in 2 A, 2 B, 2 C  -> shared:all, recurrent (2 samples, but C 40%: fails n>=3)
#   spor    in 1 A, 8 B       -> shared:A+B, sporadic
#   onlyC   in 4 C            -> unique:C
#   ab      in 3 A, 3 B       -> shared:A+B, widespread on its own pattern
mk <- function(name, smp) tibble::tibble(sample = smp, name = name, count = 50L)
A <- paste0("A", 1:10); B <- paste0("B", 1:10); C <- paste0("C", 1:5)
long <- dplyr::bind_rows(
  mk("wide", c(A[1:6], B[1:5], C[1:3])), mk("rec", c(A[1:2], B[1:2], C[1:2])),
  mk("spor", c(A[1], B[1:8])), mk("onlyC", C[1:4]), mk("ab", c(A[1:3], B[1:3])),
  mk("filler", c(A, B, C))
)
long$group <- substr(long$sample, 1, 1)
long$phylum <- "Pseudomonadota"
group_n <- c(A = 10L, B = 10L, C = 5L)
tab <- sf_species_table(long, group_n, min_samples = 3, min_prev = 0.10)
row <- function(n) tab[tab$name == n, ]

test_that("categories follow the membership pattern", {
  expect_equal(row("wide")$category, "shared:all")
  expect_equal(row("spor")$category, "shared:A+B")
  expect_equal(row("onlyC")$category, "unique:C")
  expect_equal(row("ab")$pattern, "A+B")
})

test_that("tiers use both the sample and the prevalence rule over the pattern", {
  expect_equal(row("wide")$tier, "widespread")
  expect_equal(row("rec")$tier, "recurrent")
  expect_equal(row("spor")$tier, "sporadic")
  expect_equal(row("ab")$tier, "widespread")
  expect_true(is.na(row("onlyC")$tier))
  # prevalence rule: 3 of 10 = 30% passes 10%, but a 40% cut-off fails it
  t40 <- sf_species_table(long, group_n, min_samples = 3, min_prev = 0.40)
  expect_equal(t40$tier[t40$name == "ab"], "recurrent")
})

test_that("prevalence and abundance are per group", {
  expect_equal(row("wide")$n_A, 6L)
  expect_equal(row("wide")$prev_C, 3 / 5)
  expect_equal(row("onlyC")$abund_A, 0)
  expect_gt(row("onlyC")$abund_C, 0)
})

test_that("pairwise view covers every pair and re-tiers on the pair", {
  pw <- sf_pairwise(tab, names(group_n), 3, 0.10)
  expect_setequal(unique(pw$pair), c("A vs B", "A vs C", "B vs C"))
  # 'rec' has only 2 samples per group -> recurrent on every pair
  expect_true(all(pw$tier[pw$name == "rec"] == "recurrent"))
  # 'spor' is shared only by A and B
  expect_equal(sort(unique(pw$pair[pw$name == "spor"])), "A vs B")
})

test_that("numeric linkage variables are binned; categorical are kept", {
  expect_equal(sf_bin_variable(c("x", "y", NA)), c("x", "y", NA))
  b <- sf_bin_variable(1:30, 3)
  expect_setequal(unique(b), c("low", "mid", "high"))
  expect_equal(sf_bin_variable(c(1, 2, 1), 3), c("1", "2", "1"))
})

test_that("linkage joins nodes in different groups by shared widespread species", {
  ward <- c(stats::setNames(rep(c("w1", "w2"), each = 5), A), stats::setNames(rep(c("w1", "w2"), each = 5), B),
            stats::setNames(rep("w3", 5), C))
  grp <- stats::setNames(substr(names(ward), 1, 1), names(ward))
  ws <- tab$name[!is.na(tab$tier) & tab$tier == "widespread"]
  lk <- sf_linkage(long, ward, ws, grp, min_shared = 1, scope = "between")
  expect_true(all(lk$nodes$group %in% c("A", "B", "C")))
  e <- lk$edges
  expect_true(all(substr(e$from, 1, 1) != substr(e$to, 1, 1)))      # between groups only
  expect_equal(e$shared[e$from == "A · w1" & e$to == "B · w1"], 3L) # wide + ab + filler (in every sample)
  lt <- sf_linkage_test(long, ward, ws, grp, lk, n_perm = 49, seed = 1, scope = "between")
  expect_true(all(lt$edges$p_perm > 0 & lt$edges$p_perm <= 1))
})

test_that("tree tips need min_samples detections; widespread species always kept", {
  one <- dplyr::bind_rows(long, tibble::tibble(sample = "C5", name = "single", count = 5000L, group = "C",
                                               phylum = "Pseudomonadota"))
  t1 <- sf_species_table(one, group_n, 3, 0.10)
  ws <- t1$name[!is.na(t1$tier) & t1$tier == "widespread"]
  expect_true("single" %in% sf_tree_tips(t1, names(group_n), top_n = 10, min_samples = 1))
  tips <- sf_tree_tips(t1, names(group_n), top_n = 10, min_samples = 2, widespread = ws)
  expect_false("single" %in% tips)
  expect_true(all(ws %in% tips))
})

test_that("element classes map from ABRicate databases", {
  expect_equal(sf_element_class(c("card", "vfdb", "plasmidfinder", "other")), c("ARG", "VF", "MGE", NA))
  hits <- tibble::tibble(sample = c("A1", "A1", "A2"), GENE = c("g1", "g1", "v1"),
                         DATABASE = c("card", "card", "vfdb"), name = c("wide", "wide", "wide"),
                         sequence = c("r1", "r2", "r3"))
  el <- sf_elements(hits, "wide", c("A1", "A2"), min_reads = 2)
  expect_equal(el$reads[el$gene == "g1"], 2L)
  expect_true(el$low_support[el$gene == "v1"])
})
