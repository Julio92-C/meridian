# Pairwise cluster-permutation tests must agree with the overall cluster KW
# when there are only two groups, even when samples carry different numbers
# of hits (permutations then change the group sizes).
source(file.path("..", "..", "R", "utils_cluster_perm.R"))

set.seed(1)
samples  <- paste0("s", 1:12)
patients <- rep(paste0("p", 1:6), each = 2)
grp_of   <- setNames(rep(c("A", "B"), each = 6), samples)
hits     <- c(3, 9, 2, 7, 4, 12, 5, 1, 8, 2, 6, 3)        # unequal hits per sample
unit     <- rep(samples, hits)
cl       <- rep(patients, hits)
g        <- unname(grp_of[unit])
v        <- rnorm(length(unit), mean = ifelse(g == "A", 0.4, 0))

test_that("two-group pairwise p equals the overall cluster KW p", {
  kw <- cluster_kw(v, g, cl, units = unit, nperm = 999L, seed = 7L)
  pw <- cluster_wilcox_pairs(v, g, cl, units = unit, nperm = 999L, seed = 7L)
  expect_equal(nrow(pw), 1)
  expect_equal(pw$p, kw$p.value)
})

test_that("pairwise rows are adjusted across pairs", {
  g3 <- ifelse(unit %in% c("s9", "s10", "s11", "s12"), "C", g)   # whole patients p5, p6
  pw <- cluster_wilcox_pairs(v, g3, cl, units = unit, nperm = 199L, seed = 7L)
  expect_equal(nrow(pw), 3)
  expect_true(all(pw$p.adj >= pw$p))
})
