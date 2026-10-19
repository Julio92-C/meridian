# Regression: ge_minmax_per_col() must keep matrix dims for 1-gene / 1-sample
# inputs (apply() used to collapse them to a vector -> virulome crash).
source(file.path("..", "..", "R", "utils_ge_profile.R"))

test_that("single-gene matrix keeps dims and dimnames", {
  m <- matrix(c(1, 5), nrow = 1, dimnames = list("vfG", c("s1", "s2")))
  r <- ge_minmax_per_col(m)
  expect_true(is.matrix(r)); expect_equal(dim(r), c(1L, 2L))
  expect_equal(dimnames(r), dimnames(m))
  expect_equal(unname(r[1, ]), c(0, 0))   # zero range per column
  expect_no_error(r[, c("s2", "s1"), drop = FALSE])
})

test_that("single-sample matrix keeps dims and scales to [0,1]", {
  m <- matrix(c(1, 3, 5), ncol = 1, dimnames = list(c("a", "b", "c"), "s1"))
  r <- ge_minmax_per_col(m)
  expect_equal(dim(r), c(3L, 1L)); expect_equal(unname(r[, 1]), c(0, .5, 1))
})

test_that("regular matrix and constant column behave as before", {
  m <- matrix(c(1, 2, 3, 4, 4, 4), ncol = 2, dimnames = list(letters[1:3], c("x", "y")))
  r <- ge_minmax_per_col(m)
  expect_equal(unname(r[, "x"]), c(0, .5, 1)); expect_equal(unname(r[, "y"]), c(0, 0, 0))
  expect_equal(dimnames(r), dimnames(m))
})
