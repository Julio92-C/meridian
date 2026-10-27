# Slot `composite:` overrides applied to a panel's ggplot (R/13_panels.R).
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("..", "..", "R", "13_panels.R"))

p <- ggplot2::ggplot(data.frame(x = c("a", "b"), y = 1:2, f = c("a", "b")),
                     ggplot2::aes(x, y, fill = f)) + ggplot2::geom_col()

test_that("legend_title: false blanks the legend title", {
  out <- .panels_apply_overrides(p, list(legend_title = FALSE))
  expect_s3_class(out$theme$legend.title, "element_blank")
})

test_that("legend_nrow sets the rows of the fill legend", {
  out <- .panels_apply_overrides(p, list(legend_nrow = 2))
  expect_equal(out$guides$guides$fill$params$nrow %||% out$guides$fill$nrow, 2)
})
