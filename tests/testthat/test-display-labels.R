# Tests for the display helpers in utils_panel_io.R: metadata column names
# shown as words, and room above the data for in-panel stats labels.
source(file.path("..", "..", "R", "utils_panel_io.R"))

test_that("pretty_var turns a column name into display text", {
  expect_equal(pretty_var("sample_type"), "Sample type")
  expect_equal(pretty_var("Macrolide_prophylaxis"), "Macrolide prophylaxis")
  expect_equal(pretty_var("Diet"), "Diet")
})

test_that("stats_headroom grows with the number of label lines", {
  one   <- stats_headroom(1)
  three <- stats_headroom(n_text_lines("a\nb\nc"))
  expect_equal(n_text_lines("a\nb\nc"), 3)
  expect_gt(three[3], one[3])          # top expansion (mult) grows
  expect_equal(three[1], one[1])       # bottom expansion unchanged
})
