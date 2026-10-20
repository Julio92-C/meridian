# Regression: ge_plot_beta() must skip (not crash in adonis2) when every
# hit-positive sample belongs to one group -> "contrasts can be applied only
# to factors with 2 or more levels" in the resistome stage on sparse data.
source(file.path("..", "..", "R", "utils_ge_profile.R"))

beta_hits <- function(groups) {
  data.frame(
    sample = rep(paste0("s", seq_along(groups)), each = 2),
    GENE   = rep(c("geneA", "geneB"), times = length(groups)),
    TPM    = c(5, 1, 2, 4, 3, 3, 1, 6)[seq_len(2 * length(groups))],
    Group  = rep(groups, each = 2)
  )
}

test_that("single-group hits skip beta instead of erroring", {
  logged <- character(0)
  pipeline_log <<- function(cfg, msg) logged <<- c(logged, msg)
  out <- file.path(tempdir(), c("pcoa.png", "perm.txt", "disp.txt", "scores.csv"))
  unlink(out)
  expect_no_error(r <- ge_plot_beta(
    beta_hits(c("Yes", "Yes", "Yes")), "Group", c(Yes = "red"), opts = list(),
    fig_file = out[1], ds_perm_file = out[2], ds_disp_file = out[3],
    ds_pcoa_file = out[4], log_label = "Resistome", title_prefix = "ARG",
    skip_data_label = "ARG hits", cfg = list(stats = list(permanova_permutations = 99))
  ))
  expect_null(r)
  expect_true(any(grepl("span 1 Group group", logged)))
  expect_false(any(file.exists(out)))
})
