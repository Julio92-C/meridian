# Regression: ge_plot_beta() must not crash in adonis2 when every
# hit-positive sample belongs to one group -> "contrasts can be applied only
# to factors with 2 or more levels" in the resistome stage on sparse data.
# PERMANOVA/PERMDISP are skipped; the PCoA ordination + scores CSV are still
# written (feat/cluster-aware-stats 11f20fb keeps the ordination, which
# supersedes the v1.0.2 whole-stage skip).
`%||%` <- function(a, b) if (is.null(a)) b else a
source(file.path("..", "..", "R", "00_setup.R"))        # padjust_*
source(file.path("..", "..", "R", "utils_panel_io.R"))
source(file.path("..", "..", "R", "utils_ge_profile.R"))

# Three genes so the Bray-Curtis PCoA has 2 positive axes (two genes over
# three samples is collinear and hits the separate < 2 axes skip).
beta_hits <- function(groups) {
  data.frame(
    sample = rep(paste0("s", seq_along(groups)), each = 3),
    GENE   = rep(c("geneA", "geneB", "geneC"), times = length(groups)),
    TPM    = c(5, 1, 2,  2, 4, 1,  3, 3, 6,  1, 6, 2)[seq_len(3 * length(groups))],
    Group  = rep(groups, each = 3)
  )
}

test_that("single-group hits skip PERMANOVA/PERMDISP but keep the PCoA", {
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
  expect_false(is.null(r))
  expect_true(any(grepl("PERMANOVA/PERMDISP skipped", logged)))
  expect_false(any(file.exists(out[2:3])))   # no permanova / permdisp txt
  expect_true(file.exists(out[4]))           # scores CSV still written
})
