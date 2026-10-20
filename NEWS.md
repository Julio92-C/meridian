# MERIDIAN 1.0.2

- Fix: beta diversity (`ge_plot_beta()` for resistome / virulome / mobilome,
  and R/07 `run_beta_diversity()`) now skips PERMANOVA/PCoA with a log line when
  the samples span fewer than two groups. Previously `adonis2` crashed with
  `contrasts can be applied only to factors with 2 or more levels`, e.g. when
  every ARG-positive sample in a sparse low-biomass cohort shared one group.
- Add testthat regression test (tests/testthat/test-ge-beta-groups.R).
- Default container image / SIF in `meridian.sh` bumped to v1.0.2.

# MERIDIAN 1.0.1

- Fix: `ge_minmax_per_col()` (R/utils_ge_profile.R) now keeps matrix dimensions
  for single-gene / single-sample inputs. Previously `apply()` collapsed them to a
  vector, crashing the virulome (and potentially resistome/mobilome) heatmaps with
  `incorrect number of dimensions` when only one gene survived filtering.
- Add testthat regression test (tests/testthat/test-ge-minmax.R).
- Default container image / SIF in `meridian.sh` bumped to v1.0.1.

# MERIDIAN 1.0.0

- First public release (container ghcr.io/julio92-c/meridian:v1.0.0).
