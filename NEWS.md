# MERIDIAN 1.0.1

- Fix: `ge_minmax_per_col()` (R/utils_ge_profile.R) now keeps matrix dimensions
  for single-gene / single-sample inputs. Previously `apply()` collapsed them to a
  vector, crashing the virulome (and potentially resistome/mobilome) heatmaps with
  `incorrect number of dimensions` when only one gene survived filtering.
- Add testthat regression test (tests/testthat/test-ge-minmax.R).
- Default container image / SIF in `meridian.sh` bumped to v1.0.1.

# MERIDIAN 1.0.0

- First public release (container ghcr.io/julio92-c/meridian:v1.0.0).
