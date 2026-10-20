# MERIDIAN (unreleased)

- Relative abundance (R/05): top-N taxa are now ranked by MEAN per-sample
  relative abundance (%) instead of total raw counts (which favoured the deepest
  samples), and the default `relative_abundance.top_n` is 15 (was 50). One shared
  helper (`ra_prepare_composition`) serves the per-sample, per-rank and group plots.
- Relative abundance is species rank only: R/02 now carries the kraken2 rank code
  as `tax_rank`, and R/05 keeps rank "S" so strain / subspecies rows no longer
  double-count their parent species.
- New main figure `relative_abundance_group.png` (+ `relative_abundance/relative_abundance_group.csv`):
  one stacked bar per `facet_by` level, mean of the per-sample %, n per group in
  the axis labels. The per-sample `relative_abundance.png` is kept as the
  supplementary figure; streamgraphs are unchanged in kind.
- "Others" is always light grey `#BDBDBD` in all R/05 figures; taxon colours never use it.
- Streamgraphs: title is "Species composition by <group>" (was hard-coded
  "treatment"), and the unlabelled grey band is fixed. It was a taxon in the stream's
  own top-12 (ranked by a different measure) that had no entry in the count-ranked
  palette, so it fell back to `na.value` grey and had no legend key. The stream now
  receives the same top-N and palette as the bars (`build_stream_composition(categories=)`).
- Dashboard (templates/report.qmd) shows the group figure first and the
  per-sample figure second; manifest entries added.
- Add testthat test (tests/testthat/test-relab-topn.R).

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
