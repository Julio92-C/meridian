# Multiple-testing correction in the pipeline

How `cfg$stats$padjust_method` flows through every adjustable surface, and what stays uncorrected and why.

Landed 2026-06-11 in commit `04c5363` (`feat(stats): cfg$stats$padjust_method knob`).

---

## The config knob

```yaml
# config/config_template.yaml — and any projects/<slug>/config.yaml
stats:
  padjust_method: BH    # default
```

Accepted values (validated by `padjust_method()` in `R/00_setup.R`):

| Value | Method | When to use |
|---|---|---|
| `BH` *(default)* | Benjamini-Hochberg, 1995 | Standard for genomics / microbiome compositional data. Controls FDR under independence or PRDS (positive regression dependence). Most powerful when the dependence is positive — which is the common case for shared-sample / co-regulated features. |
| `fdr` | Alias for BH | R's `p.adjust` accepts both names. |
| `BY` | Benjamini-Yekutieli, 2001 | Controls FDR under **arbitrary** dependence (including negative correlation). More conservative — adjusted p-values are roughly `ln(m) + γ` times larger than BH. Use when a reviewer asks for dependence-robust FDR, or for final-claim conservativeness. |
| `bonferroni` | Bonferroni FWER | Most conservative. Only for very small families where you want zero false-positive risk. |
| `holm`, `hochberg`, `hommel` | FWER step-down / step-up | Rare in this codebase; available for consistency with `stats::p.adjust`. |
| `none` | No correction | Reports raw p-values. Use for exploratory analysis or when reporting a single test in isolation. |

Two helpers in `R/00_setup.R` do the plumbing:

- `padjust_method(cfg)` — returns the validated method string, errors fast on typos.
- `padjust_p(p_vector, cfg)` — `stats::p.adjust(p, method = padjust_method(cfg))`, with `"none"` as a no-op pass-through. Call this at every adjustment site instead of a hard-coded `method = "BH"`.

---

## Surfaces adjusted by the knob

Every site where the pipeline currently performs `p.adjust` (or equivalent), grouped by what counts as a "family":

| Surface | R file / line | Family unit | What gets adjusted |
|---|---|---|---|
| ALDEx2 GLM | `R/08_differential_abundance.R:182` | All per-feature p-values from one GLM coefficient | Passed directly as `fdr.method` to `ALDEx2::aldex.glm`, which forwards to `p.adjust`. |
| ALDEx2 pairwise dot plot | `R/08_differential_abundance.R` Sec 6b | All per-feature p-values within one pairwise comparison | Raw `we.ep` / `wi.ep` columns re-adjusted into `we.eAdj` / `wi.eAdj` per comparison. ALDEx2's baked-in `eBH` stays available for downstream filtering. |
| VF × ARG correlation heatmap | `R/10_virulome.R:270` | All `(VF category × drug class)` cells in one go | Flatten the p-matrix → `padjust_p` → reshape → significance stars. |
| Mantel triangle (4-omics) | `R/12_network.R::.network_build_mantel_triangle` | 6 upper-triangle pairs from `{taxonomy, resistome, virulome, mobilome}` | Family-adjust the upper triangle, symmetrise. Written to `network/mantel_pairwise_padj.csv` with `(mantel_r, p_raw, p_adj, method)` per pair. Title shows `"p adj <method>"`. |
| Alpha-diversity Kruskal-Wallis | `R/06_alpha_diversity.R` | 10 metrics: Richness, Shannon, Simpson, Chao1, Pielou, Fisher, Berger-Parker, Effective species, Inverse Simpson, Total counts | Pre-pass collects all KW p-values, `padjust_p` once across the metric family, then per-metric plot titles and `alpha_stats.txt` show `"p_raw = X, p_adj (<method>) = Y"`. |
| PERMANOVA + PERMDISP (R/07) | `R/07_beta_diversity.R` | `{Bray-Curtis, Jaccard}` as one 2-test family per test type | Bray + Jaccard PERMANOVA computed up front; same for PERMDISP; family-adjust; both PCoA panels annotated with adjusted values. CSV: `beta_diversity_padj_summary.csv` with `(distance, test, r2, p_raw, p_adj, method)`. |
| GE-side alpha KW (cross-domain) | `R/utils_ge_profile.R::ge_alpha_kw` + `R/14_manifest.R` finaliser | Per metric, across `{resistome, virulome, mobilome}` | Each call appends raw p to a cross-module accumulator (`<datasets>/stats/ge_alpha_kw_raw.csv`); `padj_summary_finalise(cfg, "ge_alpha_kw")` runs at manifest time, adjusts within each metric family, writes `<datasets>/stats/ge_alpha_kw_padj_summary.csv`. Plot titles stay with raw p (rendered before all calls complete); the CSV is what the report surfaces. |

### How the cross-module accumulator works

For families that span more than one R/ module (currently just the GE-side alpha KW — resistome, virulome, and mobilome each run their own `ge_alpha_kw` calls), the pattern is:

1. `padj_summary_reset(cfg, surface)` — at pipeline startup (in `run_pipeline.R`), clear any prior accumulator so re-runs don't append to stale rows.
2. `padj_summary_record(cfg, surface, family, key, raw_p)` — each call site appends one row to `<datasets>/stats/<surface>_raw.csv`.
3. `padj_summary_finalise(cfg, surface)` — at manifest time (start of `write_manifest`), read the accumulator, group by `family`, family-adjust within each group, write `<datasets>/stats/<surface>_padj_summary.csv`.

Adding a new cross-module family means choosing a `surface` name, deciding what `family` groups it (usually a metric or comparison label), and adding one `padj_summary_record` call at the appropriate site. The reset / finalise plumbing is generic.

---

## Surfaces deliberately **not** adjusted by the knob

These are the multi-test surfaces I considered and deliberately left as raw p-values, with rationale.

| Surface | Why not adjusted |
|---|---|
| Per-organism KW in `R/09_resistome.R::run_species_deep_dive` | One test per organism in a user-specified list (default: *C. difficile*, Enterobacteriaceae). Each organism is a separate scientific claim with its own deep-dive panel; treating them as one family would conflate distinct hypotheses. |
| Per-domain abundance KW in `R/utils_ge_profile.R::ge_abundance_violin` | Each call is one log(TPM)-KW for one domain; the test isn't part of a coherent across-domain family because each domain's gene set is conceptually distinct. |
| Mobile fraction KW (`R/12::.network_build_mobile_fraction`) | Single global test per study. No multiplicity. |
| Per-(metric × domain) GE alpha KW *in plot titles* | The cross-domain family adjustment lands in the summary CSV (see above), but plot titles continue to show raw p because they render before all domains finish. This keeps the figure label legible while still surfacing adjusted values through the manifest. |

If a downstream consumer wants any of these treated as a family, the cross-module accumulator pattern is the lowest-friction extension point — see the section above.

---

## BH vs BY — when to flip

BH controls FDR under independence or positive dependence (PRDS). For compositional / abundance data this assumption almost always holds:

- Co-regulated genes ⇒ positively correlated test statistics.
- Shared samples across pairwise comparisons ⇒ positive dependence.
- Cells of a correlation matrix sharing rows/columns ⇒ positive dependence.

BY is the right call when:

- Tests can be **arbitrarily** dependent, including negative correlations (uncommon in single-study omics, more common when combining heterogeneous data sources).
- A reviewer explicitly asks for "Benjamini-Yekutieli FDR" or "dependence-robust FDR control".
- You want a conservative claim for a single load-bearing result and want to demonstrate it survives even the worst-case correction.

The price: BY p-values are roughly `(ln(m) + γ) ≈ 2.9` times larger than BH for m = 10 tests, `~5×` for m = 100. Some results that pass BH at α = 0.05 will fail BY.

Empirical example from chicken_batch1 Mantel triangle (6 pairs):

| Pair | r | p_raw | BH p_adj | BY p_adj |
|---|---|---|---|---|
| Resistome ↔ Mobilome | 0.310 | 0.0247 | **0.148** | **0.418** |
| Resistome ↔ Virulome | 0.190 | 0.0923 | 0.277 | 0.709 |
| Taxonomy ↔ Mobilome | 0.036 | 0.363 | 0.642 | 1.000 |

Same data — BH gives the most-significant pair `p_adj = 0.148`, BY gives `p_adj = 0.418`. Both are above 0.05, so neither correction calls the omics layers significantly co-correlated; the qualitative conclusion is robust. But for a tighter family the choice matters.

---

## Adding a new adjustment site

Checklist when introducing a new multi-test surface:

1. Compute raw p-values into a vector or matrix.
2. Apply `padjust_p(p, cfg)` (preferred) or `stats::p.adjust(p, method = padjust_method(cfg))`. Never hardcode `method = "BH"`.
3. Surface both raw and adjusted in the artefact:
   - For a CSV: include columns `p_raw`, `p_adj`, `method` (so the consumer can audit the choice).
   - For a plot annotation: prefer adjusted; if showing both fits, label as `p_raw = X, p_adj (<method>) = Y`.
4. If the family spans modules (rare), use the `padj_summary_record` / `padj_summary_finalise` accumulator instead of trying to gather across modules ad hoc.
5. Register any new summary CSV via `table_entry` in `R/14_manifest.R` so the report can surface it.
6. Update the surface table in this doc.

---

## Audit trail

Every manifest writes the resolved `cfg$stats$padjust_method` into `manifest.json` under `config.stats` (was already the case — `cfg$stats` is round-tripped verbatim). So any consumer of the manifest can verify which method produced any given `p_adj` column without having to re-read the project config.
