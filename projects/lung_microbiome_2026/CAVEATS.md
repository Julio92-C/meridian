# lung_microbiome_2026 — open caveats

Recorded 2026-10-05 after the first full run (`config.yaml`, Exhale vs Sputum,
outputs in `JC_2026-09_EBSP/test_run/`, log `scripts/_run_lung2026_20261005_135016.log`).
**Status: deferred — to be worked on later.** Treat all p-values from this run
as provisional until items 1–2 are resolved.

## 1. Patient random effect in the statistics — RESOLVED 2026-10-05 (opt-in)

Originally no analysis stage read `cfg$metadata$random_effect`, so all 89
samples (20 patients, 1–8 each) were treated as independent.

Fix: `R/utils_cluster_perm.R`, enabled by `stats.random_effect_tests: true`
(set in both lung_microbiome_2026 configs; off by default elsewhere so
lung_microbiome 2025 / chicken_batch2 results are unchanged).
- KW / pairwise Wilcoxon / PERMANOVA / PERMDISP keep their statistics but take
  p-values from patient-respecting permutations:
  - `sample_type` (varies within patients) → labels shuffled within each patient;
  - `immunodeficiency` (one label per patient) → whole patients' labels shuffled.
- ALDEx2 runs on counts pooled per patient (immunodeficiency) or per
  patient × sample type (sample_type).
- Validated by `scripts/validate_cluster_perm.R` (23 checks). Simulated
  between-patient null: naive KW false-positive rate 68 %, patient-aware 5 %.

Still NOT patient-aware: the Mantel tests in R/12 (taxonomy vs resistome etc.).
They don't test a grouping, but repeated samples still inflate their p-values.

## 2. PERMDISP significant for gene-level beta diversity

Resistome (PERMDISP p = 0.0016), virulome (p = 0.0005) and mobilome
(p = 0.0001) all show unequal dispersion between groups. The significant
PERMANOVA results (resistome p = 0.0029, virulome p = 0.011) may reflect
greater spread in one group, not a shift in composition.

## 3. Near-empty samples still enter alpha diversity

By decision, the 6 near-empty samples (18–80 Kraken2 reads: SPA2-1005,
SPA1-1005, SPF0-1536, SPY3-1005, SPN3-2003, EBF0-1604) stay in the metadata.
Alpha diversity covers all 89 samples (rarefaction min depth = 92 reads vs
median 12,527), so they can drive the richness / Chao1 differences.
Revisit: depth cutoff or a rarefied sensitivity analysis.

## 4. Sample coverage differs by layer

| Layer | Samples |
|---|---|
| Metadata (analysis set) | 89 |
| Re-centrifuge | 89 |
| Bracken | 84 (5 near-empty samples have no output) |
| Abricate | 82 (7 with no hits — zero-hit vs failed run still to verify on Myriad) |
| Merged Abricate × Bracken table | 76 |
| Resistome / virulome / mobilome | 56 / 65 / 38 |

Patient 1005's sputum series is the worst affected (SPA1/SPA2/SPY3 < 70 reads).

## 5. Leftover chicken-study labels (cosmetic)

- Supplementary panels are named `figS_diet_effects_{arg,vf,mge}`.
- Resistome stage emits organism panels for *C. difficile* and
  Enterobacteriaceae (0 and 1 hits here).

## 6. Taxon → ARG → MGE Sankey is nearly empty

Only 4 flows / 1 phylum / 1 MGE type: few plasmid replicons co-occur on reads
carrying ARGs. Probably not a usable figure for this cohort.

## 8. ALDEx2 Wilcoxon p = 0 with 2-patient groups — RESOLVED 2026-10-08

Fix (R/08): per pair, ALDEx2's `wi.*` columns are dropped when a group has
< 3 units (`differential_abundance.min_wilcox_group`) or the exact minimum
two-sided p, 2 / C(n1+n2, n1), cannot reach alpha; that pair falls back to
Welch. A dot plot left by an earlier run is now deleted when the plot is
skipped, so the manifest can't list it. Re-run: no taxa/genes significant in
any immunodeficiency pair; other projects unaffected (smallest pairs 6 vs 6).

Original finding:

After pooling to patients, XLA and HG have 2 units each. ALDEx2 reports
`wi.ep = 0` / `wi.eBH < 0.05` for many taxa in HG_vs_XLA (12) and SID_vs_XLA
(27), e.g. *G. haemolysans*, *H. parahaemolyticus*. An exact two-sided
Wilcoxon test cannot go below p = 0.33 (2 vs 2) or 0.056 (2 vs 7), so these
are artifacts of ALDEx2's fast Wilcoxon approximation at tiny n. Welch
(`we.eBH`) finds nothing significant in any pair. The `aldex2_dotplot_summary`
(28 features × 6 pairs) is driven by these artifacts — do not report it.
Candidate fix: ignore `wi.*` columns (or skip the pair) when any group has
< 3 units.

## 7. Bracken on long reads

nf-core/taxprofiler refuses to run Bracken on ONT data ("not evaluated for
long reads"). Bracken was run standalone with the longest available
kmer_distrib (`Bash_scripts/bracken.sh`), so Bracken counts are approximate.
