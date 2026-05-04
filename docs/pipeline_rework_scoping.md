# Pipeline rework — scoping notes

Findings from the first end-to-end smoke test against the real Chicken batch 1
data (`PC_JC_2024-11-29_Julio_Gallus`), 2026-05-01 → 2026-05-04. The
test never produced clean outputs; instead it surfaced that several modules
encode template assumptions rather than logic validated against real inputs.

This doc is the planning anchor for the rework. The next working session should
start here, not at `run_pipeline.R`.

---

## Design principle (carry into every module)

**Pipeline modules must adapt to per-project input shape.** Real studies differ
in: metadata column layout, treatment counts, presence/absence of negative
controls, Re-centrifuge column structure, and which taxonomic ranks are
populated. Anything that touches column names, sample IDs, or rank handling
must be parameterised through `cfg` or detected from the data — never
hard-indexed.

Example anti-pattern that broke this run:

```r
# R/02_clean_data.R:26
keep_cols <- c(1, 2, seq(5, ncol(rcf), by = 3))
```

Assumes a fixed column triplet structure. The real Re-centrifuge CSV repeats
sample triplets at every taxonomic rank (species/genus/family/order/class/
phylum/kingdom/domain) and ends with `Rank` and `Name` columns — none of which
this slice accounts for.

---

## Bugs fixed during the smoke test

1. **`run_pipeline.R` crashed under `Rscript`.** `sys.frame(1)$ofile` threw
   before `%||%` (also defined later) could fall back. Replaced with a
   `tryCatch` + `--file=` arg path. The README's quickstart command would
   have failed for any user before this fix.
2. **`R/02_clean_data.R` `left_join` type mismatch.** `rcf$taxid` parsed as
   character, `kraken2$taxid` as double. Coerced both join keys to character
   before the join.
3. **Per-stage progress added to `run_pipeline.R`.** Each enabled stage now
   logs `[i/n  pct%] stage_name — start` / `… — done in Xs` so the user
   can see whether a long stage is working or hung.

---

## Open issues surfaced (not fixed)

### Data-shape mismatches in `R/02_clean_data.R`

- **Re-centrifuge slicing assumption is wrong** (see anti-pattern above). The
  actual file has ~350 columns spanning 8 taxonomic ranks per sample, plus
  trailing `Rank` / `Name` columns.
- **`tidyr::separate(SEQUENCE, into = c("sequence", "taxid"), sep = "_")`**
  on line 14 emits "Additional pieces discarded in 412 rows". Any `SEQUENCE`
  value with >1 underscore loses the right-hand pieces, so `taxid` ends up
  holding a middle fragment — silent data-correctness bug.
- **Read of `Re-centrifuge .rcf.data.csv` produces duplicate column names**
  that `read_csv` auto-disambiguates as `SHARED_SUMMARY...347`, etc. The
  module never strips these, so renaming downstream is fragile.

### Config–data mismatches for chicken_batch1

| Config field | Config value | Reality on disk |
|---|---|---|
| `metadata.file` | `Datasets/chicken_metadata.csv` | This file has every-other-row blank (Excel artifact). `chicken_metadata1.csv` is the clean version (18 rows, no gaps). |
| `metadata.controls` | `[CTRL1, CTRL2]` | No control samples exist. All 18 samples (D19–D36) are treatment samples. |
| `metadata.fixed_effects` | `[Treatment, Weight]` | Real metadata has only `sample` and `Treatment` — no `Weight` column. |
| `metadata.random_effect` | `block` | No `block` column exists. |
| `inputs.bracken_report` | `Reports/Kraken2/kraken2_db2-bracken_combined_reports.csv` | Real file is `Reports/Kraken2/bracken_arranged.csv`. |
| `taxid_fixes_file` | `Metadata/taxid_fixes.csv` | No `Metadata/` folder. README "Open items" already flags that this needs to be populated. |

The placeholder values came from `config/config_template.yaml`. Before any
real run, the template needs comments explaining which fields are required
vs. optional, and `chicken_batch1/config.yaml` needs to be reconciled against
`chicken_metadata1.csv`.

### Likely issues in downstream stages (not yet hit)

Each of `R/03_normalisation.R` through `R/12_network.R` was scaffolded against
the same template assumptions and has the same risk profile. Specifically
worth checking:

- `R/06_alpha_diversity.R` and `R/07_beta_diversity.R` — assume a single
  `Treatment` grouping column; need to handle `cfg$metadata$group_cols` as a
  list of one-or-more.
- `R/08_differential_abundance.R` — pairwise ALDEx2 across treatment levels;
  the chicken data has 3 treatments → 3 pairwise comparisons. Verify the
  loop is treatment-count agnostic.
- `R/09_resistome.R` / `10_virulome.R` / `11_mobilome.R` — split by ABRicate
  database (CARD / VFDB / PlasmidFinder). The summary file has a `DATABASE`
  column; confirm the split key is read from data, not hard-coded.

---

## Plan for next session — reverse-engineer from the working scripts

The real ground truth for what each stage should produce lives in the data
folder's hand-edited scripts:

```
C:\Users\julio\OneDrive\Desktop\PC_JC_2024-11-29_Julio_Gallus\R_scripts\
```

46 scripts there, including:

- `cleanData.R` (or equivalent) — what the original `clean_data` step really does
- `aldexDA.R`, `aldex_overall_and_pairwise.R` — ALDEx2 logic
- `alpha_diversity_violinPlot.R` — alpha diversity
- `GEsNorm_PCoA.R`, `PCoA.R`, `abri_kraken2_PCoA.R` — beta diversity
- `chordDiagram*.R`, `*VennDiagram.R`, `*pHeatmap*.R` — figure scripts

### Suggested per-stage workflow

For each module `R/NN_<stage>.R`:

1. Find the corresponding hand-edited script(s) in `R_scripts/`.
2. Read them — note the actual input columns consumed and output files
   produced.
3. Identify which values are project-specific (sample IDs, treatment names,
   control IDs, paths) vs. truly fixed (statistical parameters, ABRicate
   databases). Project-specific values move into `cfg`.
4. Rewrite the module so it consumes whatever column shapes the loaded data
   actually has, driven by `cfg`.
5. Validate the module's output against the corresponding "known good" file
   already in the project's `Datasets/` (e.g. `abri_kraken2_cleaned.csv`,
   `aldex2_results.tsv`, `gephi_edges2.csv`).

### Recommended order

1. `R/02_clean_data.R` — every downstream stage depends on it. Validate
   against `Datasets/abri_kraken2_cleaned.csv` and
   `Datasets/noncontaminants_list.csv` (if present).
2. `R/01_load_inputs.R` — minor; mostly already shaped right. Add the
   `diversity_csv` loader (declared in config but unused).
3. `R/03_normalisation.R` → validate against `Datasets/genetable_normdata.csv`.
4. `R/08_differential_abundance.R` → validate against `Datasets/aldex2_results.tsv`
   and the three `*_aldex2_results.tsv` pairwise files.
5. Diversity + figures + network in any order.

### Reconcile the chicken_batch1 config along the way

Once `clean_data` is rewritten, update `projects/chicken_batch1/config.yaml`
to match the real metadata: `chicken_metadata1.csv`, no controls,
no `block`/`Weight`, fixed `bracken_report` path. Use the rewrite as the
forcing function to clean up `config/config_template.yaml` too — annotate
each field as required/optional and what data shape it implies.

---

## Test artefacts left in repo

- `projects/test_real_data/config.yaml` — kept as-is for next session.
  Outputs are redirected into `<project_root>/test_run/` so the real data
  folder stays untouched.
- `R/02_clean_data.R` — has the type-coercion patch (line 36). Will be
  rewritten in the rework.
- `run_pipeline.R` — has the Rscript fix and per-stage progress reporter.
  These should survive the rework as-is.
