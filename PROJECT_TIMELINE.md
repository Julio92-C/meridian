# Metagenomics Pipeline Automation — Timeline

**Owner:** Julio Cesar Ortega Cambara
**Started:** 2026-04-17
**Target completion:** 2026-06-12
**Status:** 🚧 in progress — Phase 1 underway

**Cadence:** ~4 working days/week (Mon–Fri, one weekday buffer)
**Validation pattern (per module):** read ground-truth script in
`PC_JC_2024-11-29_Julio_Gallus/R_scripts/` → rewrite the R/NN_*.R module
to be cfg-driven → diff output against the corresponding known-good
file in `Datasets/` (the same approach used for R/03 on 2026-05-05).

## Phases

Status legend: ✅ done · 🚧 in progress · 📋 planned

### ✅ Phase 0 — Foundations (2026-04-17 → 2026-05-05)
- [x] Repo + config-driven scaffold (R/00–R/12, run_pipeline.R, configs)
- [x] renv lockfile pinned (R 4.5.0, 180 packages)
- [x] First end-to-end smoke test on chicken_batch1 (surfaced rework scope)
- [x] R/01 load_inputs — Re-centrifuge fread fix (vroom segfault)
- [x] R/02 clean_data — Bracken merge + metadata join, cfg-driven sample detection
- [x] R/03 normalisation — TPM per (sample, GENE), validated vs `genetable_normdata.csv`
- [x] `scripts/diff_normdata.R` — reusable diff harness for future stages

### 🚧 Phase 1 — Statistical core (2026-05-06 → 2026-05-22)

Goal: every R/NN_*.R module rewritten and validated against
chicken_batch1 ground truth.

**Week of May 4–8** (3 working days remaining)
- [ ] **Wed May 6** — R/04 taxonomy (Venn + pHeatmap) — validate vs `taxa_category.csv` and the `*VennDiagram.R` / `*pHeatmap*.R` scripts
- [ ] **Thu May 7** — R/05 relative abundance (stacked bars + plotly) — validate vs `*Norm_relativeAbundance.R` outputs
- [ ] **Fri May 8** — R/09/10/11 resistome / virulome / mobilome — quick `DATABASE`-column splits over the merged table (1 day total for all three)

**Week of May 11–15** (4 working days)
- [ ] **Mon May 11** — R/06 alpha diversity (richness, Shannon, Kruskal–Wallis, violin) — validate vs `alpha_diversity_violinPlot.R` + `ARGNorm_alphaDiversity.R`
- [ ] **Tue May 12** — R/07 beta diversity (Bray–Curtis → PCoA → PERMANOVA) — read `GEsNorm_PCoA.R`, `taxaNormProfile_PCoA.R`
- [ ] **Wed May 13** — R/07 validation + buffer
- [ ] **Thu May 14** — R/08 ALDEx2 — read `aldexDA.R` + `aldex_overall_and_pairwise.R`, plan pairwise loop (treatment-count agnostic)
- [ ] _Fri May 15: buffer_

**Week of May 18–22** (4 working days)
- [ ] **Mon May 18** — R/08 ALDEx2 implement
- [ ] **Tue May 19** — R/08 validate vs `aldex2_results.tsv` + `*_aldex2_results.tsv` pairwise files
- [ ] **Wed May 20** — R/12 network — read `chordDiagramNorm_v1.2.1.R` + Gephi exports
- [ ] **Thu May 21** — R/12 validate vs `gephi_edges*.csv` / `gephi_nodes*.csv`
- [ ] _Fri May 22: buffer_

**Phase 1 deliverable (Fri May 22):** every R module validated on
chicken_batch1; row-count divergences (if any) explained as in R/03.

### 📋 Phase 2 — End-to-end report on chicken_batch1 (2026-05-25 → 2026-06-02)

Goal: full pipeline produces a publication-grade HTML report
reproducing the chicken_batch1 manuscript figures.

**Week of May 25–29** (4 working days)
- [ ] **Mon May 25** — R/13 / `templates/report.qmd` — wire every stage's outputs into the Quarto report, verify it renders
- [ ] **Tue May 26** — Reconcile `projects/chicken_batch1/config.yaml` (drop CTRL1/CTRL2, drop Weight, drop block, fix bracken path, point at `chicken_metadata1.csv`)
- [ ] **Wed May 27** — Full end-to-end run with all stages on; capture log
- [ ] **Thu May 28** — Compare rendered figures against published manuscript outputs
- [ ] _Fri May 29: buffer_

**Week of Jun 1–5** (early days)
- [ ] **Mon Jun 1** — Polish + bug fixes surfaced by the full run
- [ ] **Tue Jun 2** — Tag `v0.1.0` (chicken_batch1 reproduction baseline)

**Phase 2 deliverable (Tue Jun 2):** `v0.1.0` — tool reproduces a
peer-reviewed study end-to-end from one config file.

### 📋 Phase 3 — Multi-study validation & release (2026-06-03 → 2026-06-12)

Goal: confirm the pipeline generalises beyond chicken_batch1 by running
the three other studies in `Metagenomics_Projects_paths.txt`.

**Week of Jun 1–5** (remaining days)
- [ ] **Wed Jun 3** — chicken_batch2 (`JC_FL_2025_05_20_CK_microbiome`) — config + run
- [ ] **Thu Jun 4** — lung_microbiome — config + run

**Week of Jun 8–12** (4 working days)
- [ ] **Mon Jun 8** — hospital_microbiome — config + run
- [ ] **Tue Jun 9** — Document per-study patterns; annotate `config/config_template.yaml` (required vs optional fields)
- [ ] **Wed Jun 10** — Address any per-study issues that surfaced
- [ ] **Thu Jun 11** — Tag `v0.2.0` (multi-study validated)
- [ ] **Fri Jun 12** — Push to remote, update README status badge `alpha → beta`

**Phase 3 deliverable (Fri Jun 12):** `v0.2.0` — validated on 4 studies;
ready for outside users.

## Notes

- **Effort budget:** this plan assumes ~1 day per simple module
  (taxonomy, relative abundance, ABRicate splits) and 1.5–2 days for
  complex ones (PCoA, ALDEx2 pairwise loop, network). Today's R/03
  rewrite + load fix took ~1 focused session including a parsing-bug
  detour, which is the calibration anchor.
- **One buffer day per week is non-negotiable.** Today's session
  needed ~30% extra time on the variant-suffix taxid parsing
  investigation; without slack, every spillover pushes a deadline.
- **Don't try to make pipeline output exactly match GT row-for-row.**
  GT `cleanData.R` has a `separate(sep="_")` parsing bug that drops the
  real taxid for variant-suffix contigs. R/02 fixes this; expect ~+110
  (sample, GENE) pairs vs `genetable_normdata.csv` on chicken_batch1.
  The acceptance bar is **schema match + plausible values**, not row-
  for-row identity. (See `memory/project_chicken_batch1_taxid_parse_bug.md`.)
- **Risk — Phase 3 dependencies:** chicken_batch2 / lung / hospital may
  have data layouts that require new cfg fields the template doesn't
  yet support. If so, those days are spent extending cfg parsing rather
  than just running. Allow Phase 3 to slip a week without considering
  the project late.
- **Out of scope for v0.2.0:** publication of the tool itself,
  CI/automated tests, support for short-read data. Track separately
  if/when raised.

<!-- AUTO-DETECTED ZONE — content below this line is overwritten by `python -m agent.auto_update --write`. Edit above this line freely; it is preserved every run. -->

## Recent activity (auto-detected, last 14 days)

_Generated 2026-05-09._

- Phase 0 closed out on May 4–5 with commits `92c9502` (TPM normalisation validated against `genetable_normdata.csv`, 229 insertions across 5 files) and `8b8d18a` (`data.table::fread` fix in `R/01_load_inputs.R` to work around the vroom segfault on Re-centrifuge input) — the two largest commits in the window.
- `R/02_clean_data.R` received the Bracken merge + metadata join rewrite on May 5; on the same day `scripts/diff_normdata.R` (reusable diff harness for future stages) and `scripts/probe_taxid_parse.R` (taxid variant-suffix investigation) were added.
- All 12 recently-modified files sit inside `R/`, `scripts/`, `projects/test_real_data/`, and `config/` — no activity outside the analysis core.
- **No commits or file changes since May 5** (4 days of silence as of today, May 9); 2 Claude Code sessions logged in the window, both on or before May 5.
- The three Phase 1 tasks scheduled for the opening week — R/04 taxonomy (May 6), R/05 relative abundance (May 7), and R/09–11 resistome/virulome/mobilome (May 8) — have no corresponding files or commits; Phase 1 is 3 working days behind its opening-week plan.
