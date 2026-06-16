# Manuscript edit plan — post-audit

Generated 2026-05-29 after a deep file-by-file audit of the pipeline
(`run_pipeline.R`, all 13 `R/NN_*.R` modules, three `R/utils_*.R` files,
`templates/report.qmd`, all shipped `projects/*/config.yaml`, `renv.lock`,
`README.md`) against the current `applications_note.md`. Use this file
to drive the next pass of manuscript edits.

---

## 1. Audit findings — what's no longer accurate

### Outdated / wrong claims

1. **ResFinder is not actually used.** Only CARD is the shipped resistome
   database (`R/09_resistome.R:34, 40` — default `database = "card"`).
   ResFinder appears nowhere in `R/`, `templates/`, `config/`, or any
   project config. The manuscript and Figure 1 caption both say
   "CARD/ResFinder".
2. **R/12 produces three diagram classes, not one.** A tripartite ggraph
   network (`network.png`), a chord diagram (default CARD; overall + per
   group), and a 4-tier networkD3 Sankey (default VFDB; overall + per
   group). Plus Gephi-ready edge/node CSVs, igraph topology
   (degree / betweenness / Louvain modularity), Bray-Curtis sample
   clusters, and a degree-distribution PNG. The manuscript currently
   says "chord diagram / network".
3. **R/02 description has two errors.**
   - `R/utils_taxa.R` is NOT called from R/02. It is used by R/04, R/05,
     R/07 and R/12 (grep'd `clean_taxa_names`).
   - R/02 also performs the Bracken merge (`:119-175`) and attaches
     metadata columns declared by `cfg$metadata$fixed_effects` /
     `random_effect`. Both currently omitted.
4. **R/07 isn't just PCoA.** It also runs PERMANOVA (`vegan::adonis2`,
   9999 perms) and PERMDISP (`vegan::betadisper` + `permutest`), with
   R² / p annotated on the PCoA plot.
5. **R/08 runs DAA at two levels.** By default iterates over
   `levels = c("taxa", "gene")` (`:715, 738-755`) — pairwise ALDEx2 for
   taxa (Bracken counts) and for genes (ABRicate contig hits per
   (sample, GENE)). Plus an optional `aldex.glm` design path
   (`:114-212, 318-326`) for covariate-adjusted models when
   `cfg$differential_abundance$design` is a formula string.
6. **R/09, R/10, R/11 each produce a full PCoA + PERMANOVA + PERMDISP
   suite** via the shared `ge_plot_beta()` helper
   (`utils_ge_profile.R:505-602`). The manuscript only mentions beta
   diversity on taxa.

### Missing — substantive features the manuscript should mention

7. **Two large utility files are invisible.**
   - `R/utils_ge_profile.R` (603 lines) — the shared engine driving
     R/09–R/11 (alpha, beta, Venn, heatmap, RA, total bar, prevalence,
     palettes, renames).
   - `R/utils_prevalence.R` (139 lines) — paired (total | sample-prevalence)
     panel used by R/05, R/09, R/10, R/11.
8. **R/09 drug-class classification + optional MLS rollup**
   (`R/09_resistome.R:66-73`). Default-on: collapses the literal
   `"lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B"`
   to `MLS`, then reduces single-class vs `Multi-drug`.
9. **R/11 replicon-family classification.** Configurable regex
   (`R/11_mobilome.R:65-72`) maps each PlasmidFinder GENE to
   Col-like / IncF / IncX / Other Inc / Unknown-Other; family-level
   relative abundance, total bar, prevalence, and heatmap are emitted
   alongside the gene-level views.
10. **R/10 functional-category extraction.** VFDB categories pulled from
    the ABRicate `PRODUCT` string via the `extract_productFunction`
    helper.
11. **R/06 prefers a precomputed indices CSV.** If
    `cfg$inputs$diversity_csv` is provided (typical wf-metagenomics
    output), R/06 uses it directly; otherwise falls back to vegan via
    `compute_diversity_from_counts` (`:48-70, 110-127`). Selection
    forced via `cfg$alpha_diversity$source = "auto"|"precomputed"|"computed"`.
12. **R/05 emits extra figures** beyond the stacked bar: `species_count.png`
    (log10 horizontal bar) and `species_count_prevalence.png` (paired
    total vs sample-prevalence).
13. **R/03 (TPM) and R/09–R/11 require Bracken.**
    `R/03_normalisation.R:14-18` hard-errors if `cleaned$merged` is NULL;
    downstream gene-element stages are skipped without
    `genetable_normdata.csv`.
14. **`taxid_fixes.csv` is the structured replacement** for the legacy
    row-index taxonomy hacks the manuscript correctly contrasts against
    (`R/01_load_inputs.R:47-53` + `R/02_clean_data.R:58-61`). Columns:
    `taxid,name`.
15. **Quarto report is `dashboard` format** (`templates/report.qmd:3-9`),
    not plain HTML. Multi-page with sidebar, value boxes, per-pair
    sub-tabs, full Differential abundance tab.
16. **R/00 + R/01 always run.** Toggleable stages start at `clean_data`
    (`config_template.yaml:115-127`; `cfg$stages$*` has no `setup` /
    `load_inputs` keys).
17. **Per-stage runtime auto-instrumented** by `run_pipeline.R:43-150`
    (slowest-first breakdown logged to `cfg$outputs$log_file`). Memory
    is not.
18. **Validation aids under `scripts/`**: `diff_normdata.R` (numerical
    diff of `genetable_normdata.csv`), `probe_taxid_parse.R` (taxid-parse
    regression probe), plus `backfill_new_pngs.R`,
    `backfill_violins.R`.

### Ambiguous / caveat

- **R ≥ 4.5** is recorded in `renv.lock` and the README but is not
  enforced by any code-level version check. Current wording is
  defensible.
- **Resistome `database`** accepts a vector
  (`R/09_resistome.R:34` — `tolower(rcfg$database %||% "card")`), so
  `[card, resfinder]` would technically work — but no shipped config
  uses it and no fixture exists. Treat as "configurable but
  not the shipped default" if ResFinder is kept anywhere.
- **Re-centrifuge `by-3` column slicing** in R/02 is self-flagged TODO
  ("known to be too coarse"). Queued for rework in
  `docs/pipeline_rework_scoping.md`. May leak into §3 benchmark numbers.

### Still accurate (no change needed)

- 13 ordered modules `R/00_setup.R` through `R/12_network.R`.
- Single YAML config drives everything; module code is never
  edited per project.
- Stages toggle via `cfg$stages$*` (with the R/00 + R/01 caveat above).
- R/00 (setup), R/01 (input list), R/03 (TPM), R/04 (Venn + heatmap),
  R/05 stacked bar.
- R/02 optional negative-control subtraction (skipped when
  `cfg$metadata$controls` is empty).
- R/02 taxid-parse fix giving ~+110 (sample, GENE) rows
  (exact count still needs re-derivation from a fresh validation run).
- External CLI tools assumed run upstream.
- R ≥ 4.5; renv.

---

## 2. Per-section edit plan

### Abstract

- **Line 22** — change `"tripartite taxa-ARG-MGE network"` →
  `"tripartite sample × taxon × gene network with companion chord and
  Sankey diagrams"`.
- Optional: add a clause about per-domain PCoA/PERMANOVA
  (scope-broadening for resistome/virulome/mobilome).

### §2 Implementation — stage paragraph (line 46, main rewrite target)

- **R/02**: drop `R/utils_taxa.R` from the R/02 sentence; add Bracken
  merge + metadata-column attach + `taxid_fixes.csv` mechanism
  (the structured replacement for legacy row-index edits).
- **R/06**: `"richness, Shannon (and Simpson when an upstream
  wf-metagenomics indices CSV is supplied), with Kruskal-Wallis tests
  and violin/bar plots."`
- **R/07**: expand to `"PCoA on Hellinger-transformed Bray-Curtis
  distances using vegan, with PERMANOVA (adonis2, 9999 permutations)
  and PERMDISP overlaid."`
- **R/08**: state both taxa- and gene-level DAA; mention omnibus KW
  + pairwise Welch's t and the optional `aldex.glm` design mode.
- **R/09**: drop "ResFinder" (or downgrade to "optionally"); add
  drug-class classification + MLS rollup + per-domain PCoA suite.
- **R/10**: add VFDB functional-category extraction + per-domain
  PCoA suite.
- **R/11**: add replicon-family classification + per-domain
  PCoA suite.
- **R/12**: rewrite to include all three outputs (network + chord
  + Sankey) plus Gephi CSVs, igraph topology table (degree,
  betweenness, Louvain modularity), sample clusters, degree
  distribution.
- **New sentence** after the stage paragraph introducing the three
  utility files: `R/utils_taxa.R` (configurable name cleanup),
  `R/utils_ge_profile.R` (shared engine for R/09–R/11),
  `R/utils_prevalence.R` (paired total/prevalence panels).

### §2 Implementation — configuration paragraph (line 44)

- Add parenthetical: R/00 setup and R/01 input loading always run;
  toggleable stages start at `clean_data`.
- Add: `"TPM (R/03) and all gene-element stages (R/09–R/11) require
  Bracken outputs."`

### §2 Implementation — reporting sentence (line 48)

- Mention Quarto `dashboard` format (multi-page with sidebar +
  value boxes), not plain HTML.

### §3 Application and Validation

- Mention per-stage runtime is auto-logged by `run_pipeline.R`
  (slowest-first breakdown) so that column of the benchmark table
  is mechanically populated; memory still needs separate sourcing.
- Mention `scripts/diff_normdata.R` and `scripts/probe_taxid_parse.R`
  as the validation toolkit alongside the hand-edited script
  comparison.
- Re-flag the ~110 (sample, GENE) row delta as needing fresh-run
  confirmation before submission.
- Optional: brief disclosure of the still-coarse Re-centrifuge
  column slicing in R/02 (queued for rework).

### §4 Conclusion

- Mirror the network/chord/Sankey correction for consistency, or
  leave the existing one-line "integrated network view" phrasing
  for brevity.

### Figure 1 (`figure1_pipeline_schematic.mmd`)

- Drop the "ResFinder" label from the resistome node.
- Update the R/12 node from "Taxa × ARG × MGE network" to reflect
  network + chord + Sankey.
- Re-render the PNG via the same `mmdc` command in
  `applications_note.md:100` after editing the `.mmd` source.

### References

No new tools strictly required, but consider adding:

- **networkD3** — if Sankey gets a sentence.
- **igraph** — if topology table (degree, betweenness, Louvain
  modularity) gets a sentence.
- **ggraph** — if the tripartite network layout is named explicitly.
- **pheatmap** — listed in `pipeline_facts.md` stack but not currently
  cited; add if heatmap wording is expanded.

---

## 3. Key file:line citations (absolute paths)

For verification while editing:

- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/run_pipeline.R`
  — entry point, stage dispatcher, runtime instrumentation
  (`:20-22, 28, 43-150`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/00_setup.R`
  — `load_config()` validation (`:4-43`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/01_load_inputs.R`
  — five inputs + `taxid_fixes` CSV (`:7-45, 47-53`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/02_clean_data.R`
  — taxid regex (`:14-26`), Re-centrifuge by-3 slicing (`:31-52`),
  control subtraction (`:65-101`), ABRicate ↔ taxa link (`:111-117`),
  Bracken merge + metadata attach (`:119-175`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/03_normalisation.R`
  — TPM (`:11-63`), Bracken hard-error (`:14-18`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/04_taxonomy.R`
  — annotated heatmap (`:194-211`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/05_relative_abundance.R`
  — stacked bar + PNG/HTML (`:169-184`), count + prevalence figures
  (`:186-244`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/06_alpha_diversity.R`
  — precomputed-CSV path (`:48-70, 110-127`), vegan fallback
  (`:73-87`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/07_beta_diversity.R`
  — vegan PCoA (`:111-115`), PERMANOVA + PERMDISP (`:54-216`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/08_differential_abundance.R`
  — aldex.glm path (`:114-212, 318-326`), KW + pairwise
  (`:337-356, 357-401`), two-level loop (`:715, 738-755`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/09_resistome.R`
  — CARD default (`:34, 40`), MLS rollup (`:66-73`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/10_virulome.R`
  — VFDB default (`:32, 38`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/11_mobilome.R`
  — PlasmidFinder default (`:39, 45`), replicon families
  (`:65-72, 180-209`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/12_network.R`
  — network (`:27-71, 316-367`), chord (`:373-507`), Sankey
  (`:526-722`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/utils_taxa.R`
  — `clean_taxa_names()` (89 lines).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/utils_ge_profile.R`
  — shared gene-element engine (603 lines); beta helper (`:505-602`),
  alpha (`:88-107`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/R/utils_prevalence.R`
  — `build_count_prevalence()` / `save_count_prevalence()` (139 lines).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/templates/report.qmd`
  — Quarto dashboard (`:3-9`), DAA tab (`:628-704`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/config/config_template.yaml`
  — stage toggles (`:115-127`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/renv.lock`
  — R 4.5.0 (`:3`).
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/scripts/`
  — `diff_normdata.R`, `probe_taxid_parse.R`,
  `backfill_new_pngs.R`, `backfill_violins.R`.
- `C:/Users/julio/Desktop/Metagenomics_pipeline_automation/docs/pipeline_rework_scoping.md`
  — Re-centrifuge column-slicing rework note.
