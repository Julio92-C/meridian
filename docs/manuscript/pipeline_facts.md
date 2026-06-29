# Pipeline facts — for the Methods section

Pre-extracted so the drafting session does not need to read every
R module. If a claim here needs deeper verification, the source files
are under `../R/` (one level up from this directory).

## What the tool is

**MERIDIAN** (*Metagenomic Evaluation of Resistance, Identity & Diversity
through Integrated Analysis of Nanopore sequencing*) — a config-driven,
modular R workflow that takes Kraken2 + Bracken + ABRicate + Re-centrifuge
outputs and produces a publication-ready HTML report. Targets long-read
(nanopore / PacBio) shotgun metagenomics studies.

## Why it exists (motivation)

Running a metagenomics study typically means hand-editing 30–50 R scripts
per project: absolute paths, sample IDs pasted into `select()` calls,
manual row-index taxonomy fixes, copy-pasted chunks between studies.
This pipeline replaces that workflow with a single `config.yaml` per
study and one command (`Rscript run_pipeline.R config.yaml`) that runs
every analysis block.

## Design principles

- **One config per study.** Everything that varies (metadata file, group
  columns, controls, filter thresholds, sample IDs) lives in YAML. Module
  code never changes.
- **Modular stages.** Each analysis block is its own R module; any stage
  can be toggled on/off via `cfg$stages$*` flags.
- **Reproducible methodology.** Mirrors methods used in peer-reviewed
  long-read metagenomics studies.
- **Per-project flexibility.** Modules adapt to per-project input shape;
  nothing is hardcoded.

## Stage list (`R/NN_*.R`)

| Stage | Module | Purpose |
|-------|--------|---------|
| 00 | `00_setup.R` | Library loading, helpers, config validation |
| 01 | `01_load_inputs.R` | Read Kraken2, Bracken, ABRicate, Re-centrifuge, metadata |
| 02 | `02_clean_data.R` | Decontaminate (optional negative-control subtraction), apply taxid fixes, link ABRicate ↔ taxa |
| 03 | `03_normalisation.R` | TPM normalisation |
| 04 | `04_taxonomy.R` | Venn diagrams + heatmaps of taxa across treatment groups |
| 05 | `05_relative_abundance.R` | Stacked-bar relative abundance per sample (PNG + plotly HTML) |
| 06 | `06_alpha_diversity.R` | Shannon / Simpson / richness with violin/bar plots |
| 07 | `07_beta_diversity.R` | PCoA |
| 08 | `08_differential_abundance.R` | ALDEx2 two-level (taxa + gene): omnibus Kruskal-Wallis + pairwise Welch's t; optional `aldex.glm` design path for covariate-adjusted models |
| 09 | `09_resistome.R` | AMR gene profiling (ABRicate + CARD); drug-class classification with optional MLS rollup; per-domain alpha/Venn/heatmap/PCoA + PERMANOVA suite |
| 10 | `10_virulome.R` | Virulence factor profiling (VFDB); function category extracted from ABRicate `PRODUCT`; per-domain alpha/Venn/heatmap/PCoA + PERMANOVA suite |
| 11 | `11_mobilome.R` | Mobile genetic element profiling (PlasmidFinder); configurable replicon-family classification (Col-like, IncF, IncX, other Inc); per-domain alpha/Venn/heatmap/PCoA + PERMANOVA suite (gene- and family-level) |
| 12 | `12_network.R` | Tripartite sample × taxon × gene ggraph network + chord (overall + per-group) + 4-tier networkD3 Sankey; Gephi-compatible node/edge CSVs; igraph topology (degree, betweenness, Louvain modularity); Bray-Curtis sample clusters |
| 14 | `14_manifest.R` | Emits `manifest.json` describing every artefact, table schema and stage status (complete / skipped / failed). Stable JSON contract consumed by the downstream `metaomics-scribe` manuscript-drafting agent |

Plus three `R/utils_*` helpers:

- `utils_taxa.R` — shared taxa-name cleanup (used by R/04, R/05, R/07, R/12)
- `utils_ge_profile.R` — shared engine driving R/09–R/11 (alpha, beta, Venn, heatmap, relative abundance, total bar, prevalence, palettes, renames)
- `utils_prevalence.R` — paired total / sample-prevalence panel helper (R/05, R/09–R/11)

A Quarto template at `templates/report.qmd` is built as a `dashboard` format
(multi-page, sidebar, value boxes, per-domain panels with download links)
and is rendered at the end of the run.

`run_pipeline.R` auto-instruments per-stage wall time and logs a
slowest-first breakdown to `cfg$outputs$log_file` at the end of each run;
memory is not instrumented.

## Tool stack (cite these)

- **R** ≥ 4.5
- **Kraken2** — taxonomic classification
- **Bracken** — abundance re-estimation from Kraken2
- **ABRicate** — gene screening (CARD, VFDB, PlasmidFinder)
- **Re-centrifuge** — contamination filtering
- **ALDEx2** — differential abundance on compositional data
- **vegan** — alpha / beta diversity
- **pheatmap**, **ggplot2**, **VennDiagram**, **paletteer**, **plotly**
- **Quarto** — HTML reporting
- **renv** — reproducible R environment

## Current state

- **Phase 1 / alpha.** Scaffold runs end-to-end against real data after a
  recent smoke test. Some modules still being validated against
  ground-truth outputs.
- **Validation strategy.** Output compared against ~46 hand-edited R
  scripts from the chicken_batch1 reference dataset (PC_JC_2024-11-29
  Julio Gallus project). For each stage, the pipeline reproduces the
  known-good figures and tables to within expected differences
  (e.g. the pipeline fixes a taxid-parse bug present in the GT scripts,
  so resistome GENE counts differ by ~110 (sample, GENE) rows).

## Availability

GitHub repo (license: MIT). Install:

```bash
git clone <repo>
cd <repo>
Rscript -e 'renv::restore()'
Rscript run_pipeline.R path/to/your/config.yaml
```

No `install.packages("...")` — this is a workflow, not a package.
