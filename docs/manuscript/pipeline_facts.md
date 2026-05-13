# Pipeline facts — for the Methods section

Pre-extracted so the drafting session does not need to read every
R module. If a claim here needs deeper verification, the source files
are under `../R/` (one level up from this directory).

## What the tool is

**Metagenomics Pipeline Automation** — a config-driven, modular R workflow
that takes Kraken2 + Bracken + ABRicate + Re-centrifuge outputs and
produces a publication-ready HTML report. Targets long-read
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
| 08 | `08_differential_abundance.R` | ALDEx2 — overall and pairwise |
| 09 | `09_resistome.R` | AMR gene profiling (ABRicate + CARD/ResFinder) |
| 10 | `10_virulome.R` | Virulence factor profiling (VFDB) |
| 11 | `11_mobilome.R` | Mobile genetic element profiling |
| 12 | `12_network.R` | Tripartite chord diagram / network of taxa × ARGs × MGEs |

Plus `R/utils_taxa.R` — shared taxa-name cleanup helper.

A Quarto template at `templates/report.qmd` renders all outputs into a
single HTML report at the end.

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
