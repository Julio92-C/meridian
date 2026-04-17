# Metagenomics Pipeline Automation

Automation of the long-read metagenomics analysis pipeline used across microbiome
studies by Julio Cesar Ortega Cambara (UWL). Each study currently requires
hand-editing of ~30-50 R scripts with hardcoded paths, sample lists and
metadata column names. This repo replaces that workflow with a single
config-driven runner plus a Quarto HTML report.

## Reference studies

The pipeline is derived from two completed projects (papers under review):

- **Chicken batch 1** — *Long-reads metagenomics elucidates the effects of Dulse
  supplementation on the poultry caecal bacteriome and its associated genetic
  repertoire.* (Frontiers submission draft, November 2024.)
- **Lung microbiome** — EBSP/UCLH cohort, PhD thesis chapter.

Additional studies that will reuse the pipeline: Chicken batch 2
(`JC_FL_2025_05_20_CK_microbiome`) and Hospital microbiome.

## Pipeline stages

Ordered by the runner, grouped as the user requested:

| # | Stage               | Tool / method                                                 |
|---|---------------------|---------------------------------------------------------------|
| 0 | Load inputs         | Kraken2, Bracken, ABRicate, Re-centrifuge (raw reports)       |
| 1 | Clean & decontaminate | Merge reports, remove controls/crossover taxa, taxid fixes   |
| 2 | Normalisation       | TPM for gene elements (GEs)                                   |
| 3 | Taxonomy            | Venn diagrams, pHeatmap, relative-abundance table             |
| 4 | Relative abundance  | Stacked bar plots at phylum / class / genus / species         |
| 5 | Alpha diversity     | Richness + Shannon; Kruskal–Wallis (fixed/random effects)     |
| 6 | Beta diversity      | Bray–Curtis → PCoA; PERMANOVA (`adonis2`, 9999 permutations)  |
| 7 | Differential abundance | ALDEx2 (128 Monte-Carlo Dirichlet instances, prevalence ≥ 2)|
| 8 | Resistome           | ABRicate + CARD → ARG profile, Venn, pHeatmap, KW stats       |
| 9 | Virulome            | ABRicate + VFDB → VF profile, Sankey, pHeatmap                |
|10 | Mobilome            | ABRicate + PlasmidFinder → MGE profile, pHeatmap              |
|11 | Network analysis    | Tripartite (Treatment × taxa × GE) + chord diagrams (Gephi)   |
|12 | Report              | Quarto → single HTML per study                                |

The methodology matches §2 of the Chicken batch 1 manuscript. See
`reference_methodology/chicken_batch1_methods.md` for direct quotes and tool
versions.

## Repo layout

```
.
├── R/                  function library (sourced by run_pipeline.R)
├── config/             config schema + template
├── projects/           per-study config files + taxid fix tables
├── templates/          report.qmd + scaffolding for new studies
├── reference_methodology/   extracted methods from the reference papers
├── docs/               conversation log, decisions, diagrams
└── run_pipeline.R      master entry point
```

## Running the pipeline for a study

```bash
Rscript run_pipeline.R projects/chicken_batch1/config.yaml
```

This reads the config, sources `R/`, runs every stage, writes cleaned CSVs
into the study's `Datasets/`, writes figures into `Figures/`, then renders
`templates/report.qmd` to `<study>/Reports/pipeline_report.html`.

## Status

See `docs/conversation_log.md` for ongoing decisions. The scaffold is being
built incrementally from the two reference studies. Current state: directory
skeleton, config schema, methodology notes, and empty module stubs.
