# Conversation log

Chronological record of the design decisions made while building this
automation. Each entry is timestamped. Future conversations append here.

---

## 2026-04-17 — Initial scoping

**User goal**
Automate the long-read metagenomics pipeline that Julio runs for multiple
microbiome studies. The pipeline is currently ~30–50 hand-edited R scripts
per study, with hardcoded paths, sample IDs and metadata column names copied
between projects. The only thing that really changes between studies is the
metadata file.

**Studies surveyed (`Metagenomics_Projects_paths.txt`)**
- Chicken gut microbiome, batch 1 — `PC_JC_2024-11-29_Julio_Gallus` *(completed, paper under review — reference project)*
- Chicken gut microbiome, batch 2 — `JC_FL_2025_05_20_CK_microbiome`
- Hospital microbiome — `Hospital_microbiome`
- Lung microbiome — `Lung_microbiome` *(completed, thesis chapter — secondary reference)*

**Main sources of manual work observed**
1. Hardcoded absolute paths in every R script (some Lung paths were even
   copy-pasted into the Chicken scripts).
2. Sample IDs hardcoded inside `select(...)` calls, e.g. chicken batch 2
   `cleanData.R` lines 70–74.
3. Manual taxid → name fixes by row index
   (`noncontaminants_list[26, 4] <- "Escherichia coli O157"`, repeated ×20+).
4. Metadata column names (`Treatment`, `Weight`) hardcoded into every
   plotting / stats script.
5. Scripts are re-run one at a time in RStudio for each new study.

**Methodology extracted from Chicken batch 1 paper**
(see `reference_methodology/chicken_batch1_methods.md` for the long version)

| Block | Tool / method |
|---|---|
| Host removal       | Samtools v1.17 |
| Taxonomy           | Kraken2 v2.1.2 (PlusPF-8) + Bracken v2.7.0 |
| AMR / VF / MGE     | ABRicate v1.0.1 — CARD / VFDB / PlasmidFinder |
| Decontamination    | Re-centrifuge |
| Normalisation      | TPM per gene element |
| Alpha diversity    | Richness + Shannon, Kruskal–Wallis (block as random effect) |
| Beta diversity     | Bray–Curtis → PCoA, PERMANOVA `adonis2` 9999 perms |
| Differential abund.| ALDEx2 — 128 MC Dirichlet, prevalence ≥ 2 |
| Network            | Tripartite Treatment × taxa × GE; igraph / Gephi |

---

## 2026-04-17 — Decisions

- **Reference study:** Chicken batch 1 (requested: "you can use chicken
  batch 1 or lung microbiome" — both completed). Lung config stubbed so the
  config format can be validated against a cohort-style metadata schema
  alongside the randomised-block chicken design.
- **Config format:** single `config.yaml` per study, placed under
  `projects/<study_id>/config.yaml`. Schema documented in
  `config/config_template.yaml`.
- **Master runner:** `run_pipeline.R` sources every `R/NN_*.R` module in
  order, reads the config, and dispatches by the `stages:` map. One command
  per study.
- **Manual taxid edits replaced** by `Metadata/taxid_fixes.csv` (taxid → name)
  joined into the noncontaminants table — removes the brittle row-index fixes.
- **Report format:** Quarto HTML (`templates/report.qmd`, rendered per study
  to `Reports/pipeline_report.html`). Chosen over plain R Markdown because
  Quarto bundles pandoc and is the path Posit recommends going forward.
- **Output policy:** study `Datasets/` / `Figures/` / `Reports/` are
  git-ignored (large, regenerable). Only configs, code, and small reference
  files (methodology notes, taxid fix tables) are versioned.

---

## 2026-04-17 — Repo scaffold written

Files created this session:

```
README.md
.gitignore
run_pipeline.R
config/config_template.yaml
reference_methodology/chicken_batch1_methods.md
templates/report.qmd
projects/chicken_batch1/config.yaml
projects/lung_microbiome/config.yaml
R/00_setup.R
R/01_load_inputs.R
R/02_clean_data.R
R/03_normalisation.R
R/04_taxonomy.R
R/05_relative_abundance.R
R/06_alpha_diversity.R
R/07_beta_diversity.R
R/08_differential_abundance.R
R/09_resistome.R
R/10_virulome.R
R/11_mobilome.R
R/12_network.R
docs/conversation_log.md
```

---

## 2026-04-17 — User preferences captured

- **GitHub visibility:** Private.
- **HTML reporting:** Quarto — user will install from quarto.org.
  `run_pipeline.R` already detects Quarto on PATH and skips the report
  stage gracefully if it's missing.

## 2026-04-17 — Remote pushed

- Remote: <https://github.com/Julio92-C/Metagenomics_pipeline_automation.git>
  (name kept as `Metagenomics_pipeline_automation`, matching the local dir).
- User created an empty repo with an auto-generated README and LICENSE.
  Local history was rebased onto `origin/main`; the README conflict was
  resolved by keeping the detailed local README. The GitHub-supplied
  LICENSE is now in the repo.
- Branch `main` tracks `origin/main`; initial scaffold pushed successfully.

## Open items

- User to install Quarto, then the report stage will run automatically on
  the next `Rscript run_pipeline.R projects/chicken_batch1/config.yaml`.
- First validation run against Chicken batch 1 — compare regenerated figures
  against the published figures in `Figures/` before declaring parity.
- Populate `Metadata/taxid_fixes.csv` for Chicken batch 1 by extracting the
  hardcoded fixes from the original `relativeAbundance.R`.
- Install R packages `{yaml, vegan, ape, igraph, pheatmap, VennDiagram,
  paletteer, ALDEx2, plotly, htmlwidgets, rprojroot}` via `renv`.
