<h1 align="center">🧬 Metagenomics Pipeline Automation</h1>

<p align="center">
  <b>A config-driven, reproducible analysis pipeline for long-read metagenomics studies.</b><br/>
  From raw Kraken2 / Bracken / ABRicate reports to a publication-ready HTML report — in one command.
</p>

<p align="center">
  <img alt="R" src="https://img.shields.io/badge/R-%E2%89%A5%204.5-276DC3?logo=r&logoColor=white">
  <img alt="Quarto" src="https://img.shields.io/badge/Quarto-HTML%20report-74AADB?logo=quarto&logoColor=white">
  <img alt="Platform" src="https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey">
  <img alt="License" src="https://img.shields.io/badge/license-MIT-blue">
  <img alt="Status" src="https://img.shields.io/badge/status-alpha-orange">
</p>

---

## ✨ Why this exists

Running a metagenomics study shouldn't mean hand-editing 30–50 R scripts
every time the metadata changes. Absolute paths, sample IDs pasted into
`select()` calls, manual row-index taxonomy fixes, copy-pasted chunks
between projects — it adds up to days of error-prone work per study.

**Metagenomics Pipeline Automation** replaces that workflow with a single
`config.yaml` per study and one command to run every analysis block, from
decontamination to the final HTML report.

---

## 🚀 Features

- 🧾 **One config per study.** Everything that varies between studies
  (metadata file, grouping columns, controls, filter thresholds) lives in a
  single YAML file. The code never changes.
- 🧩 **Modular stages.** Each analysis block is its own R module; toggle
  any stage on/off with a flag.
- 🛡️ **Reproducible methodology.** Mirrors the methods used in peer-reviewed
  long-read metagenomics studies (Kraken2 + Bracken + ABRicate +
  Re-centrifuge + ALDEx2 + vegan).
- 📊 **Publication-ready figures.** Violin plots, PCoA ordinations, Venn
  diagrams, pHeatmaps, Sankey and chord diagrams, tripartite networks.
- 🖥️ **Single-command HTML report.** Quarto renders everything into a
  self-contained `pipeline_report.html`.
- 🔍 **Transparent logs.** Every run writes a timestamped log with the
  statistics produced at each stage.

---

## 🧪 Analysis blocks

| # | Block                      | Icon | Tool / method                                                     |
|---|----------------------------|:----:|-------------------------------------------------------------------|
| 0 | Load inputs                | 📥  | Kraken2, Bracken, ABRicate, Re-centrifuge raw reports             |
| 1 | Clean & decontaminate      | 🧼  | Merge reports, control subtraction, taxid fixes                   |
| 2 | Normalisation              | ⚖️  | TPM per gene element (GE)                                         |
| 3 | Taxonomy                   | 🌿  | Venn + pHeatmap across treatment groups                           |
| 4 | Relative abundance         | 📊  | Stacked bar plots (phylum → species) + interactive plotly HTML    |
| 5 | Alpha diversity            | 🎻  | Richness, Shannon; Kruskal–Wallis, violin plots                   |
| 6 | Beta diversity             | 🧭  | Bray–Curtis → PCoA; PERMANOVA (`adonis2`, 9999 permutations)      |
| 7 | Differential abundance     | 🧮  | ALDEx2 — 128 MC Dirichlet instances, CLR pairwise tests           |
| 8 | Resistome                  | 💊  | ABRicate + CARD → ARG profile                                     |
| 9 | Virulome                   | 🦠  | ABRicate + VFDB → virulence factor profile                        |
| 10| Mobilome                   | 🧬  | ABRicate + PlasmidFinder → MGE profile                            |
| 11| Network analysis           | 🕸️  | Tripartite Treatment × taxa × GE; igraph + Gephi-ready edge lists |
| 12| Report                     | 📄  | Quarto → single self-contained HTML per study                     |

> Methodology reproduced from a peer-review-stage long-read metagenomics
> manuscript. See [`reference_methodology/chicken_batch1_methods.md`](reference_methodology/chicken_batch1_methods.md)
> for tool versions, parameters and direct quotes.

---

## 📦 Installation

**Prerequisites**

- R ≥ 4.5
- [Quarto](https://quarto.org/docs/get-started/) ≥ 1.4
- Git

**Clone and install R dependencies**

```bash
git clone https://github.com/Julio92-C/Metagenomics_pipeline_automation.git
cd Metagenomics_pipeline_automation

Rscript -e 'install.packages(c(
  "yaml", "readr", "dplyr", "tidyr", "purrr", "stringr",
  "ggplot2", "ggpubr", "paletteer", "plotly", "htmlwidgets",
  "vegan", "ape", "igraph", "pheatmap", "VennDiagram",
  "tibble", "rprojroot"
))'

# ALDEx2 is on Bioconductor
Rscript -e 'if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install("ALDEx2")'
```

---

## ⚡ Quickstart

1. **Copy the config template** and point it at your study:

   ```bash
   cp config/config_template.yaml projects/my_study/config.yaml
   ```

2. **Edit `projects/my_study/config.yaml`** — set `project_root`, the
   metadata file path, `sample_id_col`, `group_cols` and your control
   sample IDs. Nothing else needs to change.

3. **Run the pipeline:**

   ```bash
   Rscript run_pipeline.R projects/my_study/config.yaml
   ```

4. **Open the report** at `<project_root>/Reports/pipeline_report.html`. 🎉

Each stage can be skipped by flipping a flag under `stages:` in the config.

---

## 🗂️ Repository layout

```
.
├── 📁 R/                       function library (sourced by run_pipeline.R)
│   ├── 00_setup.R              paths, logging, config loader
│   ├── 01_load_inputs.R        Kraken2 / Bracken / ABRicate / Re-centrifuge
│   ├── 02_clean_data.R         merge + decontamination
│   ├── 03_normalisation.R      TPM
│   ├── 04_taxonomy.R           Venn + heatmap
│   ├── 05_relative_abundance.R stacked bars + plotly
│   ├── 06_alpha_diversity.R    richness, Shannon, Kruskal–Wallis
│   ├── 07_beta_diversity.R     Bray–Curtis, PCoA, PERMANOVA
│   ├── 08_differential_abundance.R   ALDEx2
│   ├── 09_resistome.R          CARD / ARGs
│   ├── 10_virulome.R           VFDB / VFs
│   ├── 11_mobilome.R           PlasmidFinder / MGEs
│   └── 12_network.R            tripartite + topology metrics
├── 📁 config/                  schema + annotated template
├── 📁 projects/                one subfolder per study, each with config.yaml
├── 📁 templates/               Quarto report template
├── 📁 reference_methodology/   methods extracted from reference papers
├── 📁 docs/                    conversation log, design decisions
├── 📝 run_pipeline.R           master entry point
└── 📄 README.md
```

---

## 🧠 How it replaces manual work

| Manual pattern in old scripts                           | Replaced by                                          |
|---------------------------------------------------------|------------------------------------------------------|
| `"C:/Users/.../Study_X/Reports/..."` hardcoded paths    | `project_root` + relative paths in `config.yaml`     |
| `select("PO4B4", "PO5B4", ...)` sample ID lists         | Sample IDs read dynamically from metadata            |
| `noncontaminants_list[26, 4] <- "Escherichia coli O157"`| `Metadata/taxid_fixes.csv` lookup table              |
| Re-running 30+ scripts by hand                          | One `Rscript run_pipeline.R <config>` command        |
| Paper-ready figure assembly in PowerPoint               | Quarto HTML report with every figure embedded        |

---

## 🤝 Contributing

Contributions are very welcome — bug reports, methodology improvements,
support for new databases, extra analysis blocks.

1. Fork the repo and create a feature branch (`git checkout -b feature/awesome`).
2. Keep each R module self-contained and driven by `cfg`.
3. Add a short section to `docs/conversation_log.md` describing the design
   choice (the repo keeps a running design journal, not just code).
4. Open a pull request with a short demo config that shows the change.

For larger changes (new pipeline stage, breaking config schema change),
please open an issue first so we can discuss the design.

---

## 📚 Citation

If you use this pipeline in your research, please cite the upstream
methodology papers:

- Wood, D. E., Lu, J. & Langmead, B. (2019). *Improved metagenomic analysis
  with Kraken 2.* **Genome Biology**, 20, 257.
- Lu, J. *et al.* (2022). *Metagenome analysis using the Kraken software
  suite.* **Nature Protocols**, 17, 2815–2839.
- Fernandes, A. D. *et al.* (2014). *Unifying the analysis of high-throughput
  sequencing datasets… ALDEx2.* **Microbiome**, 2, 15.
- Oksanen, J. *et al.* **vegan: Community Ecology Package** (R package).

Citation details for the reference manuscripts will be added once the
papers are accepted.

---

## 👤 Author & contact

**Julio Cesar Ortega Cambara**
PhD Researcher — University of West London
<juliocesar921016@gmail.com> · [GitHub @Julio92-C](https://github.com/Julio92-C)

---

## 📜 License

MIT — see [`LICENSE`](LICENSE).

---

<p align="center">
  Made with ☕, 🧬, and a lot of <code>Rscript</code>.
</p>
