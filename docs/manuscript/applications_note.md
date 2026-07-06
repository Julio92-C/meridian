# MERIDIAN: a config-driven R workflow for long-read shotgun metagenomics analysis

**Authors / Affiliations**

Julio C. Ortega Cambara<sup>1</sup>, Piotr Cuber<sup>4</sup>, Pedro Humberto Lebre<sup>1</sup>, Raju Misra<sup>2</sup>, Timothy D. McHugh<sup>3</sup>, Sylvia Rofael<sup>3</sup>, David Lowe<sup>3</sup>, Hermine Mkrtchyan<sup>1*</sup>

<sup>1</sup>School of Biomedical Science, University of West London, London, United Kingdom
<sup>2</sup>UK Health Security Agency, London, United Kingdom
<sup>3</sup>Centre for Clinical Microbiology, University College London, London, United Kingdom
<sup>4</sup>Natural History Museum, London, United Kingdom

<sup>*</sup>Correspondence: Hermine Mkrtchyan — Hermine.Mkrtchyan@uwl.ac.uk

{{TODO: ORCIDs for all authors}}

---

## Abstract

**Motivation.** Long-read shotgun metagenomics studies are typically analysed by hand-editing dozens of R scripts per project — pasting sample identifiers, patching taxonomy strings and copying chunks between studies. Existing R workflows largely target short-read amplicon data or a single analysis (e.g. resistome risk scoring), leaving long-read shotgun studies without an integrated, configuration-driven option.

**Results.** We present **MERIDIAN** (*Metagenomic Evaluation of Resistance, Identity & Diversity through Integrated Analysis of Nanopore sequencing*), a modular R workflow that consumes bioinformatics tools outputs and produces a single Quarto HTML dashboard, publication-ready composite figures and a supplementary tables workbook. The workflow covers decontamination, TPM normalisation, taxonomy, alpha/beta diversity with PERMANOVA, two-level ALDEx2 differential abundance, resistome, virulome, and mobilome profiling, and a tripartite sample × taxon × gene network reconstruction. A study is parameterised by a single YAML file; module code is not edited per project, and every run emits a `manifest.json` artefact catalogue for downstream integration.

**Availability and Implementation.** Source at https://github.com/Julio92-C/Metagenomics_pipeline_automation (MIT). This note describes release **v0.1.0** (https://github.com/Julio92-C/Metagenomics_pipeline_automation/releases/tag/v0.1.0). Run with `git clone --branch v0.1.0 <repo> && Rscript run_pipeline.R config.yaml` after `renv::restore()`.

**Contact.** Hermine.Mkrtchyan@uwl.ac.uk

---

## 1 Introduction

Long-read shotgun metagenomics (Oxford Nanopore, PacBio) is increasingly central to One Health surveillance, where threats intersect humans, animals and ecosystems [1]. By profiling microbial communities and their functional genetic elements in a single experiment, it links taxonomic identity to determinants of resistance, virulence, and mobility across environmental, clinical, and livestock matrices. The downstream analysis stack is well established — Kraken2 [2] and Bracken [3] for taxonomic assignment and abundance estimation, ABRicate [4] for screening antimicrobial resistance genes (ARGs), virulence factors (VFs) and mobile genetic elements (MGEs) against CARD [5], VFDB [6] and PlasmidFinder [7], Re-centrifuge [8] for contamination assessment, and R packages such as vegan [9] and ALDEx2 [10] for ecological and differential analyses — but the integration between them is not.

In practice, each new study is analysed by adapting a set of bespoke R scripts: hard-coded input paths, sample-identifier filters embedded as literal strings, per-study negative-control renames, and taxonomy fixes patched by row index. Consequently, individual script copies diverge, and the lack of a shared interface prevents automated output logging. Comparing different studies requires manual reconciliation of both code and data. The reference cohort for the present work (`Chicken Batch 1`; ~46 hand-edited scripts) typifies this pattern. The L-ARRAP pipeline [11] addresses one slice of this problem for long-read resistome risk assessment, and the rANOMALY workflow [12] demonstrates the value of a config-driven R workflow with auto-generated reports for amplicon data. To our knowledge, no comparable workflow exists for end-to-end long-read shotgun metagenomics that covers taxonomy, diversity, differential abundance, and full resistome–virulome–mobilome profiling in a single configurable run.

We present an alpha-stage R workflow that fills this gap. The tool is in active development; this note describes its design, current scope, and the real-data validation strategy we are using to harden each stage.

## 2 Implementation

The workflow is a stand-alone R project, not an installable package. A single entry point (`run_pipeline.R`) consumes one YAML configuration file per study and dispatches to fifteen modules (`R/00_setup.R` through `R/13_panels.R`; plus the manifest stage `R/14` which runs last); a Quarto dashboard template (`templates/report.qmd`) renders all outputs into a single multi-page HTML deliverable, a publication-panels stage assembles main and supplementary composite figures plus a supplementary tables workbook, and a final manifest stage emits `manifest.json` describing every artefact produced (**Figure 1**).

![**Figure 1.** MERIDIAN workflow. A single `config.yaml` drives fifteen R modules (R/00–R/14) consuming Kraken2/Bracken, ABRicate, Re-centrifuge and sample metadata. After setup and decontamination (R/00–R/02), community analyses (R/04–R/08) and functional profiling against CARD, VFDB and PlasmidFinder (TPM in R/03, per-domain analyses in R/09–R/11) converge in a tripartite sample × taxon × gene network (R/12). R/12 fans out to a Quarto dashboard and a publication-panels stage (R/13); R/14 finally writes `manifest.json` (contract v1.2), a kind-tagged catalogue of every artefact. Source: `generate_figure1.py`.](figure1_pipeline_schematic.png)

**Configuration-driven.** Everything that varies between studies — metadata file path, grouping columns, negative-control identifiers, abundance thresholds, sample selection — is declared in `config.yaml`. Module code is not edited per-project; dataset-specific adjustments live in `config.yaml` or in the structured fix tables it references. R/00 (setup) and R/01 (input loading) run by default; every subsequent stage can be toggled via `cfg$stages$*` flags, e.g, users can re-run only the resistome modules after revising thresholds. TPM normalisation (R/03) and all gene-element stages (R/09–R/11) require ABRicate outputs and are skipped otherwise. Per-study NCBI taxonomy fixes, when needed, are declared in a structured `taxid_fixes.csv` rather than as row-index edits. A single `cfg$stats$padjust_method` knob (Benjamini-Hochberg is used by default; full list in `docs/STATISTICS.md`) routes every multiple-testing correction across the pipeline.

**Stages.** The pipeline organises its work into five conceptual blocks.

(1) *Ingest and cleaning* (R/00–R/02) load upstream tables — Kraken2/Bracken, ABRicate, Re-centrifuge contamination tags, study metadata — accommodating per-project input shapes rather than assuming a fixed column layout, apply optional negative-control subtraction, merge Bracken counts with the metadata covariates declared in `cfg$metadata$fixed_effects` / `random_effect`, and link ABRicate calls to taxonomic assignments. Per-study taxonomy patches are applied from `taxid_fixes.csv` when supplied.

(2) *Gene-element normalisation* (R/03) computes per-(sample, GENE) TPM from ABRicate hits and gene lengths.

(3) *Community analyses* (R/04–R/08) cover the standard taxonomic toolkit. R/04–R/05 produce Venn diagrams and annotated heatmaps after a shared taxa-name cleanup, plus top-N stacked-bar relative-abundance plots (static PNG and interactive plotly HTML) with paired total-count and sample-prevalence panels. Alpha diversity (R/06) reports Shannon, Simpson, richness, Pielou, and Fisher metrics with Kruskal-Wallis tests and violin/bar plots, preferring a pre-computed indices CSV (e.g., Oxford Nanopore wf-metagenomics output) and falling back to the vegan R package. Beta diversity (R/07) computes PCoA on Hellinger-transformed Bray-Curtis distances with PERMANOVA (`adonis2`, 9999 permutations) and PERMDISP (`betadisper`) overlaid on the ordination. ALDEx2 differential abundance (R/08) runs at two levels — taxa (Bracken counts) and genes (ABRicate contig hits per (sample, GENE)), combining an omnibus Kruskal-Wallis with pairwise Welch's t tests, with an optional `aldex.glm` design path for covariate-adjusted models.

(4) *Functional profiling* (R/09–R/11) runs a shared analysis engine three times against domain-specific databases: the resistome (ABRicate/CARD, with drug-class classification), the virulome (VFDB, with functional categories extracted from the ABRicate `PRODUCT` field) and the mobilome (PlasmidFinder, with configurable replicon-family classification — Col-like, IncF, IncX, other Inc). Each domain emits its own alpha diversity, Venn diagram, gene- and category-level heatmaps, relative-abundance and total-count panels, and a PCoA + PERMANOVA + PERMDISP suite.

(5) *Integration, reporting and manifest* (R/12–R/14). R/12 reports igraph topology metrics (degree, betweenness, Louvain modularity) alongside its network, chord and Sankey outputs, and exports Gephi-compatible node/edge tables. R/13 assembles main and supplementary composite figures (PNG and TIFF independently toggleable) driven by a slot-YAML layout, plus `supplementary_tables.xlsx`. R/14 writes `manifest.json` (contract v1.2), a kind-tagged catalogue of every figure, table and intermediate artefact the pipeline produced, providing a stable contract for downstream integration.

**Reporting.** The Quarto template renders all outputs as a multi-page HTML dashboard with per-domain panels (taxonomy, diversity, differential abundance, resistome, virulome, mobilome, network) and embedded download links. Dependencies are pinned via `renv`, and the project requires R ≥ 4.5; external command-line tools (Kraken2, Bracken, ABRicate, Re-centrifuge) are assumed to have been run upstream.

## 3 Application and Validation

We are validating the pipeline against two real long-read chicken caecum metagenomics studies. The reference cohort (`Chicken batch 1`, PC_JC_2024-11-29; 18 samples; 3 dietary groups — Dulse, Reference diet, Soyabean meal) has approximately 46 hand-edited R scripts and known-good figures and tables already in place, so each stage can be compared output-for-output against the corresponding ground-truth artefact. The scaling and feature-coverage cohort (`Chicken batch 2`, JC_FL_2025-05-20; 32 samples; 4-level composite `Treatment_Bird` grouping — Control_W4, Control_W5, Dulse_W4, Dulse_W5) exercises a random effect (`block`, pens 1–8), two declared negative controls (CTRL1, CTRL2), and a larger per-stage feature load. Two purpose-built validation scripts (`scripts/diff_normdata.R` and `scripts/probe_taxid_parse.R`) provide numerical-diff and regression checks on the most error-prone intermediates, and `run_pipeline.R` auto-instruments per-stage wall time and peak R-process memory (via a `gc()`-bracketed wrapper equivalent to the `peakRAM` methodology) and prints a slowest-first breakdown to the run log.

**Table 1.** Per-stage primary output on the two validation cohorts (B1 = Chicken batch 1, n=18; B2 = Chicken batch 2, n=32). Same statistical settings (PERMANOVA 9999 permutations, ALDEx2 `mc.samples` = 128); 2026-06-23 sequential run.

| Stage | B1 primary output | B2 primary output |
| --- | --- | --- |
| load_inputs | 5 upstream tables parsed | 5 upstream tables parsed |
| clean_data | `taxid_fixes.csv` applied; no controls declared | no `taxid_fixes.csv`; 2 controls (CTRL1, CTRL2) subtracted |
| normalisation | gene-element TPM table written | gene-element TPM table written |
| taxonomy | 2,971 / 11,196 (taxon, sample) rows kept (count ≥ 5) | 5,570 / 6,549 (taxon, sample) rows kept |
| relative_abundance | 1,700 / 72,152 rows kept (count > 17); top-50 + Others | 7,557 / 18,073 rows kept; top-50 + Others |
| alpha_diversity | 9 metrics × 18 samples; KW per metric | 9 metrics × 32 samples; KW per metric |
| beta_diversity | PERMANOVA R² = 0.145, p = 0.22; PERMDISP p = 0.027 | PERMANOVA R² = 0.121, p = 0.15; PERMDISP p = 0.053 |
| differential_abundance | 31 taxa + 157 gene features; 3 pairwise contrasts | 53 taxa + 221 gene features; 6 pairwise contrasts |
| resistome (CARD) | 348 rows / 85 genes / 18 drug classes | 648 rows / 99 genes / 20 drug classes |
| virulome (VFDB) | 320 rows / 140 genes / 6 VF functions (17 / 18 samples) | 492 rows / 156 genes / 8 VF functions (26 / 32 samples) |
| mobilome (PlasmidFinder) | 65 rows / 10 replicons / 3 families | 137 rows / 25 replicons / 5 families |
| network | 247 edges, 100 nodes (18 + 42 + 40); modularity = 0.31 | 576 edges, 210 nodes (32 + 77 + 101); modularity = 0.40 |
| panels | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` |
| manifest | `manifest.json` (contract v1.2; 85 kB) | `manifest.json` (contract v1.2; 80 kB) |
| report | Quarto dashboard HTML (~14 MB) | Quarto dashboard HTML (~18 MB) |

Across both cohorts (**Table 1**), every enabled stage produced its expected primary artefact under identical settings. Feature counts scale with cohort size and grouping complexity: B2 (n=32, 4-level `Treatment_Bird`) yields ~1.5–2× the resistome (648 vs 348 rows), virulome (492 vs 320) and mobilome (137 vs 65) rows of B1, with 53 vs 31 differentially abundant taxa across 6 vs 3 pairwise contrasts and a denser integrated network (576 edges/210 nodes/modularity 0.40 vs 247/100/0.31). The pipeline closes cleanly on both runs, producing a contract-v1.2 `manifest.json` (~80–85 kB) and a Quarto HTML dashboard (~14–18 MB). The performance profile (**Figure 2**) illustrated that the end-to-end runtime grew from 6m 41s on B1 to 8m 19s on B2 (+24 %, sublinear in both sample count and pairwise-contrast count), dominated on both cohorts by ALDEx2 differential abundance (32 % of B1 total, 40 % of B2). Peak R-process memory tracks the same monotone profile — climbing through the gene-element domains as the joint long-form `(sample × gene)` tables grow (~1 GB plateau) and tops out during R/13 panel assembly (1.33 GB on B1, 1.51 GB on B2), where every domain's PNGs and tables are loaded simultaneously to lay out; the `report` stage's modest peak (~0.6 GB) reflects work delegated to the external Quarto process.

![**Figure 2.** Pipeline performance on the two validation cohorts. (A) Per-stage wall time (log scale) and (B) peak R-process heap, sorted by B2 runtime. (C) Cumulative wall time and (D) running peak memory in pipeline order. R 4.5.0 single-thread analysis, including PERMANOVA 9999 permutations, ALDEx2 `mc.samples` = 128. B1 = Chicken batch 1 (n=18, 3 dietary groups); B2 = Chicken batch 2 (n=32, 4-level `Treatment_Bird`, `block` random effect, two negative controls).](figure2_pipeline_performance.png)

One expected and deliberate divergence concerns the resistome stage: the ground-truth `cleanData.R` misparses taxids carrying variant suffixes, dropping the corresponding ABRicate rows. R/02 fixes this parsing bug, so the pipeline retains 110 additional (sample, GENE) rows relative to the hand-edited output (`scripts/diff_normdata.R` on the 2026-06-23 B1 run: pipeline = 733 rows, ground truth = 623; 121 rows gained by the parsing fix, 11 rows present in the ground truth that the pipeline currently drops — concentrated in a handful of plasmid replicon and accessory resistance calls, queued for inspection rather than treated as a regression). Jaccard agreement on the `(sample, GENE)` key set is 0.82. Across the other stages with directly comparable ground-truth artefacts, sample-set agreement is exact for alpha diversity (18/18), and Jaccard for the network Gephi exports against the most recent hand-edited snapshot is 0.85 (nodes) and 0.87 (edges). Per-stage comparisons against the remaining hand-edited artefacts — taxonomy, beta diversity, and differential abundance — are deferred to a separate validation effort. The Re-centrifuge column-slicing heuristic in R/02 is acknowledged as still coarse and is queued for rework.

## 4 Conclusion

MERIDIAN consolidates community analyses, functional profiling, and an integrated network reconstruction into one configuration-driven R project for long-read shotgun metagenomics. Single-thread end-to-end execution completes in minutes with a low memory footprint, scaling sublinearly with cohort size and grouping complexity. It is in active development; we are releasing this alpha-stage workflow (v0.1.0) now to invite feedback from groups running comparable analyses.

---

**Acknowledgements.** The authors thank the **Centre for Innovation in Genomics and Microbiome Sciences (CIGMiS)**, which generated the chicken caecum metagenomics datasets used to develop and validate the pipeline. We are grateful to **colleagues** for feedback on the configuration interface and Quarto report layout, and to the School of Biomedical Science at the University of West London for hosting the project. **[Optional: compute providers — e.g. UWL / UCL HPC, cloud credits.]**

**Funding.** This study was funded by a Vice-Chancellor scholarship provided by the University of West London.

**Conflict of Interest.** The authors declare that they have no competing interests.

---

## References

1. Al-Khalaifah H, Rahman MH, Al-Surrayai T, Al-Dhumair A, Al-Hasan M. A One-Health Perspective of Antimicrobial Resistance (AMR): Human, Animals and Environmental Health. *Life (Basel)* 2025;15(10):1598. doi:10.3390/life15101598.
2. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol* 2019;20:257. doi:10.1186/s13059-019-1891-0.
3. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci* 2017;3:e104. doi:10.7717/peerj-cs.104.
4. Seemann T. ABRicate: mass screening of contigs for antimicrobial and virulence genes. https://github.com/tseemann/abricate (accessed 2026-06-23).
5. Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res* 2023;51:D690–D699. doi:10.1093/nar/gkac920.
6. Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Res* 2022;50:D912–D917. doi:10.1093/nar/gkab1107.
7. Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014;58:3895–3903. doi:10.1128/AAC.02412-14.
8. Martí JM. Recentrifuge: robust comparative analysis and contamination removal for metagenomics. *PLoS Comput Biol* 2019;15:e1006967. doi:10.1371/journal.pcbi.1006967.
9. Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, et al. vegan: Community Ecology Package. R package version 2.7-3, 2026. doi:10.32614/CRAN.package.vegan.
10. Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome* 2014;2:15. doi:10.1186/2049-2618-2-15.
11. Li Y, Gao Y, Liu X, Mao Y, Wang M, Qin Y, et al. Quantifying antibiotic resistome risks across environmental niches: the L-ARRAP for long-read metagenomic profiling. *Brief Bioinform* 2025;26(5):bbaf535. doi:10.1093/bib/bbaf535.
12. Theil S, Rifa E. rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis. *F1000Res* 2021;10:7. doi:10.12688/f1000research.27268.1.
13. Allaire JJ, et al. Quarto. https://quarto.org (accessed 2026-06-23).
14. R Core Team. R: A Language and Environment for Statistical Computing (version 4.5.0). R Foundation for Statistical Computing, Vienna, 2025. https://www.R-project.org/.
15. Pedersen TL. ggraph: An Implementation of Grammar of Graphics for Graphs and Networks. R package version 2.2.2, 2025. doi:10.32614/CRAN.package.ggraph.
16. Csárdi G, Nepusz T. The igraph software package for complex network research. *InterJournal* 2006;Complex Systems:1695. R package version 2.3.1; doi:10.5281/zenodo.7682609. https://igraph.org.
17. Allaire JJ, Gandrud C, Russell K, Yetman C. networkD3: D3 JavaScript Network Graphs from R. R package version 0.4.1, 2025. doi:10.32614/CRAN.package.networkD3.
