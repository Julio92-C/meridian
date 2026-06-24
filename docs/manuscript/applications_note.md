# MERIDIAN: a config-driven R workflow for long-read shotgun metagenomics analysis

**Authors / Affiliations**

Julio C. Ortega Cambara<sup>1*</sup>, Hermine Mkrtchyan<sup>1</sup>, Raju Misra<sup>2</sup>, Timothy D. McHugh<sup>3</sup>, Sylvia Rofael<sup>3</sup>, David Lowe<sup>3</sup>, Piotr Cuber<sup>4</sup>, Pedro Humberto Lebre<sup>1</sup>

<sup>1</sup>School of Biomedical Science, University of West London, London, United Kingdom
<sup>2</sup>UK Health Security Agency, London, United Kingdom
<sup>3</sup>Centre for Clinical Microbiology, University College London, London, United Kingdom
<sup>4</sup>Natural History Museum, London, United Kingdom

<sup>*</sup>Correspondence: Hermine Mkrtchyan — Hermine.Mkrtchyan@uwl.ac.uk

{{TODO: ORCIDs for all authors}}

---

## Abstract

**Motivation.** Long-read shotgun metagenomics studies are typically analysed by hand-editing dozens of R scripts per project — pasting sample identifiers, patching taxonomy strings and copying chunks between studies. Existing R workflows largely target short-read amplicon data or a single analysis (e.g. resistome risk scoring), leaving long-read shotgun studies without an integrated, configuration-driven option.

**Results.** We present **MERIDIAN** (*Metagenomic Evaluation of Resistance, Identity & Diversity through Integrated Analysis of Nanopore sequencing*), an alpha-stage modular R workflow that consumes Kraken2, Bracken, ABRicate and Re-centrifuge outputs and produces a single Quarto HTML dashboard together with a set of publication-ready composite figures (PNG and TIFF) and a supplementary tables workbook. The workflow covers decontamination, TPM normalisation, taxonomy, alpha/beta diversity with PERMANOVA, two-level ALDEx2 differential abundance, and resistome/virulome/mobilome profiling joined by a tripartite sample × taxon × gene network with companion chord and Sankey diagrams. A study is parameterised by one YAML file; module code is not edited per-project, and every run emits a `manifest.json` (contract v1.2) artefact catalogue tagging every figure and table by kind for downstream tooling.

**Availability and Implementation.** Source at https://github.com/Julio92-C/Metagenomics_pipeline_automation (MIT). Run with `git clone <repo> && Rscript run_pipeline.R config.yaml` after `renv::restore()`.

**Contact.** Hermine.Mkrtchyan@uwl.ac.uk

---

## 1 Introduction

By simultaneously profiling microbial communities and their functional genetic elements across environmental, clinical, and veterinary or livestock matrices, long-read shotgun metagenomics (Oxford Nanopore, PacBio) serves as a powerful One Health tool for tracking health threats that intersect humans, animals, and ecosystems [1]. The downstream analysis stack is well established — Kraken2 [2] and Bracken [3] for taxonomic assignment and abundance estimation, ABRicate [4] against CARD [5], VFDB [6] and PlasmidFinder [7] for antimicrobial resistance genes (ARGs), virulence factors, and mobile genetic elements (MGEs) gene screening, Re-centrifuge [8] for contamination filtering, and R packages such as vegan [9] and ALDEx2 [10] for ecological and differential analyses — but the integration between them is not.

In practice, each new study is analysed by adapting a set of bespoke R scripts: hard-coded input paths, sample-identifier filters embedded as literal strings, per-study negative-control renames, and taxonomy fixes patched by row index. The resulting per-study script copies drift apart, no shared interface emits a machine-readable record of which outputs each run produced, and cross-study comparisons require reconciling code as well as data. The reference cohort for the present work (`chicken_batch1`; ~46 hand-edited scripts) typifies this pattern. The L-ARRAP pipeline [11] addresses one slice of this problem for long-read resistome risk assessment, and the rANOMALY workflow [12] demonstrates the value of a config-driven R workflow with auto-generated reports for amplicon data. To our knowledge, no comparable workflow exists for end-to-end long-read shotgun metagenomics covering taxonomy, diversity, differential abundance, and the full resistome–virulome–mobilome profiling in one configurable run.

We present an early-stage R workflow that fills this gap. The tool is in active development; this note describes its design, current scope, and the real-data validation strategy we are using to harden each stage.

## 2 Implementation

The workflow is a stand-alone R project, not an installable package. A single entry point (`run_pipeline.R`) consumes one YAML configuration file per study and dispatches to fifteen ordered modules (`R/00_setup.R` through `R/14_panels.R`); a Quarto dashboard template (`templates/report.qmd`) renders all outputs into a single multi-page HTML deliverable, a publication-panels stage assembles main and supplementary composite figures plus a supplementary tables workbook, and a final manifest stage emits `manifest.json` describing every artefact produced (Figure 1).

![**Figure 1.** MERIDIAN workflow. A single `config.yaml` parameterises a study and drives fifteen ordered R modules (R/00–R/14) consuming upstream outputs (Kraken2/Bracken, ABRicate, Re-centrifuge) plus sample metadata. After setup, input loading and decontamination (R/00–R/02), community analyses (R/04–R/08) and functional profiling against CARD, VFDB and PlasmidFinder (TPM normalisation in R/03 and per-domain analyses in R/09–R/11; alpha/beta diversity and ALDEx2 DAA shared between the two columns) feed a tripartite sample × taxon × gene network with companion chord and Sankey diagrams (R/12). R/12 then fans out to two output columns: a Quarto dashboard (`templates/report.qmd`, HTML with per-domain panels) and a publication-panels stage (R/14, composite figures plus `supplementary_tables.xlsx`), whose PNG/TIFF panels and tables workbook are the publication-ready outputs. R/13 finally writes `manifest.json` (contract v1.2; kind-tagged figure and table catalogue), cataloguing the dashboard, the R/14 panels and every intermediate artefact. Source: `generate_figure1.py`.](figure1_pipeline_schematic.png)

**Configuration-driven.** Everything that varies between studies — metadata file path, grouping columns, negative-control identifiers, abundance thresholds, sample selection — is declared in `config.yaml`. Module code is not edited per-project; dataset-specific adjustments live in `config.yaml` or in the structured fix tables it references. R/00 (setup) and R/01 (input loading) always run; every subsequent stage can be toggled via `cfg$stages$*` flags so a user can, for example, re-run only the resistome modules after revising thresholds. TPM normalisation (R/03) and all gene-element stages (R/09–R/11) require ABRicate outputs and are skipped otherwise. Per-study NCBI taxonomy fixes, when needed, are declared in a structured `taxid_fixes.csv` rather than as row-index edits. A single `cfg$stats$padjust_method` knob (Benjamini-Hochberg by default; full list in `docs/STATISTICS.md`) routes every multiple-testing correction across the pipeline.

**Stages.** The pipeline organises its work into five conceptual blocks.

(1) *Ingest and cleaning* (R/00–R/02) load upstream tables — Kraken2/Bracken, ABRicate, Re-centrifuge contamination tags, study metadata — accommodating per-project input shapes rather than assuming a fixed column layout, apply optional negative-control subtraction, merge Bracken counts with the metadata covariates declared in `cfg$metadata$fixed_effects` / `random_effect`, and link ABRicate calls to taxonomic assignments. Per-study taxonomy patches are applied from `taxid_fixes.csv` when supplied.

(2) *Gene-element normalisation* (R/03) computes per-(sample, GENE) TPM from ABRicate hits and gene lengths.

(3) *Community analyses* (R/04–R/08) cover the standard taxonomic toolkit: Venn diagrams and annotated heatmaps after a shared taxa-name cleanup; top-N stacked-bar relative-abundance plots (static PNG and interactive plotly HTML) with paired total-count and sample-prevalence panels; alpha diversity reporting Shannon, Simpson, richness, Pielou and Fisher metrics with Kruskal-Wallis tests and violin/bar plots, preferring a pre-computed indices CSV (e.g. Oxford Nanopore wf-metagenomics output) and falling back to vegan; PCoA on Hellinger-transformed Bray-Curtis distances with PERMANOVA (`adonis2`, 9999 permutations) and PERMDISP (`betadisper`) overlaid on the ordination; and ALDEx2 differential abundance at two levels — taxa (Bracken counts) and genes (ABRicate contig hits per (sample, GENE)) — combining an omnibus Kruskal-Wallis with pairwise Welch's t tests, with an optional `aldex.glm` design path for covariate-adjusted models.

(4) *Functional profiling* (R/09–R/11) runs a shared analysis engine three times against domain-specific databases: the resistome (ABRicate/CARD, with drug-class classification and a macrolide-lincosamide-streptogramin rollup), the virulome (VFDB, with functional categories extracted from the ABRicate `PRODUCT` field) and the mobilome (PlasmidFinder, with configurable replicon-family classification — Col-like, IncF, IncX, other Inc). Each domain emits its own alpha diversity, Venn diagram, gene- and category-level heatmaps, relative-abundance and total-count panels, and a PCoA + PERMANOVA + PERMDISP suite.

(5) *Integration, reporting and manifest* (R/12–R/14). R/12 assembles a tripartite ggraph network of sample × taxon × gene plus chord diagrams (overall and per group) and four-tier networkD3 Sankey diagrams, exports Gephi-compatible node/edge tables, and reports igraph topology metrics (degree, betweenness, Louvain modularity) and Bray-Curtis sample clusters. R/14 assembles publication-ready main and supplementary composite figures (PNG and TIFF outputs are independently toggleable) driven by a Frontiers slot YAML, alongside `supplementary_tables.xlsx`. R/13 finally writes `manifest.json` (contract v1.2), a kind-tagged catalogue of every figure, table and intermediate artefact the run produced, providing a stable JSON contract that decouples the pipeline from downstream report consumers and manuscript-drafting tools.

**Reporting.** The Quarto template renders all outputs as a multi-page HTML dashboard with per-domain panels (taxonomy, diversity, differential abundance, resistome, virulome, mobilome, network) and embedded download links. Dependencies are pinned via `renv` and the project requires R ≥ 4.5; external command-line tools (Kraken2, Bracken, ABRicate, Re-centrifuge) are assumed to have been run upstream.

## 3 Application and Validation

We are validating the pipeline against two real long-read chicken caecum metagenomics studies. `chicken_batch1` (PC_JC_2024-11-29; 18 samples; 3 dietary groups — Dulce, Reference diet, Soyabean meal) is the reference cohort: approximately 46 hand-edited R scripts and known-good figures and tables already exist for it, so each stage can be compared output-for-output against the corresponding ground-truth artefact. `chicken_batch2` (JC_FL_2025-05-20; 32 samples; 4-level composite `Treatment_Bird` grouping — Control_W4, Control_W5, Dulse_W4, Dulse_W5) is a scaling and feature-coverage cohort: it exercises a random effect (`block`, pens 1–8), two declared negative controls (CTRL1, CTRL2) subtracted in R/02, and a larger per-stage feature load. Two purpose-built validation scripts (`scripts/diff_normdata.R` and `scripts/probe_taxid_parse.R`) provide numerical-diff and regression checks on the most error-prone intermediates, and `run_pipeline.R` auto-instruments per-stage wall time and peak R-process memory (via a `gc()`-bracketed wrapper equivalent to the `peakRAM` methodology) and prints a slowest-first breakdown to the run log, mechanically populating both the runtime and peak-MB columns of Table 1.

**Table 1.** Per-stage wall time, peak R-process memory and primary output for both validation cohorts on the same workstation (Windows 11, R 4.5, single thread); identical pipeline commit and identical statistical settings (PERMANOVA 9999 permutations; ALDEx2 `mc.samples` = 128). **B1** = chicken_batch1 (PC_JC_2024-11-29; n=18; 3 dietary groups). **B2** = chicken_batch2 (JC_FL_2025-05-20; n=32; 4-level `Treatment_Bird` grouping, `block` random effect, two declared negative controls). Sample sizes shown in some output cells refer to samples with non-zero hits at that stage. Runtimes and peak-MB are from the same 2026-06-23 sequential run (B1 = 6m 41s, B2 = 8m 19s; 15/15 stages each) on the post-instrumentation build; the peak-MB column reports the gc()-tracked maximum used heap for each stage (equivalent to `peakRAM` methodology) and the Total row shows the overall maximum across the run.

| Stage | B1 runtime | B1 peak MB | B1 primary output | B2 runtime | B2 peak MB | B2 primary output |
| --- | --- | --- | --- | --- | --- | --- |
| load_inputs | 0.7 s | 168 | 5 upstream tables parsed | 1.3 s | 205 | 5 upstream tables parsed |
| clean_data | 2.1 s | 245 | `taxid_fixes.csv` applied; no controls declared | 6.4 s | 522 | no `taxid_fixes.csv`; 2 controls (CTRL1, CTRL2) subtracted |
| normalisation | 0.4 s | 160 | gene-element TPM table written | 0.5 s | 205 | gene-element TPM table written |
| taxonomy | 2.2 s | 230 | 2,971 / 11,196 (taxon, sample) rows kept (count ≥ 5) | 2.3 s | 308 | 5,570 / 6,549 (taxon, sample) rows kept |
| relative_abundance | 27.6 s | 362 | 1,700 / 72,152 rows kept (count > 17); top-50 + Others | 29.6 s | 574 | 7,557 / 18,073 rows kept; top-50 + Others |
| alpha_diversity | 32.6 s | 472 | 9 metrics × 18 samples; KW per metric | 31.4 s | 587 | 9 metrics × 32 samples; KW per metric |
| beta_diversity | 5.9 s | 472 | PERMANOVA R² = 0.145, p = 0.22; PERMDISP p = 0.027 | 5.2 s | 482 | PERMANOVA R² = 0.121, p = 0.15; PERMDISP p = 0.053 |
| differential_abundance | 2m 7s | 1028 | 31 taxa + 157 gene features; 3 pairwise contrasts | 3m 18s | 1133 | 53 taxa + 221 gene features; 6 pairwise contrasts |
| resistome (CARD) | 29.3 s | 1028 | 348 rows / 85 genes / 18 drug classes | 25.9 s | 1102 | 648 rows / 99 genes / 20 drug classes |
| virulome (VFDB) | 16.0 s | 1028 | 320 rows / 140 genes / 6 VF functions (17 / 18 samples) | 17.3 s | 1132 | 492 rows / 156 genes / 8 VF functions (26 / 32 samples) |
| mobilome (PlasmidFinder) | 15.3 s | 1028 | 65 rows / 10 replicons / 3 families | 16.1 s | 1100 | 137 rows / 25 replicons / 5 families |
| network | 53.3 s | 1183 | 247 edges, 100 nodes (18 + 42 + 40); modularity = 0.31 | 1m 26s | 1418 | 576 edges, 210 nodes (32 + 77 + 101); modularity = 0.40 |
| panels (R/14) | 49.0 s | 1330 | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` | 49.0 s | 1505 | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` |
| manifest | 2.8 s | 922 | `manifest.json` (contract v1.2; 85 kB) | 3.1 s | 1040 | `manifest.json` (contract v1.2; 80 kB) |
| report | 36.9 s | 611 | Quarto dashboard HTML (~14 MB) | 27.4 s | 647 | Quarto dashboard HTML (~18 MB) |
| **Total (15 stages)** | **6m 41s** | **1330** | — | **8m 19s** | **1505** | — |

Differential abundance dominates wall time on both cohorts (32 % of total on B1, 40 % on B2), driven by the ALDEx2 Monte-Carlo loop scaling with feature count, sample count and pairwise-contrast count; the network stage is next-most expensive (~13–17 %) because Sankey assembly iterates over every gene–category mapping. Across the recorded stages, increasing the sample count from 18 to 32 and the pairwise-contrast count from 3 to 6 (B1 → B2) raised end-to-end runtime from 6m 41s to 8m 19s (+24 %), markedly sublinear in either axis and consistent with ALDEx2 still being the dominant path rather than per-sample fixed overhead. Peak R-process memory tracks the same monotone profile — climbing through the gene-element domains as the joint long-form `(sample × gene)` tables grow (~1.0 GB plateau) and topping out during R/14 panel assembly (1.33 GB on B1, 1.51 GB on B2), where every domain's PNGs and tables are loaded simultaneously to lay out the composite figures and write `supplementary_tables.xlsx`. The `report` stage's peak is comparatively modest (~0.6 GB) because the heavy lifting is delegated to the external Quarto process and is therefore not captured by the in-process gc() instrumentation. Figure 2 visualises the same numbers: panels A and B give the per-stage wall time and peak heap as paired bars (stages sorted by B2 runtime), and panels C and D show the cumulative trajectories in pipeline order, making the ALDEx2 inflection and the R/14 memory peak readable at a glance. Tabulated PERMANOVA and PERMDISP statistics for the per-domain stages (resistome, virulome, mobilome), per-metric KW p-values for alpha diversity, and all DAA per-pair feature tables are surfaced in the Quarto dashboard rather than reproduced here.

![**Figure 2.** Pipeline performance on the two validation cohorts. (A) Per-stage wall time (log scale) and (B) peak R-process heap, both with stages sorted by B2 runtime. (C) Cumulative wall time and (D) running peak memory across the 15 stages in pipeline order. Numbers reproduce Table 1; both runs used identical pipeline commit `eee122c`, R 4.5.0 single thread, PERMANOVA 9999 permutations, ALDEx2 `mc.samples` = 128. B1 = chicken_batch1 (PC_JC_2024-11-29; n=18; 3 dietary groups). B2 = chicken_batch2 (JC_FL_2025-05-20; n=32; 4-level `Treatment_Bird`, `block` random effect, two declared negative controls). Source: `generate_figure2.R`.](figure2_pipeline_performance.png)

One expected and deliberate divergence concerns the resistome stage: the ground-truth `cleanData.R` mis-parses taxids carrying variant suffixes, dropping the corresponding ABRicate rows. R/02 fixes this parsing bug, so the pipeline retains 110 additional (sample, GENE) rows relative to the hand-edited output (`scripts/diff_normdata.R` on the 2026-06-23 B1 run: pipeline = 733 rows, ground truth = 623; 121 rows gained by the parsing fix, 11 rows present in the ground truth that the pipeline currently drops — concentrated in a handful of plasmid replicon and accessory resistance calls, queued for inspection rather than treated as a regression). Jaccard agreement on the `(sample, GENE)` key set is 0.82. Across the other stages with directly comparable ground-truth artefacts, sample-set agreement is exact on alpha diversity (18 / 18), and Jaccard on the network gephi exports against the most-recent hand-edited snapshot is 0.85 (nodes) and 0.87 (edges); per-stage agreement against the remaining hand-edited artefacts (taxonomy, beta diversity, differential abundance) is deferred to a separate validation effort. The Re-centrifuge column-slicing heuristic in R/02 is acknowledged as still coarse and is queued for rework; Table 1 will be re-derived once that change lands.

## 4 Conclusion

MERIDIAN is an alpha-stage workflow that consolidates a long-read shotgun metagenomics analysis — taxonomy, diversity, differential abundance, resistome, virulome, mobilome, and an integrated network view — into one configuration-driven R project with a single Quarto dashboard, a publication-panels stage and a downstream-readable `manifest.json` contract. It is in active development and validation against real reference data is ongoing; we release it now to invite feedback from groups running comparable analyses and to encourage convergence on a shared, scriptable layout for long-read metagenomics studies.

---

**Acknowledgements.** The authors thank the **[clinical / experimental collaborators]** who generated the chicken caecum metagenomics datasets used to develop and validate the pipeline. We are grateful to **[early-tester names / labs]** for feedback on the configuration interface and Quarto report layout, and to the School of Biomedical Science at the University of West London for hosting the project. **[Optional: compute providers — e.g. UWL / UCL HPC, cloud credits.]**

**Funding.** This study was funded by a Vice-Chancellor scholarship provided by the University of West London.

**Conflict of Interest.** The authors declare that they have no competing interests.

---

## References

1. Al-Khalaifah H, Rahman MH, Al-Surrayai T, Al-Dhumair A, Al-Hasan M. A One-Health Perspective of Antimicrobial Resistance (AMR): Human, Animals and Environmental Health. *Life (Basel)* 2025;15(10):1598. doi:10.3390/life15101598.
2. Wood DE, Salzberg SL. Kraken: ultrafast metagenomic sequence classification using exact alignments. *Genome Biol* 2014;15:R46. doi:10.1186/gb-2014-15-3-r46. Kraken2: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol* 2019;20:257. doi:10.1186/s13059-019-1891-0.
3. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci* 2017;3:e104. doi:10.7717/peerj-cs.104.
4. Seemann T. ABRicate: mass screening of contigs for antimicrobial and virulence genes. https://github.com/tseemann/abricate (accessed 2026-06-23).
5. Alcock BP, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res* 2023;51:D690–D699. doi:10.1093/nar/gkac920.
6. Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Res* 2022;50:D912–D917. doi:10.1093/nar/gkab1107.
7. Carattoli A, Zankari E, García-Fernández A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014;58:3895–3903. doi:10.1128/AAC.02412-14.
8. Martí JM. Recentrifuge: robust comparative analysis and contamination removal for metagenomics. *PLoS Comput Biol* 2019;15:e1006967. doi:10.1371/journal.pcbi.1006967.
9. Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, et al. vegan: Community Ecology Package. R package version 2.7-3, 2026. doi:10.32614/CRAN.package.vegan.
10. Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome* 2014;2:15. doi:10.1186/2049-2618-2-15.
11. Li Y, Gao Y, Liu X, Mao Y, Wang M, Qin Y, Zhang C, Chen Q, Ning K, Wang Z, Han M. Quantifying antibiotic resistome risks across environmental niches: the L-ARRAP for long-read metagenomic profiling. *Brief Bioinform* 2025;26(5):bbaf535. doi:10.1093/bib/bbaf535.
12. Theil S, Rifa E. rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis. *F1000Res* 2021;10:7. doi:10.12688/f1000research.27268.1.
13. Allaire JJ, et al. Quarto. https://quarto.org (accessed 2026-06-23).
14. R Core Team. R: A Language and Environment for Statistical Computing (version 4.5.0). R Foundation for Statistical Computing, Vienna, 2025. https://www.R-project.org/.
15. Pedersen TL. ggraph: An Implementation of Grammar of Graphics for Graphs and Networks. R package version 2.2.2, 2025. doi:10.32614/CRAN.package.ggraph.
16. Csárdi G, Nepusz T. The igraph software package for complex network research. *InterJournal* 2006;Complex Systems:1695. R package version 2.3.1; doi:10.5281/zenodo.7682609. https://igraph.org.
17. Allaire JJ, Gandrud C, Russell K, Yetman C. networkD3: D3 JavaScript Network Graphs from R. R package version 0.4.1, 2025. doi:10.32614/CRAN.package.networkD3.
