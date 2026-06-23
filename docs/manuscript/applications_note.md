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

We are validating the pipeline against two real long-read chicken caecum metagenomics studies. `chicken_batch1` (PC_JC_2024-11-29; 18 samples; 3 dietary groups — Dulce, Reference diet, Soyabean meal) is the reference cohort: approximately 46 hand-edited R scripts and known-good figures and tables already exist for it, so each stage can be compared output-for-output against the corresponding ground-truth artefact. `chicken_batch2` (JC_FL_2025-05-20; 32 samples; 4-level composite `Treatment_Bird` grouping — Control_W4, Control_W5, Dulse_W4, Dulse_W5) is a scaling and feature-coverage cohort: it exercises a random effect (`block`, pens 1–8), two declared negative controls (CTRL1, CTRL2) subtracted in R/02, and a larger per-stage feature load. Two purpose-built validation scripts (`scripts/diff_normdata.R` and `scripts/probe_taxid_parse.R`) provide numerical-diff and regression checks on the most error-prone intermediates, and `run_pipeline.R` auto-instruments per-stage wall time and prints a slowest-first breakdown to the run log, mechanically populating the runtime columns of Table 1.

**Table 1.** Per-stage wall time and primary output for both validation cohorts on the same workstation (Windows 11, R 4.5, single thread); identical pipeline commit and identical statistical settings (PERMANOVA 9999 permutations; ALDEx2 `mc.samples` = 128). **B1** = chicken_batch1 (PC_JC_2024-11-29; n=18; 3 dietary groups). **B2** = chicken_batch2 (JC_FL_2025-05-20; n=32; 4-level `Treatment_Bird` grouping, `block` random effect, two declared negative controls). Sample sizes shown in some output cells refer to samples with non-zero hits at that stage. Runtimes are from the 2026-05-31 pipeline.log and pre-date the R/14 panels stage; a fresh end-to-end re-run is queued before submission (see Drafting notes).

| Stage | B1 runtime | B1 primary output | B2 runtime | B2 primary output |
| --- | --- | --- | --- | --- |
| load_inputs | 0.8 s | 5 upstream tables parsed | 1.7 s | 5 upstream tables parsed |
| clean_data | 2.0 s | `taxid_fixes.csv` applied; no controls declared | 6.0 s | no `taxid_fixes.csv`; 2 controls (CTRL1, CTRL2) subtracted |
| normalisation | 0.2 s | gene-element TPM table written | 0.3 s | gene-element TPM table written |
| taxonomy | 2.7 s | 2,971 / 11,196 (taxon, sample) rows kept (count ≥ 5) | 3.1 s | 5,570 / 6,549 (taxon, sample) rows kept |
| relative_abundance | 6.5 s | 1,700 / 72,152 rows kept (count > 17); top-50 + Others | 10.9 s | 7,557 / 18,073 rows kept; top-50 + Others |
| alpha_diversity | 12.1 s | 9 metrics × 18 samples; KW per metric | 19.8 s | 9 metrics × 32 samples; KW per metric |
| beta_diversity | 3.5 s | PERMANOVA R² = 0.145, p = 0.22; PERMDISP p = 0.028 | 5.1 s | PERMANOVA R² = 0.121, p = 0.15; PERMDISP p = 0.056 |
| differential_abundance | 1m 54s | 31 taxa + 157 gene features; 3 pairwise contrasts | 3m 27s | 53 taxa + 221 gene features; 6 pairwise contrasts |
| resistome (CARD) | 9.0 s | 348 rows / 85 genes / 18 drug classes | 10.1 s | 648 rows / 99 genes / 20 drug classes |
| virulome (VFDB) | 8.0 s | 320 rows / 140 genes / 6 VF functions (17 / 18 samples) | 9.7 s | 492 rows / 156 genes / 8 VF functions (26 / 32 samples) |
| mobilome (PlasmidFinder) | 7.7 s | 65 rows / 10 replicons / 3 families | 10.1 s | 137 rows / 25 replicons / 5 families |
| network | 33.4 s | 247 edges, 100 nodes (18 + 42 + 40); modularity = 0.30 | 1m 7s | 576 edges, 210 nodes (32 + 77 + 101); modularity = 0.40 |
| panels (R/14) | {{TODO: re-run}} | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` | {{TODO: re-run}} | main + supplementary composite figures (PNG / TIFF) + `supplementary_tables.xlsx` |
| manifest | 0.5 s | `manifest.json` (contract v1.2; 47 kB) | 0.4 s | `manifest.json` (contract v1.2; 57 kB) |
| report | 25.1 s | Quarto dashboard HTML (~14 MB) | 40.5 s | Quarto dashboard HTML (~18 MB) |
| **Total (15 stages)** | **{{TODO: re-run}}** | — | **{{TODO: re-run}}** | — |

Differential abundance dominates wall time on both cohorts (≈50 % of the pre-R/14 total), driven by the ALDEx2 Monte-Carlo loop scaling with feature count, sample count and pairwise-contrast count; the network stage is next-most expensive because Sankey assembly iterates over every gene–category mapping. Across the recorded stages, doubling the sample count and pairwise-contrast count (B1 → B2) raised end-to-end runtime by ~73 %, suggesting near-linear scaling on the dominant ALDEx2 path. Tabulated PERMANOVA and PERMDISP statistics for the per-domain stages (resistome, virulome, mobilome), per-metric KW p-values for alpha diversity, and all DAA per-pair feature tables are surfaced in the Quarto dashboard rather than reproduced here. Peak memory is auto-instrumented in `run_pipeline.R` (gc()-based, equivalent to peakRAM's methodology) and emitted per stage in the run log; the Table 1 peak-MB column will be populated from the next end-to-end re-run on both cohorts.

One expected and deliberate divergence concerns the resistome stage: the ground-truth `cleanData.R` mis-parses taxids carrying variant suffixes, dropping the corresponding ABRicate rows. R/02 fixes this parsing bug, so the pipeline retains approximately {{TODO: confirm exact count from a fresh validation run}} additional (sample, GENE) rows relative to the hand-edited output. We treat this as a correctness improvement rather than a failure to reproduce, and the validation log flags it accordingly. The Re-centrifuge column-slicing heuristic in R/02 is acknowledged as still coarse and is queued for rework; Table 1 will be re-derived once that change lands.

## 4 Conclusion

MERIDIAN is an alpha-stage workflow that consolidates a long-read shotgun metagenomics analysis — taxonomy, diversity, differential abundance, resistome, virulome, mobilome, and an integrated network view — into one configuration-driven R project with a single Quarto dashboard, a publication-panels stage and a downstream-readable `manifest.json` contract. It is in active development and validation against real reference data is ongoing; we release it now to invite feedback from groups running comparable analyses and to encourage convergence on a shared, scriptable layout for long-read metagenomics studies.

---

**Acknowledgements.** {{TODO: acknowledgements — collaborators, reviewers, compute providers}}

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

---

## Drafting notes

Placeholders inserted, and the data needed to fill each one:

Resolved since first draft:

- **Pipeline rename to MERIDIAN** (commit 8f5fb3f, 2026-06-23) — manuscript title now leads with "MERIDIAN:"; first mention in Abstract Results introduces the expanded form (*Metagenomic Evaluation of Resistance, Identity & Diversity through Integrated Analysis of Nanopore sequencing*); Figure 1 caption and §4 Conclusion follow. `generate_figure1.py` docstring banner and the orphaned `figure1_pipeline_schematic.mmd` header updated for consistency; `pipeline_facts.md` "What the tool is" entry renamed. Figure 1 PNG/SVG do not need re-rendering — the pipeline name does not appear on the rendered figure itself.
- Author block, affiliations, corresponding author and email — filled.
- GitHub URL — filled (https://github.com/Julio92-C/Metagenomics_pipeline_automation).
- Figure 1 — canonical source is `generate_figure1.py` (matplotlib, 300 DPI, custom colour-coded layout), added in the "polish Figure 1" pass to replace the original Mermaid render. Regenerate via `& "C:\Users\julio\AppData\Local\Programs\Python\Python313\python.exe" generate_figure1.py` (plain `python` resolves to Inkscape's interpreter and lacks matplotlib) — writes both `figure1_pipeline_schematic.png` and `.svg`. The old `figure1_pipeline_schematic.mmd` is legacy/orphaned: do not re-render it with `mmdc`, that overwrites the polished figure with the plain Mermaid style. Layout refreshed on 2026-06-19: R/12 now fans out to a two-column output row (Quarto dashboard + R/13 manifest on the left, R/14 panels + Publication outputs on the right), with a dotted cross-feed from R/14 panels into R/13 manifest. Caption updated in lockstep.
- R/11 mobilome database — confirmed as PlasmidFinder (reference 7).
- Module count refreshed to fifteen (R/00–R/14); R/13 manifest stage + new R/14 publication-panels stage + `manifest.json` (contract v1.2) added across Abstract, §1, §2, Figure 1.
- §1 ¶1 rewritten with One Health framing; new ref [1] Al-Khalaifah et al. 2025 inserted; refs 2–17 renumbered accordingly.
- §1 ¶2 limitations made concrete (comment #33 / #34 from 2026-06-15 DOCX review): removed the vague "30–50 scripts" anchor, named the symptoms (drift across script copies; no machine-readable artefact catalogue; cross-study comparisons reconcile code as well as data), and anchored the script-count claim on the concrete `chicken_batch1` ~46-script baseline.
- §2 "Module code is never edited" softened (comment #38) → "Module code is not edited per-project; dataset-specific adjustments live in `config.yaml` or in the structured fix tables it references."
- §2 `taxid_fixes.csv` description reframed (comment #41) → optional capability ("when needed") rather than a core stage feature.
- §2 Stages paragraph rewritten (comment #42) into five conceptual blocks — ingest & cleaning / gene-element normalisation / community analyses / functional profiling / integration, reporting and manifest — with module IDs in parentheses rather than as the spine.
- §2 dependency line corrected: TPM normalisation (R/03) and gene-element stages (R/09–R/11) now stated as requiring **ABRicate** outputs (not Bracken). R/03 operates on per-(sample, GENE) ABRicate hits with `GENE`/`START`/`END` columns; Bracken contributes via the R/02 merge but is not the carrier for the TPM math.
- §2 "resistome–virulome–mobilome triad" → "profiling".
- §2 utilities sentence deleted; "Shared utilities and Reporting" → "Reporting".
- §2 new sentence on `cfg$stats$padjust_method` knob (Benjamini-Hochberg default; full list in `docs/STATISTICS.md`).
- Quarto deliverable corrected from "HTML report" to "dashboard with per-domain panels".
- §3 validation: auto-instrumented per-stage runtime + `scripts/diff_normdata.R` + `scripts/probe_taxid_parse.R` now cited; Re-centrifuge column-slicing caveat noted.
- §3 Table 1 populated from the 2026-05-31 pipeline.log on both validation cohorts (chicken_batch1 PC_JC_2024-11-29 = 3m 46s; chicken_batch2 JC_FL_2025-05-20 = 6m 32s). Two-batch scaling table chosen over ground-truth-agreement table; agreement column deferred. Memory column intentionally omitted — instrumentation pending.
- §3 Table 1 R/14 panels row added as `{{TODO: re-run}}`; pre-R/14 caveat added to caption.
- Funding (Vice-Chancellor scholarship, University of West London) and Conflict of Interest (none declared) filled.
- Software references (4, 9, 13, 14, 15, 16, 17) pinned and verified on 2026-06-23 via `scripts/_collect_citations.R` (runs `citation()` under the project renv). vegan 2.7-3, ggraph 2.2.2, igraph 2.3.1, networkD3 0.4.1, R 4.5.0. ABRicate (4) and Quarto (13) carry `(accessed 2026-06-23)` since both are software with no DOI. Re-run the helper before submission to refresh access dates and catch any post-update version drift.
- Peak-memory instrumentation added to `run_pipeline.R` on 2026-06-23 (gc()-based — same methodology peakRAM uses internally, no extra dep). Each stage now logs `peak %.0f MB` on completion and the slowest-first breakdown gains a peak-MB column. Caveat: only R-managed allocation is captured, so the `report` stage (which shells out to Quarto) under-reports. Filling the Table 1 peak-memory column happens automatically on the Table 1 B1+B2 re-run.

Still open:

- `{{TODO: ORCIDs for all authors}}` — eight ORCIDs to collect.
- **Re-run Table 1 (B1 + B2)** on the post-rebuild commit so panels (R/14) runtime + the updated 15-stage totals replace the `{{TODO: re-run}}` placeholders. Other rows should also be refreshed — drift since 2026-05-31 (palette resolver, kraken2 ancestry surface, Sankey `top_samples` filter, fig08_sankey_taxon_arg_mge drop + slot renumbering, heatmap/violin/Venn polish round) may have shifted earlier stages.
- Ground-truth agreement column for chicken_batch1 — Table 1 currently reports per-batch runtime + primary output but does not yet quantify per-stage agreement (row overlap, abundance correlation, figure-by-figure identity) against the ~46 hand-edited scripts. Needs `scripts/diff_normdata.R` (and per-stage equivalents) run against the hand-edited baseline before submission.
- `{{TODO: confirm exact count}}` (resistome row delta) — memory says "~+110 (sample, GENE) rows". Confirm the precise number from the validation run before publication.
- `{{TODO: acknowledgements}}`.
- Ref 6b (VFDB foundational citation): optional — if the 2005 founding paper is also desired, add as 6b alongside the current 2022 update.

Factual claims I was not fully certain about:

- The exact upstream tool list cited above (Kraken2 + Bracken + ABRicate + Re-centrifuge + ALDEx2 + vegan + Quarto + R) matches the brief verbatim; CARD, VFDB and PlasmidFinder are pulled from `pipeline_facts.md` and cited because they are the databases the ABRicate stage targets.
- The pipeline description states "fifteen ordered modules (R/00 through R/14)" — verified against `R/` listing (00–14 inclusive plus six `utils_*` helpers).
- R/11 mobilome stage uses **PlasmidFinder** (confirmed by author). Cited as reference 7 (Carattoli et al. 2014, AAC).
- Runtime / memory characteristics for R/14 panels and the post-2026-05-31 drift are not yet measured; runtimes are routed through the benchmark-table placeholder ({{TODO: re-run}}).
- The TPM-normalisation choice is taken from `pipeline_facts.md` Stage 03 and confirmed against `R/03_normalisation.R`: TPM is computed per (sample, GENE) over ABRicate hits, so the dependency is ABRicate (not Bracken as the first draft asserted).
- The phrase "approximately 46 hand-edited R scripts" follows the memory snapshot ("~46"); confirm the exact figure before submission.
