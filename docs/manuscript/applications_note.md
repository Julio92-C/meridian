# A config-driven R workflow for long-read shotgun metagenomics: taxonomy, resistome, virulome and mobilome in one report

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

**Results.** We present Metagenomics Pipeline Automation, an alpha-stage modular R workflow that consumes Kraken2, Bracken, ABRicate and Re-centrifuge outputs and produces a single Quarto HTML dashboard covering decontamination, TPM normalisation, taxonomy, alpha/beta diversity with PERMANOVA, two-level ALDEx2 differential abundance, and resistome/virulome/mobilome profiling joined by a tripartite sample × taxon × gene network with companion chord and Sankey diagrams. A study is parameterised by one YAML file; module code is never edited per project, and every run emits a `manifest.json` artefact catalogue for downstream tooling.

**Availability and Implementation.** Source at https://github.com/Julio92-C/Metagenomics_pipeline_automation (MIT). Run with `git clone <repo> && Rscript run_pipeline.R config.yaml` after `renv::restore()`.

**Contact.** Hermine.Mkrtchyan@uwl.ac.uk

---

## 1 Introduction

Long-read shotgun metagenomics (Oxford Nanopore, PacBio) is now routinely used to profile microbial communities together with antimicrobial resistance genes (ARGs), virulence factors and mobile genetic elements (MGEs) in environmental, clinical and agricultural samples. The downstream analysis stack is well established — Kraken2 [1] and Bracken [2] for taxonomic assignment and abundance estimation, ABRicate [3] against CARD [4], VFDB [5] and PlasmidFinder [6] for gene screening, Re-centrifuge [7] for contamination filtering, and R packages such as vegan [8] and ALDEx2 [9] for ecological and differential analyses — but the integration between them is not.

In practice, each new study is analysed by adapting a set of bespoke R scripts: absolute paths are rewritten, sample identifiers are pasted into filtering expressions, control samples are renamed inside `subset()` calls, and taxonomy strings are patched by row index. We routinely observe 30–50 such scripts per study. The L-ARRAP pipeline [10] addresses one slice of this problem for long-read resistome risk assessment, and the rANOMALY workflow [11] demonstrates the value of a config-driven R workflow with auto-generated reports for amplicon data. To our knowledge, no comparable workflow exists for end-to-end long-read shotgun metagenomics covering taxonomy, diversity, differential abundance, and the full resistome–virulome–mobilome triad in one configurable run.

We present an early-stage R workflow that fills this gap. The tool is in active development; this note describes its design, current scope, and the real-data validation strategy we are using to harden each stage.

## 2 Implementation

The workflow is a stand-alone R project, not an installable package. A single entry point (`run_pipeline.R`) consumes one YAML configuration file per study and dispatches to fourteen ordered modules (`R/00_setup.R` through `R/13_manifest.R`); a Quarto dashboard template (`templates/report.qmd`) renders all outputs into a single multi-page HTML deliverable, and a final manifest stage emits `manifest.json` describing every artefact produced (Figure 1).

![**Figure 1.** Metagenomics Pipeline Automation workflow. A single `config.yaml` parameterises a study and drives fourteen ordered R modules (R/00–R/13) consuming upstream outputs (Kraken2/Bracken, ABRicate, Re-centrifuge) plus sample metadata. Community analyses (R/04–R/08) and functional profiling against CARD, VFDB and PlasmidFinder (R/09–R/11) feed a tripartite sample × taxon × gene network with companion chord and Sankey diagrams (R/12); a Quarto dashboard (R/report) and a `manifest.json` contract for downstream tooling (R/13) are emitted at the end. Source: `figure1_pipeline_schematic.mmd`.](figure1_pipeline_schematic.png)

**Configuration-driven.** Everything that varies between studies — metadata file path, grouping columns, negative-control identifiers, abundance thresholds, sample selection — is declared in `config.yaml`. Module code is never edited. R/00 (setup) and R/01 (input loading) always run; every subsequent stage can be toggled via `cfg$stages$*` flags so a user can, for example, re-run only the resistome modules after revising thresholds. TPM normalisation (R/03) and all gene-element stages (R/09–R/11) require Bracken outputs and are skipped otherwise. Taxonomy fixes that used to be row-index edits in the hand-edited reference scripts are now declared in a structured `taxid_fixes.csv` keyed by NCBI taxid.

**Stages.** R/00 handles library loading, helper sourcing and configuration validation. R/01 reads Kraken2/Bracken tables, ABRicate gene calls, Re-centrifuge contamination tags and the study metadata, accommodating per-project input shapes rather than assuming a fixed column layout. R/02 applies the structured `taxid_fixes.csv` table, performs optional negative-control subtraction, merges Bracken counts with the metadata covariates declared in `cfg$metadata$fixed_effects` / `random_effect`, and links ABRicate calls to taxonomic assignments. R/03 normalises read counts to transcripts-per-million (TPM). R/04 builds Venn diagrams and annotated heatmaps after a shared taxa-name cleanup (`R/utils_taxa.R`). R/05 produces stacked-bar relative-abundance plots (static PNG and interactive plotly HTML) with top-N collapse, plus paired total-count and sample-prevalence panels. R/06 prefers a pre-computed indices CSV (e.g. Oxford Nanopore wf-metagenomics output) and falls back to vegan otherwise, reporting Shannon/Simpson/richness/Pielou/Fisher metrics with Kruskal-Wallis tests and violin/bar plots. R/07 runs PCoA on Hellinger-transformed Bray-Curtis distances using vegan, with PERMANOVA (`adonis2`, 9999 permutations) and PERMDISP (`betadisper`) overlaid on the ordination. R/08 runs ALDEx2 differential abundance at two levels — taxa (Bracken counts) and genes (ABRicate contig hits per (sample, GENE)) — combining an omnibus Kruskal-Wallis with pairwise Welch's t tests, and supports an optional `aldex.glm` design path for covariate-adjusted models. R/09–R/11 profile the resistome (ABRicate/CARD, with drug-class classification and an MLS rollup), virulome (VFDB, with functional categories extracted from the ABRicate `PRODUCT` field) and mobilome (PlasmidFinder, with configurable replicon-family classification — Col-like, IncF, IncX, other Inc); each emits its own alpha diversity, Venn, gene- and category-level heatmaps, relative-abundance and total-count panels, and a PCoA + PERMANOVA + PERMDISP suite. R/12 assembles a tripartite ggraph network of sample × taxon × gene plus chord diagrams (overall and per group) and four-tier networkD3 Sankey diagrams, exports Gephi-compatible node/edge tables, and reports igraph topology metrics (degree, betweenness, Louvain modularity) and Bray-Curtis sample clusters. R/13 finally writes a `manifest.json` describing every artefact, schema and stage status, providing a stable JSON contract that decouples the pipeline from downstream report consumers and manuscript-drafting tools.

**Shared utilities and reporting.** Three helper files keep stage code thin: `R/utils_taxa.R` (configurable name cleanup), `R/utils_ge_profile.R` (the shared alpha/beta/Venn/heatmap/relative-abundance engine driving R/09–R/11) and `R/utils_prevalence.R` (paired total/sample-prevalence panels used by R/05 and R/09–R/11). The Quarto template renders all outputs as a multi-page HTML dashboard with per-domain panels (taxonomy, diversity, differential abundance, resistome, virulome, mobilome, network) and embedded download links. Dependencies are pinned via `renv` and the project requires R ≥ 4.5; external command-line tools (Kraken2, Bracken, ABRicate, Re-centrifuge) are assumed to have been run upstream.

## 3 Application and Validation

We are validating the pipeline against two real long-read chicken caecum metagenomics studies. `chicken_batch1` (PC_JC_2024-11-29; 18 samples; 3 dietary groups — Dulce, Reference diet, Soyabean meal) is the reference cohort: approximately 46 hand-edited R scripts and known-good figures and tables already exist for it, so each stage can be compared output-for-output against the corresponding ground-truth artefact. `chicken_batch2` (JC_FL_2025-05-20; 32 samples; 4-level composite `Treatment_Bird` grouping — Control_W4, Control_W5, Dulse_W4, Dulse_W5) is a scaling and feature-coverage cohort: it exercises a random effect (`block`, pens 1–8), two declared negative controls (CTRL1, CTRL2) subtracted in R/02, and a larger per-stage feature load. Two purpose-built validation scripts (`scripts/diff_normdata.R` and `scripts/probe_taxid_parse.R`) provide numerical-diff and regression checks on the most error-prone intermediates, and `run_pipeline.R` auto-instruments per-stage wall time and prints a slowest-first breakdown to the run log, mechanically populating the runtime columns of Table 1.

**Table 1.** Per-stage wall time and primary output for both validation cohorts on the same workstation (Windows 11, R 4.5, single thread); identical pipeline commit and identical statistical settings (PERMANOVA 9999 permutations; ALDEx2 `mc.samples` = 128). **B1** = chicken_batch1 (PC_JC_2024-11-29; n=18; 3 dietary groups). **B2** = chicken_batch2 (JC_FL_2025-05-20; n=32; 4-level `Treatment_Bird` grouping, `block` random effect, two declared negative controls). Sample sizes shown in some output cells refer to samples with non-zero hits at that stage.

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
| report | 25.1 s | Quarto dashboard HTML (~14 MB) | 40.5 s | Quarto dashboard HTML (~18 MB) |
| manifest | 0.5 s | `manifest.json` (47 kB) | 0.4 s | `manifest.json` (57 kB) |
| **Total (14 stages)** | **3m 46s** | — | **6m 32s** | — |

Differential abundance dominates wall time on both cohorts (≈50 %), driven by the ALDEx2 Monte-Carlo loop scaling with feature count, sample count and pairwise-contrast count; the network stage is next-most expensive because Sankey assembly iterates over every gene–category mapping. Across the 14 stages, doubling the sample count and pairwise-contrast count (B1 → B2) raises end-to-end runtime by ~73 %, suggesting near-linear scaling on the dominant ALDEx2 path. Tabulated PERMANOVA and PERMDISP statistics for the per-domain stages (resistome, virulome, mobilome), per-metric KW p-values for alpha diversity, and all DAA per-pair feature tables are surfaced in the Quarto dashboard rather than reproduced here. Peak memory is not yet auto-instrumented and is therefore omitted from Table 1; a separate `peakRAM`-based profiling pass is planned before submission.

One expected and deliberate divergence concerns the resistome stage: the ground-truth `cleanData.R` mis-parses taxids carrying variant suffixes, dropping the corresponding ABRicate rows. R/02 fixes this parsing bug, so the pipeline retains approximately {{TODO: confirm exact count from a fresh validation run}} additional (sample, GENE) rows relative to the hand-edited output. We treat this as a correctness improvement rather than a failure to reproduce, and the validation log flags it accordingly. The Re-centrifuge column-slicing heuristic in R/02 is acknowledged as still coarse and is queued for rework; Table 1 will be re-derived once that change lands.

## 4 Conclusion

Metagenomics Pipeline Automation is an alpha-stage workflow that consolidates a long-read shotgun metagenomics analysis — taxonomy, diversity, differential abundance, resistome, virulome, mobilome, and an integrated network view — into one configuration-driven R project with a single Quarto dashboard and a downstream-readable `manifest.json` contract. It is in active development and validation against real reference data is ongoing; we release it now to invite feedback from groups running comparable analyses and to encourage convergence on a shared, scriptable layout for long-read metagenomics studies.

---

**Acknowledgements.** {{TODO: acknowledgements — collaborators, reviewers, compute providers}}

**Funding.** {{TODO: funding sources and grant numbers}}

**Conflict of Interest.** {{TODO: declare conflicts (or state "none declared")}}

---

## References

1. Wood DE, Salzberg SL. Kraken: ultrafast metagenomic sequence classification using exact alignments. *Genome Biol* 2014;15:R46. doi:10.1186/gb-2014-15-3-r46. Kraken2: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol* 2019;20:257. doi:10.1186/s13059-019-1891-0.
2. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci* 2017;3:e104. doi:10.7717/peerj-cs.104.
3. Seemann T. ABRicate: mass screening of contigs for antimicrobial and virulence genes. https://github.com/tseemann/abricate. {{TODO: verify ref — software, no DOI}}
4. Alcock BP, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res* 2023;51:D690–D699. doi:10.1093/nar/gkac920.
5. Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Res* 2022;50:D912–D917. doi:10.1093/nar/gkab1107.
6. Carattoli A, Zankari E, García-Fernández A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014;58:3895–3903. doi:10.1128/AAC.02412-14.
7. Martí JM. Recentrifuge: robust comparative analysis and contamination removal for metagenomics. *PLoS Comput Biol* 2019;15:e1006967. doi:10.1371/journal.pcbi.1006967.
8. Oksanen J, et al. vegan: Community Ecology Package. R package. {{TODO: verify ref — pin version cited}}
9. Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome* 2014;2:15. doi:10.1186/2049-2618-2-15.
10. Li Y, Gao Y, Liu X, Mao Y, Wang M, Qin Y, Zhang C, Chen Q, Ning K, Wang Z, Han M. Quantifying antibiotic resistome risks across environmental niches: the L-ARRAP for long-read metagenomic profiling. *Brief Bioinform* 2025;26(5):bbaf535. doi:10.1093/bib/bbaf535.
11. Theil S, Rifa E. rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis. *F1000Res* 2021;10:7. doi:10.12688/f1000research.27268.1.
12. Allaire JJ, et al. Quarto. https://quarto.org. {{TODO: verify ref — software, no DOI}}
13. R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/. {{TODO: verify ref — pin version cited (≥4.5)}}
14. Pedersen TL. ggraph: An Implementation of Grammar of Graphics for Graphs and Networks. R package. {{TODO: verify ref — pin version via `citation("ggraph")`}}
15. Csárdi G, Nepusz T. The igraph software package for complex network research. *InterJournal* 2006; Complex Systems:1695. https://igraph.org. {{TODO: verify ref — confirm citation format from `citation("igraph")`}}
16. Allaire JJ, Gandrud C, Russell K, Yetman CJ. networkD3: D3 JavaScript Network Graphs from R. R package. {{TODO: verify ref — pin version via `citation("networkD3")`}}

---

## Drafting notes

Placeholders inserted, and the data needed to fill each one:

Resolved since first draft:

- Author block, affiliations, corresponding author and email — filled.
- GitHub URL — filled (https://github.com/Julio92-C/Metagenomics_pipeline_automation).
- Figure 1 — schematic created at `figure1_pipeline_schematic.mmd` (Mermaid source). Rendered to PNG (`figure1_pipeline_schematic.png`) via `npx --yes -p @mermaid-js/mermaid-cli mmdc -i figure1_pipeline_schematic.mmd -o figure1_pipeline_schematic.png -w 1800 -b white`. Regenerate (or render to SVG) from the `.mmd` source as needed.
- R/11 mobilome database — confirmed as PlasmidFinder (reference 6).
- Module count refreshed to fourteen (R/00–R/13); new R/13 manifest stage and `manifest.json` contract added across Abstract, §2, Figure 1.
- §2 stage paragraph rewritten against the 2026-05-29 audit (`manuscript_edit_plan.md`): R/02 (taxid_fixes + Bracken merge), R/06 (precomputed indices priority), R/07 (Hellinger-Bray-Curtis + PERMANOVA + PERMDISP), R/08 (two-level + omnibus KW + pairwise + aldex.glm option), R/09 (drug-class + MLS rollup, ResFinder dropped), R/10 (VFDB function extraction), R/11 (replicon families), R/12 (network + chord + Sankey + Gephi + igraph topology).
- Quarto deliverable corrected from "HTML report" to "dashboard with per-domain panels".
- Utility files surfaced (`utils_taxa.R`, `utils_ge_profile.R`, `utils_prevalence.R`).
- §3 validation: auto-instrumented per-stage runtime + `scripts/diff_normdata.R` + `scripts/probe_taxid_parse.R` now cited; Re-centrifuge column-slicing caveat noted.
- §3 Table 1 populated from the 2026-05-31 pipeline.log on both validation cohorts (chicken_batch1 PC_JC_2024-11-29 = 3m 46s; chicken_batch2 JC_FL_2025-05-20 = 6m 32s). Two-batch scaling table chosen over ground-truth-agreement table; agreement column deferred. Memory column intentionally omitted — instrumentation pending.

Still open:

- `{{TODO: ORCIDs for all authors}}` — eight ORCIDs to collect.
- Ground-truth agreement column for chicken_batch1 — Table 1 currently reports per-batch runtime + primary output but does not yet quantify per-stage agreement (row overlap, abundance correlation, figure-by-figure identity) against the ~46 hand-edited scripts. Needs `scripts/diff_normdata.R` (and per-stage equivalents) run against the hand-edited baseline before submission.
- Peak memory column — `run_pipeline.R` does not currently auto-instrument memory; needs a `peakRAM`-based pass on both cohorts before submission.
- `{{TODO: confirm exact count}}` (resistome row delta) — memory says "~+110 (sample, GENE) rows". Confirm the precise number from the validation run before publication.
- `{{TODO: acknowledgements}}`, `{{TODO: funding sources and grant numbers}}`, `{{TODO: declare conflicts}}`.
- `{{TODO: verify ref}}` — remaining open items after PubMed verification on 2026-05-29:
  - Ref 3 (ABRicate): no PubMed entry — software-only. Current GitHub URL citation is correct; consider adding `(accessed YYYY-MM-DD)`.
  - Ref 5 (VFDB): **resolved** — now cites Liu et al. 2022 (NAR 50:D912–D917, doi:10.1093/nar/gkab1107), the latest VFDB update (PMID 34850947). If the 2005 foundational citation is also desired, add it as 5b.
  - Ref 8 (vegan): no PubMed entry — R package. Use the canonical `citation("vegan")` output and pin the version actually used (`renv.lock` records this).
  - Ref 9 (ALDEx2): **resolved** — Fernandes et al. 2014, *Microbiome* 2:15 (doi:10.1186/2049-2618-2-15, PMID 24910773) confirmed as the canonical ALDEx2 citation by PubMed.
  - Ref 10 (L-ARRAP): **resolved** — author list and full title corrected from PubMed (PMID 41066697).
  - Ref 12 (Quarto): software-only — use `citation("quarto")` from R, or cite the website with access date.
  - Ref 13 (R Core Team): use the output of `citation()` in R; pin to the version in `renv.lock`.
  - Refs 14–16 (ggraph / igraph / networkD3): added so the network/Sankey/topology vocabulary in §2 has citations. Pin versions via `citation("ggraph")`, `citation("igraph")`, `citation("networkD3")` and confirm formatting before submission.

Factual claims I was not fully certain about:

- The exact upstream tool list cited above (Kraken2 + Bracken + ABRicate + Re-centrifuge + ALDEx2 + vegan + Quarto + R) matches the brief verbatim; CARD, VFDB and PlasmidFinder are pulled from `pipeline_facts.md` and cited because they are the databases the ABRicate stage targets.
- The pipeline description states "fourteen ordered modules (R/00 through R/13)" — verified against `R/` listing (00–13 inclusive plus three `utils_*` helpers).
- R/11 mobilome stage uses **PlasmidFinder** (confirmed by author). Cited as reference 6 (Carattoli et al. 2014, AAC).
- Runtime / memory characteristics are not in the facts file, so I made no claim about them in the body text and routed them through the benchmark-table placeholder.
- The TPM-normalisation choice is taken from `pipeline_facts.md` Stage 03; no further reference is cited for TPM in metagenomics — the user may want to add one if a specific method paper is being followed.
- The phrase "approximately 46 hand-edited R scripts" follows the memory snapshot ("~46"); confirm the exact figure before submission.
