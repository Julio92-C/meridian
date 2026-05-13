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

**Results.** We present Metagenomics Pipeline Automation, an alpha-stage modular R workflow that consumes Kraken2, Bracken, ABRicate and Re-centrifuge outputs and produces a single Quarto HTML report covering decontamination, TPM normalisation, taxonomy, alpha/beta diversity, ALDEx2 differential abundance, and resistome/virulome/mobilome profiling with a tripartite taxa-ARG-MGE network. A study is parameterised by one YAML file; module code is never edited per project.

**Availability and Implementation.** Source at https://github.com/Julio92-C/Metagenomics_pipeline_automation (MIT). Run with `git clone <repo> && Rscript run_pipeline.R config.yaml` after `renv::restore()`.

**Contact.** Hermine.Mkrtchyan@uwl.ac.uk

---

## 1 Introduction

Long-read shotgun metagenomics (Oxford Nanopore, PacBio) is now routinely used to profile microbial communities together with antimicrobial resistance genes (ARGs), virulence factors and mobile genetic elements (MGEs) in environmental, clinical and agricultural samples. The downstream analysis stack is well established — Kraken2 [1] and Bracken [2] for taxonomic assignment and abundance estimation, ABRicate [3] against CARD [4], VFDB [5] and PlasmidFinder [6] for gene screening, Re-centrifuge [7] for contamination filtering, and R packages such as vegan [8] and ALDEx2 [9] for ecological and differential analyses — but the integration between them is not.

In practice, each new study is analysed by adapting a set of bespoke R scripts: absolute paths are rewritten, sample identifiers are pasted into filtering expressions, control samples are renamed inside `subset()` calls, and taxonomy strings are patched by row index. We routinely observe 30–50 such scripts per study. The L-ARRAP pipeline [10] addresses one slice of this problem for long-read resistome risk assessment, and the rANOMALY workflow [11] demonstrates the value of a config-driven R workflow with auto-generated reports for amplicon data. To our knowledge, no comparable workflow exists for end-to-end long-read shotgun metagenomics covering taxonomy, diversity, differential abundance, and the full resistome–virulome–mobilome triad in one configurable run.

We present an early-stage R workflow that fills this gap. The tool is in active development; this note describes its design, current scope, and the real-data validation strategy we are using to harden each stage.

## 2 Implementation

The workflow is a stand-alone R project, not an installable package. A single entry point (`run_pipeline.R`) consumes one YAML configuration file per study and dispatches to thirteen ordered modules (`R/00_setup.R` through `R/12_network.R`); a Quarto template (`templates/report.qmd`) renders all outputs into a single HTML report (Figure 1).

![**Figure 1.** Metagenomics Pipeline Automation workflow. A single `config.yaml` parameterises a study and drives thirteen ordered R modules (R/00–R/12) consuming upstream outputs (Kraken2/Bracken, ABRicate, Re-centrifuge) plus sample metadata. Community analyses (R/04–R/08) and functional profiling against CARD/ResFinder, VFDB and PlasmidFinder (R/09–R/11) feed a tripartite taxa × ARG × MGE network (R/12); all outputs are rendered into a single Quarto HTML report. Source: `figure1_pipeline_schematic.mmd`.](figure1_pipeline_schematic.png)

**Configuration-driven.** Everything that varies between studies — metadata file path, grouping columns, negative-control identifiers, abundance thresholds, sample selection — is declared in `config.yaml`. Module code is never edited. Each stage can be toggled via `cfg$stages$*` flags so a user can, for example, re-run only the resistome modules after revising thresholds. This contrasts sharply with the hand-edited reference scripts (Section 3), where per-study sample identifiers and control names are embedded throughout the analysis code.

**Stages.** R/00 handles library loading, helper sourcing and configuration validation. R/01 reads Kraken2/Bracken tables, ABRicate gene calls, Re-centrifuge contamination tags and the study metadata, accommodating per-project input shapes rather than assuming a fixed column layout. R/02 performs optional negative-control subtraction, applies a shared taxid-cleanup helper (`R/utils_taxa.R`) and links ABRicate calls to taxonomic assignments. R/03 normalises read counts to transcripts-per-million (TPM). R/04–R/07 produce taxonomy Venn diagrams and heatmaps, stacked-bar relative-abundance plots (static PNG and interactive plotly HTML), Shannon/Simpson/richness diversity with violin and bar plots, and PCoA ordinations using vegan. R/08 runs ALDEx2 for overall and pairwise differential abundance on compositional data. R/09–R/11 profile the resistome (ABRicate/CARD/ResFinder), virulome (VFDB) and mobilome (PlasmidFinder [6]) respectively. R/12 assembles a tripartite chord diagram / network linking taxa, ARGs and MGEs.

**Reporting and reproducibility.** The Quarto report template stitches every figure and table into a single HTML deliverable per study. Dependencies are pinned via `renv` and the project requires R ≥ 4.5; external command-line tools (Kraken2, Bracken, ABRicate, Re-centrifuge) are assumed to have been run upstream.

## 3 Application and Validation

We are validating each module against a real long-read chicken gut microbiome dataset (`chicken_batch1`, internal project code PC_JC_2024-11-29) for which approximately 46 hand-edited R scripts and known-good figures and tables already exist. Each pipeline stage is compared, output-for-output, against the corresponding ground-truth artefact, with the configuration file set to reproduce the original analysis decisions.

{{TODO: insert benchmark table — stage × ground-truth artefact × pipeline output × agreement metric (e.g. row overlap, correlation of abundance values, identical figure rendering). Numbers to be populated from validation runs.}}

One expected and deliberate divergence concerns the resistome stage: the ground-truth `cleanData.R` mis-parses taxids carrying variant suffixes, dropping the corresponding ABRicate rows. R/02 fixes this parsing bug, so the pipeline retains approximately {{TODO: confirm exact count}} additional (sample, GENE) rows relative to the hand-edited output. We treat this as a correctness improvement rather than a failure to reproduce, and the validation table flags it accordingly.

End-to-end runtime, peak memory and the per-stage outputs that have so far reached parity with ground truth are reported in the benchmark table.

## 4 Conclusion

Metagenomics Pipeline Automation is an alpha-stage workflow that consolidates a long-read shotgun metagenomics analysis — taxonomy, diversity, differential abundance, resistome, virulome, mobilome, and an integrated network view — into one configuration-driven R project with a single Quarto report. It is in active development and validation against real reference data is ongoing; we release it now to invite feedback from groups running comparable analyses and to encourage convergence on a shared, scriptable layout for long-read metagenomics studies.

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
5. Chen L, Yang J, Yu J, et al. VFDB: a reference database for bacterial virulence factors. *Nucleic Acids Res* 2005;33:D325–D328. doi:10.1093/nar/gki008. {{TODO: verify ref — confirm latest edition cited}}
6. Carattoli A, Zankari E, García-Fernández A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014;58:3895–3903. doi:10.1128/AAC.02412-14.
7. Martí JM. Recentrifuge: robust comparative analysis and contamination removal for metagenomics. *PLoS Comput Biol* 2019;15:e1006967. doi:10.1371/journal.pcbi.1006967.
8. Oksanen J, et al. vegan: Community Ecology Package. R package. {{TODO: verify ref — pin version cited}}
9. Fernandes AD, Reid JN, Macklaim JM, et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome* 2014;2:15. doi:10.1186/2049-2618-2-15. {{TODO: verify ref — confirm canonical ALDEx2 citation}}
10. Li X, et al. L-ARRAP: Long-read Antibiotic Resistome Risk Assessment Pipeline. *Brief Bioinform* 2025;26(5):bbaf535. doi:10.1093/bib/bbaf535.
11. Theil S, Rifa E. rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis. *F1000Res* 2021;10:7. doi:10.12688/f1000research.27268.1.
12. Allaire JJ, et al. Quarto. https://quarto.org. {{TODO: verify ref — software, no DOI}}
13. R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/. {{TODO: verify ref — pin version cited (≥4.5)}}

---

## Drafting notes

Placeholders inserted, and the data needed to fill each one:

Resolved since first draft:

- Author block, affiliations, corresponding author and email — filled.
- GitHub URL — filled (https://github.com/Julio92-C/Metagenomics_pipeline_automation).
- Figure 1 — schematic created at `figure1_pipeline_schematic.mmd` (Mermaid source). Rendered to PNG (`figure1_pipeline_schematic.png`) via `npx @mermaid-js/mermaid-cli mmdc -i figure1_pipeline_schematic.mmd -o figure1_pipeline_schematic.png -w 1800 -b white`. Regenerate (or render to SVG) from the `.mmd` source as needed.
- R/11 mobilome database — confirmed as PlasmidFinder (reference 6).

Still open:

- `{{TODO: ORCIDs for all authors}}` — eight ORCIDs to collect.
- `{{TODO: insert benchmark table}}` — quantitative comparison of pipeline output vs hand-edited chicken_batch1 ground truth, per stage. Suggested columns: stage, ground-truth artefact, pipeline output, agreement metric, notes. Needs end-to-end runtime, peak memory, per-stage agreement, and the exact resistome row-count delta.
- `{{TODO: confirm exact count}}` (resistome row delta) — memory says "~+110 (sample, GENE) rows". Confirm the precise number from the validation run before publication.
- `{{TODO: acknowledgements}}`, `{{TODO: funding sources and grant numbers}}`, `{{TODO: declare conflicts}}`.
- `{{TODO: verify ref}}` — references 3 (ABRicate, software, no formal paper), 5 (VFDB — confirm whether to cite the 2005 original or a more recent update), 8 (vegan — pin the cited version), 9 (ALDEx2 — confirm the canonical citation, possibly Fernandes 2013 *PLoS ONE* instead of 2014 *Microbiome*), 12 (Quarto — software citation), 13 (R — version-pinned citation).

Factual claims I was not fully certain about:

- The exact upstream tool list cited above (Kraken2 + Bracken + ABRicate + Re-centrifuge + ALDEx2 + vegan + Quarto + R) matches the brief verbatim; CARD, VFDB and PlasmidFinder are pulled from `pipeline_facts.md` and cited because they are the databases the ABRicate stage targets.
- The pipeline description states "thirteen ordered modules (R/00 through R/12)" — this matches the stage table in `pipeline_facts.md` (00 through 12 inclusive = 13 modules).
- R/11 mobilome stage uses **PlasmidFinder** (confirmed by author). Cited as reference 6 (Carattoli et al. 2014, AAC).
- Runtime / memory characteristics are not in the facts file, so I made no claim about them in the body text and routed them through the benchmark-table placeholder.
- The TPM-normalisation choice is taken from `pipeline_facts.md` Stage 03; no further reference is cited for TPM in metagenomics — the user may want to add one if a specific method paper is being followed.
- The phrase "approximately 46 hand-edited R scripts" follows the memory snapshot ("~46"); confirm the exact figure before submission.
