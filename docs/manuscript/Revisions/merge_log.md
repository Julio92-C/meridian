# MERIDIAN Applications Note — revision merge log

**Date:** 2026-09-01
**Baseline:** `applications_note.md` @ b7b1391 (release v0.1.0 pin)
**Sources merged:** `Revisions/applications_note_HM.docx` (HM, 2026-08-31 — 46 tracked insertions + 25 comments) and `Revisions/applications_note_PL.docx` (PL, 2026-07-24 — 28 ins / 36 del + 10 comments)
**Strategy (agreed):** adopt HM's full section rewrites as the new prose, layer PL's non-conflicting edits, reconcile conflicts case-by-case. Target journal: **Briefings in Bioinformatics**.

---

## Decisions locked in

| # | Fork | Decision |
|---|------|----------|
| 1 | Merge base | HM rewrite as base + PL layered |
| 2 | Author list | HM's 4-author list: Ortega Cambara, Cuber, Lebre, Mkrtchyan ("Hermine V Mkrtchyan"). Affiliations renumbered to 1 = new Centre (CIGMiS / School of Medicine and Bioscience, UWL), 2 = Natural History Museum (Cuber). **Julio to confirm the four removed authors (Misra, McHugh, Rofael, Lowe) have agreed to be dropped / moved to acknowledgements.** |
| 3 | Abstract lead | HM's Background/Results/Availability/Corresponding-author structure, with PL's "resource-heavy pipelines requiring HPC access" point folded into Background. |
| 4 | Journal | Briefings in Bioinformatics (HM). Abstract structure already matches this Applications-Note style; "Implementation" heading retained (valid for BiB). |

## Changes applied (mechanical / editorial — done)

- **Author block + affiliation** per decision 2.
- **Abstract** rewritten to HM structure; "Motivation."→"Background:", "Contact."→"Corresponding author:"; PL HPC point + "sustainable" blended in; typos cleaned ("a a single", "generated"→"generates", "characterization"→"characterisation", "two-levels"→"two-level", "networks"→"network", "mobile genetic elements repertoires"→"…element repertoires").
- **Availability** de-duplicated (HM had renv::restore twice); kept the single renv sentence; "license"→"licence".
- **Introduction** paras 1–3 to HM prose; "reproducible, and reusable"→"reproducible and reusable".
- **Novelty claim narrowed** (Intro para 2): HM wanted "this is the first…"; PL objected that broad-coverage pipelines exist (nf-core/funcscan). Reconciled to "no integrated end-to-end **R** workflow covers … in a single configurable run" — true and defensible (funcscan is a Nextflow pipeline, not an R workflow, and does not span taxonomy+diversity+DA+RVM in one run). See open item O2.
- **Implementation** para: "single entry script that reads…", "complies"→"compiles", "publication-panels workflow"→"stage".
- **Configuration-driven** para: fixed garbled HM merge ("through live in config.yaml"→"through config.yaml settings"); kept accurate "gene-element stages (R/09–R/11)" rather than HM's "gene-mobile genetic element integration stages"; removed stray trailing ")".
- **Stages**: labels + connective prose to HM; "(1) Data input and decontamination", tense fixed (imports/applies/merges/links); "pipieline"→"pipeline"; stage-3 clarified what KW vs ALDEx2 do (HM comment); stage-4 "(sample, GENE))"→"(sample, GENE)", stray ")" removed; stage-5 added "integrates results across all domains" opener.
- **Reporting** para to HM.
- **Validation** para 1 to HM prose; **removed the unverified "publically available" claim** → replaced with a {{TODO}} flag (see O1); tense normalised.
- **Results** para: "cohorts"→"datasets", em-dash→comma per HM.
- **Divergence** para: HM synonyms applied (divergence→discrepancy, fixes→corrects, gained→recovered, Across→For); "The R/02"→"R/02".
- **Conclusion** to HM prose + HM's preferred closing sentence ("…to solicit feedback from researchers conducting comparable long-read metagenomic studies…"); kept "(v0.1.0)".
- **Acknowledgements**: "School of Biomedical Science"→"School of Medicine and Bioscience" for consistency with the new affiliation.
- **PL formatting fixes** (renv::restore backticks) already subsumed by the rewritten Availability sentence.

## Open items — need Julio's decision / real work (NOT auto-merged)

These are substantive scientific/authorship points raised in comments. They cannot be resolved by merging text and are deliberately left for you.

- **O1 — Public data / reproducibility (BOTH reviewers, high priority). RESOLVED 2026-09-01.** Batch 1 is publicly available at NCBI BioProject **PRJNA1406192** (already cited by a published paper); Batch 2 release is pending publication of its associated study. §3 now states the accession; the `{{TODO}}` is removed; a **Data Availability** statement was added to the back matter. Mock/curated-community validation (PL's stronger suggestion) remains a future-work item — see O3.
- **O2 — Count transformation for differential abundance (PL, high). RESOLVED 2026-09-02.** Added a justification sentence at the end of §2 stage (3): ALDEx2 uses a centred log-ratio (CLR) transform (avoiding raw/total-count artefacts), and gene features are length-normalised as TPM (R/03). Cites the existing Fernandes 2014 [10] reference; no new citation added. This states the pipeline already applies the compositionally-correct transform PL was asking for.
- **O3 — Validation vs the same scripts (HM + PL).** PL: running data through the pipeline's own descendant scripts is not independent validation; suggests cross-checking against an independent curated pipeline. HM echoes "test on public datasets". Scope decision for you.
- **O4 — Local-run feasibility / Kraken2 memory (PL). RESOLVED 2026-09-02.** Added a sentence to the §2 Reporting paragraph: MERIDIAN consumes pre-computed outputs, so its footprint is modest/single-threaded (~1.5 GB peak; §3) and decoupled from the heavy memory of upstream classification (e.g. a full Kraken2 DB), which can run on separate high-memory/cloud infrastructure. Directly answers the low-spec-PC question.
- **O5 — Pipeline name uniqueness / CR (PL).** Confirm "MERIDIAN" is distinct and searchable; check for trademark/name collisions before submission.
- **O6 — Author removals (decision 2).** Confirm Misra, McHugh, Rofael, Lowe agree to removal / acknowledgement move; update Acknowledgements and Funding/COI accordingly.
- **O7 — ORCIDs** still outstanding ({{TODO}} in author block).

## Comments already addressed by HM's own rewrite (no action)

HM's inserted prose already resolved these HM comments: define AMR abbreviation (now spelled out in abstract); opening sentence for community-analyses stage; clarify KW vs statistical comparisons (stage 3); "network analysis" phrasing (stage 5); MIT-licence note; renv reproducibility sentence.
