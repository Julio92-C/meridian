---
name: Review_log
description: Reviewer-style review of applications_note.md (MERIDIAN Applications Note)
date: 2026-06-24
journal: Bioinformatics (OUP) — Applications Note
english: British
scope: whole manuscript
---

# MERIDIAN Applications Note — Reviewer-style review

**File reviewed:** `docs/manuscript/applications_note.md` (canonical source; `.docx` is rendered)
**Date:** 2026-06-24
**English variant:** British
**Target journal:** OUP *Bioinformatics* — Applications Note (≤2 pages; ~1,200–1,500 words excl. refs; abstract ≤150 words)
**Scope:** Whole manuscript

## Summary

The draft is well organised, technically dense and the upstream/downstream logic (config → modules → manifest → Quarto + panels) reads coherently. The dominant issue for an Applications Note is **length**: the body runs ~2,116 words against the ≤1,500-word budget (+41 %) and the Abstract runs 205 words against the ≤150-word cap (+37 %). The Conclusion and Implementation each carry one paragraph that could be cut without information loss. Secondary issues: a recurring "Dulce" → "Dulse" typo, a confusing R/13 vs R/14 narrative order, an overloaded Table 1 caption (which doubles as Figure-2 prose), and a malformed Kraken / Kraken2 combined reference. No fabricated-numbers issues found — every figure traces to Table 1 or `manifest.json`.

## Findings

| # | Section | Issue | Severity | Suggested action |
|---|---|---|---|---|
| 1 | Whole | Body ~2,116 words vs ≤1,500 budget; Abstract ~205 vs ≤150 | High | Set explicit cut targets per section and prune in passes 2–7 |
| 2 | §3 Application | "Dulce" should be "Dulse" (*Palmaria palmata*) | High | Find/replace within §3 |
| 3 | §2 / Fig 1 caption | Module-order confusion: numerical R/13 → R/14 vs execution order R/14 panels → R/13 manifest | High | State once, up-front, that the manifest stage runs last by design |
| 4 | Abstract | Results sub-section overlong; technical bolt-ons ("contract v1.2", "tagging every figure and table by kind for downstream tooling") | High | Rewrite to ~95 words on Results; drop contract version and tooling clause |
| 5 | §1 Introduction | Opening sentence 53 words, multi-clause; reads as front-loaded | Medium | Split; lead with the One Health framing, follow with the technology |
| 6 | §1 Introduction | "for antimicrobial resistance genes (ARGs), virulence factors, and mobile genetic elements (MGEs) gene screening" — "gene screening" trails | Medium | Reword: "for screening antimicrobial resistance genes (ARGs), virulence factors and mobile genetic elements (MGEs)" |
| 7 | §2 Implementation | Stage (3) is a single sentence of ~140 words covering five distinct analyses | Medium | Split into 3–4 sentences or convert to a tight list |
| 8 | §2 Implementation | Stage (5) repeats per-stage details already in stages (1)–(4); also re-describes R/12, R/13, R/14 prose already given in the opening paragraph and Fig 1 caption | Medium | Compress §2 stage (5) to two sentences; let Fig 1 carry the architecture |
| 9 | §3 Application | Table 1 caption (~135 words) carries Figure-2 description, methodology footnotes and cohort descriptions | Medium | Move Figure-2 description into Fig 2 caption only; trim Table 1 caption to ≤60 words |
| 10 | §3 Application | Post-table paragraph re-explains peak-memory tracking already in the Table 1 caption ("equivalent to peakRAM methodology" appears twice) | Medium | Remove the duplicate from the post-table paragraph |
| 11 | §3 Application | Long compound sentence on Jaccard agreement / sample-set comparison mixes resistome, alpha, network and "deferred" stages | Medium | Split into a per-stage micro-table or two short sentences |
| 12 | §4 Conclusion | 175 words for an Applications Note (typical ≤80); ends aspirationally ("encourage convergence on a shared, scriptable layout…") | Medium | Trim to one tight paragraph; keep release intent, drop community-building phrasing |
| 13 | References | Ref 2 combines Kraken (2014) and Kraken2 (2019) into one entry | High | Split into two numbered references; only Kraken2 needs citation if Kraken1 is not cited elsewhere |
| 14 | References | Style inconsistency: refs 5, 9 use "et al." vs full author list; ref 9 lists 11 authors then "et al." | Low | Pick one rule (Bioinformatics: ≤6 authors list all; >6 list first 6 + et al.) |
| 15 | §3 Application | `chicken_batch1`, `chicken_batch2` slugs read as internal project IDs | Low | Replace body mentions with "the reference cohort" / "the scaling cohort"; keep slugs in Table 1 footer only |
| 16 | §3 Application | "Figure 2 visualises the same numbers" — light filler | Low | "Figure 2 plots the same data." or drop |
| 17 | §1 / §2 | "long-read shotgun metagenomics" vs "long-read metagenomics" alternates | Low | Pick "long-read shotgun metagenomics" throughout; drop the bare form |
| 18 | Whole | Phrase "alpha-stage" appears in Abstract, §1, §2, §4 | Low | Keep in §1 + §4 only; drop from Abstract Results and §2 opener |

## Detailed findings

### Finding 1 — Length overrun
**Section / paragraph:** Whole manuscript
**Why it matters:** OUP *Bioinformatics* desk-rejects Applications Notes that materially exceed the page budget. Current body excl. refs/tables/captions is ~2,116 words (target ≤1,500); Abstract is 205 words (target ≤150).
**Current section totals:**
- Abstract: 205 words → cut **~55 words**
- §1 Introduction: 273 words → cut **~30 words** (down to ~240)
- §2 Implementation: 647 words → cut **~150 words** (down to ~500)
- §3 Application: 704 words → cut **~250 words** (down to ~450)
- §4 Conclusion: 175 words → cut **~95 words** (down to ~80)
**Proposed action:** Use findings 4, 7, 8, 9, 10, 11, 12 as the trim levers; do not attempt a uniform shave.

### Finding 2 — "Dulce" / "Dulse" typo
**Section / paragraph:** §3 ¶1
**Reviewer-style comment:** The cohort description lists dietary groups as "Dulce, Reference diet, Soyabean meal", but the test ingredient is dulse (*Palmaria palmata*), and "Dulse" is used correctly later (chicken_batch2: "Dulse_W4, Dulse_W5").
**Proposed action:** Replace "Dulce" with "Dulse" in §3 ¶1.
**Status:** Done — 2026-06-24 (line 62).

### Finding 3 — Module-order narrative confusion
**Section / paragraph:** §2 ¶1, §2 stage (5), Fig 1 caption
**Reviewer-style comment:** §2 ¶1 says "fifteen ordered modules (R/00_setup.R through R/14_panels.R)" — implying R/14 runs last. Stage (5) then says "R/14 assembles publication-ready… R/13 finally writes manifest.json" — manifest runs last. Fig 1 caption agrees with the latter. A reviewer encountering "ordered R/00–R/14" and then a R/13-after-R/14 description will read it as an error.
**Proposed action:** In §2 ¶1, replace "fifteen ordered modules (R/00_setup.R through R/14_panels.R)" with a phrase that makes the execution order explicit, e.g. "fifteen modules (R/00_setup.R through R/14_panels.R; R/13 writes the run manifest last, after R/14 panel assembly)". Then drop the redundant clarification in stage (5).
**Status:** Done — 2026-06-24. Two-edit commit: §2 ¶1 reworded with parenthetical execution-order note; §2 stage (5) opening reworded from "(R/12–R/14)" to "(R/12, then R/14, then R/13)". Fig 1 caption already correct.

**Follow-up — 2026-06-29:** The wording fix above was a workaround for a file-numbering mismatch that predated it (run order was swapped in commit `84299a4` but the files were never renamed to match). That mismatch is now resolved at the source: `R/13_manifest.R` → `R/14_manifest.R` and `R/14_panels.R` → `R/13_panels.R` (see `RENAME_R13_R14_HANDOFF.md`), so file numbering now matches execution order directly and the manuscript's own R/00–R/14 module-count language is no longer self-contradictory. The quoted reviewer comment and proposed action above describe the pre-rename numbering and are left as-is as the historical record of the finding; no further wording change is needed.

### Finding 4 — Abstract Results overlong + technical bolt-ons
**Section / paragraph:** Abstract → Results
**Reviewer-style comment:** The Results sub-section is ~165 words and ends with two technical bolt-ons ("`manifest.json` (contract v1.2) artefact catalogue tagging every figure and table by kind for downstream tooling") that belong in §2 or §3, not the Abstract.
**Proposed action:** Rewrite Results to ~95 words. Drop contract version, drop "tagging every figure and table by kind for downstream tooling", soften "publication-ready composite figures (PNG and TIFF) and a supplementary tables workbook" to "publication-ready composite figures and a supplementary tables workbook".
**Status:** Done — 2026-06-24. Cuts applied: dropped "alpha-stage" + "(PNG and TIFF)" + "(contract v1.2)" + downstream-tooling clause + "companion". Expansion retained at user request. Abstract now 183 words (was 205); still ~33 over cap — Motivation/Availability remain candidates for a second pass.

### Finding 5 — Intro opening sentence overloaded
**Section / paragraph:** §1 ¶1, sentence 1
**Reviewer-style comment:** "By simultaneously profiling microbial communities and their functional genetic elements across environmental, clinical, and veterinary or livestock matrices, long-read shotgun metagenomics (Oxford Nanopore, PacBio) serves as a powerful One Health tool for tracking health threats that intersect humans, animals, and ecosystems [1]." — 53 words, two subordinate clauses before the main verb.
**Proposed action:** Split into two sentences.
**Draft replacement text:** "Long-read shotgun metagenomics (Oxford Nanopore, PacBio) is increasingly central to One Health surveillance, where threats intersect humans, animals and ecosystems [1]. By profiling microbial communities and their functional genetic elements in a single experiment, it ties taxonomic identity to resistance, virulence and mobility determinants across environmental, clinical and livestock matrices."
**Status:** Done — 2026-06-24 (§1 ¶1, sentence 1 split into two).

### Finding 6 — "gene screening" trailing modifier
**Section / paragraph:** §1 ¶1
**Reviewer-style comment:** "ABRicate [4] against CARD [5], VFDB [6] and PlasmidFinder [7] for antimicrobial resistance genes (ARGs), virulence factors, and mobile genetic elements (MGEs) gene screening" — "gene screening" reads as a noun-cluster afterthought.
**Proposed action:** Reword.
**Draft replacement text:** "ABRicate [4] for screening antimicrobial resistance genes (ARGs), virulence factors and mobile genetic elements (MGEs) against CARD [5], VFDB [6] and PlasmidFinder [7]"
**Status:** Done — 2026-06-24 (§1 ¶1).

### Finding 7 — §2 stage (3) is one ~140-word sentence
**Section / paragraph:** §2 stage (3)
**Reviewer-style comment:** Five distinct analyses (Venn/heatmap, stacked-bar, alpha, beta, ALDEx2) are chained with semicolons into one paragraph-as-sentence. Readers skim Applications Notes; long sentences hide content.
**Proposed action:** Split into 3 or 4 sentences. Keep methodological specifics (Hellinger BC + 9999 perms + PERMDISP) intact.
**Status:** Done — 2026-06-24. Split into one lead sentence + 4 module-tagged sentences (R/04–R/05, R/06, R/07, R/08). All methodological specifics retained.

### Finding 8 — §2 stage (5) repeats Fig 1 caption material
**Section / paragraph:** §2 stage (5)
**Reviewer-style comment:** Stage (5) re-describes R/12 (network + chord + Sankey + Gephi + igraph topology) at the same level of detail as the Fig 1 caption, then re-describes R/14 and R/13. With the figure as the authoritative reference, the prose can be much terser.
**Proposed action:** Compress to ~2 sentences naming the three deliverables (network outputs + composite figures/tables + manifest) and pointing to Fig 1 for layout.
**Status:** Done — 2026-06-24. Compressed to lead + 3 sentences (~85 words; was ~135). Dropped duplicated architecture; kept Gephi/igraph/Bray-Curtis/slot-YAML technical detail.

### Finding 9 — Table 1 caption overloaded
**Section / paragraph:** §3 — Table 1 caption (~135 words)
**Reviewer-style comment:** The current caption mixes (a) what the table shows, (b) cohort descriptions, (c) machine and software settings, (d) `peakRAM` methodology footnote, and (e) implicitly references Figure 2 in the surrounding paragraph. Table captions should be ≤60 words and cover only (a) + the legend.
**Proposed action:** Move cohort descriptions and methodology footnotes into a short "Notes" line below the table or into the post-table paragraph; trim the caption to "Per-stage wall time, peak R-process memory and primary output for the two validation cohorts on the same workstation; runtimes from a 2026-06-23 sequential run."
**Status:** Done — 2026-06-24. Caption trimmed to ~55 words (was ~135). Cohort descriptors retained in §3 ¶1; runtime/peak-MB methodology kept only in §3 ¶1 line 62 (resolves Finding 10 — no edit needed there).

### Finding 10 — Duplicated peakRAM methodology note
**Section / paragraph:** §3 — Table 1 caption + paragraph following the table
**Reviewer-style comment:** "equivalent to the `peakRAM` methodology" appears in both the caption and the prose immediately following the table.
**Proposed action:** Keep in caption, drop from prose (after caption trim per finding 9, may need to land in a single short Methods note).
**Status:** Done — 2026-06-24 (resolved by Finding 9 caption trim; `peakRAM` now only at §3 ¶1 line 62 prose).

### Finding 11 — Jaccard / sample-set agreement compound sentence
**Section / paragraph:** §3 ¶3, sentence beginning "Jaccard agreement on the (sample, GENE)…"
**Reviewer-style comment:** The sentence reports a resistome Jaccard, alpha sample-set match, network node/edge Jaccards, and explicitly defers three more stages. The defer-clause weakens the otherwise concrete numbers around it.
**Proposed action:** Split: one sentence reports the three agreement numbers (resistome 0.82; alpha 18/18; network nodes 0.85 / edges 0.87), one sentence acknowledges the deferred comparisons.
**Status:** Done — 2026-06-24. Split at the semicolon; `gephi` → `Gephi`; deferred-list parens → em-dashes; "agreement" → "comparison" in the deferred clause.

### Finding 12 — Conclusion length and tone
**Section / paragraph:** §4
**Reviewer-style comment:** 175 words for an Applications Note Conclusion is generous; phrasing such as "to encourage convergence on a shared, scriptable layout for long-read metagenomics studies" reads as advocacy rather than method-paper register.
**Proposed action:** Trim to ~80 words. Retain (a) what MERIDIAN consolidates, (b) alpha state + active development, (c) invitation to feedback. Drop the convergence clause.
**Status:** Done — 2026-06-24. Conclusion paragraph now 56 words (was ~92). Correction: original 175-word count was an artefact of including Acks/Funding/COI under the §4 heading; Conclusion proper was ~92 words to start.

### Finding 13 — Kraken / Kraken2 combined reference
**Section / paragraph:** References — entry 2
**Reviewer-style comment:** Entry 2 fuses Wood & Salzberg 2014 (Kraken) and Wood, Lu & Langmead 2019 (Kraken2) into a single numbered reference. Only Kraken2 is used in the pipeline; the 2014 paper is not directly relevant.
**Proposed action:** Replace entry 2 with the Kraken2 paper alone. Renumber if needed (no other text-citation impact).
**Draft replacement text:** "2. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol* 2019;20:257. doi:10.1186/s13059-019-1891-0."
**Status:** Done — 2026-06-24 (reference list entry 2).

### Finding 14 — Author-list inconsistency in references
**Section / paragraph:** References — entries 5, 9 (and check others)
**Reviewer-style comment:** Ref 5 uses "Alcock BP, et al."; ref 9 lists 11 authors then "et al."; refs 1, 3, 6, 7, 10, 11 list all authors. OUP *Bioinformatics* standard is up to 6 authors, then "et al."
**Proposed action:** Apply the 6-author rule uniformly across the reference list.
**Status:** Done — 2026-06-24. Refs 9 + 11 truncated to first 6 + et al. Refs 5 + 7 expanded to first 6 + et al. (CARD authors from academic.oup.com; PlasmidFinder authors from PubMed PMID 24777092 — DOI [10.1128/AAC.02412-14](https://doi.org/10.1128/AAC.02412-14)). Ref 13 (Quarto) left as "Allaire JJ, et al." — no canonical citation exists on quarto.org.

### Finding 15 — Internal project slugs in body text
**Section / paragraph:** §3 throughout
**Reviewer-style comment:** `chicken_batch1` (PC_JC_2024-11-29) and `chicken_batch2` (JC_FL_2025-05-20) are file-system slugs. In body text, "the reference cohort" and "the scaling cohort" read more cleanly; the slugs can stay in Table 1 / Figure 2 captions for reproducibility.
**Proposed action:** Replace body mentions; keep slugs in captions and the methods note.
**Status:** Done — 2026-06-24 (§3 ¶1). Cohort role now subject; slug parenthetical. Net cut ~10 words.

### Finding 16 — "Figure 2 visualises the same numbers"
**Section / paragraph:** §3 — closing sentence before Figure 2
**Proposed action:** "Figure 2 plots the same data: …" or drop the lead-in and let the caption do the work.
**Status:** Done — 2026-06-24 (§3 ¶2). Lead-in tightened; trailing meta-observation dropped (~17 words cut).

### Finding 17 — Inconsistent term: long-read shotgun vs long-read metagenomics
**Section / paragraph:** §1 / §2 / §4
**Proposed action:** Use "long-read shotgun metagenomics" on first mention in each section; "the workflow" or "MERIDIAN" thereafter. Audit and unify.
**Status:** Closed — 2026-06-24, no edit needed. Audit found canonical "long-read shotgun metagenomics" in title, abstract, §1, §4. §3 ¶1 line 62 "long-read chicken caecum metagenomics studies" omits "shotgun" but is acceptable in context. Ref 11 title kept verbatim.

### Finding 18 — Repetition of "alpha-stage"
**Section / paragraph:** Abstract, §1 ¶3, §2 ¶1 (implied), §4 ¶1
**Proposed action:** Keep "alpha-stage" in §1 ¶3 (development-state framing) and §4 ¶1 (release intent). Drop from Abstract Results ("an alpha-stage modular R workflow" → "a modular R workflow"); not in §2 by name.
**Status:** Done — 2026-06-24. After Findings 4 + 12 the term lived in §1 as "early-stage" and §4 as "alpha-stage" — unified on "alpha-stage" at §1 line 36 (matches pipeline_facts.md "Phase 1 / alpha").

## Workflow

Per skill default: walk these one finding at a time, propose draft text, wait for sign-off, then commit to `applications_note.md`. After each commit, re-read the touched paragraph for orphan punctuation / connectors. Log status as `Done — <date>` against each finding.

Order I propose: 2 → 13 → 3 → 4 → 6 → 5 → 7 → 8 → 9 → 10 → 11 → 12 → 14 → 15 → 16 → 17 → 18.
Rationale: cheap wins and citation fixes first (2, 13), structural narrative fixes next (3), abstract (4), then incremental prose tightening that cumulatively delivers the word-count cut required by finding 1. Finding 1 itself is closed by the cumulative effect of 4 + 7 + 8 + 9 + 10 + 11 + 12.

## Post-walk addendum — Option 2 (Table 1 / Figure 2 overlap)

After completing the 18-finding walk, user asked whether Table 1 should move to supplementary given its overlap with Figure 2. Decision: keep Table 1 in main text but strip the duplicated columns (runtime + peak MB) — Figure 2 already plots those. Table 1 retained as "Stage | B1 primary output | B2 primary output", preserving the PERMANOVA / ALDEx2 / resistome / network validation evidence in the main text. Post-table analytical paragraph collapsed from ~300 words (performance accounting) to ~110 words (Figure 2 reading + dashboard pointer). §3 ¶1 closing sentence redirected to Figure 2; Figure 2 caption "Numbers reproduce Table 1" clause dropped.

Net effect (Option 2):
- §3 went from 621 → 462 words (−159)
- Body excl refs/tables/captions: 1,983 → 1,824 words

## Final word-count state (after all 18 findings + Option 2)

| Section | Original | Final | Δ | Target |
|---|---|---|---|---|
| Abstract | 205 | 183 | −22 | ≤150 (still +33) |
| §1 Introduction | 273 | 280 | +7 | ~240 |
| §2 Implementation | 647 | 637 | −10 | ~500 |
| §3 Application | 704 | 462 | −242 | ~450 (✓ hit) |
| §4 Conclusion (incl. Acks/Funding/COI) | 175 | 150 | −25 | — (Conclusion paragraph proper: 56 words ✓) |
| **Body total** | **2,116** | **1,824** | **−292 (−14 %)** | ≤1,500 (still +324) |

Closing notes: §3 now within budget; §2 and Abstract still over. Remaining trim candidates if a second pass is wanted: §2 ¶1 (introductory architecture paragraph), §2 stage descriptions, Abstract Motivation.
