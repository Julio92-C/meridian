# Drafting brief — Bioinformatics Applications Note

## Goal

Produce a complete first-draft manuscript at `applications_note.md`
formatted as a **Bioinformatics (OUP) Applications Note**.

## Length and structure

- ≤2 published pages (~1,200–1,500 words excluding refs)
- Sections, in order:

1. **Title** (≤120 chars, no jargon)
2. **Authors / Affiliations** — `{{TODO: author block}}` placeholder
3. **Abstract** — exactly 4 labelled mini-sections:
   *Motivation*, *Results*, *Availability and Implementation*, *Contact*.
   Total ≤150 words.
4. **1 Introduction** — ~250 words; problem framing + gap in existing tools.
5. **2 Implementation** — config-driven architecture, modular stages
   (R/00–R/12), tool stack, Quarto report. Include a figure callout:
   `{{Figure 1: pipeline schematic}}`.
6. **3 Application / Validation** — describe the chicken_batch1 real-data
   benchmark and how pipeline output was compared against hand-edited
   ground-truth scripts. **Do not fabricate numerical results.** Use
   `{{TODO: insert benchmark table}}` where specific numbers belong.
7. **4 Conclusion** — one short paragraph.
8. **Acknowledgements / Funding / Conflict of Interest** — placeholders.
9. **References** — numbered style; cite the actual upstream tools
   (Kraken2, Bracken, ABRicate, Re-centrifuge, ALDEx2, vegan, Quarto, R)
   with real DOIs/URLs where known, otherwise `{{TODO: verify ref}}`.

## Template to mirror (structurally, not verbatim)

**L-ARRAP** — Li et al. 2025, *Briefings in Bioinformatics*,
DOI 10.1093/bib/bbaf535. Fetch via WebFetch if permitted, otherwise
follow the structural arc from memory: "fragmented existing tools →
unified pipeline → validation on real data → availability."

**Secondary**: **rANOMALY** — Theil & Rifa 2021, *F1000Research*,
DOI 10.12688/f1000research.27268.1, for the R-workflow + auto-report
framing.

## Hard constraints

- **No fabricated results.** Every number/benchmark you would normally
  cite becomes a `{{TODO}}` placeholder. End the draft with a
  `## Drafting notes` section listing every placeholder and what data
  is needed to fill it.
- **Workflow framing, not R-package framing.** This is published as
  a workflow on GitHub. Availability says "git clone + Rscript
  run_pipeline.R", not `install.packages()`.
- **Tone**: concise methods-paper register. No marketing language.
- **Honest about state**: the pipeline is alpha / Phase-1. Frame as
  "we present" but do not overclaim production-readiness.

## Deliverables

1. `applications_note.md` — the draft.
2. A short summary back to the user when complete: word count,
   list of `{{TODO}}` placeholders, structural decisions taken, and
   any factual claims you were unsure about.
