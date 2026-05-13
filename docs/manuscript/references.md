# References — target journal + template papers

## Target journal

**Bioinformatics** (Oxford University Press), Applications Note format.

Author guidelines: https://academic.oup.com/bioinformatics/pages/instructions_for_authors

Key constraints:
- ≤2 published pages (~1,500 words body, excluding refs/figures)
- Single figure max (typically)
- Numbered references, OUP citation style
- Abstract structured as: Motivation / Results / Availability and
  Implementation / Contact

## Backup journal options (in order of fit)

1. **Bioinformatics Advances** — OUP sister journal, more flexible format
2. **F1000Research** — rolling open review, lower barrier (rANOMALY's venue)
3. **GigaScience** — rewards reproducibility, FAIR data, containers
4. **BMC Bioinformatics** — Software section
5. **Briefings in Bioinformatics** — if scope expands beyond Applications Note (L-ARRAP's venue)

## Template papers

### Primary — closest scope match

**L-ARRAP: Long-read Antibiotic Resistome Risk Assessment Pipeline**
Li et al., 2025. *Briefings in Bioinformatics* 26(5).
DOI: [10.1093/bib/bbaf535](https://doi.org/10.1093/bib/bbaf535)
PMID: 41066697, PMC: PMC12510455

Why it's the template:
- Same input modality: long-read (nanopore / PacBio) metagenomics
- Same outputs: ARGs + MGEs + bacterial pathogens, integrated risk score
- Same framing: pipeline paper, validates against existing tools
- Narrative arc to mirror: "existing tools are short-read-only /
  fragmented → we built an integrated long-read pipeline → here is
  validation on three datasets → availability"

### Secondary — for R-workflow + auto-report angle

**rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis**
Theil & Rifa, 2021. *F1000Research* 10:7.
DOI: [10.12688/f1000research.27268.1](https://doi.org/10.12688/f1000research.27268.1)
PMID: 33537122, PMC: PMC7836088

Why useful:
- Same form factor: R-based reproducible workflow with auto-generated
  Rmarkdown reports
- Useful for the "config-driven, reproducible, auto-reported" framing
- Different domain (amplicon vs shotgun) so do not copy phrasing —
  use only as structural template

### Honourable mention — R protocol reference

**microeco: comprehensive R protocol for microbiome data**
Liu et al., 2025. *Nature Protocols*.
DOI: [10.1038/s41596-025-01239-4](https://doi.org/10.1038/s41596-025-01239-4)
PMID: 40770112

Step-by-step R protocol covering preprocessing, alpha/beta diversity,
differential abundance, machine learning on amplicon + metagenomic +
metabolomic data. Useful as a reference for how to present each analysis
stage rigorously.

## Citation attribution

Above DOIs were retrieved via PubMed. Any factual claims about these
papers should preserve the PubMed attribution.
