# wetlands_flyway — independent validation cohort

**One Health axis:** wildlife + natural/anthropogenic environment (East Atlantic Flyway)

## Source
Perlas A, Reska T, Sánchez-Cano A, et al. *Real-time genomic pathogen,
resistance, and host range characterization from passive water sampling of
wetland ecosystems.* Appl Environ Microbiol 2026;92(5):e0254325.
doi:10.1128/aem.02543-25 (PMID 42029155).

- **Platform:** Oxford Nanopore MinION (R10.4.1), long-read shotgun metagenomics
- **Samples:** 24 DNA metagenomes — 12 sites × 2 biological replicates
- **Data:** ENA BioProject **PRJEB96272**
- **Study code:** https://github.com/ttmgr/Wetland_Health

## Metadata decode (`Metadata/wetlands_flyway_metadata.csv`)
Sample-title pattern `SQK-RBK114-24_barcodeNN.<SITE>.<REP>` where `<SITE>` =
`[country][land-use][n]`:
- **country:** G = Germany, F = France, S = Spain
- **land_use:** A = anthropogenic (impacted), N = natural
- This split yields exactly **6 anthropogenic** (GA1-3, FA1-2, SA1) and **6 natural**
  (GN1-3, FN1, SN1-2) sites — matching the paper's "6 impacted vs 6 natural".
  `land_use` is the primary contrast; `country` is a 3-level covariate; `site` is the
  natural random-effect candidate (2 reps per site).

> ⚠️ **Inferred decode.** The country/land-use letters are inferred from the sample
> titles + the paper's design, not from an explicit deposit field. Sanity-check
> against the paper's site table before final reporting.

## Download the reads (ENA)
The `fastq_ftp` field is empty for this project — reads were submitted as
`*.fastq.tar.gz` archives. Pull them via the ENA `submitted_ftp` field or SRA:
```bash
# Get the submitted-file URLs:
# https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB96272&result=read_run&fields=run_accession,submitted_ftp&format=tsv
# Only the 24 RBK114 runs are DNA shotgun (see run_accession list in the metadata CSV).
# Exclude the RPB114 runs (RNA virome / cDNA).
```
Untar each archive and rename the FASTQ to its `sample` id (e.g. `GA1_1`) so
upstream outputs carry matching sample names.

## Then (you run these upstream, per MERIDIAN's contract)
1. Kraken2 + Bracken → `Reports/Kraken2/`
2. ABRicate (CARD / VFDB / PlasmidFinder) → `Reports/Abricate/summary_reportClean.csv`
3. Re-centrifuge → `Reports/Re-centrifuge/*.rcf.data.csv`
4. Set `project_root` in `config.yaml`, then run:
   `Rscript run_pipeline.R projects/wetlands_flyway/config.yaml`

## Caveats
- **RPB114 runs excluded.** The 12 `SQK-RPB114-24_*` runs are the RNA-virome/cDNA
  subset, not DNA shotgun — do not add them to the metadata.
- Environmental water metagenomes are host-poor and diverse; expect a different
  taxonomy/resistome shape from gut cohorts. That contrast is the point of the test.
- No hand-edited ground truth — validates **portability/robustness**, not
  output-for-output correctness.
