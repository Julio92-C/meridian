# hospital_wastewater_eg — independent validation cohort

**One Health axis:** human / clinical + built environment (hospital effluent)

## Source
Radwan HM, El Menofy NG, Tharwat EK, Mysara M, Radwan SMR. *Metagenomic
profiling of microbial communities and the resistome within Egyptian hospital
wastewater and tap water.* Sci Rep 2026;16(1). doi:10.1038/s41598-026-49481-4
(PMID 42062403).

- **Platform:** Oxford Nanopore MinION (R9.4.1), whole-genome shotgun metagenomics
- **Samples:** 20 — 10 hospital wastewater (HWW) + 10 tap water; 5 hospitals × 2 seasons
- **Data:** NCBI BioProject **PRJNA1266479** (20 runs, single-end Nanopore FASTQ)
- **Databases in source study:** CARD, ResFinder, PlasmidFinder (maps to MERIDIAN resistome/mobilome)

## Metadata decode (`Metadata/hospital_wastewater_eg_metadata.csv`)
- SRA alias `NSW` = hospital wastewater; `NP` = tap ("potable") water. Encoded in
  `water_type` (wastewater | tap_water) — the primary PERMANOVA/ALDEx2 contrast.
- `hospital` and `season` are left **NA**: the aliases encode hospital × season but
  the deposit alone doesn't reveal the mapping. Fill from the paper Methods/Table
  before using either as a covariate. **Do not** treat tap water as a control — it
  is a biological comparison group (`controls: []`).

## Download the reads (SRA/ENA)
ENA direct FASTQ (one URL per run — see the `fastq_ftp` column of the ENA
filereport). Example for the first two:
```bash
# All 20 fastq URLs: https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1266479&result=read_run&fields=run_accession,fastq_ftp&format=tsv
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR336/097/SRR33664697/SRR33664697_1.fastq.gz  # HWW1 (1SW)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR336/095/SRR33664695/SRR33664695_1.fastq.gz  # TAP1 (1P)
```
Or via SRA toolkit: `prefetch SRR33664697 && fasterq-dump SRR33664697`.
Rename each FASTQ to its `sample` id (HWW1..HWW10, TAP1..TAP10) so upstream
outputs carry matching sample names.

## Then (you run these upstream, per MERIDIAN's contract)
1. Kraken2 + Bracken → `Reports/Kraken2/`
2. ABRicate (CARD / VFDB / PlasmidFinder) → `Reports/Abricate/summary_reportClean.csv`
3. Re-centrifuge → `Reports/Re-centrifuge/*.rcf.data.csv`
4. Set `project_root` in `config.yaml`, then run:
   `Rscript run_pipeline.R projects/hospital_wastewater_eg/config.yaml`

## Caveats
- Read depth is uneven (1.6k–158k reads/run); very low-depth runs may thin the
  resistome/mobilome tables. Inspect `min_count_per_sample` if stages come back sparse.
- No hand-edited ground truth exists — this validates **portability/robustness**, not
  output-for-output correctness.
