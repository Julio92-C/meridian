# Related papers in scope — shortlist for citation / framing

Compiled 2026-05-29 via PubMed (legal notice: information from PubMed; all DOIs included as links).

## A. Direct templates (already in `references.md`)

### 1. L-ARRAP — closest scope match (long-read AMR pipeline)

- **Title:** Quantifying antibiotic resistome risks across environmental niches: the L-ARRAP for long-read metagenomic profiling.
- **Authors:** Li Y, Gao Y, Liu X, Mao Y, Wang M, Qin Y, Zhang C, Chen Q, Ning K, Wang Z, Han M.
- **Journal:** *Briefings in Bioinformatics* 2025;26(5):bbaf535.
- **DOI:** [10.1093/bib/bbaf535](https://doi.org/10.1093/bib/bbaf535)
- **PMID:** 41066697, **PMC:** PMC12510455
- **Full text saved:** `templates/L-ARRAP_Li_2025_BriefBioinform.md`
- **Why useful:** same input modality (long-read), same scope (ARGs + MGEs + pathogens), narrative arc: "fragmented existing tools → integrated long-read pipeline → multi-niche validation → availability". This is the closest peer paper.

### 2. rANOMALY — R-workflow + auto-report angle

- **Title:** rANOMALY: AmplicoN wOrkflow for Microbial community AnaLYsis.
- **Authors:** Theil S, Rifa E.
- **Journal:** *F1000Research* 2021;10:7.
- **DOI:** [10.12688/f1000research.27268.1](https://doi.org/10.12688/f1000research.27268.1)
- **PMID:** 33537122, **PMC:** PMC7836088
- **Full text saved:** `templates/rANOMALY_Theil_Rifa_2021_F1000Research.md`
- **Why useful:** same form factor (R package + Rmarkdown auto-report); structural template for "config-driven, reproducible, auto-reported" framing. Different domain (amplicon vs shotgun) — do not copy phrasing, use as structural template only.

### 3. microeco — R protocol reference

- **Title:** A workflow for statistical analysis and visualization of microbiome omics data using the R microeco package.
- **Authors:** Liu C, Mansoldo FRP, Li H, Vermelho AB, Zeng RJ, Li X, Yao M.
- **Journal:** *Nature Protocols* 2025;21(4):1300–1324.
- **DOI:** [10.1038/s41596-025-01239-4](https://doi.org/10.1038/s41596-025-01239-4)
- **PMID:** 40770112
- **Full text:** not in PMC OA (Nature Protocols subscription) — abstract metadata retrieved only.
- **Why useful:** comprehensive R protocol covering preprocessing, alpha/beta diversity, differential abundance, ML on amplicon + metagenomic + metabolomic data. Reference for how to present each analysis stage rigorously, not for direct citation.

## B. Additional candidates surfaced 2026-05-29

### 4. BALROG-MON — long-read AMR Nextflow pipeline

- **Title:** BALROG-MON: a high-throughput pipeline for Bacterial AntimicrobiaL Resistance annOtation of Genomes-Metagenomic Oxford Nanopore.
- **Authors:** Bird E, Pickens V, Molik D, Silver K, Nayduch D.
- **Journal:** *microPublication Biology* 2025.
- **DOI:** [10.17912/micropub.biology.001427](https://doi.org/10.17912/micropub.biology.001427)
- **PMID:** 40046038, **PMC:** PMC11880933
- **Relevance:** Nextflow-based, not R-based; same problem space (long-read AMR + pathogen + ARG annotation) as our tool. Worth citing in the Introduction when listing existing long-read pipelines alongside L-ARRAP — strengthens the "no comparable R workflow exists" claim.

### 5. Reska et al. — long-read air-microbiome pipeline

- **Title:** Air monitoring by nanopore sequencing.
- **Authors:** Reska T, Pozdniakova S, Borràs S, et al.
- **Journal:** *ISME Communications* 2024;4(1):ycae099.
- **DOI:** [10.1093/ismeco/ycae099](https://doi.org/10.1093/ismeco/ycae099)
- **PMID:** 39081363, **PMC:** PMC11287864
- **Relevance:** another long-read metagenomics workflow (different niche). Useful only if the Introduction needs more "long-read shotgun pipelines exist" precedent. Probably skip unless you want to broaden citation breadth.

### 6. Urrutia-Angulo et al. — short vs long-read resistome comparison

- **Title:** Resistome and microbiome profiling of bovine milk following antimicrobial dry cow therapy: insights from short- and long-read metagenomic sequencing.
- **Authors:** Urrutia-Angulo L, Lavín JL, Oporto B, Aduriz G, Hurtado A, Ocejo M.
- **Journal:** *Frontiers in Microbiomes* 2025;4:1672438.
- **DOI:** [10.3389/frmbi.2025.1672438](https://doi.org/10.3389/frmbi.2025.1672438)
- **PMID:** 41852383, **PMC:** PMC12993662
- **Relevance:** explicitly highlights that "sequencing platform and bioinformatic pipeline choices substantially influence resistome profiling outcomes". A useful supporting citation for the Introduction's motivation paragraph ("no standardised analytical approach in resistome studies").

### 7. Fang et al. — mNGS review

- **Title:** Emerging role of metagenomic next-generation sequencing in infectious disease diagnostics: Clinical integration and future directions.
- **Authors:** Fang T, Yuan F, Chen Y, et al.
- **Journal:** *mLife* 2026;5(2):148-163.
- **DOI:** [10.1002/mlf2.70078](https://doi.org/10.1002/mlf2.70078)
- **PMID:** 42079444, **PMC:** PMC13131327
- **Relevance:** broad review of mNGS in infectious disease diagnostics, including long-read platforms. Useful citation if a clinical-application sentence is added to the Introduction.

### 8. Amoutzias et al. — bacterial pathogen genomics review

- **Title:** The Notable Achievements and the Prospects of Bacterial Pathogen Genomics.
- **Authors:** Amoutzias GD, Nikolaidis M, Hesketh A.
- **Journal:** *Microorganisms* 2022;10(5):1040.
- **DOI:** [10.3390/microorganisms10051040](https://doi.org/10.3390/microorganisms10051040)
- **PMID:** 35630482, **PMC:** PMC9148168
- **Relevance:** review of pathogen genomics including long-read assemblies. Probably skip — too broad for a 2-page applications note's reference budget.

## C. Citation hygiene findings (now applied)

| Manuscript ref | Status after PubMed verification | Action taken |
|---|---|---|
| 1 — Kraken/Kraken2 (Wood 2014, 2019) | Correct (PMID 24580807 / 31779668) | No change |
| 2 — Bracken (Lu 2017) | Correct | No change |
| 3 — ABRicate (Seemann) | Software-only, no PubMed | TODO note clarified |
| 4 — CARD (Alcock 2023) | Correct | No change |
| 5 — VFDB (Chen 2005) | Stale — superseded by VFDB 2022 update | **Replaced** with Liu B et al. 2022 *Nucleic Acids Res* 50:D912–D917 ([10.1093/nar/gkab1107](https://doi.org/10.1093/nar/gkab1107), PMID 34850947) |
| 6 — PlasmidFinder (Carattoli 2014) | Correct | No change |
| 7 — Re-centrifuge (Martí 2019) | Correct | No change |
| 8 — vegan | R package, no PubMed | TODO note clarified |
| 9 — ALDEx2 (Fernandes 2014) | Correct (PMID 24910773) | **Author list expanded** to full six authors; TODO cleared |
| 10 — L-ARRAP (Li 2025) | Title/authors incomplete | **Replaced** with verified PubMed metadata: Li Y, Gao Y, Liu X et al. + full title |
| 11 — rANOMALY (Theil & Rifa 2021) | Correct | No change |
| 12 — Quarto | Software-only, no PubMed | TODO note clarified |
| 13 — R Core Team | Software-only, no PubMed | TODO note clarified |

Information from PubMed (per the PubMed MCP attribution requirement).
