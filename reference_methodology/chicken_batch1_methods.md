# Chicken batch 1 — reference methodology

Extracted from the Frontiers submission draft
`Paper 1Impact_of_Dulse_on_Poultry_Gut_Microbiome_reviewed_PL_JC_Frontiers_Template_PL.docx`
(Chicken batch 1 project directory). These are the exact methods and tool
versions the automated pipeline must reproduce.

## Sequencing & pre-processing

- Platform: Oxford Nanopore long reads.
- Host-read removal: **Samtools v1.17** against
  `Gallus gallus` reference (NCBI Annotation Release 106, `GCF_016699485.2`).

## Taxonomic classification

- **Kraken2 v2.1.2** — Wood *et al.* 2019.
- Database: PlusPF-8 (accessed 26 Apr 2025); RefSeq archaea, bacteria,
  protozoa, fungi, viruses + human genome + UniVec_Core.
- Abundance estimation: **Bracken v2.7.0**.

## ARGs / VFs / MGEs

- **ABRicate v1.0.1** on Kraken2-classified reads.
- Databases:
  - CARD → ARGs (resistome).
  - VFDB (Chen *et al.* 2016) → virulence factors (virulome).
  - PlasmidFinder → MGEs (mobilome).

## Decontamination

- **Re-centrifuge**: removes negative controls and crossover taxa.
- Contaminant read counts are retrieved from the metagenomics dataset; the
  Kraken2 report is then joined to the ABRicate summary via a custom R step
  ("LR-GEs contigs").

## Normalisation

- Raw read counts normalised as **TPM** (transcripts per million) per gene
  element.

## Alpha diversity

- Metrics: **richness** (count of unique taxa / GEs), **Shannon index**.
- Statistical test: **Kruskal–Wallis** across treatment groups.
- Block number used as a **random factor** in the models (per metadata
  design — randomised block design, chicken batch 1/2).

## Beta diversity

- Dissimilarity: **Bray–Curtis**.
- Ordination: **PCoA**.
- Test: **PERMANOVA** via `adonis2` (vegan), **9999 permutations**, α = 0.05.

## Differential abundance

- **ALDEx2** with **128 Monte-Carlo Dirichlet instances**.
- Low-prevalence filter: species present in **< 2 samples** removed before DA.
- Pairwise CLR comparisons across dietary treatments.

## Network analysis

- Tripartite network: **Treatment × taxa × genetic elements** (ARG/VF/MGE).
- Force-directed layout in Gephi (iterated to minimal energy state).
- Topology metrics: degree centrality, betweenness, modularity.
- Node colour/size encode node type and experimental group affiliation.

## Additional plots produced in the paper

- Venn diagrams of shared ARGs / VFs / taxa across treatment groups.
- pHeatmap of ARG / VF / MGE abundance.
- Chord diagrams for taxa–GE relationships.
- Sankey diagrams for virulence factor profiles.
- Violin plots of alpha diversity per treatment.
- Stacked bar plots of relative abundance at multiple taxonomic levels.

## Metadata design (Chicken batches)

- Randomised block design.
- Random effect: **Block**.
- Fixed effects: **Treatment** and **Weight** (and possibly their interaction).
- Treatment levels (batch 2): Control (0% seaweed) → Seaweed (6% inclusion).
- Treatment levels (batch 1): Control, Reference diet, Soyabean meal, Dulse.

## Lung microbiome variant

Same core stages, but the metadata shape differs: subject-level clinical
variables (`Immunodeficiency diagnosis`, `Treatment` as antibiotic, `Sex`,
`Age`, `Immunocompromised type`) replace the block/treatment design. The
pipeline must accept any set of grouping columns named in the config.
