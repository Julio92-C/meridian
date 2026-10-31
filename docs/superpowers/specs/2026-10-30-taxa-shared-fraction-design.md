# Taxa shared fraction, metadata linkage and taxonomic tree — design

**Status:** approved 2026-10-30 (defaults: ≥3 samples and ≥10% in every group; no change to composite panels) · **Branch:** `feat/taxa-shared-fraction` (off `fix/hospital-species-aldex2`, PR #12) · **Date:** 2026-10-30

## Why

The taxonomy section produces less than the resistome, virulome and mobilome sections: a Venn diagram (2–5 groups only) and a top-N heatmap. The hospital study showed that the overlap between groups needs more than a count to interpret:

- 590 of 2,255 species were shared between air and surfaces.
- But 496 of the 590 occurred in only one sample of at least one compartment.
- Only 13 species were widespread in both; most were hospital-associated Gram-negative bacteria.

The prototype is Python, in hosMicro `manuscript/scripts/shared_fraction/`. This stage brings that analysis into meridian for any number of groups, and adds two pieces:

- links between the levels of any metadata variable (ward, room, clinical status);
- a taxonomic tree.

## Scope

**New files**
- `R/04_shared_fraction.R` (stage `shared_fraction`, run after `taxonomy`; file names must match `^\d+_` to be sourced)
- `R/utils_shared_fraction.R` (pure functions, unit-tested)
- `R/04_taxa_tree.R` (stage `taxa_tree`)

**Changed files**
- `R/14_manifest.R`: new `kind`s.
- `run_pipeline.R`
- `config/config_template.yaml`

**Out of scope for v1**
- A molecular phylogeny of reads. See "Taxonomic tree".
- Automatic filtering of implausible gene–taxon links. These are flagged, not removed.

## Inputs

Inputs are the same as the taxonomy stage, so the numbers agree with the Venn diagram from the same run:

- `community_filter()` → `clean_taxa_names()` applied to `cleaned$noncontaminants`, at species rank, with the >N-read floor and the shallow-sample drop;
- the group from `cfg$metadata$group_cols[[1]]`;
- the element hits (`abri_kraken2_cleaned`), for genes carried per species;
- the linkage variables from the metadata (see config).

## Definitions (any number of groups G ≥ 2)

**Per species, per group g:**
- **n_g:** the number of samples where the species is detected.
- **prev_g:** n_g as a fraction of the samples in g.
- **abund_g:** the mean relative abundance over all samples in g.

**Membership pattern:** the set of groups with n_g > 0.

**Categories:**
- `unique:<g>` — the species is present in one group only;
- `shared:<g1+g2…>` — present in some groups but not all (only when G > 2);
- `shared:all` — present in every group.

**Tiers:** these apply to any shared species and are evaluated over the groups in its pattern.
- **widespread:** in every group of the pattern, `n_g ≥ widespread_min_samples` **and** `prev_g ≥ widespread_min_prev`;
- **recurrent:** not widespread, but `n_g ≥ 2` in every group of the pattern;
- **sporadic:** everything else, i.e. one sample in at least one group.

**Pairwise view:** for every pair of groups (i, j), the species present in both, tiered on those two groups only. With G = 2 there is one pair, and it equals `shared:all`.

## Outputs

**Tables** (in `Datasets/shared_fraction/`):

| File | `kind` |
|---|---|
| `shared_fraction_species.csv`: one row per species, with n/prev/abund per group, pattern, category and tier | `shared_fraction_species` |
| `shared_fraction_pairwise.csv`: one row per pair and species | `shared_fraction_pairwise` |
| `shared_fraction_widespread_elements.csv`: widespread species × gene, with the number of reads per link and a `low_support` flag (reads < `elements_min_reads`) | `shared_fraction_elements` |
| `shared_fraction_linkage_edges.csv`: one row per variable, node pair and species count | `shared_fraction_linkage_edges` |
| `taxa_tree_nodes.csv` | `taxa_tree_nodes` |

**Figures** (in `Figures/shared_fraction/`, PNG):

| Figure | `kind` | Content |
|---|---|---|
| Prevalence–abundance | `shared_fraction_prevalence_abundance` | All species; colour = category; widespread species circled |
| Pairwise prevalence | `shared_fraction_pairwise` | prev_i against prev_j for every pair (faceted when G > 2); size = abundance; colour = tier |
| Membership | `shared_fraction_upset` | UpSet plot of membership patterns with bars stacked by tier. Replaces the 5-group Venn limit |
| Tiers × phylum | `shared_fraction_tiers` | Stacked bars |
| Widespread species | `shared_fraction_widespread` | Prevalence per group, with ARG/VF/MGE counts |
| Heatmap | `shared_fraction_heatmap` | Top-N shared species, ranked by minimum prevalence across groups; samples ordered by group, then by the first linkage variable |
| Linkage | `shared_fraction_linkage` | One per linkage variable (see below) |
| Taxonomic tree | `taxa_tree` | See below |

**Composite (`R/13`):** not changed for now (decision 2026-10-30). The figures are standalone and registered in the manifest; a composite can be added later.

## Linkage by metadata variable

`shared_fraction.linkage_vars` takes a list of metadata columns, such as `[ward]`, `[room]`, `[patient_status, antibiotic_exposure]`.

**Nodes** are (group, level) pairs, e.g. `Air · ward B` or `Surface · SR 16`.
- Categorical variables are used as they are.
- Numeric variables are binned into quantiles (`linkage_bins`, default 3), e.g. low, mid and high CRP.
- Samples with a missing value go to an `NA` node, which is dropped unless `linkage_keep_na: true`.

**Edges** join two nodes that share at least `linkage_min_shared` widespread species (default 1); the edge weight is that count.
- By default only nodes in different groups are joined (`linkage_scope: between`); `within` and `all` are also allowed.

**Layout** is a k-partite layout: one column per group, with nodes ordered by level.
- This is schematic, not spatial. A true floor-plan map would need coordinates, which a later `linkage_coords` CSV could supply.
- The figure is paired with a presence matrix (species × node).

**Statistics** (optional, `linkage_test: true`) use a permutation test: sample labels are shuffled within each group, to ask whether edge weights exceed chance. The seed comes from `cluster_seed(cfg)`.

**Hospital example:** `linkage_vars: [ward]` with a derived `room` column for the surfaces reproduces the note's Figure 2.

## Taxonomic tree

There is no molecular phylogeny for read-level profiles. The tree is built from the NCBI lineage the pipeline already carries (phylum → class → order → family → genus → species), as GraPhlAn and MetaPhlAn cladograms do.

**Built with:** `igraph` and `ggraph` (circular dendrogram), both already in `renv.lock`, so no new dependencies. The edge list is keyed by full lineage path, so repeated labels such as "unclassified class" under different parents stay separate.

**Tips:** the top-N species (`taxa_tree.top_n`, default 60), ranked by mean abundance among species detected in at least `taxa_tree.min_samples` samples (default 2), so a species seen once in a low-depth sample cannot enter on its relative abundance alone. The widespread shared species are always included.

**Encoding:**
- clade background = phylum;
- tip colour = membership category;
- tip size = mean abundance;
- outer ring of tiles = prevalence in each group;
- a marker on tips that carry an ARG, VF or MGE.

**Caption wording:** "Taxonomic tree based on NCBI lineage; branch lengths are not evolutionary distances".

**Later option:** a genome phylogeny of MAGs from the GTDB-Tk tree. That would need `ape`/`ggtree` added to renv and is outside v1.

## Config (template defaults)

```yaml
shared_fraction:
  enabled: true                 # auto-skipped when < 2 groups
  widespread_min_samples: 3
  widespread_min_prev: 0.10     # both conditions must hold in every group of the pattern
  heatmap_top_n: 30
  elements_min_reads: 2
  linkage_vars: []              # e.g. [ward]
  linkage_bins: 3
  linkage_min_shared: 1
  linkage_scope: between
  linkage_keep_na: false
  linkage_test: false
taxa_tree:
  enabled: true
  top_n: 60
  min_samples: 2              # tips detected in at least this many samples
```

**Hospital config** (to reproduce the note):
- `widespread_min_samples: 5`
- `widespread_min_prev: 0`
- `linkage_vars: [ward]`, plus the surface room column

## Validation

1. **Unit tests** (`tests/testthat/test-shared-fraction.R`), using a synthetic 3-group table with known patterns, check:
   - categories, tiers and pairwise counts;
   - binning of a numeric linkage variable;
   - that G = 1 skips cleanly;
   - that G = 6 works where the Venn would skip.
2. **Hospital regression:**
   - The air-only / shared / surface-only counts equal the taxonomy Venn of the same run (690/586/966 with pipeline name-cleaning). The 4–6 species differences from the paper's 696/590/969 are the known placeholder-name drops.
   - The widespread set equals the note's 13 species, minus any species lost to name-cleaning, which is to be listed.
3. **Panels:** the composite renders and the manifest validates.

## Manifest contract

The new `kind`s need rows in metaomics-scribe `docs/MANIFEST_SCHEMA.md` §`kind` vocabulary. That is a **MINOR** bump in the manifest version. Existing kinds are unchanged, so older agents still read the manifest.

## Decisions (2026-10-30)

1. Template default for widespread: at least 3 samples **and** at least 10% prevalence in every group of the pattern.
2. No change to the composite figure panels (`R/13`) for now.
3. The existing taxonomy heatmap annotation (Core / Accessory / Unique) is left unchanged.
