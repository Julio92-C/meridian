# Pipeline v2 — manifest and figure gaps to close

Literature-grounded punch-list for [`Metagenomics_pipeline_automation`](https://github.com/Julio92-C/Metagenomics_pipeline_automation).
Derived from systematic review of 8 published shotgun-metagenomics papers on poultry gut microbiome (2022–2025, primarily Frontiers in Microbiology) using the `metaomics-scribe` literature agent (session 2026-06-10).

> **Related**: see [`docs/STATISTICS.md`](STATISTICS.md) for the family-wide multiple-testing correction system (`cfg$stats$padjust_method`, default BH, optional BY) and a per-surface map of where it applies.

**Comparator papers consulted:**
- Yang et al. 2023, PMC10222538 — broiler resistome + virulome + correlations
- Dierikx et al. 2022, fmicb.2022.833790 — chicken gut resistome, Frontiers Microbiology
- Zhou et al. 2024, fmicb.2024.1509461 — organic vs conventional feeding, Frontiers Microbiology
- Nasser et al. 2024, fmicb.2024.1487595 — village chicken faecal microbiota, Frontiers Microbiology
- Xu et al. 2023, PMC9751072 — resistomes for interconnected humans, soil, livestock (ISME J)
- Liu et al. 2023, PMC10511805 — retail chicken resistome, mobilome, composition (Poultry Science)
- Cai et al. 2022, fmicb.2022.930289 — indigenous Chinese chicken microbiome + ARGs, Frontiers Microbiology
- Laying Hens 2025, PMC12838260 — fecal microbiome + resistome, laying hens

---

## Scope and version bump

This file covers three categories of change:

| Category | What changes | Version impact |
|---|---|---|
| **A. R script fixes** | Fix existing rendered panels (no new kinds, no manifest change) | None — manifest unchanged |
| **B. v1.1 carry-forward** | Kinds documented in `PIPELINE_V1.1_GAPS.md` not yet emitted | `1.0 → 1.1` (MINOR, additive) |
| **C. v1.2 new kinds** | New kinds required by the literature-grounded figure restructure | `1.1 → 1.2` (MINOR, additive) |

All changes are **additive only** — old manifests continue to validate; the agent reflows or skips slots whose new kinds are absent.

Bump `manifest_version` in `R/14_manifest.R` to `"1.2"` once all Category B + C items are emitted.

---

## Section A — R script fixes (no manifest change needed)

These are rendering defects in existing pipeline output. Fix the R plot code; re-run the pipeline; the agent picks up the corrected PNGs automatically.

### A1 — Fig 1: add trend-line overlay to alpha-bar plots

**File to fix:** whichever R script produces `alpha_diversity/richness_bar.png` and `alpha_diversity/shannon_bar.png`

**Fix:** Add a `geom_line(aes(group = 1))` connecting the per-treatment group means (or a smoothing line), matching the style used in comparator papers 2 and 3. The dashed horizontal line (current global mean) is useful but does not replace the trend line that reviewers expect to show within-treatment progression.

```r
# Add after geom_bar():
stat_summary(aes(group = treatment), fun = mean, geom = "line",
             colour = "black", linewidth = 0.8, linetype = "dashed")
```

---

### A2 — Fig 7: equal-radius circles for nested Venns

**File to fix:** `R/utils_ge_profile.R::ge_plot_category_venn` (called by R/09 resistome, R/10 virulome, R/11 mobilome, R/12 network)

**Reality check (2026-06-10):** Original spec diagnosed a "factor-level ordering issue" — that's not the bug. `ge_plot_category_venn` already passes groups in sorted alphabetic order via `split()`, which gives `Dulce / Reference diet / Soyabean meal`. The visual nesting the spec described is correct data (one treatment genuinely has fewer unique VFs); the issue is `VennDiagram::venn.diagram`'s area-scaling makes the small set invisible.

**Fix:** Pass `scaled = FALSE` to the existing `do.call(VennDiagram::venn.diagram, ...)` call so circles are equal-radius regardless of set sizes. Optionally surface as a `cfg$<domain>$venn_scaled` knob (default FALSE).

---

### A3 — Fig 10: switch network to force-directed layout

**File to fix:** `R/12_network.R::.network_plot_ggraph`

**Reality check (2026-06-10):** ✅ Already done. `cfg$network$layout %||% "fr"` is Fruchterman-Reingold; the star layout the spec described was from an earlier pipeline version. Optional follow-up: edge colouring by sign (positive/negative correlation) is not yet wired — edge weights are positive-only (Σ sampleCount or binary 1). To get signed edges we'd need a different network model (e.g. Spearman/Pearson on the gene/taxa profile matrix), which is a v1.2/v2 scope item, not an A-level fix.

---

## Section B — v1.1 carry-forward (kinds documented but not yet emitted)

These are copied verbatim from `PIPELINE_V1.1_GAPS.md`. Mark each as ✅ complete once the pipeline emits it.

| # | `kind` | Slot | Stage | Status |
|---|---|---|---|---|
| B1 | `composition_unique_species` | figS1_unique_sp | relative_abundance | ✅ shipped (fa8d651, 2026-06-10) |
| B2 | `composition_shared_species` | figS2_shared_sp | relative_abundance | ✅ shipped (fa8d651) |
| B3 | `taxa_daa_heatmap` | figS3_top30_daa | differential_abundance | ✅ shipped (fa8d651) |
| B4 | `ge_alpha_shannon_bar` (+ `_violin`) | figS4 B, figS7 B, figS12 B | resistome / virulome / mobilome | ✅ shipped (fa8d651) |
| B5 | `ge_venn_args` | fig4_resistome panel D | resistome | ✅ shipped (fa8d651) |
| B6 | `species_count_genmap` | figS6 (C. diff), figS11 (Entero) | resistome / virulome | ✅ shipped (4422862, 2026-06-10) |
| B7 | `connectivity_venn_taxa` | figS13 panel A | network | ✅ shipped (fa8d651) |
| B8 | `connectivity_venn_genesets` | figS13 panel B | network | ✅ shipped (fa8d651) |
| B9 | `sankey_png` | figS9_vf_sankey | network | ✅ shipped (fa8d651) + height fix (938eec1) |
| B10 | `taxa_ges_composition` | (no slot) | — | 🗑 retired (2026-06-11) — no `fig9_taxa_ges` slot in v2 layout; cross-domain composition is covered by `fig9_network` + `fig10_multiomics_integration` (Sankey + Mantel + mobile-fraction). |

Full emit specs for B1–B9 are in `docs/PIPELINE_V1.1_GAPS.md`. B10 has been retired — see status note above.

---

## Section C — v1.2 new kinds (literature-grounded restructure)

The comparator-paper survey recommends a restructured main-figure layout.  Each item below is a new `kind` entry needed to support that layout.

---

### C1 — `rarefaction_curves`

**Literature basis:** Papers 3 (Zhou et al.), 4 (Nasser et al.) include rarefaction curves in main or supplementary to demonstrate adequate sequencing depth. All 8 papers cite it as a standard QC check.

- **Used by:** `figS1_rarefaction` (new supplementary slot)
- **Stage:** `alpha_diversity` (or `taxonomy` / QC stage)
- **Render hint:** One line per sample, x = number of reads subsampled, y = observed species, lines coloured by treatment group, vertical dashed line at actual sequencing depth. Use `vegan::rarecurve()` output fed to ggplot2.

```jsonc
{
  "path":         "test_run/Figures/alpha_diversity/rarefaction_curves.png",
  "kind":         "rarefaction_curves",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Rarefaction curves for all 18 samples coloured by treatment group, demonstrating adequate sequencing depth."
}
```

---

### C2 — `alpha_violin` — per-treatment-group variant for Richness, Chao1, Simpson

**Literature basis:** All 8 comparator papers use **per-treatment-group** violin + jitter plots (not per-sample bar charts) for alpha diversity. 2×2 layout with Richness, Shannon, Chao1, and Simpson is the dominant convention (Papers 3, 8, Eimeria paper).

**Reality check (2026-06-10):** R/06_alpha_diversity.R already emits violin + bar for *nine* metrics including Shannon, Richness, and Simpson — `richness_violin.png`, `shannon_diversity_index_violin.png`, `simpson_s_index_violin.png` are all present in chicken_batch1's manifest. **Only Chao1 is missing.** The other three are covered by existing `alpha_violin` entries with their `metric` field set.

Scope: add Chao1 estimator computation to R/06 + the violin/bar emissions, then the manifest auto-classifies it as another `alpha_violin` / `alpha_bar` entry.

- **Used by:** New `fig2_alpha_diversity` slot (replaces current `fig1_taxa_alpha` bar-chart layout)
- **Stage:** `alpha_diversity`
- **Render hint:** Per-treatment violin with overlaid jitter (individual sample points), box plot IQR inside violin, Kruskal-Wallis p-value annotated, Dunn post-hoc significance brackets between groups. Consistent treatment colour palette across all four panels.

```jsonc
// Emit one entry per metric — example for Richness:
{
  "path":         "test_run/Figures/alpha_diversity/richness_violin.png",
  "kind":         "alpha_violin",
  "metric":       "Richness",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Observed species richness by dietary treatment group."
}
// Repeat for metric = "Chao1" and metric = "Simpson index"
// Shannon diversity index violin likely already emitted — verify and add if missing
```

---

### C3 — `pcoa_jaccard`

**Literature basis:** All 8 papers show at least two beta-diversity metrics. Bray-Curtis + Jaccard is the most common pairing (Papers 2, 5, 8). Jaccard (presence/absence) captures the rare biosphere and is often the only ordination that shows treatment separation when abundant taxa dominate Bray-Curtis.

- **Used by:** New `fig3_beta_diversity` slot, panel B (alongside existing `pcoa_scatter` Bray-Curtis)
- **Stage:** `beta_diversity`
- **Render hint:** Same styling as the existing Bray-Curtis PCoA: treatment ellipses, PERMANOVA R² + p annotated, PC variance % on axes. Distance = Jaccard binary.

```jsonc
{
  "path":         "test_run/Figures/beta_diversity/jaccard_pcoa.png",
  "kind":         "pcoa_jaccard",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "PCoA of Jaccard (presence/absence) dissimilarity. PERMANOVA R² and p-value annotated."
}
```

---

### C4 — `genus_heatmap`

**Literature basis:** Papers 7 (Fig 2D), 1 (Fig 2) both present a genus-level abundance heatmap with hierarchical clustering as a complement to stacked-bar composition. The heatmap reveals within-group clustering structure that bars cannot show.

- **Used by:** `fig4_taxa_composition` panel C (new combined composition slot)
- **Stage:** `relative_abundance` or `taxonomy`
- **Render hint:** Top 30 genera by mean relative abundance. Rows = genera, columns = 18 samples (grouped by treatment with annotation bar). Hierarchical clustering on both axes (complete linkage, Bray-Curtis distance). Colour scale = log10-transformed relative abundance, white at 0. `pheatmap` or `ComplexHeatmap`.

```jsonc
{
  "path":         "test_run/Figures/relative_abundance/genus_heatmap_top30.png",
  "kind":         "genus_heatmap",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Heatmap of the top 30 genera by mean relative abundance, hierarchically clustered on both axes."
}
```

---

### C5 — `aldex2_maplot`

**Literature basis:** ALDEx2 is the differential-abundance method already used in this pipeline. MA plots (log-ratio vs. mean abundance) are the standard ALDEx2 output visualisation (Papers 2, 7 use analogous per-comparison DA plots). Dedicated DA figures are present as main figures in 6 of the 8 papers.

- **Used by:** New `fig5_differential_abundance` slot, panels A–C (one per pairwise comparison)
- **Stage:** `differential_abundance`
- **Render hint:** x-axis = mean log2 CLR-transformed abundance (A), y-axis = between-condition difference (M). Significant features (BH-adjusted expected overlap < 0.05 AND |effect| > 0.5) coloured by direction (up = red, down = blue). Top 10 significant taxa labelled. One panel per comparison: Dulce vs. Reference, Dulce vs. Soyabean, Reference vs. Soyabean. Use `aldex.plot()` output or rebuild in ggplot2.

Emit one entry per comparison using the `pair` field:

```jsonc
{
  "path":         "test_run/Figures/differential_abundance/aldex2_maplot_dulce_vs_reference.png",
  "kind":         "aldex2_maplot",
  "format":       "png",
  "pair":         ["Dulce", "Reference diet"],
  "caption_seed": "ALDEx2 MA plot: Dulce seaweed vs. Reference diet. Significant features (BH p < 0.05, |effect| > 0.5) coloured by direction."
}
```

```jsonc
{
  "path":         "test_run/Figures/differential_abundance/aldex2_maplot_dulce_vs_soyabean.png",
  "kind":         "aldex2_maplot",
  "pair":         ["Dulce", "Soyabean meal"],
  ...
}
```

```jsonc
{
  "path":         "test_run/Figures/differential_abundance/aldex2_maplot_reference_vs_soyabean.png",
  "kind":         "aldex2_maplot",
  "pair":         ["Reference diet", "Soyabean meal"],
  ...
}
```

---

### C6 — `aldex2_dotplot`

**Literature basis:** A cross-comparison effect-size summary (dot plot or forest plot) synthesising the three pairwise ALDEx2 results into one panel is used in Papers 1 (Fig 7) and 2 (Fig 6). It makes the manuscript text tractable: "taxon X was significantly enriched in Dulce vs. both other treatments."

- **Used by:** `fig5_differential_abundance` slot, panel D
- **Stage:** `differential_abundance`
- **Render hint:** Rows = taxa significant in at least one comparison. Columns = three pairwise comparisons (faceted or grouped). x-axis = ALDEx2 effect size (with 95% CI). Dot colour = direction. Annotate BH-adjusted significance as * / ** / ***. Requires combining output from all three ALDEx2 `aldex.ttest()` calls.

```jsonc
{
  "path":         "test_run/Figures/differential_abundance/aldex2_dotplot_summary.png",
  "kind":         "aldex2_dotplot",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Effect-size summary across all three pairwise ALDEx2 comparisons. Taxa significant in at least one comparison shown."
}
```

---

### C7 — `arg_upset`

**Literature basis:** Paper 8 (Fig 4D) uses an UpSet plot for the ARG landscape to show treatment-specific and shared gene families. UpSet scales better than a 3-set Venn when many intersections are present and the intersection sizes vary widely, which is common in resistome studies.

- **Used by:** New `fig7_resistome_composition` slot, panel B
- **Stage:** `resistome`
- **Render hint:** Sets = three treatment groups. Elements = ARG gene family identifiers. Bar chart of intersection sizes at top. `UpSetR::upset()` or `ComplexHeatmap::UpSet()`. Minimum intersection size = 1 (show all). Colour bars by dominant treatment group.

```jsonc
{
  "path":         "test_run/Figures/resistome/arg_upset_treatments.png",
  "kind":         "arg_upset",
  "format":       "png",
  "domain":       "resistome",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "UpSet plot of ARG gene families shared and unique across dietary treatment groups."
}
```

---

### C8 — `arg_circos`

**Literature basis:** Papers 5 (Fig 3B) and 7 (Fig 5B) use Circos plots to show proportional contribution of drug classes to total ARG load and/or their distribution across taxa or treatment groups. More compact than stacked bars when there are 15+ drug classes.

- **Used by:** New `fig7_resistome_composition` slot, panel C
- **Stage:** `resistome`
- **Render hint:** Outer tracks = drug classes (arc width proportional to total TPM). Inner chords = links to treatment group or sample cluster. Colours per drug class consistent with `fig5_drug_heatmap`. Use `circlize` R package. Not the same as `chord` (which links ARGs to taxa — that is Fig 9 / Fig 6).

```jsonc
{
  "path":         "test_run/Figures/resistome/arg_circos_drugclass.png",
  "kind":         "arg_circos",
  "format":       "png",
  "domain":       "resistome",
  "caption_seed": "Circos plot showing proportional contribution of antibiotic drug classes to total ARG load across dietary treatment groups."
}
```

---

### C9 — `vf_arg_corr_heatmap`

**Literature basis:** Paper 1 (Fig 7) presents a Spearman correlation matrix between VF categories and ARG drug classes per treatment. Paper 8 includes this in supplementary. Condensing to a single cross-category heatmap is the standard approach when the study spans both virulome and resistome.

- **Used by:** New `fig8_virulome_mobilome` slot, panel D
- **Stage:** `resistome` + `virulome` combined (post-merge step)
- **Render hint:** Rows = VF functional categories (adherence, invasion, toxin, etc.). Columns = ARG drug classes. Values = Spearman r. Colour = red-white-blue divergent scale (red = positive, blue = negative). Significance stars (* p < 0.05, ** p < 0.01, *** p < 0.001 after BH correction) in each cell. `pheatmap` or `corrplot`.

```jsonc
{
  "path":         "test_run/Figures/resistome/vf_arg_correlation_heatmap.png",
  "kind":         "vf_arg_corr_heatmap",
  "format":       "png",
  "caption_seed": "Spearman correlation heatmap between virulence-factor functional categories and antibiotic drug classes."
}
```

---

### C10 — `mantel_triangle`

**Literature basis:** Paper 8 (Supplementary Fig S1C) uses a Mantel test to quantify pairwise community-level correlations across omics layers. Bringing it to the main integration figure is justified given this study spans four layers (taxonomy, resistome, virulome, mobilome).

- **Used by:** New `fig10_multiomics_integration` slot, panel B
- **Stage:** `network` or a dedicated `integration` stage
- **Render hint:** Symmetric matrix (triangle) with pairs: taxonomy ↔ resistome, taxonomy ↔ virulome, taxonomy ↔ mobilome, resistome ↔ virulome, resistome ↔ mobilome, virulome ↔ mobilome. Each cell = Mantel r (font) + permutation p-value. Colour intensity = |r|. Use `vegan::mantel()` with Bray-Curtis distance matrices and 9999 permutations.

```jsonc
{
  "path":         "test_run/Figures/network/mantel_correlation_triangle.png",
  "kind":         "mantel_triangle",
  "format":       "png",
  "caption_seed": "Mantel test correlation triangle showing pairwise community-level correlations between taxonomic, resistome, virulome, and mobilome distance matrices."
}
```

---

### C11 — `sankey_taxon_arg_mge`

> **Status (2026-06-11):** ✅ Shipped. Phylum surface delivered via a kraken2 pre-order tree walk in `R/02_clean_data.R::build_taxid_phylum`, sidestepping the rcf Rank/Name path — no R/02 rank-aware rewrite required (that rewrite stays queued in `docs/pipeline_rework_scoping.md` for the other rcf-slicing issues). C11 itself lives in `R/12_network.R::.network_build_sankey_taxon_arg_mge`; same CARD + PlasmidFinder co-occurrence definition as C12 `mobile_fraction_bar`. Skips cleanly on chicken_batch1 (no CARD+PF co-occurrence); render path verified with a synthetic three-contig fixture.


**Literature basis:** Paper 8 (Fig 5 area) uses a Sankey/alluvial diagram for taxon → ARG class → MGE type flow. Paper 5 uses Circos for similar multi-omic flow. This is distinct from `sankey_png` (VF taxa-to-function) and `arg_circos` (drug classes only).

- **Used by:** New `fig10_multiomics_integration` slot, panel A
- **Stage:** `network` or `integration`
- **Render hint:** Three columns of nodes: taxonomic phylum (left), ARG drug class (centre), MGE type (right). Ribbon width = relative contribution by mean TPM. Use `ggalluvial` in R or `networkD3::sankeyNetwork()` rendered to PNG. Consistent colour palette with rest of manuscript.

```jsonc
{
  "path":         "test_run/Figures/network/sankey_taxon_arg_mge.png",
  "kind":         "sankey_taxon_arg_mge",
  "format":       "png",
  "caption_seed": "Sankey diagram showing relative abundance flow from bacterial phylum to ARG drug class to MGE type across the chicken caecal metagenome."
}
```

---

### C12 — `mobile_fraction_bar`

**Literature basis:** Papers 5 and 8 explicitly quantify the mobilisable fraction of the resistome — what fraction of ARGs (and VFs) are found on predicted mobile contigs. This is a key biosafety conclusion in dietary intervention studies and is increasingly expected by Frontiers reviewers.

- **Used by:** New `fig10_multiomics_integration` slot, panel C
- **Stage:** `network` or `mobilome` (requires contig-level linkage of ARG + MGE)
- **Render hint:** Stacked bar chart, one bar per treatment group (and one total). y-axis = percentage of total ARG TPM. Stacks = mobile (carried on a contig also harbouring an IS element, integron, or plasmid replicon) vs. non-mobile. Repeat for VFs if data allows. Statistical annotation for treatment differences.

```jsonc
{
  "path":         "test_run/Figures/network/mobile_arg_fraction_bar.png",
  "kind":         "mobile_fraction_bar",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Proportion of ARG abundance carried on predicted mobile genetic elements (contigs co-harbouring MGE markers) by treatment group."
}
```

---

## Proposed restructured figure slot layout (v2)

This is the updated `journals/frontiers_microbiome.yaml` slot layout recommended by the literature survey. Update the YAML once the pipeline emits the required kinds.

| Slot | Title | Panels | Layout | New kinds needed |
|---|---|---|---|---|
| `fig1_taxa_alpha` *(rename → `fig1_study_overview`)* | Study overview and QC | A: study design schematic, B: rarefaction curves | 1×2 | `rarefaction_curves` |
| `fig2_alpha_diversity` *(was fig1)* | Alpha diversity across treatments | A: Richness violin, B: Shannon violin, C: Chao1 violin, D: Simpson violin | 2×2 | `alpha_violin` (Chao1, Simpson metrics) |
| `fig3_beta_diversity` *(was panel in fig2)* | Beta diversity and community structure | A: Bray-Curtis PCoA, B: Jaccard PCoA | 1×2 | `pcoa_jaccard` |
| `fig4_taxa_composition` *(was fig2+fig3)* | Taxonomic composition | A: phylum stacked bar, B: genus stacked bar, C: genus heatmap (top 30) | 2+1 row | `genus_heatmap` |
| `fig5_differential_abundance` *(new)* | Differential abundance (ALDEx2) | A–C: MA plots × 3 comparisons, D: effect-size summary | 2×2 | `aldex2_maplot`, `aldex2_dotplot` |
| `fig6_resistome_diversity` *(was part of fig4)* | Resistome diversity and ordination | A: ARG richness violin, B: ARG Shannon violin, C: drug-class grouped bar, D: ARG PCoA | 2×2 | — (existing kinds) |
| `fig7_resistome_composition` *(was fig4+5)* | Resistome composition | A: top-40 ARG heatmap, B: ARG UpSet, C: drug-class Circos | 2+1 row | `arg_upset`, `arg_circos` |
| `fig8_virulome_mobilome` *(was fig7+fig8 merged)* | Virulome and mobilome profiles | A: VF category bar, B: MGE type pie, C: MGE family bar, D: VF×ARG correlation heatmap | 2×2 | `vf_arg_corr_heatmap` |
| `fig9_network` *(was fig10)* | ARG–taxa co-occurrence network | A: force-directed network, B: top-10 degree bar | 2:1 split | — (existing `network_graph`) |
| `fig10_multiomics_integration` *(was fig9 area)* | Multi-omic integration | A: taxon→ARG→MGE Sankey, B: Mantel triangle, C: mobile fraction bar | 1+2 row | `sankey_taxon_arg_mge`, `mantel_triangle`, `mobile_fraction_bar` |

### Supplementary slots (v2)

| Slot | Title | Panels | New kinds needed |
|---|---|---|---|
| `figS1_rarefaction` *(new)* | Rarefaction curves | A: all samples per treatment | `rarefaction_curves` |
| `figS2_unique_sp` *(was S1)* | Unique-species composition | A | `composition_unique_species` (B-carry) |
| `figS3_shared_sp` *(was S2)* | Shared-species composition | A | `composition_shared_species` (B-carry) |
| `figS4_taxa_levels` *(new)* | Full taxonomic breakdown | A: class, B: order, C: family, D: species | existing `composition_stacked_bar` at each level |
| `figS5_top30_daa` *(was S3)* | Top-30 DAA heatmap | A | `taxa_daa_heatmap` (B-carry) |
| `figS6_arg_alpha` *(was S4)* | ARG alpha diversity | A: richness bar, B: Shannon bar | `ge_alpha_shannon_bar` (B-carry) |
| `figS7_arg_pheatmap` *(was S5)* | Full ARG gene heatmap | A | existing `ge_pheatmap_genes` |
| `figS8_arg_beta` *(new)* | ARG beta diversity (extended) | A: Jaccard PCoA, B: NMDS | `pcoa_jaccard` (domain: resistome), `nmds_scatter` |
| `figS9_vf_alpha` *(was S7)* | VF alpha diversity | A: richness bar, B: Shannon bar | `ge_alpha_shannon_bar` (B-carry) |
| `figS10_vf_pheatmap` *(was S8)* | VF function heatmap | A | existing `ge_pheatmap_functions` |
| `figS11_vf_all_pheatmap` *(was S10)* | All-VF gene heatmap (A: adhesion, B: secretion) | A, B | existing `ge_pheatmap_genes` |
| `figS12_vf_sankey` *(was S9)* | VF taxa-to-function sankey | A | `sankey_png` (B-carry) |
| `figS13_cdiff_genmap` *(was S6)* | *C. difficile* count + gene map | A, B | `species_count_genmap` (B-carry) |
| `figS14_entero_genmap` *(was S11)* | Enterobacteriaceae count + gene map | A, B | `species_count_genmap` (B-carry) |
| `figS15_mge_alpha` *(was S12)* | MGE alpha diversity | A: richness bar, B: Shannon bar | `ge_alpha_shannon_bar` (B-carry) |
| `figS16_connectivity` *(was S13)* | Cross-domain connectivity | A: taxa Venn, B: geneset Venn, C: Shannon violin, D: PCoA | existing (B-carry) |
| `figS17_corr_extended` *(new)* | Extended correlation matrices | A: genus × ARG heatmap, B: MGE × ARG heatmap | `vf_arg_corr_heatmap` or separate `extended_corr_heatmap` |

---

## Summary checklist

### Immediate (A-level, no manifest change)
- [x] A1 — Add trend-line overlay to alpha-bar plots in R (per-group mean line, fa8d651 → 25a8438)
- [x] A2 — Reality check: not a group-ordering bug; real fix is `scaled=FALSE` (still pending render-side change)
- [x] A3 — Reality check: already Fruchterman-Reingold; no change needed

### Short-term (B-level, manifest v1.1, already specified)
- [x] B1–B5, B7–B9 shipped fa8d651
- [x] B6 species_count_genmap shipped 4422862
- [x] B10 taxa_ges_composition — retired 2026-06-11 (no v2 slot; coverage already provided by fig9_network + fig10_multiomics_integration)

### Medium-term (C-level, manifest v1.2)
- [x] C1 `rarefaction_curves` (25a8438)
- [x] C2 `alpha_violin` for Chao1 (Richness + Simpson already emitted; Shannon already there)
- [x] C3 `pcoa_jaccard` (25a8438)
- [x] C4 `genus_heatmap` (13ba39e)
- [x] C5 `aldex2_maplot` (13ba39e, 3 pairwise comparisons)
- [x] C6 `aldex2_dotplot` (13ba39e; fires only when ≥1 feature significant)
- [x] C7 `arg_upset` (25a8438) + `vf_upset` + `mge_upset` (087b5ab)
- [x] C8 `arg_circos` (13ba39e + circos label fix 087b5ab)
- [x] C9 `vf_arg_corr_heatmap` (13ba39e)
- [x] C10 `mantel_triangle` (4422862)
- [x] C11 `sankey_taxon_arg_mge` (2026-06-11; kraken2 pre-order phylum walk in R/02 + ggalluvial render in R/12; skips cleanly when no CARD+PF co-occurrence)
- [x] C12 `mobile_fraction_bar` (4422862; fires only when contigs co-harbour CARD+PlasmidFinder)

### After C-level is emitted
- [x] Bump `manifest_version` to `"1.2"` in R pipeline writer
- [ ] Update `journals/frontiers_microbiome.yaml` figure slots to v2 layout (agent repo)
- [ ] Re-run `scratch/demo_real.py` — expect 10/10 main + 17/17 supplementary (agent repo)

---

## Notes on backward compatibility

- All new kinds are **optional** from the manifest schema's perspective — the agent reflows or skips slots that don't resolve.
- Renaming existing slots in `frontiers_microbiome.yaml` (e.g. `fig1_taxa_alpha` → `fig1_study_overview`) is a journal config change, not a manifest contract change — no version bump required for that.
- The `SUPPORTED_MAJOR` constant in the agent remains `1`; it accepts `1.0`, `1.1`, and `1.2` without modification.
