# Pipeline v1.1 — manifest gaps to close

Punch-list for [`Metagenomics_pipeline_automation`](https://github.com/Julio92-C/Metagenomics_pipeline_automation) so a future run emits everything `metaomics-scribe` needs to compose the full chicken_batch1 Frontiers manuscript (10 main + 13 supplementary figures) without ad-hoc filename parsing.

Each entry below lists:

- **`kind`** — the controlled vocabulary string the agent matches on
- **Used by** — which manuscript figure slot in `journals/frontiers_microbiome.yaml` needs it
- **Emit as** — a manifest JSON snippet to add to the named stage's `figures` array, plus the file path convention to write the PNG to (relative to `outputs.project_root`, per `MANIFEST_SCHEMA.md` §`outputs`)
- **Render hint** — short description of what should actually be plotted; matches what was hand-rendered for the manuscript

This is a **MINOR** bump (`1.0` → `1.1`): purely additive. Old manifests continue to validate against the agent; the agent reflows or skips slots whose new kinds aren't yet emitted.

The full v1.1 contract is in `docs/MANIFEST_SCHEMA.md` (changelog at the bottom + *(1.1)*-tagged rows in the `kind` vocabulary table). This file is the actionable subset.

---

## 1. `composition_unique_species`

- **Used by:** `figS1_unique_sp` (Fig S1)
- **Stage:** `relative_abundance`
- **Render hint:** the existing stacked-bar composition plot, filtered to species that are unique to one treatment group (not shared with any other group). Same layout, axes, and legend style as the existing `composition_stacked_bar`.

```jsonc
{
  "path":         "test_run/Figures/relative_abundance/unique_species_relative_abundance.png",
  "kind":         "composition_unique_species",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Relative abundance of species unique to each treatment group."
}
```

## 2. `composition_shared_species`

- **Used by:** `figS2_shared_sp` (Fig S2)
- **Stage:** `relative_abundance`
- **Render hint:** same as above but filtered to species shared by **all** treatment groups.

```jsonc
{
  "path":         "test_run/Figures/relative_abundance/shared_species_relative_abundance.png",
  "kind":         "composition_shared_species",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Relative abundance of species shared across all treatment groups."
}
```

## 3. `taxa_daa_heatmap`

- **Used by:** `figS3_top30_daa` (Fig S3)
- **Stage:** `differential_abundance`
- **Render hint:** sample-level heatmap of the top 30 differentially-abundant taxa (rows = taxa, columns = samples grouped by treatment); the existing `top_candidates.pdf` content but as a publication-grade PNG with the colour bar and treatment header strip on top.

```jsonc
{
  "path":         "test_run/Figures/differential_abundance/taxa/top30_daa_heatmap.png",
  "kind":         "taxa_daa_heatmap",
  "format":       "png",
  "domain":       "taxa",
  "caption_seed": "Sample-level heatmap of the top 30 differentially-abundant taxa across treatment groups."
}
```

## 4. `ge_alpha_shannon_bar`

- **Used by:** `figS4_arg_alpha`, `figS7_vf_alpha`, `figS12_mge_alpha` (Figs S4, S7, S12)
- **Stages:** `resistome`, `virulome`, `mobilome` (one entry per stage)
- **Render hint:** matches the existing `ge_alpha_richness_bar` style — per-sample bar coloured by treatment, with the Kruskal-Wallis p-value annotated — but plotting the Shannon diversity index instead of richness. Set `domain` to the gene-set domain.

```jsonc
{
  "path":         "test_run/Figures/resistome/alpha_shannon_bar.png",
  "kind":         "ge_alpha_shannon_bar",
  "format":       "png",
  "domain":       "resistome",
  "caption_seed": "Shannon diversity of the resistome by treatment."
}
```

Emit one per domain (`resistome`, `virulome`, `mobilome`), changing the `path`, `domain`, and `caption_seed` accordingly.

Optionally also emit `ge_alpha_shannon_violin` for symmetry with the existing `ge_alpha_richness_violin`; not required by any current Frontiers slot but cheap to add while you're there.

## 5. `ge_venn_args`

- **Used by:** `fig4_resistome` panel D (Fig 4)
- **Stage:** `resistome`
- **Render hint:** 3-set Venn diagram of ARG genes (not drug classes) shared and unique across treatment groups. The current manifest already has `ge_venn_drug_classes`; this is the gene-level counterpart.

```jsonc
{
  "path":         "test_run/Figures/resistome/venn_args.png",
  "kind":         "ge_venn_args",
  "format":       "png",
  "domain":       "resistome",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Venn diagram of ARG genes shared and unique across treatment groups."
}
```

## 6. `species_count_genmap`

- **Used by:** `figS6_cdiff_genmap` (Fig S6 — *C. difficile*), `figS11_entero_genmap` (Fig S11 — Enterobacteriaceae)
- **Stage:** `resistome` (and/or `virulome` depending on the organism)
- **Render hint:** a 2-panel figure for a per-organism deep-dive: panel A = per-sample count of the organism by treatment with Kruskal-Wallis p-value; panel B = gene neighbourhood map (gggenes-style arrows) per sample. The pipeline currently emits this as a hand-made TIFF; we need each panel as its own PNG so the journal slot can stitch them.

Emit **two** entries with the same `kind` but distinguished by the new `panel` annotation field — or, simpler, by adding `count` and `genmap` to the path so order is deterministic:

```jsonc
{
  "path":         "test_run/Figures/resistome/c_difficile_count_per_treatment.png",
  "kind":         "species_count_genmap",
  "format":       "png",
  "domain":       "resistome",
  "organism":     "Clostridioides difficile",
  "caption_seed": "Clostridioides difficile per-sample count grouped by treatment."
}
```

```jsonc
{
  "path":         "test_run/Figures/resistome/c_difficile_gene_map.png",
  "kind":         "species_count_genmap",
  "format":       "png",
  "domain":       "resistome",
  "organism":     "Clostridioides difficile",
  "caption_seed": "Per-sample ARG neighbourhood map for Clostridioides difficile."
}
```

Repeat the pair for Enterobacteriaceae (Fig S11). The `organism` annotation is a new optional field — see also `Figure.organism` in the v1.1 schema doc bump.

## 7. `connectivity_venn_taxa`

- **Used by:** `figS13_connectivity` panel A (Fig S13)
- **Stage:** `network`
- **Render hint:** 3-set Venn of **taxa** shared and unique across treatment groups, derived from the co-occurrence-network member list (not the raw taxa from `taxonomy`).

```jsonc
{
  "path":         "test_run/Figures/network/connectivity_venn_taxa.png",
  "kind":         "connectivity_venn_taxa",
  "format":       "png",
  "groups":       ["Dulce", "Reference diet", "Soyabean meal"],
  "caption_seed": "Network-connected taxa shared and unique across treatment groups."
}
```

## 8. `connectivity_venn_genesets`

- **Used by:** `figS13_connectivity` panel B (Fig S13)
- **Stage:** `network`
- **Render hint:** 3-set Venn of **gene elements** (ARGs ∩ VFs ∩ MGEs) participating in the co-occurrence network across the three treatments.

```jsonc
{
  "path":         "test_run/Figures/network/connectivity_venn_genesets.png",
  "kind":         "connectivity_venn_genesets",
  "format":       "png",
  "groups":       ["ARGs", "VFs", "MGEs"],
  "caption_seed": "Network-connected gene elements (ARGs / VFs / MGEs)."
}
```

## 9. `sankey_png`

- **Used by:** `figS9_vf_sankey` (Fig S9)
- **Stage:** `network`
- **Render hint:** static PNG render of the existing interactive HTML sankey (`sankey_overall.html`). The HTML version is fine for the supplementary website; the manuscript needs an image. Same data, no styling changes.

```jsonc
{
  "path":         "test_run/Figures/network/sankey_overall.png",
  "kind":         "sankey_png",
  "format":       "png",
  "caption_seed": "VF taxa-to-function sankey (PNG render of the interactive supplementary)."
}
```

Keep emitting the HTML `sankey` entry too — the website needs it. They coexist as separate entries pointing at different files.

---

## Updates to existing entries (no new `kind`)

These don't need new kinds, just additional optional annotations on existing figure entries:

- **`Figure.domain`** *(already in v1.1 schema, just start populating it consistently)* — set `domain: "resistome"`, `"virulome"`, or `"mobilome"` on every `ge_*` figure so the agent can disambiguate `ge_pcoa` (etc.) between domains. Currently the chicken_batch1 manifest sets this on most `ge_*` entries; just make sure new ones do too.
- **`Figure.organism`** *(new optional field, v1.1)* — see §6 above. Lets the agent target a specific per-organism deep-dive plot without parsing filenames.

---

## Validation checklist

When the pipeline ships v1.1, re-run the agent demo against the new project:

```powershell
uv run python scratch/demo_real.py
```

A successful v1.1 run will show **10/10 main + 13/13 supplementary** under `runs/chicken_batch1/figures/` instead of the current 10/10 + 8/13.

Bump `manifest_version` to `"1.1"` in the pipeline's writer module (R/14_manifest.R or wherever it's emitted) once these are in. The agent's `SUPPORTED_MAJOR` is `1`, so it accepts any `1.x` automatically — no agent-side change required.
