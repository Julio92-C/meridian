# Hospital alpha-diversity panel polish — design spec

**Date:** 2026-09-23
**Status:** Approved (design), pending implementation
**Scope:** project-specific (hospital_microbiome only). **Does NOT touch the
pipeline** (`R/NN_*.R`, `R/utils_*.R`, `templates/`, or any config schema).

## Goal

Produce a publication-tidied version of the custom richness + shannon
per-sample bar panel for the hospital paper
(`Figures/alpha_diversity/alpha_richness_shannon_panel.png`). The current
image is a one-off composite with no generating script in the repo; this spec
introduces a small, self-contained script that rebuilds it with cleaner
styling while keeping the **same layout**.

Requirement (from the user): **keep the existing A/B faceted bar layout — just
tidy the styling**. No collision-repel, no re-ordering, no reference-line
redesign.

## Non-goals

- No changes to any pipeline module, util, or config schema.
- No new metrics, no layout redesign, no sample re-ordering.
- No TIFF output (the original is PNG-only).

## Data source (read-only pipeline output)

`<project_root>/<datasets_dir>/alpha_diversity.csv`, already emitted by
`R/06_alpha_diversity.R`. Columns: `sample, richness, shannon, simpson, chao1,
ID, type, description, ward, location`. The script consumes this CSV directly —
it does not recompute diversity and does not source any pipeline R file.

## Component

**New file:** `projects/hospital_microbiome/polish_alpha_panel.R` — standalone
`Rscript` entry point.

Dependencies (all already in `renv`): `yaml`, `readr`, `dplyr`, `tidyr`,
`ggplot2`, `patchwork` (cowplot fallback for the A/B tags + stacking).

### Flow

1. Read `projects/hospital_microbiome/config.yaml` to resolve `project_root`,
   `outputs$datasets_dir`, `outputs$figures_dir`, and the `colors$type` map
   (Air `#828385` / Surface `#ed8cb6` — reused as the single colour source so
   the polish stays in sync with the pipeline recolor).
2. Read `alpha_diversity.csv`; keep `sample`, `type`, `richness`, `shannon`.
3. For each metric (`richness`, `shannon`) build one faceted bar panel:
   - `geom_bar(stat = "identity")`, fill = `type`, `scale_fill_manual` grey/pink.
   - Red dashed **overall-mean** hline; black dashed **per-group-mean** hlines
     (per facet) — identical semantics to R/06.
   - Per-bar value labels **removed** — with 64 samples they crowded/clipped
     even rotated; the bars + y-axis carry the values. (Iteration history:
     retained → 2x-scaled → vertical → dropped for the cleanest read.)
   - `facet_wrap(~ type, scales = "free_x", nrow = 1)`.
   - KW p + BH p_adj computed from the CSV (`kruskal.test`, `p.adjust(method
     = "BH")` — mirrors the pipeline) and placed in a **light subtitle**, not
     the bold title.
4. Combine richness (top) + shannon (bottom) with `patchwork` (`plot_annotation
   (tag_levels = "A")` for bold **A**/**B** tags).
5. Write the composite to
   `<project_root>/<figures_dir>/alpha_diversity/alpha_richness_shannon_panel.png`,
   **overwriting** the original, at 300 dpi, ~13 × 7 in.

### Tidy styling knobs (the actual "polish")

- **Word-readability:** the 13x7in export is shrunk ~50% in Word, halving
  apparent font size. A single `word_scale` (=2) multiplies an explicit,
  BALANCED `sz` hierarchy (title / subtitle / axis title / axis text / sample
  IDs / strip / tag) so text stays legible on the page *and* the title stays
  only slightly larger than the body — avoiding the runaway-heading bug that a
  blanket `base_size * 2` caused.
- bold centred metric title, light grey subtitle carrying KW / p_adj.
- `axis.text.x`: 45°, `hjust = 1`.
- tuned `geom_col(width = 0.8)`, clean grey facet strips, `legend.position =
  "none"` (facets already name the groups).

## Testing / validation

No unit tests (single figure script). Validation is visual:

1. Run `Rscript projects/hospital_microbiome/polish_alpha_panel.R`; expect exit
   0 and the overwritten PNG.
2. Read the PNG back and confirm: same A/B faceted layout, grey/pink bars, both
   dashed mean lines present, KW/BH stats now in the subtitle, cleaner
   typography, no regressions in the numbers vs the source CSV.
3. User sign-off on the rendered figure.

Safety: the target lives in the gitignored `test_run/` sandbox and is fully
regenerable from the CSV, so overwriting is low-risk.

## Files touched

| File | Change |
| --- | --- |
| `projects/hospital_microbiome/polish_alpha_panel.R` | new — standalone polish script |
| `docs/superpowers/specs/2026-09-23-hospital-alpha-panel-polish-design.md` | new — this spec |

## Open items

None. Colours, layout, output path, and format are all fixed above.
