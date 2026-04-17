# 09_resistome.R — ARG profile analysis. Uses CARD rows from ABRicate output.
# Produces Venn, pHeatmap, relative-abundance bar plot, and group-level
# Kruskal-Wallis stats per §3 "Resistome" of the Chicken batch 1 paper.

run_resistome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Resistome (CARD)")
  .run_ge_block(cleaned, cfg, db = "card", label = "ARG",
                out_subdir = "resistome")
}
