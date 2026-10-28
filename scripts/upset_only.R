# Re-draw ONLY the resistome / virulome / mobilome UpSet PNGs from the
# already-written genetable_normdata.csv, without re-running the stages
# (so PERMANOVA / KW outputs on disk stay untouched). Use after changing
# cfg$panels$upset_text_scale, then run scripts/panels_only.R.
#
# Usage:
#   Rscript scripts/upset_only.R projects/<slug>/config.yaml

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript scripts/upset_only.R <config.yaml>")
`%||%` <- function(a, b) if (is.null(a)) b else a

mods  <- sort(list.files("R", pattern = "^\\d+_.*\\.R$",  full.names = TRUE))
utils <- sort(list.files("R", pattern = "^utils_.*\\.R$", full.names = TRUE))
invisible(lapply(c(mods, utils), source))

cfg   <- load_config(args[[1]])
norm  <- ge_load_norm_table(cfg, "UpSet")
stopifnot(!is.null(norm))
group <- cfg$metadata$group_cols[[1]]

domains <- list(
  list(key = "resistome", db = "card",          label = "Resistome",
       file = "arg_upset_treatments.png", ylab = "ARGs count"),
  list(key = "virulome",  db = "vfdb",          label = "Virulome",
       file = "vf_upset_treatments.png",  ylab = "VFs count"),
  list(key = "mobilome",  db = "plasmidfinder", label = "Mobilome",
       file = "mge_upset_treatments.png", ylab = "MGEs count")
)

for (d in domains) {
  dcfg <- cfg[[d$key]] %||% list()
  dbs  <- tolower(dcfg$database %||% d$db)
  df   <- norm[tolower(norm$DATABASE) %in% dbs, , drop = FALSE]
  if (nrow(df) == 0) { message(d$label, ": no rows — skipped"); next }
  df <- ge_ensure_group_column(df, cfg, group, d$label)
  lv <- sort(unique(as.character(df[[group]])))
  pal <- resolve_top_level_colors(group, lv, cfg) %||%
         ge_named_palette(lv, dcfg$group_colors,
                          dcfg$group_palette %||% "ggsci::default_nejm")
  out <- file.path(cfg$project_root, cfg$outputs$figures_dir, d$key, d$file)
  ge_plot_category_upset(df, category_col = "GENE", group = group,
                         file = out, log_label = d$label, cfg = cfg,
                         pal_group = pal, value_label = d$ylab)
  message("wrote ", out)
}
