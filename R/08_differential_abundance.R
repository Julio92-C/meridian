# 08_differential_abundance.R — ALDEx2 DA across group pairs.
# Matches the paper: 128 MC Dirichlet instances, low-prevalence filter n<2,
# CLR-transformed pairwise comparisons.

run_aldex <- function(cleaned, cfg) {
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    pipeline_log(cfg, "ALDEx2 not installed — skipping DA")
    return(invisible(NULL))
  }
  pipeline_log(cfg, "Differential abundance (ALDEx2)")

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid   <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]
  mc    <- cfg$stats$aldex_mc_samples %||% 128
  minp  <- cfg$filters$min_prevalence_aldex %||% 2

  counts <- cleaned$noncontaminants |>
    dplyr::group_by(sample, name) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0)

  gene_names <- counts$name
  mat <- as.matrix(counts[, -1])
  rownames(mat) <- gene_names
  keep <- rowSums(mat > 0) >= minp
  mat  <- mat[keep, , drop = FALSE]

  sample_meta <- meta[match(colnames(mat), meta[[sid]]), , drop = FALSE]
  conds <- as.character(sample_meta[[group]])

  ds_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "aldex2")
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  pairs <- utils::combn(unique(conds), 2, simplify = FALSE)
  for (pr in pairs) {
    keep_s <- conds %in% pr
    x <- ALDEx2::aldex(
      mat[, keep_s, drop = FALSE], conds[keep_s],
      mc.samples = mc, test = "t", effect = TRUE
    )
    out_name <- paste0(pr[1], "_vs_", pr[2], "_aldex2.tsv")
    readr::write_tsv(tibble::rownames_to_column(x, "feature"),
                     file.path(ds_dir, out_name))
    pipeline_log(cfg, sprintf("ALDEx2 %s vs %s → %s", pr[1], pr[2], out_name))
  }
  invisible(NULL)
}
