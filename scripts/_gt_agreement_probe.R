#!/usr/bin/env Rscript
# One-shot probe: for the few stages whose GT artefact lives next to ours,
# report headline agreement (row counts, key overlap, Jaccard on the key sets).
# Not a full per-stage agreement column — just enough to scope what a real
# column could carry.

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

root <- "C:/Users/julio/Desktop/PC_JC_2024-11-29_Julio_Gallus"

compare_keyset <- function(label, ours_path, gt_path, key_cols) {
  cat("=== ", label, " ===\n", sep = "")
  if (!file.exists(ours_path)) { cat("  ours MISSING: ", ours_path, "\n", sep = ""); return(invisible()) }
  if (!file.exists(gt_path))   { cat("  gt   MISSING: ", gt_path,   "\n", sep = ""); return(invisible()) }
  ours <- suppressWarnings(read_csv(ours_path, show_col_types = FALSE))
  gt   <- suppressWarnings(read_csv(gt_path,   show_col_types = FALSE))
  cat(sprintf("  ours: %d rows x %d cols (%s)\n", nrow(ours), ncol(ours),
              paste(colnames(ours), collapse = ", ")))
  cat(sprintf("  gt  : %d rows x %d cols (%s)\n", nrow(gt),   ncol(gt),
              paste(colnames(gt),   collapse = ", ")))
  ok_keys <- all(key_cols %in% colnames(ours)) && all(key_cols %in% colnames(gt))
  if (!ok_keys) {
    cat("  KEY COLUMNS MISSING — can't compute Jaccard\n\n"); return(invisible())
  }
  ok <- function(df) do.call(paste, c(df[, key_cols, drop = FALSE], sep = "|"))
  ours_k <- ok(ours); gt_k <- ok(gt)
  inter  <- length(intersect(ours_k, gt_k))
  uni    <- length(union(ours_k, gt_k))
  cat(sprintf("  intersect : %d  | union : %d  | Jaccard = %.3f\n", inter, uni, inter/uni))
  cat(sprintf("  only ours : %d  | only gt : %d\n",
              length(setdiff(ours_k, gt_k)), length(setdiff(gt_k, ours_k))))
  cat("\n")
}

# clean_data — R/02 intermediates
compare_keyset("clean_data: abri_kraken2_cleaned",
               file.path(root, "test_run/Datasets/abri_kraken2_cleaned.csv"),
               file.path(root, "Datasets/abri_kraken2_cleaned.csv"),
               c("sample", "GENE"))

compare_keyset("clean_data: abri_kraken2Bracken_merged",
               file.path(root, "test_run/Datasets/abri_kraken2Bracken_merged.csv"),
               file.path(root, "Datasets/abri_kraken2Bracken_merged.csv"),
               c("sample", "GENE"))

# resistome — normalised gene table (the well-validated case)
compare_keyset("resistome: genetable_normdata",
               file.path(root, "test_run/Datasets/genetable_normdata.csv"),
               file.path(root, "Datasets/genetable_normdata.csv"),
               c("sample", "GENE"))

# network — gephi exports. The GT has multiple snapshots (gephi_nodes2/3, gephi_edges2/3);
# compare against the latest-suffixed one as the most-evolved baseline.
for (suffix in c("2", "3")) {
  compare_keyset(sprintf("network: gephi_nodes (vs gephi_nodes%s)", suffix),
                 file.path(root, "test_run/Datasets/network/gephi_nodes.csv"),
                 file.path(root, sprintf("Datasets/gephi_nodes%s.csv", suffix)),
                 c("Id"))
  compare_keyset(sprintf("network: gephi_edges (vs gephi_edges%s)", suffix),
                 file.path(root, "test_run/Datasets/network/gephi_edges.csv"),
                 file.path(root, sprintf("Datasets/gephi_edges%s.csv", suffix)),
                 c("Source", "Target"))
}

# alpha_diversity — schema differs but row count + sample overlap is meaningful.
compare_keyset("alpha_diversity: per-sample table",
               file.path(root, "test_run/Datasets/alpha_diversity.csv"),
               file.path(root, "Datasets/alphaDiversity_DataFile2.csv"),
               c("sample"))
