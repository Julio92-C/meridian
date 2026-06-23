#!/usr/bin/env Rscript
# Diff pipeline-produced genetable_normdata.csv against the chicken_batch1
# ground-truth file from PC_JC_2024-11-29_Julio_Gallus.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

root <- "C:/Users/julio/Desktop/PC_JC_2024-11-29_Julio_Gallus"
ours <- read_csv(file.path(root, "test_run/Datasets/genetable_normdata.csv"),
                 show_col_types = FALSE)
gt   <- read_csv(file.path(root, "Datasets/genetable_normdata.csv"),
                 show_col_types = FALSE)

cat(sprintf("ours: %d rows × %d cols\n", nrow(ours), ncol(ours)))
cat(sprintf("gt  : %d rows × %d cols\n", nrow(gt),   ncol(gt)))
cat("\nours cols:", paste(colnames(ours), collapse = ", "), "\n")
cat("gt   cols:", paste(colnames(gt),   collapse = ", "), "\n\n")

key <- c("sample", "GENE")
ours <- ours[do.call(order, ours[, key]), ]
gt   <- gt  [do.call(order, gt[,   key]), ]

cat("--- key alignment ---\n")
cat(sprintf("identical (sample,GENE) sets : %s\n",
            identical(as.character(ours$sample), as.character(gt$sample)) &&
              identical(as.character(ours$GENE), as.character(gt$GENE))))
cat(sprintf("ours unique pairs           : %d\n",
            nrow(dplyr::distinct(ours[, key]))))
cat(sprintf("gt   unique pairs           : %d\n",
            nrow(dplyr::distinct(gt[,   key]))))

cat("\n--- anti-join (ours vs gt) ---\n")
ours_keys <- ours[, key]
gt_keys   <- gt  [, key]
only_ours <- dplyr::anti_join(ours_keys, gt_keys, by = key)
only_gt   <- dplyr::anti_join(gt_keys,   ours_keys, by = key)
cat(sprintf("only in ours: %d (sample,GENE) pairs\n", nrow(only_ours)))
cat(sprintf("only in gt  : %d (sample,GENE) pairs\n", nrow(only_gt)))
cat("\nper-sample row counts (top 20 by absolute diff):\n")
ours_per <- as.data.frame(table(ours$sample))
gt_per   <- as.data.frame(table(gt$sample))
names(ours_per) <- c("sample", "ours_n")
names(gt_per)   <- c("sample", "gt_n")
per <- merge(ours_per, gt_per, by = "sample", all = TRUE)
per$diff <- per$ours_n - per$gt_n
per <- per[order(-abs(per$diff)), ]
print(head(per, 20), row.names = FALSE)
cat("\nfirst 15 (sample,GENE) only in ours:\n")
print(head(only_ours, 15), row.names = FALSE)
cat("\nfirst 15 (sample,GENE) only in gt:\n")
print(head(only_gt, 15), row.names = FALSE)

if (nrow(ours) == nrow(gt) &&
    identical(as.character(ours$sample), as.character(gt$sample)) &&
    identical(as.character(ours$GENE),   as.character(gt$GENE))) {
  cat("\n--- per-column diff (rows aligned) ---\n")
  for (cc in intersect(colnames(ours), colnames(gt))) {
    if (cc %in% key) next
    v1 <- ours[[cc]]; v2 <- gt[[cc]]
    if (is.numeric(v1) && is.numeric(v2)) {
      d <- max(abs(v1 - v2), na.rm = TRUE)
      cat(sprintf("  [num] %-12s max abs diff = %g\n", cc, d))
    } else {
      eq <- identical(as.character(v1), as.character(v2))
      cat(sprintf("  [chr] %-12s identical    = %s\n", cc, eq))
      if (!eq) {
        diffs <- which(as.character(v1) != as.character(v2))
        cat(sprintf("        first diff at row %d: ours=%s | gt=%s\n",
                    diffs[1], v1[diffs[1]], v2[diffs[1]]))
      }
    }
  }
}
