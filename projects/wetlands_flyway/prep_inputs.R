#!/usr/bin/env Rscript
# prep_inputs.R — Materialise pipeline-ready inputs for the wetlands_flyway
# validation study from the upstream files deposited under Reports/. Run ONCE
# before invoking run_pipeline.R against projects/wetlands_flyway/config.yaml.
#
# Why a pre-step is necessary (mirrors hospital_wastewater_eg, but SIMPLER —
# sample names already match the metadata, so no alias reconciliation):
#   1. Kraken2 (db1) and Bracken (db2) arrive as multi-sample *combined reports*
#      with a positional `#S<N>\t<sample>_se_dbX...` preamble and `<N>_all/<N>_lvl`
#      columns. The preamble sample names ARE the metadata `sample` ids
#      (GA1_1, FA2_1, ...), so we map position -> sample directly (no alias map).
#        - Kraken2: strip the preamble only (keep the full report in NCBI DFS
#          order with lvl_type — needed by R/02 build_taxid_ancestry()).
#        - Bracken: drop the per-sample `_lvl` columns and rename `<N>_all` ->
#          sample to produce the arranged wide bracken_arranged.csv R/02 needs.
#   2. No Re-centrifuge .rcf.data.csv exists. R/01 hard-stops without one, so we
#      synthesise it in R/02's expected 3-header-row (cnt/una/sco triplet) shape
#      from the Kraken2 db1 `_all` counts. No negative controls (config
#      controls: []) -> R/02's control_max = -Inf, nothing subtracted.
#   3. Abricate summary_reportClean.csv needs NO reshaping: `#FILE` is
#      "<sample>_kraken2ID_PFP.fasta" and R/02 now strips the _kraken2ID... tail
#      (preserving the underscore in GA1_1). It is read directly from the deposit.
#   4. The wf-metagenomics-diversity.csv is used directly (config diversity_csv);
#      its sample column headers already match the metadata.
#
# Outputs (created/overwritten):
#   Reports/Kraken2/kraken2_db1_combined_reports.txt   (preamble-stripped)
#   Reports/Bracken/bracken_arranged.csv
#   Reports/Re-centrifuge/wetlands_flyway.rcf.data.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT <- "C:/Users/julio/Desktop/Pipeline_Validation/wetlands_flyway"
REP  <- file.path(ROOT, "Reports")
MET  <- file.path(ROOT, "Metadata")

cat("=== prep_inputs.R (wetlands_flyway) ===\n")

meta <- read_csv(file.path(MET, "wetlands_flyway_metadata.csv"),
                 show_col_types = FALSE, progress = FALSE)
sample_cols <- as.character(meta$sample)               # GA1_1 .. SN2_2 (24)
cat(sprintf("[meta] %d samples (e.g. %s)\n", length(sample_cols), sample_cols[1]))

# Parse a kraken2/bracken multi-sample combined report. The `#S<N>` preamble
# maps positional index -> sample id directly (no alias indirection here).
read_combined <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pre   <- grep("^#S\\d+\\t", lines, value = TRUE)
  idx   <- as.integer(sub("^#S(\\d+)\\t.*$", "\\1", pre))
  samp  <- sub("_se_db\\d.*$", "", sub("^#S\\d+\\t", "", pre))   # GA1_1_se_db1... -> GA1_1
  miss  <- setdiff(samp, sample_cols)
  if (length(miss) > 0) {
    stop(sprintf("%s: preamble samples not in metadata: %s",
                 basename(path), paste(miss, collapse = ", ")))
  }
  pos_to_sample <- setNames(samp, as.character(idx))

  hdr_idx <- which(grepl("^#perc\\ttot_all", lines))[1]
  stopifnot("could not find data header (#perc\\ttot_all)" = !is.na(hdr_idx))
  body <- paste(lines[hdr_idx:length(lines)], collapse = "\n")
  d <- read_delim(I(body), delim = "\t", show_col_types = FALSE, progress = FALSE)
  list(d = d, pos_to_sample = pos_to_sample,
       header_lines = lines[hdr_idx:length(lines)])
}

rename_all_cols <- function(d, pos_to_sample) {
  d <- d[, !grepl("^\\d+_lvl$", colnames(d)), drop = FALSE]   # drop per-sample lvl (keep tot_lvl)
  nm <- colnames(d)
  hit <- grepl("^\\d+_all$", nm)
  nm[hit] <- unname(pos_to_sample[sub("_all$", "", nm[hit])])
  colnames(d) <- nm
  d
}

# === 1. Kraken2 db1: strip preamble, keep full DFS report ==================
k2 <- read_combined(file.path(REP, "Kraken", "kraken2_db1_combined_reports.txt"))
dir.create(file.path(REP, "Kraken2"), recursive = TRUE, showWarnings = FALSE)
writeLines(k2$header_lines,
           file.path(REP, "Kraken2", "kraken2_db1_combined_reports.txt"))
cat(sprintf("[out] Reports/Kraken2/kraken2_db1_combined_reports.txt  (%d rows, preamble stripped)\n",
            length(k2$header_lines) - 1L))

k2_wide <- rename_all_cols(k2$d, k2$pos_to_sample)
stopifnot("kraken samples missing after rename" =
            all(sample_cols %in% colnames(k2_wide)))
cat(sprintf("[k2] %d taxa x %d samples\n", nrow(k2_wide), length(sample_cols)))

# === 2. Bracken db2: arrange into wide CSV keyed by sample =================
brk <- read_combined(file.path(REP, "Bracken",
                               "kraken2_db2-bracken_combined_reports.txt"))
brk_wide <- rename_all_cols(brk$d, brk$pos_to_sample)
stopifnot("bracken samples missing after rename" =
            all(sample_cols %in% colnames(brk_wide)))
brk_out <- brk_wide[, c("#perc", "tot_all", "tot_lvl",
                        sample_cols, "taxid", "name"), drop = FALSE]
write_csv(brk_out, file.path(REP, "Bracken", "bracken_arranged.csv"))
cat(sprintf("[out] Reports/Bracken/bracken_arranged.csv  (%d taxa, %d samples)\n",
            nrow(brk_out), length(sample_cols)))

# === 3. Synthetic Re-centrifuge .rcf.data.csv ==============================
write_combined_rcf <- function(wide, samples, out_path) {
  N <- length(samples)
  header_row <- c("Samples", as.vector(rbind(
    paste0("Kraken_outputs/", samples),
    paste0("Kraken_outputs/", samples),
    paste0("Kraken_outputs/", samples))))
  stats_row <- c("Stats", rep(c("cnt", "una", "sco"), N))
  id_row    <- c("Id",    rep("", 3L * N))

  mat <- matrix("", nrow = nrow(wide), ncol = 1L + 3L * N)
  mat[, 1L] <- as.character(wide$taxid)
  for (i in seq_along(samples)) {
    mat[, 1L + 3L * (i - 1L) + 1L] <- as.character(wide[[samples[i]]])
    mat[, 1L + 3L * (i - 1L) + 2L] <- "0"
    mat[, 1L + 3L * (i - 1L) + 3L] <- "0"
  }
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  con <- file(out_path, "w"); on.exit(close(con))
  writeLines(paste(header_row, collapse = ","), con)
  writeLines(paste(stats_row,  collapse = ","), con)
  writeLines(paste(id_row,     collapse = ","), con)
  write.table(mat, con, sep = ",", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}
write_combined_rcf(k2_wide, sample_cols,
                   file.path(REP, "Re-centrifuge", "wetlands_flyway.rcf.data.csv"))
cat(sprintf("[out] Reports/Re-centrifuge/wetlands_flyway.rcf.data.csv  (%d taxa, %d samples, no controls)\n",
            nrow(k2_wide), length(sample_cols)))

cat("=== Done. Run: Rscript run_pipeline.R projects/wetlands_flyway/config.yaml ===\n")
