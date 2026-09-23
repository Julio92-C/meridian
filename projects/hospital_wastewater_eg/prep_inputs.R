#!/usr/bin/env Rscript
# prep_inputs.R — Materialise pipeline-ready inputs for the hospital_wastewater_eg
# validation study from the upstream files deposited under metaOmics/. Run ONCE
# before invoking run_pipeline.R against projects/hospital_wastewater_eg/config.yaml.
#
# Why a pre-step is necessary (this public-data cohort deviates from baseline):
#   1. Sample names in every upstream output are SRA aliases (1SW..10SW, 1P..10P),
#      not the human-readable metadata `sample` ids (HWW1..HWW10, TAP1..TAP10).
#      We reconcile alias -> sample via the metadata `sra_alias` column so all
#      emitted files are keyed by HWW/TAP and downstream metadata joins succeed.
#   2. Kraken2 (db1) and Bracken (db2) arrive as multi-sample *combined reports*
#      with a positional `#S<N>\t<alias>_se_dbX...` preamble and `<N>_all/<N>_lvl`
#      columns. Kraken's and Bracken's position->alias orderings DIFFER, so each
#      file is mapped via its OWN preamble.
#        - Kraken2: R/02 uses it for taxid->name AND build_taxid_ancestry(), which
#          needs the FULL report in NCBI DFS pre-order with `lvl_type`. So we only
#          strip the preamble (keep every data row + column) — never collapse it.
#        - Bracken: R/02 hard-requires the arranged wide CSV. We drop `_lvl`
#          columns and rename `<N>_all` -> sample to produce bracken_arranged.csv.
#   3. No Re-centrifuge .rcf.data.csv exists. R/01 hard-stops without one, so we
#      synthesise it in R/02's expected 3-header-row (cnt/una/sco triplet) shape
#      from the Kraken2 db1 `_all` counts. This study has NO negative controls
#      (config controls: []), so R/02's control_max = -Inf and nothing is
#      subtracted — the rcf is purely the per-sample count matrix.
#   4. Abricate summary_reportClean.csv `#FILE` prefixes are aliases
#      (`1SW_kraken2ID_PFP.fasta`). R/02 derives sample = sub("_.*","") = "1SW",
#      which would never match HWW1. We re-key each `#FILE` prefix alias->sample
#      and write a NEW summary_reportReconciled.csv (raw file left untouched).
#
# Outputs (created/overwritten):
#   metaOmics/Reports/Kraken2/kraken2_db1_combined_reports.txt   (preamble-stripped)
#   metaOmics/Reports/Bracken/bracken_arranged.csv
#   metaOmics/Reports/Re-centrifuge/hww_eg.rcf.data.csv
#   metaOmics/Reports/Abricate/summary_reportReconciled.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

ROOT <- "C:/Users/julio/Desktop/Pipeline_Validation/hospital_wastewater_eg/metaOmics"
REP  <- file.path(ROOT, "Reports")
MET  <- file.path(ROOT, "Metadata")

cat("=== prep_inputs.R (hospital_wastewater_eg) ===\n")

# === 0. Metadata alias -> sample map (authoritative) =======================
meta <- read_csv(file.path(MET, "hospital_wastewater_eg_metadata.csv"),
                 show_col_types = FALSE, progress = FALSE)
stopifnot(all(c("sample", "sra_alias") %in% colnames(meta)))
alias2sample <- setNames(as.character(meta$sample), as.character(meta$sra_alias))
cat(sprintf("[map] %d alias->sample pairs (e.g. %s -> %s)\n",
            length(alias2sample), names(alias2sample)[1], alias2sample[[1]]))

# Parse a kraken2/bracken multi-sample combined report. Returns the data frame
# (from the `#perc` header onward) plus the positional index -> sample map built
# from THIS file's own `#S<N>` preamble, and the raw lines (for verbatim re-emit).
read_combined <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pre   <- grep("^#S\\d+\\t", lines, value = TRUE)
  idx   <- as.integer(sub("^#S(\\d+)\\t.*$", "\\1", pre))
  raw   <- sub("^#S\\d+\\t", "", pre)
  alias <- sub("_se_db\\d.*$", "", raw)             # 8SW_se_db1.kraken2... -> 8SW
  samp  <- unname(alias2sample[alias])
  if (any(is.na(samp))) {
    stop(sprintf("unmapped aliases in %s: %s", basename(path),
                 paste(alias[is.na(samp)], collapse = ", ")))
  }
  pos_to_sample <- setNames(samp, as.character(idx))

  hdr_idx <- which(grepl("^#perc\\ttot_all", lines))[1]
  stopifnot("could not find data header (#perc\\ttot_all)" = !is.na(hdr_idx))
  body <- paste(lines[hdr_idx:length(lines)], collapse = "\n")
  d <- read_delim(I(body), delim = "\t", show_col_types = FALSE, progress = FALSE)

  list(d = d, pos_to_sample = pos_to_sample,
       header_lines = lines[hdr_idx:length(lines)])
}

# === 1. Kraken2 db1: strip preamble, keep full DFS report ==================
# R/02 select(taxid,name) + build_taxid_ancestry() need the whole report in
# order with lvl_type. We re-emit from the `#perc` header onward verbatim.
k2 <- read_combined(file.path(REP, "Kraken", "kraken2_db1_combined_reports.txt"))
dir.create(file.path(REP, "Kraken2"), recursive = TRUE, showWarnings = FALSE)
writeLines(k2$header_lines,
           file.path(REP, "Kraken2", "kraken2_db1_combined_reports.txt"))
cat(sprintf("[out] Reports/Kraken2/kraken2_db1_combined_reports.txt  (%d rows, preamble stripped)\n",
            length(k2$header_lines) - 1L))

# Kraken2 wide `_all` counts keyed by sample — feeds the synthetic rcf.
rename_all_cols <- function(d, pos_to_sample) {
  d <- d[, !grepl("^\\d+_lvl$", colnames(d)), drop = FALSE]  # drop per-SAMPLE lvl cols (keep tot_lvl)
  nm <- colnames(d)
  hit <- grepl("^\\d+_all$", nm)
  nm[hit] <- unname(pos_to_sample[sub("_all$", "", nm[hit])])
  colnames(d) <- nm
  d
}
k2_wide <- rename_all_cols(k2$d, k2$pos_to_sample)
sample_cols <- unname(alias2sample)                        # HWW1..TAP10, canonical order
stopifnot("kraken samples missing after rename" =
            all(sample_cols %in% colnames(k2_wide)))
cat(sprintf("[k2] %d taxa x %d samples\n", nrow(k2_wide), length(sample_cols)))

# === 2. Bracken db2: arrange into wide CSV keyed by sample =================
brk <- read_combined(file.path(REP, "Bracken",
                               "kraken2_db2-bracken_combined_reports.txt"))
brk_wide <- rename_all_cols(brk$d, brk$pos_to_sample)
stopifnot("bracken samples missing after rename" =
            all(sample_cols %in% colnames(brk_wide)))
# Column order mirrors the chicken_batch1 baseline: meta, samples, taxid, name.
brk_out <- brk_wide[, c("#perc", "tot_all", "tot_lvl",
                        sample_cols, "taxid", "name"), drop = FALSE]
write_csv(brk_out, file.path(REP, "Bracken", "bracken_arranged.csv"))
cat(sprintf("[out] Reports/Bracken/bracken_arranged.csv  (%d taxa, %d samples)\n",
            nrow(brk_out), length(sample_cols)))

# === 3. Synthetic Re-centrifuge .rcf.data.csv ==============================
# R/02 layout: row1 = Samples,Kraken_outputs/<s>,<s>,<s>,...; row2 = Stats,cnt,una,sco,...
# row3 = Id,,,...; then taxid + per-sample [cnt, una=0, sco=0] triplets. R/02
# keeps cols c(1,2,seq(5,.,3)) = taxid + each sample's cnt. Counts = kraken _all.
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
                   file.path(REP, "Re-centrifuge", "hww_eg.rcf.data.csv"))
cat(sprintf("[out] Reports/Re-centrifuge/hww_eg.rcf.data.csv  (%d taxa, %d samples, no controls)\n",
            nrow(k2_wide), length(sample_cols)))

# === 4. Abricate: re-key `#FILE` prefix alias -> sample ====================
abri <- read_csv(file.path(REP, "Abricate", "summary_reportClean.csv"),
                 show_col_types = FALSE, progress = FALSE)
file_alias <- sub("_.*", "", abri$`#FILE`)                # 1SW_kraken2ID_PFP.fasta -> 1SW
new_sample <- unname(alias2sample[file_alias])
if (any(is.na(new_sample))) {
  stop(sprintf("unmapped Abricate aliases: %s",
               paste(unique(file_alias[is.na(new_sample)]), collapse = ", ")))
}
abri$`#FILE` <- sub("^[^_]+", "", abri$`#FILE`)           # drop old alias prefix
abri$`#FILE` <- paste0(new_sample, abri$`#FILE`)          # prepend sample id
write_csv(abri, file.path(REP, "Abricate", "summary_reportReconciled.csv"))
cat(sprintf("[out] Reports/Abricate/summary_reportReconciled.csv  (%d rows, %d samples)\n",
            nrow(abri), length(unique(new_sample))))

cat("=== Done. Run: Rscript run_pipeline.R projects/hospital_wastewater_eg/config.yaml ===\n")
