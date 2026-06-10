#!/usr/bin/env Rscript
# prep_inputs.R — Materialise pipeline-ready inputs for the Hospital_microbiome
# study from the four upstream files the user placed under metaOmics/. Run this
# ONCE before invoking run_pipeline.R against projects/hospital_microbiome/config.yaml.
#
# Why a pre-step is necessary (study deviates from the chicken_batch1 baseline):
#   1. Two Kraken2 multi-sample reports (H1A = 52 Air, H1S = 36 Surface) with
#      positional column headers (1_all, 2_all, ...). The original sample names
#      live in the file preamble (#S<N>\t<name>_se_db1.kraken2...). We strip
#      the preamble, rename columns to real sample IDs, then concatenate Air +
#      Surface into one combined report.
#   2. No Bracken output exists. R/03_normalisation.R hard-requires the
#      abri_kraken2Bracken_merged table, so we synthesise a Bracken-shaped file
#      from the merged Kraken2 counts (per-sample columns + taxid + name +
#      tot_all + #perc). Substitutes Bracken counts with Kraken2 counts.
#   3. No standard Re-centrifuge .rcf.data.csv exists. Instead the project has
#      a flat NX_AH4_contam.csv (taxid + NX, N3 controls + 11 Air samples). We
#      build a synthetic rcf in R/02's expected 3-header-row + triplet-column
#      shape, using the merged Kraken2 counts as the per-sample values and the
#      NX/N3 contam columns as the negative-control baseline R/02 will subtract.
#   4. Abricate Summary.csv is in pre-cleanup form. The Air A37R replicate is
#      collapsed into A37 (sum) and A25R is dropped, matching the original
#      abri_kraken2_cleanData.R script's sample reconciliation.
#   5. Metadata MetadataLocations2.csv (64 samples) is extended with two
#      Control rows (NX, N3) so downstream stages that iterate metadata see
#      every sample column present in the rcf.
#
# Outputs (overwrites any existing files at these paths):
#   metaOmics/Reports/Kraken2/kraken2_db1_combined_reports.txt
#   metaOmics/Reports/Bracken/bracken_arranged.csv
#   metaOmics/Reports/Re-centrifuge/H1AS_NX_N3.rcf.data.csv
#   metaOmics/Reports/Abricate/summary_reportClean.csv
#   metaOmics/Metadata/hospital_metadata.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

ROOT <- "C:/Users/julio/Desktop/Hospital_microbiome/metaOmics"
DAT  <- file.path(ROOT, "Datasets")
REP  <- file.path(ROOT, "Reports")
MET  <- file.path(ROOT, "Metadata")

CONTROLS <- c("NX", "N3")

cat("=== prep_inputs.R (hospital_microbiome) ===\n")

# === 1. Read both Kraken2 reports, rename positional cols to sample IDs ====

read_kraken2_named <- function(path, batch_label) {
  cat(sprintf("[k2] reading %s\n", basename(path)))
  lines <- readLines(path, warn = FALSE)
  # Preamble: "#S<N>\t<sample>_se_db1.kraken2.kraken2.report.txt"
  pre_lines <- grep("^#S\\d+\\t", lines, value = TRUE)
  idx_pos   <- as.integer(sub("^#S(\\d+)\\t.*$", "\\1", pre_lines))
  raw_name  <- sub("^#S\\d+\\t", "", pre_lines)
  samp_name <- sub("_se_db1\\.kraken2\\.kraken2\\.report\\.txt$", "", raw_name)
  pos_to_sample <- setNames(samp_name, as.character(idx_pos))
  cat(sprintf("  [%s] %d positional sample names parsed\n",
              batch_label, length(pos_to_sample)))

  hdr_idx <- which(grepl("^#perc\\ttot_all", lines))[1]
  stopifnot("could not find data header" = !is.na(hdr_idx))
  body_tsv <- paste(lines[hdr_idx:length(lines)], collapse = "\n")
  d <- read_delim(I(body_tsv), delim = "\t",
                  show_col_types = FALSE, progress = FALSE)

  # Drop *_lvl (per-level counts; we only want *_all = cumulative counts)
  d <- d[, !grepl("_lvl$", colnames(d)), drop = FALSE]
  # Rename "<N>_all" -> pos_to_sample[<N>]
  nm <- colnames(d)
  rename_idx <- grepl("^\\d+_all$", nm)
  pos_keys   <- sub("_all$", "", nm[rename_idx])
  nm[rename_idx] <- unname(pos_to_sample[pos_keys])
  colnames(d) <- nm
  d
}

h1a <- read_kraken2_named(file.path(DAT, "H1A_kraken2_db1_combined_reports.txt"),
                          "Air")
h1s <- read_kraken2_named(file.path(DAT, "H1S_kraken2_db1_combined_reports.txt"),
                          "Surface")

# === 2. Long-format counts, sample-name reconciliation ====================

META_COLS <- c("#perc", "tot_all", "tot_lvl", "lvl_type", "taxid", "name")

to_long_k2 <- function(d) {
  scols <- setdiff(colnames(d), META_COLS)
  d |>
    select(taxid, name, all_of(scols)) |>
    mutate(across(all_of(scols), as.character)) |>   # readr may infer mixed double/character across sample cols; unify before pivot
    pivot_longer(all_of(scols), names_to = "sample", values_to = "count") |>
    mutate(count = suppressWarnings(as.numeric(count)),
           taxid = as.character(taxid))
}

k2_long <- bind_rows(to_long_k2(h1a), to_long_k2(h1s)) |>
  mutate(count = tidyr::replace_na(count, 0))

# Sample reconciliation: A37R -> A37 (sum); drop A25R; A60B-A63B not present
k2_long <- k2_long |>
  filter(sample != "A25R") |>
  mutate(sample = ifelse(sample == "A37R", "A37", sample)) |>
  group_by(taxid, name, sample) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

cat(sprintf("[k2] %d taxa x %d samples (combined Air+Surface, reconciled)\n",
            length(unique(k2_long$taxid)), length(unique(k2_long$sample))))

# Wide table once — reused for both bracken and rcf synthesis
k2_wide <- pivot_wider(k2_long, names_from = sample,
                       values_from = count, values_fill = 0)
sample_cols <- setdiff(colnames(k2_wide), c("taxid", "name"))

# === 3. Combined Kraken2 (minimal taxid->name lookup) ======================
# R/02 only uses inputs$kraken2 for the taxid->name join (line 29 in 02_clean_data.R),
# so we just need a 2-column TSV. Keep #perc + tot_all + tot_lvl as placeholder
# zeros so the file looks like a standard kraken2 multi-report (defensive against
# any future stage that grepls headers).

k2_lookup <- k2_wide |>
  mutate(`#perc` = 0, tot_all = 0, tot_lvl = 0, lvl_type = "") |>
  select(`#perc`, tot_all, tot_lvl, lvl_type, taxid, name) |>
  distinct(taxid, .keep_all = TRUE)

write_tsv(k2_lookup,
          file.path(REP, "Kraken2", "kraken2_db1_combined_reports.txt"))
cat(sprintf("[out] Reports/Kraken2/kraken2_db1_combined_reports.txt  (%d taxids)\n",
            nrow(k2_lookup)))

# === 4. Synthetic Bracken (per-sample counts wide; R/03 needs this) ========
# R/02 (lines 142-179) pivots bracken on sample columns and joins on
# (taxid, sample, name). Optional rename: TotalBCount = tot_all, PercB = #perc.

brk_out <- k2_wide |>
  mutate(tot_all = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) |>
  mutate(`#perc` = round(100 * tot_all / sum(tot_all), 4))
brk_out <- brk_out[, c("#perc", "tot_all", sample_cols, "taxid", "name")]
write_csv(brk_out, file.path(REP, "Bracken", "bracken_arranged.csv"))
cat(sprintf("[out] Reports/Bracken/bracken_arranged.csv  (%d taxa, %d samples)\n",
            nrow(brk_out), length(sample_cols)))

# === 5. Synthetic Re-centrifuge .rcf.data.csv ==============================
# Reshape NX_AH4_contam.csv (taxid + NX + N3 + 11 Air samples) into R/02's
# expected layout: 3 header rows + per-sample triplet columns (cnt/una/sco).
# We inject NX/N3 control columns (from contam) alongside ALL Kraken2 sample
# columns (real counts), so R/02 can subtract control_max = max(NX, N3) per
# row from each sample. Samples not in metadata get filtered downstream.

contam <- read_csv(file.path(DAT, "NX_AH4_contam.csv"),
                   show_col_types = FALSE, progress = FALSE) |>
  mutate(taxid = as.character(taxid))

# Pull the NX, N3 control vectors (named by taxid) for fast lookup.
contam_lookup <- function(ctrl_name) {
  v <- setNames(as.numeric(contam[[ctrl_name]]), contam$taxid)
  function(tids) {
    out <- unname(v[as.character(tids)])
    ifelse(is.na(out), 0, out)
  }
}
lookup_NX <- contam_lookup("NX")
lookup_N3 <- contam_lookup("N3")

rcf_wide <- k2_wide
rcf_wide$NX <- lookup_NX(rcf_wide$taxid)
rcf_wide$N3 <- lookup_N3(rcf_wide$taxid)
rcf_wide$name <- NULL   # rcf has no name column; R/02 left-joins it via taxid_name

# Column order: taxid, then all sample cols + controls last (R/02 doesn't care)
rcf_sample_cols <- c(sample_cols, CONTROLS)
rcf_wide <- rcf_wide[, c("taxid", rcf_sample_cols), drop = FALSE]

write_combined_rcf <- function(wide, out_path) {
  samples <- setdiff(colnames(wide), "taxid")
  N <- length(samples)
  header_row <- c("Samples",
                  as.vector(rbind(
                    paste0("Kraken_outputs/", samples),
                    paste0("Kraken_outputs/", samples),
                    paste0("Kraken_outputs/", samples)
                  )))
  stats_row  <- c("Stats", rep(c("cnt", "una", "sco"), N))
  id_row     <- c("Id",    rep("", 3L * N))

  data_mat <- matrix("", nrow = nrow(wide), ncol = 1L + 3L * N)
  data_mat[, 1L] <- as.character(wide$taxid)
  for (i in seq_along(samples)) {
    data_mat[, 1L + 3L*(i-1L) + 1L] <- as.character(wide[[samples[i]]])
    data_mat[, 1L + 3L*(i-1L) + 2L] <- "0"
    data_mat[, 1L + 3L*(i-1L) + 3L] <- "0"
  }

  con <- file(out_path, "w")
  on.exit(close(con))
  writeLines(paste(header_row, collapse = ","), con)
  writeLines(paste(stats_row,  collapse = ","), con)
  writeLines(paste(id_row,     collapse = ","), con)
  write.table(data_mat, con, sep = ",", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}

write_combined_rcf(rcf_wide,
                   file.path(REP, "Re-centrifuge", "H1AS_NX_N3.rcf.data.csv"))
cat(sprintf("[out] Reports/Re-centrifuge/H1AS_NX_N3.rcf.data.csv  (%d taxa, %d samples + 2 controls)\n",
            nrow(rcf_wide), length(sample_cols)))

# === 6. Cleaned Abricate Summary ===========================================
# Source has sample = "<id>.fasta". Drop A25R rows; rename A37R -> A37.
# A60B-A63B aren't present anyway. R/02 will further normalise #FILE via its
# own gsub/sub on lines 18-19.

abri <- read_csv(file.path(DAT, "Summary.csv"),
                 show_col_types = FALSE, progress = FALSE)
abri_sample <- sub("\\.fasta$", "", abri$`#FILE`)
abri <- abri[abri_sample != "A25R", , drop = FALSE]
abri$`#FILE` <- ifelse(sub("\\.fasta$", "", abri$`#FILE`) == "A37R",
                       "A37.fasta", abri$`#FILE`)
write_csv(abri, file.path(REP, "Abricate", "summary_reportClean.csv"))
cat(sprintf("[out] Reports/Abricate/summary_reportClean.csv  (%d rows)\n",
            nrow(abri)))

# === 7. Metadata with NX/N3 control rows ===================================

meta <- read_csv(file.path(MET, "MetadataLocations2.csv"),
                 show_col_types = FALSE, progress = FALSE)
ctrl_rows <- tibble(
  ID          = CONTROLS,
  sample      = CONTROLS,
  type        = "Control",
  description = c("Negative extraction control (NX)", "Negative control (N3)"),
  ward        = NA_character_,
  location    = NA_character_
)
meta_out <- bind_rows(meta, ctrl_rows)
write_csv(meta_out, file.path(MET, "hospital_metadata.csv"))
cat(sprintf("[out] Metadata/hospital_metadata.csv  (%d rows: %d samples + %d controls)\n",
            nrow(meta_out), nrow(meta), length(CONTROLS)))

cat("=== Done. Run pipeline with projects/hospital_microbiome/config.yaml ===\n")
