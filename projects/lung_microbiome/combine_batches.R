#!/usr/bin/env Rscript
# combine_batches.R — Merge the 3 ONT sequencing batches (April / May / June
# 2025) into single pipeline-ready inputs under metaOmics/Reports/. Run this
# once before invoking run_pipeline.R against projects/lung_microbiome/config.yaml.
#
# Per-batch sample IDs are reconciled via:
#   - April rcf: CTRL1 -> NC-19
#   - May rcf:   CTRL1 -> NC-10, CTRL2 -> NC-13
#   - June rcf/abricate/diversity: mapping CSV at Metadata/june_sample_mapping.csv
#     (wf-metagenomics IDs -> metadata sample IDs, derived from sequence_info.xlsx
#     via patient-ID join with sample_metadata.csv).
#
# Overlapping samples (same metadata ID appearing in >1 batch) get their counts
# summed across batches — user-confirmed: batches are treated as independent
# re-preps so additive combination is valid. Affected samples:
#   SP-16, SP-15, SP-03, SP-10, EB-23 (April + May)
#   EB-17                            (April + June)
#
# Outputs (overwrites any existing files at these paths):
#   metaOmics/Reports/Abricate/summary_reportClean.csv
#   metaOmics/Reports/Bracken/bracken_arranged.csv
#   metaOmics/Reports/Kraken/kraken2_db1_combined_reports.txt
#   metaOmics/Reports/Re-centrifuge/EBSP_samples.rcf.data.csv
#   metaOmics/Reports/wf-metagenomics/wf-metagenomics-diversity.csv
#     (NA for the 6 overlap samples — precomputed metrics cannot be summed)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

ROOT <- "C:/Users/julio/Desktop/Lung_microbiome/metaOmics"
REP  <- file.path(ROOT, "Reports")

BATCHES <- list(
  apr = "JC_PC_2025-04-11_EBSP_UCLH",
  may = "JC_PC_2025-05-13_EBSP_UCLH",
  jun = "JC_2025-06-18_EBSP_UCLH"
)

# Per-batch rcf control rename (CTRL labels -> real NC-* IDs)
CTRL_RENAME <- list(
  apr = c("CTRL1" = "NC-19"),
  may = c("CTRL1" = "NC-10", "CTRL2" = "NC-13"),
  jun = c("CTRL1" = "NC-22", "CTRL2" = "NC-30")
)

# Per-batch negative-control sample IDs (post-rename). Used for per-batch
# baseline subtraction (equivalent of R/02's pooled control_max filter, but
# scoped to the batch a sample was sequenced in — see project memory for why).
BATCH_CONTROLS <- list(
  apr = c("NC-19"),
  may = c("NC-10", "NC-13"),
  jun = c("NC-22", "NC-30")
)
ALL_CONTROLS <- unique(unlist(BATCH_CONTROLS))

# June wf-metagenomics naming -> metadata naming
june_map_df <- read.csv(file.path(ROOT, "Metadata/june_sample_mapping.csv"),
                        stringsAsFactors = FALSE)
JUN_MAP <- setNames(june_map_df$metadata_id, june_map_df$wf_id)
JUN_MAP <- JUN_MAP[nchar(JUN_MAP) > 0 & !is.na(JUN_MAP)]

apply_map <- function(x, m) {
  ifelse(x %in% names(m), unname(m[x]), x)
}

# Per-batch contaminant filter: for each taxid, take max count across THIS
# batch's controls; zero out any (taxid, sample) cell whose count is not
# strictly greater than that max; then drop the control columns. Mirrors the
# semantics of R/02's filter but scoped to a single batch so April's NC-19
# doesn't suppress reads in June samples (and vice versa).
subtract_controls_rcf <- function(rcf_wide, controls, batch_name) {
  ctrl_present <- intersect(controls, colnames(rcf_wide))
  ctrl_missing <- setdiff(controls, colnames(rcf_wide))
  if (length(ctrl_missing) > 0) {
    warning(sprintf("[%s] declared controls not found in rcf columns: %s",
                    batch_name, paste(ctrl_missing, collapse = ", ")))
  }
  if (length(ctrl_present) == 0) {
    cat(sprintf("  [%s] no controls found in rcf -- subtraction skipped\n",
                batch_name))
    return(rcf_wide)
  }
  ctrl_mat <- as.matrix(rcf_wide[, ctrl_present, drop = FALSE])
  ctrl_mat[is.na(ctrl_mat)] <- 0
  ctrl_max <- apply(ctrl_mat, 1L, max)
  sample_cols <- setdiff(colnames(rcf_wide), c("taxid", ctrl_present))
  n_zeroed <- 0L; total_before <- 0; total_after <- 0
  for (s in sample_cols) {
    vals <- rcf_wide[[s]]
    drops <- !is.na(vals) & vals <= ctrl_max
    n_zeroed <- n_zeroed + sum(drops)
    total_before <- total_before + sum(vals, na.rm = TRUE)
    rcf_wide[[s]][drops] <- 0
    total_after <- total_after + sum(rcf_wide[[s]], na.rm = TRUE)
  }
  cat(sprintf("  [%s] controls {%s}: zeroed %d (taxon,sample) cells; reads %.0f -> %.0f (-%.1f%%)\n",
              batch_name, paste(ctrl_present, collapse = ","),
              n_zeroed, total_before, total_after,
              100 * (total_before - total_after) / max(total_before, 1)))
  rcf_wide[, c("taxid", sample_cols), drop = FALSE]
}

cat("=== combine_batches.R ===\n")
cat("Loaded", length(JUN_MAP), "June wf->metadata mappings\n\n")

# === 1. ABRICATE ===========================================================
cat("[1/5] Combining Abricate summaries...\n")

read_abri <- function(batch_suffix, rename_map = NULL) {
  f <- file.path(REP, "Abricate",
                 paste0("summary_reportClean_", batch_suffix, ".csv"))
  d <- read_csv(f, show_col_types = FALSE, progress = FALSE)
  if (!is.null(rename_map)) {
    base <- sub("_.*", "", sub("\\.fasta$", "", d$`#FILE`))
    new_base <- apply_map(base, rename_map)
    # Drop rows whose new_base is "" (orphans not in mapping)
    keep <- nchar(new_base) > 0
    d <- d[keep, , drop = FALSE]
    new_base <- new_base[keep]
    d$`#FILE` <- sub("^[^_]+", "", d$`#FILE`)  # strip old prefix
    d$`#FILE` <- paste0(new_base, d$`#FILE`)
  }
  d
}

abri <- bind_rows(
  read_abri(BATCHES$apr),
  read_abri(BATCHES$may),
  read_abri(BATCHES$jun, rename_map = JUN_MAP)
)
# Drop control rows: cells from NC-* samples (mostly empty anyway, since
# Abricate runs per-sample on contigs). They'd be dropped at R/02's
# inner_join with noncontaminants too, but stripping here keeps the
# combined file consistent with rcf/bracken.
abri_sample <- sub("_.*", "", sub("\\.fasta$", "", abri$`#FILE`))
abri <- abri[!abri_sample %in% ALL_CONTROLS, , drop = FALSE]
write_csv(abri, file.path(REP, "Abricate", "summary_reportClean.csv"))
cat("  written:", nrow(abri), "rows (controls stripped)\n\n")

# === 2. BRACKEN ============================================================
cat("[2/5] Combining Bracken arranged tables...\n")

read_brk <- function(batch_suffix) {
  f <- file.path(REP, "Bracken",
                 paste0("bracken_arranged_", batch_suffix, ".csv"))
  read_csv(f, show_col_types = FALSE, progress = FALSE)
}

brk_apr <- read_brk(BATCHES$apr)
brk_may <- read_brk(BATCHES$may)
brk_jun <- read_brk(BATCHES$jun)

# NOTE: June bracken is already correctly labeled — the user pre-renamed
# wf-EB-04 -> EB-16 (patient 1604) and wf-EB-16 -> EB-38 (patient 1536)
# during manual relabeling, so both columns hold real, distinct samples.
# Verified by raw-data inspection: EB-16 has 21339 reads / 516 nonzero taxa,
# EB-38 has 34501 reads / 521 nonzero taxa. No rename or drop is applied.

# Drop each batch's controls so the combined bracken contains only the 46
# analysis samples. Bracken doesn't carry a contaminant-floor concept of its
# own (R/02's filter is rcf-driven), so we just remove the columns. The
# downstream R/02 merge would inner-join them away anyway, but stripping
# keeps the combined file self-consistent and prevents NC-* rows from being
# pivoted in stages that pivot bracken before the join.
drop_ctrl_cols <- function(brk, controls, label) {
  to_drop <- intersect(controls, colnames(brk))
  if (length(to_drop) > 0) {
    cat(sprintf("  [%s] dropping bracken control cols: %s\n",
                label, paste(to_drop, collapse = ",")))
    brk[, !colnames(brk) %in% to_drop, drop = FALSE]
  } else {
    brk
  }
}
brk_apr <- drop_ctrl_cols(brk_apr, BATCH_CONTROLS$apr, "apr")
brk_may <- drop_ctrl_cols(brk_may, BATCH_CONTROLS$may, "may")
brk_jun <- drop_ctrl_cols(brk_jun, BATCH_CONTROLS$jun, "jun")

META_COLS <- c("#perc", "tot_all", "tot_lvl", "taxid", "name")

to_long_brk <- function(d) {
  scols <- setdiff(colnames(d), META_COLS)
  d |>
    select(taxid, name, all_of(scols)) |>
    pivot_longer(all_of(scols), names_to = "sample", values_to = "count") |>
    mutate(count = suppressWarnings(as.numeric(count)))
}

brk_long <- bind_rows(
  to_long_brk(brk_apr),
  to_long_brk(brk_may),
  to_long_brk(brk_jun)
) |>
  mutate(count = tidyr::replace_na(count, 0))

brk_summed <- brk_long |>
  group_by(taxid, name, sample) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

brk_wide <- pivot_wider(brk_summed, names_from = sample,
                        values_from = count, values_fill = 0)

sample_cols <- setdiff(colnames(brk_wide), c("taxid", "name"))
tot_all_v <- rowSums(brk_wide[, sample_cols, drop = FALSE], na.rm = TRUE)
brk_out <- tibble(
  `#perc`  = round(100 * tot_all_v / sum(tot_all_v), 4),
  tot_all  = tot_all_v
)
brk_out <- bind_cols(brk_out, brk_wide[, sample_cols], brk_wide[, c("taxid", "name")])
write_csv(brk_out, file.path(REP, "Bracken", "bracken_arranged.csv"))
cat("  written:", nrow(brk_out), "rows,", length(sample_cols), "sample cols\n\n")

# === 3. KRAKEN2 (taxid -> name lookup) =====================================
cat("[3/5] Combining Kraken2 taxid->name lookup...\n")

read_k2 <- function(batch_suffix) {
  f <- file.path(REP, "Kraken",
                 paste0("kraken2_db1_combined_reports_", batch_suffix, ".txt"))
  read_delim(f, delim = "\t", show_col_types = FALSE, progress = FALSE)
}

k2_all <- bind_rows(
  read_k2(BATCHES$apr) |> select(taxid, name),
  read_k2(BATCHES$may) |> select(taxid, name),
  read_k2(BATCHES$jun) |> select(taxid, name)
) |> distinct()

write_tsv(k2_all, file.path(REP, "Kraken", "kraken2_db1_combined_reports.txt"))
cat("  written:", nrow(k2_all), "(taxid, name) pairs\n\n")

# === 4. RE-CENTRIFUGE ======================================================
cat("[4/5] Combining Re-centrifuge data...\n")

# Read one batch rcf, strip fake-sample columns, return clean (taxid, samples...)
# data.frame. Mirrors R/02_clean_data.R's parsing logic.
read_rcf_clean <- function(batch_suffix, rename_map = NULL) {
  f <- file.path(REP, "Re-centrifuge",
                 paste0("recentrifuge_", batch_suffix, ".rcf.data.csv"))
  rcf <- fread(f, header = TRUE, colClasses = "character", data.table = FALSE,
               showProgress = FALSE)
  # Drop the Stats and Id data rows (data row 1 and 2)
  rcf <- rcf[-c(1, 2), , drop = FALSE]
  # Drop trailing Details columns (Rank, Name)
  rcf <- rcf[, !grepl("^Details($|\\.)", colnames(rcf)), drop = FALSE]
  # Keep col 1 (Samples = taxid) + every 3rd column starting at col 2 (cnt cols)
  keep_cols <- c(1, 2, seq(5, ncol(rcf), by = 3))
  keep_cols <- keep_cols[keep_cols <= ncol(rcf)]
  rcf <- rcf[, keep_cols, drop = FALSE]
  # Clean column names: strip Kraken_outputs/ prefix and disambig suffix
  colnames(rcf) <- vapply(colnames(rcf), function(n) {
    n <- gsub(".*?/", "", n)
    gsub("\\..*", "", n)
  }, character(1))
  names(rcf)[1] <- "taxid"
  # Drop pseudo-sample columns (EXCLUSIVE_*, CTRL_*, SHARED_* etc.)
  is_pseudo <- grepl("(_EXCLUSIVE|_CTRL|SHARED)", colnames(rcf)) &
               colnames(rcf) != "taxid"
  rcf <- rcf[, !is_pseudo, drop = FALSE]
  # Apply rename map
  if (!is.null(rename_map)) {
    nm <- colnames(rcf)
    new_nm <- ifelse(nm %in% names(rename_map),
                     unname(rename_map[nm]), nm)
    valid <- nchar(new_nm) > 0
    rcf <- rcf[, valid, drop = FALSE]
    colnames(rcf) <- new_nm[valid]
  }
  # Convert sample columns to numeric (taxid stays character)
  for (cn in setdiff(colnames(rcf), "taxid")) {
    rcf[[cn]] <- suppressWarnings(as.numeric(rcf[[cn]]))
  }
  rcf
}

# For each batch: control rename (CTRL1->NC-19 etc) THEN June wf-id remap.
# April/May: only the CTRL rename. June: CTRL rename + wf-id remap. We compose
# the maps so a single sweep applies both for June.
apr_rename <- CTRL_RENAME$apr
may_rename <- CTRL_RENAME$may
jun_rename <- c(JUN_MAP, CTRL_RENAME$jun)
# Resolve duplicates: CTRL labels override wf-id mapping if both keys present
# (CN-22 etc. handled via JUN_MAP; CTRL1/CTRL2 handled via CTRL_RENAME$jun).
jun_rename <- jun_rename[!duplicated(names(jun_rename))]

rcf_apr <- read_rcf_clean(BATCHES$apr, apr_rename)
rcf_may <- read_rcf_clean(BATCHES$may, may_rename)
rcf_jun <- read_rcf_clean(BATCHES$jun, jun_rename)

# Per-batch contaminant subtraction (replaces R/02's pooled-controls filter
# for this study). Applied BEFORE bind_rows so each batch's controls only
# affect that batch's samples. Control columns are dropped on the way out,
# so the combined rcf is control-free; config sets `controls: []` so R/02
# sees nothing to validate and skips its own subtraction.
rcf_apr <- subtract_controls_rcf(rcf_apr, BATCH_CONTROLS$apr, "apr")
rcf_may <- subtract_controls_rcf(rcf_may, BATCH_CONTROLS$may, "may")
rcf_jun <- subtract_controls_rcf(rcf_jun, BATCH_CONTROLS$jun, "jun")

to_long_rcf <- function(d) {
  scols <- setdiff(colnames(d), "taxid")
  d |>
    pivot_longer(all_of(scols), names_to = "sample", values_to = "count") |>
    mutate(count = tidyr::replace_na(count, 0))
}

rcf_long <- bind_rows(
  to_long_rcf(rcf_apr),
  to_long_rcf(rcf_may),
  to_long_rcf(rcf_jun)
)
rcf_summed <- rcf_long |>
  group_by(taxid, sample) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
rcf_wide <- pivot_wider(rcf_summed, names_from = sample,
                        values_from = count, values_fill = 0)

# Write the combined rcf preserving R/02's expected 2-header-row layout:
#   row 1: Samples, Kraken_outputs/<s1>, Kraken_outputs/<s1>, Kraken_outputs/<s1>, ...
#   row 2: Stats,   cnt, una, sco, cnt, una, sco, ...
#   row 3: Id,      (empty cells)
#   row 4+: <taxid>, cnt, 0, 0, cnt, 0, 0, ... (una/sco set to 0; R/02 only reads cnt)
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

  # Build data matrix
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
                   file.path(REP, "Re-centrifuge", "EBSP_samples.rcf.data.csv"))
cat("  written:", nrow(rcf_wide), "taxa,",
    length(setdiff(colnames(rcf_wide), "taxid")), "sample cols\n\n")

# === 5. WF-METAGENOMICS DIVERSITY ==========================================
cat("[5/5] Combining wf-metagenomics diversity tables...\n")

read_div <- function(batch_suffix, rename_map = NULL) {
  f <- file.path(REP, "wf-metagenomics",
                 paste0("wf-metagenomics-diversity_", batch_suffix, ".csv"))
  d <- read_csv(f, show_col_types = FALSE, progress = FALSE)
  # Drop the trailing "total" column (per-batch aggregate, not a sample)
  if ("total" %in% colnames(d)) d$total <- NULL
  if (!is.null(rename_map)) {
    nm <- colnames(d)
    new_nm <- ifelse(nm %in% names(rename_map),
                     unname(rename_map[nm]), nm)
    valid <- nchar(new_nm) > 0
    d <- d[, valid, drop = FALSE]
    colnames(d) <- new_nm[valid]
  }
  d
}

div_apr <- read_div(BATCHES$apr)
div_may <- read_div(BATCHES$may)
# For June: use ONLY the wf->metadata map (no CTRL remap; June div uses CN-22/CN-30
# which JUN_MAP already handles).
div_jun <- read_div(BATCHES$jun, rename_map = JUN_MAP)

# Identify samples that appear in ≥2 batches — precomputed metrics for those
# cannot be summed, so they get NA in the combined table. The pipeline's R/06
# will see NA rows and either drop them at the metadata join or fall back to
# computing from the combined counts.
samples_per_batch <- list(
  setdiff(colnames(div_apr), "Indices"),
  setdiff(colnames(div_may), "Indices"),
  setdiff(colnames(div_jun), "Indices")
)
counts <- table(unlist(samples_per_batch))
overlap_samples <- names(counts)[counts >= 2]
cat("  overlap samples (set to NA in combined):",
    paste(overlap_samples, collapse = ", "), "\n")

# Outer-join by Indices
div_all <- div_apr |>
  full_join(div_may, by = "Indices") |>
  full_join(div_jun, by = "Indices")

# Resolve duplicated column names (introduced by full_join when same sample
# appears in multiple batches): collapse to a single column per sample,
# valued NA if the sample is in overlap_samples.
all_samples <- unique(unlist(samples_per_batch))
res <- tibble(Indices = div_all$Indices)
for (s in all_samples) {
  cols_with_s <- grep(paste0("^", make.names(s), "($|\\.)|^", s, "($|\\.)"),
                      colnames(div_all), value = TRUE)
  cols_with_s <- intersect(cols_with_s,
                           grep(paste0("^", s), colnames(div_all), value = TRUE))
  if (s %in% overlap_samples) {
    res[[s]] <- NA_character_
  } else {
    # Single batch — find the column whose base name matches s exactly
    match_col <- cols_with_s[1]
    res[[s]] <- div_all[[match_col]]
  }
}

write_csv(res, file.path(REP, "wf-metagenomics", "wf-metagenomics-diversity.csv"))
cat("  written:", nrow(res), "metrics,",
    length(all_samples), "samples (",
    length(overlap_samples), "NA'd)\n\n")

cat("=== Done. ===\n")
