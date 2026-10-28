#!/usr/bin/env Rscript
# prep_inputs.R — Materialise pipeline-ready inputs for lung_microbiome_2026
# (JC_2026-09_EBSP: RFH 2026 EBSP batch, 4 ONT runs EBSP0203/0304/0405R/0506R).
# Run ONCE (from the repo root) before run_pipeline.R, and again whenever the
# upstream reports or Metadata/excluded_samples.csv change:
#   Rscript projects/lung_microbiome_2026/prep_inputs.R
#
# Why a pre-step is necessary:
#   1. Upstream is nf-core/taxprofiler 2.0.1, not wf-metagenomics, so file
#      layouts differ from the 2025 lung study:
#        - Kraken2 combined report carries a `#Number of...` / `#S<N>` preamble
#          that readr::read_delim in R/01 cannot parse -> strip it (keep the
#          full report in DFS order with lvl_type for build_taxid_ancestry()).
#        - taxprofiler skips Bracken for OXFORD_NANOPORE, so Bracken was run
#          standalone (Bash_scripts/bracken.sh) and combined with KrakenTools
#          combine_kreports.py (`#<sample>` preamble, `<sample>_all/_lvl`
#          columns) -> reshape to the arranged wide `#perc,tot_all,<samples>,
#          taxid,name` layout R/02 expects.
#        - Re-centrifuge was run with `-c 1`, which names the control CTRL1
#          and appends derived *_CTRL_* / *_EXCLUSIVE_* / SHARED_SUMMARY /
#          Details columns. R/02's by-3 slicing assumes plain per-sample
#          cnt/una/sco triplets -> keep only the analysis samples' triplets.
#   2. Analysis-set decisions live in Metadata/excluded_samples.csv (built by
#      Metadata/build_sample_metadata.py), `action` column:
#        exclude          -> dropped (dropped patients 2413/2454/2470, SPN0-2025)
#        qc               -> dropped (CNN0-0609 / rcf CTRL1: QC record only, no subtraction)
#        merge:<target>   -> counts summed into <target> (EBR0-2466 -> EBF0-2466)
#      Abricate hits of a merged sample are relabelled to the target (hit-level
#      union = summed read evidence).
#   3. Low-depth samples (5 with 18-80 reads have no Bracken output) are KEPT
#      in the metadata by decision (2026-10-05); they simply drop out of the
#      Bracken-joined tables. Abricate has no rows for 7 samples (zero hits at
#      80/80 or failed runs — verify on Myriad).
#
# Outputs (created/overwritten; originals in Reports/ are never touched):
#   Reports/pipeline_inputs/kraken2_db1_combined_reports.txt
#   Reports/pipeline_inputs/bracken_arranged.csv
#   Reports/pipeline_inputs/EBSP_samples.rcf.data.csv
#   Reports/pipeline_inputs/summary_reportClean.csv
# Separate folder because Windows paths are case-insensitive (Reports/kraken2
# == Reports/Kraken2), so writing to the template's default paths would
# overwrite the upstream originals.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(data.table)
})

ROOT <- "C:/Users/julio/Desktop/Lung_microbiome/ONT_Long_read_seq/JC_2026-09_EBSP"
REP  <- file.path(ROOT, "Reports")
MET  <- file.path(ROOT, "Metadata")
OUT  <- file.path(REP, "pipeline_inputs")
dir.create(OUT, showWarnings = FALSE)

cat("=== prep_inputs.R (lung_microbiome_2026) ===\n")

# --- Analysis set + merge map ------------------------------------------------
meta     <- read_csv(file.path(MET, "sample_metadata.csv"), show_col_types = FALSE)
excluded <- read_csv(file.path(MET, "excluded_samples.csv"), show_col_types = FALSE)
analysis <- as.character(meta$sample)
merges   <- excluded |>
  filter(startsWith(action, "merge:")) |>
  transmute(source = sample, target = sub("^merge:", "", action))
stopifnot("merge target not in analysis set" = all(merges$target %in% analysis))
merge_to <- setNames(merges$target, merges$source)          # EBR0-2466 -> EBF0-2466
cat(sprintf("[meta] %d analysis samples; %d excluded/QC; merges: %s\n",
            length(analysis), nrow(excluded) - nrow(merges),
            paste(sprintf("%s->%s", merges$source, merges$target), collapse = ", ")))

report_coverage <- function(label, present) {
  miss <- setdiff(analysis, present)
  cat(sprintf("[%s] %d/%d analysis samples present%s\n", label,
              length(intersect(analysis, present)), length(analysis),
              if (length(miss)) paste0(" — missing: ", paste(miss, collapse = ", ")) else ""))
}

# Lines from the `#perc` header onward (drops combine_kreports preambles).
read_from_perc <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^#perc\t", lines)
  stopifnot(length(start) == 1)
  lines[start:length(lines)]
}

# === 1. Kraken2: strip preamble =============================================
k2_lines <- read_from_perc(file.path(REP, "kraken2", "kraken2_db1_combined_reports.txt"))
writeLines(k2_lines, file.path(OUT, "kraken2_db1_combined_reports.txt"))
cat(sprintf("[out] kraken2_db1_combined_reports.txt  (%d taxa rows)\n", length(k2_lines) - 1L))

# === 2. Bracken: combined kreport -> arranged wide ==========================
brk <- fread(text = read_from_perc(file.path(REP, "bracken", "bracken_db1_combined_reports.txt")),
             sep = "\t", header = TRUE, colClasses = "character", data.table = FALSE)
brk$name <- trimws(brk$name)
all_cols <- grep("_all$", names(brk), value = TRUE)
all_cols <- setdiff(all_cols, "tot_all")
brk_samples <- sub("_all$", "", all_cols)
counts <- as.data.frame(lapply(brk[all_cols], as.numeric), check.names = FALSE)
names(counts) <- brk_samples

for (src in names(merge_to)) {
  tgt <- merge_to[[src]]
  if (src %in% names(counts) && tgt %in% names(counts)) {
    counts[[tgt]] <- counts[[tgt]] + counts[[src]]
    cat(sprintf("[bracken] summed %s into %s\n", src, tgt))
  } else if (src %in% names(counts)) {
    stop(sprintf("bracken: merge target %s absent while source %s present", tgt, src))
  }
}
keep_brk <- intersect(analysis, names(counts))           # metadata order
report_coverage("bracken", keep_brk)
counts   <- counts[keep_brk]
tot_all  <- rowSums(counts)
n_total  <- sum(tot_all[brk$name %in% c("root", "unclassified")])
brk_out  <- data.frame(`#perc` = round(100 * tot_all / n_total, 4), tot_all = tot_all,
                       counts, taxid = brk$taxid, name = brk$name,
                       check.names = FALSE)
write_csv(brk_out, file.path(OUT, "bracken_arranged.csv"))
cat(sprintf("[out] bracken_arranged.csv  (%d taxa, %d samples)\n", nrow(brk_out), length(keep_brk)))

# === 3. Re-centrifuge: keep analysis samples' cnt/una/sco triplets ==========
rcf <- fread(file.path(REP, "Re-centrifuge", "EBSP_samples.rcf.data.csv"),
             header = FALSE, colClasses = "character", data.table = FALSE)
hdr   <- sub("^Kraken_outputs/", "", unlist(rcf[1, ]))
stats <- unlist(rcf[2, ])
stopifnot(hdr[1] == "Samples", stats[1] == "Stats")

triplet <- function(s) {
  idx <- which(hdr == s)
  if (length(idx) != 3 || !identical(unname(stats[idx]), c("cnt", "una", "sco"))) {
    stop(sprintf("rcf: sample %s does not map to one cnt/una/sco triplet", s))
  }
  idx
}
num  <- function(x) { v <- suppressWarnings(as.numeric(x)); v[is.na(v)] <- 0; v }
body <- 3:nrow(rcf)                                       # rows after Samples/Stats (Id + taxa)

for (src in names(merge_to)) {
  tgt <- merge_to[[src]]
  if (!src %in% hdr) next
  s_idx <- triplet(src); t_idx <- triplet(tgt)
  for (k in 1:2) {                                         # cnt, una: sum
    a <- rcf[body, t_idx[k]]; b <- rcf[body, s_idx[k]]
    both_blank <- a == "" & b == ""
    summed <- as.character(num(a) + num(b))
    summed[both_blank] <- ""
    rcf[body, t_idx[k]] <- summed
  }
  fill <- rcf[body, t_idx[3]] == "" & rcf[body, s_idx[3]] != ""   # sco: keep target's, else source's
  rcf[body[fill], t_idx[3]] <- rcf[body[fill], s_idx[3]]
  cat(sprintf("[rcf] summed %s into %s (cnt/una)\n", src, tgt))
}

keep_rcf <- intersect(analysis, hdr)
report_coverage("rcf", keep_rcf)
cols    <- c(1L, unlist(lapply(keep_rcf, triplet)))
rcf_out <- rcf[, cols]                                     # row 1 keeps the original Kraken_outputs/<s> names
fwrite(rcf_out, file.path(OUT, "EBSP_samples.rcf.data.csv"), col.names = FALSE, quote = FALSE)
cat(sprintf("[out] EBSP_samples.rcf.data.csv  (%d rows, %d samples; dropped CTRL1 + derived columns)\n",
            nrow(rcf_out) - 3L, length(keep_rcf)))

# === 4. Abricate: analysis samples only, merged sample relabelled ===========
abr <- read_csv(file.path(REP, "abricate", "summary_reportClean.csv"),
                col_types = cols(.default = col_character()))
abr <- abr[abr$`#FILE` != "#FILE", ]                      # defensive: stray repeated headers
abr_sample <- sub("_kraken2ID.*$", "", abr$`#FILE`)
for (src in names(merge_to)) {
  hit <- abr_sample == src
  abr$`#FILE`[hit] <- sub(paste0("^", src), merge_to[[src]], abr$`#FILE`[hit])
  abr_sample[hit]  <- merge_to[[src]]
  if (any(hit)) cat(sprintf("[abricate] relabelled %d hits %s -> %s\n", sum(hit), src, merge_to[[src]]))
}
keep_abr <- abr_sample %in% analysis
report_coverage("abricate", unique(abr_sample[keep_abr]))
write_csv(abr[keep_abr, ], file.path(OUT, "summary_reportClean.csv"), na = "")
cat(sprintf("[out] summary_reportClean.csv  (%d hits; dropped %d hits from excluded/QC samples)\n",
            sum(keep_abr), sum(!keep_abr)))

cat("=== Done. Run: Rscript run_pipeline.R projects/lung_microbiome_2026/config.yaml ===\n")
