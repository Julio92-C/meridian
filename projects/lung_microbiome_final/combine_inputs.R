#!/usr/bin/env Rscript
# combine_inputs.R — Preprocessing step 4/4: merge the two RFH lung microbiome
# batches into one pipeline-ready dataset under Lung_microbiome/Final_Analysis/.
#
#   RFH_2025 : 46 sequenced analysis samples (Oct 2024 - Apr 2025), 3 ONT runs,
#              wf-metagenomics upstream, merged by
#              projects/lung_microbiome/combine_batches.R (per-run negative
#              control subtraction applied; re-sequenced samples summed).
#              Inputs: Lung_microbiome/metaOmics/{Reports,Metadata}
#   RFH_2026 : 89 analysis samples (Feb - Jun 2026), 4 ONT runs,
#              nf-core/taxprofiler 2.0.1 + standalone Bracken, prepared by
#              projects/lung_microbiome_2026/prep_inputs.R.
#              Inputs: ONT_Long_read_seq/JC_2026-09_EBSP/{Reports/pipeline_inputs,Metadata}
#
# Preceded by (Final_Analysis/Metadata/, Python):
#   step 1  build_batch1_2025_info.py -> intermediate/batch1_2025_sample_info.csv
#   step 2  build_clinical_info.py    -> intermediate/{peak_flow_by_sample,respiratory_diagnosis}.csv
#   step 3  build_sample_name_map.py  -> sample_name_map.csv (old -> new names)
# This step renames every sample to the unified code (<T><Y><M><V>-<PPPP>),
# sums technical replicates that share a new name (EB-42D -> E5A0-1242,
# SP-41 -> S4D0-2340), and writes the cleaned metadata. Decisions and counts
# go to Final_Analysis/Logs/04_combine_inputs.log, not into the metadata.
#
# Run from the repo root:  Rscript projects/lung_microbiome_final/combine_inputs.R
#
# Outputs (Final_Analysis/):
#   Metadata/combined_sample_metadata.csv            (all analysis samples)
#   Metadata/combined_sample_metadata_diagnosis.csv  (samples with an immunodeficiency diagnosis)
#   Reports/pipeline_inputs/kraken2_combined_taxonomy.txt   (merged NCBI tree, DFS)
#   Reports/pipeline_inputs/bracken_arranged.csv
#   Reports/pipeline_inputs/EBSP_combined.rcf.data.csv
#   Reports/pipeline_inputs/summary_reportClean.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(data.table)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

LUNG <- "C:/Users/julio/Desktop/Lung_microbiome"
S25  <- list(rep = file.path(LUNG, "metaOmics/Reports"),
             met = file.path(LUNG, "metaOmics/Metadata"))
S26  <- list(rep = file.path(LUNG, "ONT_Long_read_seq/JC_2026-09_EBSP/Reports/pipeline_inputs"),
             met = file.path(LUNG, "ONT_Long_read_seq/JC_2026-09_EBSP/Metadata"))
ROOT <- file.path(LUNG, "Final_Analysis")
OUT  <- file.path(ROOT, "Reports/pipeline_inputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, "Logs"), showWarnings = FALSE)

log_con <- file(file.path(ROOT, "Logs/04_combine_inputs.log"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)
cat(sprintf("=== combine_inputs.R — %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M")))

# === 1. Name map ==============================================================
nm <- read_csv(file.path(ROOT, "Metadata/sample_name_map.csv"),
               col_types = cols(.default = col_character()))
use <- nm |>
  filter(status == "analysis" | startsWith(status, "merged_into:")) |>
  mutate(target = ifelse(status == "analysis", new_sample, sub("^merged_into:", "", status)))
to_new <- setNames(use$target, use$old_sample)        # every input name -> final name
cat(sprintf("[map] %d input names -> %d analysis samples (%d replicate(s) summed)\n",
            nrow(use), n_distinct(use$target), sum(use$status != "analysis")))
for (i in which(use$status != "analysis")) {
  cat(sprintf("      %s -> %s (summed)\n", use$old_sample[i], use$target[i]))
}

# === 2. Metadata ==============================================================
pid <- function(x) sprintf("ID-%04d", as.integer(sub("^ID-", "", as.character(x))))

b1 <- read_csv(file.path(ROOT, "Metadata/intermediate/batch1_2025_sample_info.csv"),
               col_types = cols(.default = col_character())) |>
  transmute(old_sample = sample, seq_runs, gDNA = as.numeric(gDNA_ng_ul))
m26 <- read_csv(file.path(S26$met, "sample_metadata.csv"),
                col_types = cols(.default = col_character())) |>
  transmute(old_sample = sample, seq_runs = seq_run, gDNA = as.numeric(gDNA_ng_ul))

per_input <- use |>
  left_join(bind_rows(b1, m26), by = "old_sample") |>
  filter(batch == "RFH_2025" | status == "analysis")  # 2026 replicate already summed upstream
stopifnot("inputs without run/gDNA info" = !any(is.na(per_input$seq_runs)))

# Sequencing runs, chronological -> categorical seq_run R1..R7 (decision 2026-10-08).
# A sample sequenced in more than one run (reads summed) takes its first run.
RUNS <- c(R1 = "JC_PC_2025-04-11", R2 = "JC_PC_2025-05-13", R3 = "JC_2025-06-18",
          R4 = "EBSP0203", R5 = "EBSP0304", R6 = "EBSP0405R", R7 = "EBSP0506R")
run_code <- setNames(names(RUNS), RUNS)
cat("[runs] seq_run codes: ", paste(sprintf("%s=%s", names(RUNS), RUNS), collapse = ", "), "\n", sep = "")
seq_info <- per_input |>
  group_by(sample = target) |>
  summarise(runs = list(sort(unique(unname(run_code[unlist(strsplit(seq_runs, "+", fixed = TRUE))])))),
            gDNA_ng_ul = round(mean(gDNA), 2), .groups = "drop") |>
  mutate(seq_run = vapply(runs, function(r) r[1], ""))
stopifnot("unknown sequencing run name" = !any(is.na(unlist(seq_info$runs))))
multi <- seq_info[lengths(seq_info$runs) > 1, ]
for (i in seq_len(nrow(multi))) {
  cat(sprintf("[runs] %s sequenced in %s -> seq_run %s (first run)\n",
              multi$sample[i], paste(multi$runs[[i]], collapse = " + "), multi$seq_run[i]))
}
seq_info$runs <- NULL
for (s in unique(per_input$target[duplicated(per_input$target)])) {
  g <- per_input[per_input$target == s, ]
  cat(sprintf("[gDNA] %s: mean of %s = %.2f ng/ul\n", s,
              paste(sprintf("%s %.2f", g$old_sample, g$gDNA), collapse = " and "), mean(g$gDNA)))
}

pat <- read_csv(file.path(LUNG, "Metadata/patient_metadata_clean.csv"), show_col_types = FALSE) |>
  transmute(ID = pid(ID), Sex, Age_2026 = Age,
            # HG and PHGG merged into one hypogammaglobulinaemia group (decision 2026-10-05)
            Immunodeficiency_diagnosis = recode(Immunodeficiency_diagnosis, PHGG = "HG"),
            Treatment, Immunocompromised_type)

# Peak flow (mean of the 3 readings per visit) and respiratory diagnosis (batch-1 'Chest' notes + nurse returns;
# verbatim) from step 2 (build_clinical_info.py).
pf <- read_csv(file.path(ROOT, "Metadata/intermediate/peak_flow_by_sample.csv"),
               col_types = cols(.default = col_character())) |>
  transmute(ID, sample_type, collection_date = as.Date(sample_date),
            Peak_flow_mean = as.numeric(Peak_flow_mean)) |>
  distinct(ID, sample_type, collection_date, .keep_all = TRUE)
# Bronchiectasis (Yes/No; empty = not known) is step 2's category of that text; its basis
# (explicit/inferred/...) stays in the step-2 log and intermediate file.
rdx <- read_csv(file.path(ROOT, "Metadata/intermediate/respiratory_diagnosis.csv"),
                col_types = cols(.default = col_character())) |>
  mutate(Bronchiectasis = na_if(Bronchiectasis, "")) |>
  select(ID, Respiratory_diagnosis, Bronchiectasis)

# batch: sampling batch, categorical B1 = RFH_2025 (Oct 2024 - Apr 2025),
# B2 = RFH_2026 (Feb - Jun 2026). Letter prefixes keep batch/seq_run categorical
# (bare 1/2 would be read as numeric and break the pipeline's discrete palettes).
meta <- nm |>
  filter(status == "analysis") |>
  transmute(sample = new_sample, ID = pid(ID), sample_type,
            batch = c(RFH_2025 = "B1", RFH_2026 = "B2")[batch],
            collection_date = as.Date(collection_date), visit = as.integer(visit)) |>
  left_join(seq_info, by = "sample") |>
  left_join(pf, by = c("ID", "sample_type", "collection_date")) |>
  left_join(pat, by = "ID") |>
  # Grouped patient-level variables for the between-patient comparisons (Methods Table 3.4):
  #   Macrolide_prophylaxis  any macrolide (azithromycin or clarithromycin; the raw sheet spells
  #                          it 'Claritromycin') in the long-term treatment, alone or combined,
  #                          = Yes; every other regimen, including None, = No (Julio, 2026-10-13)
  #   Age_band               Age_2026 < 50 / 50-69 / >= 70 years. Labels sort in age order
  #                          (the pipeline orders group levels alphabetically); all patients are adults.
  mutate(Macrolide_prophylaxis = ifelse(is.na(Treatment), NA,
                                        ifelse(grepl("Azithromycin|Clarith?romycin", Treatment), "Yes", "No")),
         Age_band = cut(Age_2026, c(-Inf, 50, 70, Inf), right = FALSE,
                        labels = c("18-49", "50-69", "70+")) |> as.character()) |>
  left_join(rdx, by = "ID") |>
  group_by(ID) |>
  mutate(in_both_cohorts = n_distinct(batch) == 2) |>
  group_by(ID, batch) |>
  mutate(days_since_first_sample = ifelse(
    batch == "B2", as.integer(collection_date - min(collection_date)), NA_integer_)) |>
  ungroup() |>
  # patient x batch cluster for the patient-aware EB-vs-SP test (kept by request)
  mutate(ID_cohort = paste(ID, batch, sep = "|")) |>
  select(sample, ID, sample_type, batch, seq_run, collection_date, visit, gDNA_ng_ul,
         Peak_flow_mean, Sex, Age_2026, Immunodeficiency_diagnosis, Respiratory_diagnosis,
         Bronchiectasis, Treatment, Macrolide_prophylaxis, Immunocompromised_type, Age_band,
         in_both_cohorts, days_since_first_sample,
         ID_cohort) |>
  arrange(batch, ID, collection_date, sample)

# Spelling corrections (Metadata/spelling_corrections.csv; column,from,to,reason).
# Treatment values are replaced whole; Respiratory_diagnosis free text is patched as
# a fixed substring. Source spreadsheets are left untouched.
sp <- read_csv(file.path(ROOT, "Metadata/spelling_corrections.csv"),
               col_types = cols(.default = col_character()))
for (i in seq_len(nrow(sp))) {
  col <- sp$column[i]
  hit <- if (col == "Respiratory_diagnosis") grepl(sp$from[i], meta[[col]], fixed = TRUE)
         else meta[[col]] %in% sp$from[i]
  hit[is.na(hit)] <- FALSE
  if (any(hit)) {
    meta[[col]][hit] <- if (col == "Respiratory_diagnosis")
      sub(sp$from[i], sp$to[i], meta[[col]][hit], fixed = TRUE) else sp$to[i]
  }
  cat(sprintf("[spelling] %s: '%s' -> '%s' (%s): %d sample(s)
",
              col, sp$from[i], sp$to[i], sp$reason[i], sum(hit)))
}

stopifnot(!anyDuplicated(meta$sample), !any(is.na(meta$seq_run)), !any(is.na(meta$batch)))
# each run must belong to exactly one batch
stopifnot(all(tapply(meta$batch, meta$seq_run, function(b) length(unique(b))) == 1))
no_cov <- meta$sample[is.na(meta$Immunodeficiency_diagnosis)]
if (length(no_cov)) {
  cat(sprintf("[meta] %d sample(s) without patient covariates (recovered 2026-10-08; kept for EB-vs-SP and batch analyses, left out of the diagnosis file): %s
",
              length(no_cov), paste(no_cov, collapse = ", ")))
}
miss_pf <- meta$sample[is.na(meta$Peak_flow_mean)]
cat(sprintf("[meta] Peak_flow_mean (mean of 3 readings) for %d/%d samples%s
",
            sum(!is.na(meta$Peak_flow_mean)), nrow(meta),
            if (length(miss_pf)) paste0("; missing: ", paste(miss_pf, collapse = ", ")) else ""))
cat(sprintf("[meta] Respiratory_diagnosis for %d/%d patients (batch-1 'Chest' notes + research-nurse returns 2026-10-13)
",
            n_distinct(meta$ID[!is.na(meta$Respiratory_diagnosis)]), n_distinct(meta$ID)))
write_csv(meta, file.path(ROOT, "Metadata/combined_sample_metadata.csv"), na = "")
write_csv(meta[!is.na(meta$Immunodeficiency_diagnosis), ],
          file.path(ROOT, "Metadata/combined_sample_metadata_diagnosis.csv"), na = "")
# Subsets for the patient-level contrasts whose grouping is missing or single-patient for
# some samples (the pipeline does not drop samples with an empty group):
#   Bronchiectasis          Yes vs No; patients not known (not assessed / not documented) left out
#   Immunocompromised_type  P vs S; the single P/S patient (ID-1761) left out
subsets <- list(
  bronchiectasis         = !is.na(meta$Bronchiectasis),
  immunocompromised_type = meta$Immunocompromised_type %in% c("P", "S"))
for (nm_s in names(subsets)) {
  sub <- meta[subsets[[nm_s]], ]
  write_csv(sub, file.path(ROOT, sprintf("Metadata/combined_sample_metadata_%s.csv", nm_s)), na = "")
  cat(sprintf("[meta] %s subset: %d samples / %d patients (left out: %s)
", nm_s, nrow(sub),
              n_distinct(sub$ID), paste(sort(unique(meta$ID[!subsets[[nm_s]]])), collapse = ", ")))
}
for (v in c("Sex", "Immunodeficiency_diagnosis", "Immunocompromised_type", "Bronchiectasis",
            "Macrolide_prophylaxis", "Age_band")) {
  tb <- meta |> filter(!is.na(.data[[v]])) |> group_by(level = .data[[v]]) |>
    summarise(n = n(), p = n_distinct(ID), .groups = "drop")
  cat(sprintf("[meta] %s: %s
", v, paste(sprintf("%s %d samples/%d patients", tb$level, tb$n, tb$p),
                                            collapse = "; ")))
}
cat(sprintf("[meta] %d samples (%s) / %d patients (%d in both batches); columns: %s\n",
            nrow(meta), paste(names(table(meta$batch)), table(meta$batch), sep = "=", collapse = ", "),
            n_distinct(meta$ID), n_distinct(meta$ID[meta$in_both_cohorts]),
            paste(names(meta), collapse = ", ")))
samples <- meta$sample

# === 3. Kraken2 taxonomy: merge the two NCBI trees ==========================
# R/02 needs one DFS-ordered report (taxid, lvl_type, name) with unique taxids
# for name + phylum..genus ancestry lookups. Parse both indented trees, union
# nodes by taxid (2026 names/placement win: newer NCBI taxonomy, e.g.
# Bacillota vs Firmicutes), then write a single DFS walk. Not sample-specific.
read_tree <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#(Number|Total|S[0-9]+\t)", lines)]
  hdr   <- strsplit(lines[1], "\t")[[1]]
  rows  <- strsplit(lines[-1], "\t")
  ti <- match("taxid", hdr); li <- match("lvl_type", hdr); ni <- match("name", hdr)
  data.frame(taxid = vapply(rows, `[`, "", ti),
             lvl_type = vapply(rows, `[`, "", li),
             raw_name = vapply(rows, function(r) r[ni], ""),
             stringsAsFactors = FALSE)
}
tree_parents <- function(tr) {
  depth  <- (nchar(tr$raw_name) - nchar(sub("^ +", "", tr$raw_name))) / 2
  # The merged 2025 report has a few indentation jumps of >1 level; attach
  # such nodes to the nearest listed ancestor (the deepest one on the stack).
  parent <- character(nrow(tr)); stack <- character(0)
  for (i in seq_len(nrow(tr))) {
    d <- depth[i]
    stack <- stack[seq_len(min(length(stack), d))]
    parent[i] <- if (d == 0 || length(stack) == 0) "" else stack[length(stack)]
    stack <- c(stack, tr$taxid[i])
  }
  tr$parent <- parent; tr$name <- trimws(tr$raw_name); tr
}
t26 <- tree_parents(read_tree(file.path(S26$rep, "kraken2_db1_combined_reports.txt")))
t25 <- tree_parents(read_tree(file.path(S25$rep, "Kraken/kraken2_db1_combined_reports.txt")))
nodes <- rbind(t26[, c("taxid", "lvl_type", "name", "parent")],
               t25[!t25$taxid %in% t26$taxid, c("taxid", "lvl_type", "name", "parent")])
orphan <- nodes$parent != "" & !nodes$parent %in% nodes$taxid
nodes$parent[orphan] <- "1"                   # 2025-only node with unknown parent -> root
nodes$parent[nodes$parent == ""] <- ".TOP"    # [[""]] cannot index a list
kids <- split(nodes$taxid, factor(nodes$parent, levels = unique(nodes$parent)))
idx  <- setNames(seq_len(nrow(nodes)), nodes$taxid)
out_rows <- character(0)
walk <- function(tx, depth) {
  i <- idx[[tx]]
  out_rows[length(out_rows) + 1L] <<- paste(nodes$taxid[i], nodes$lvl_type[i],
                                           paste0(strrep("  ", depth), nodes$name[i]), sep = "\t")
  for (k in kids[[tx]] %||% character(0)) walk(k, depth + 1L)
}
for (top in kids[[".TOP"]]) walk(top, 0L)
stopifnot(length(out_rows) == nrow(nodes))
writeLines(c("taxid\tlvl_type\tname", out_rows), file.path(OUT, "kraken2_combined_taxonomy.txt"))
cat(sprintf("[kraken2] merged taxonomy: %d taxa (2026 %d + 2025-only %d; %d re-rooted orphans)\n",
            nrow(nodes), nrow(t26), sum(!t25$taxid %in% t26$taxid), sum(orphan)))

# === 4. Bracken: rename, sum replicates, outer join by taxid ================
rd_brk <- function(p) {
  b <- read_csv(p, show_col_types = FALSE, col_types = cols(.default = col_character()))
  b$name <- trimws(b$name); b
}
b25 <- rd_brk(file.path(S25$rep, "Bracken/bracken_arranged.csv"))
b26 <- rd_brk(file.path(S26$rep, "bracken_arranged.csv"))
i25 <- intersect(names(to_new), names(b25)); i26 <- intersect(names(to_new), names(b26))
long <- bind_rows(
  tidyr::pivot_longer(b25[, c("taxid", i25)], all_of(i25), names_to = "old", values_to = "count"),
  tidyr::pivot_longer(b26[, c("taxid", i26)], all_of(i26), names_to = "old", values_to = "count")
) |> mutate(sample = unname(to_new[old]), count = as.numeric(count))
tax_names <- bind_rows(b26[, c("taxid", "name")], b25[, c("taxid", "name")]) |>
  distinct(taxid, .keep_all = TRUE)                       # 2026 name wins
wide <- long |>
  group_by(taxid, sample) |> summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0)
keep_b <- intersect(samples, names(wide))
wide <- wide[, c("taxid", keep_b)] |> left_join(tax_names, by = "taxid")
tot  <- rowSums(wide[keep_b])
n_tot <- sum(tot[wide$name %in% c("root", "unclassified")])
brk_out <- data.frame(`#perc` = round(100 * tot / n_tot, 4), tot_all = tot,
                      wide[keep_b], taxid = wide$taxid, name = wide$name,
                      check.names = FALSE)
write_csv(brk_out, file.path(OUT, "bracken_arranged.csv"))
cat(sprintf("[bracken] %d taxa x %d samples; no Bracken output for: %s\n",
            nrow(brk_out), length(keep_b), paste(setdiff(samples, keep_b), collapse = ", ")))

# === 5. Re-centrifuge: rename, sum replicate triplets, outer join on taxid ==
rd_rcf <- function(p) fread(p, header = FALSE, colClasses = "character", data.table = FALSE)
long_rcf <- function(rcf) {
  hdr <- sub("^Kraken_outputs/", "", trimws(unlist(rcf[1, ])))
  st  <- trimws(unlist(rcf[2, ]))
  ins <- intersect(names(to_new), hdr)
  rbindlist(lapply(ins, function(s) {
    ix <- which(hdr == s)
    if (length(ix) != 3 || !identical(unname(st[ix]), c("cnt", "una", "sco")))
      stop("rcf: ", s, " is not one cnt/una/sco triplet")
    data.table(taxid = rcf[-(1:3), 1], sample = to_new[[s]],
               cnt = rcf[-(1:3), ix[1]], una = rcf[-(1:3), ix[2]], sco = rcf[-(1:3), ix[3]])
  }))
}
num <- function(x) { v <- suppressWarnings(as.numeric(x)); v[is.na(v)] <- 0; v }
lr <- rbind(long_rcf(rd_rcf(file.path(S25$rep, "Re-centrifuge/EBSP_samples.rcf.data.csv"))),
            long_rcf(rd_rcf(file.path(S26$rep, "EBSP_samples.rcf.data.csv"))))
lr <- lr[cnt != "" | una != "" | sco != ""]
agg <- lr[, .(cnt = if (all(cnt == "")) "" else as.character(sum(num(cnt))),
              una = if (all(una == "")) "" else as.character(sum(num(una))),
              sco = if (any(sco != "")) sco[sco != ""][1] else ""),   # score: first replicate's
          by = .(taxid, sample)]
keep_r <- intersect(samples, unique(agg$sample))
taxa   <- unique(lr$taxid)
mat <- matrix("", nrow = length(taxa), ncol = 1 + 3 * length(keep_r))
mat[, 1] <- taxa
ti <- match(agg$taxid, taxa); si <- match(agg$sample, keep_r)
ok <- !is.na(si)
mat[cbind(ti[ok], 1 + 3 * (si[ok] - 1) + 1)] <- agg$cnt[ok]
mat[cbind(ti[ok], 1 + 3 * (si[ok] - 1) + 2)] <- agg$una[ok]
mat[cbind(ti[ok], 1 + 3 * (si[ok] - 1) + 3)] <- agg$sco[ok]
header <- rbind(c("Samples", paste0("Kraken_outputs/", rep(keep_r, each = 3))),
                c("Stats",   rep(c("cnt", "una", "sco"), length(keep_r))),
                c("Id",      rep("", 3 * length(keep_r))))
fwrite(as.data.frame(rbind(header, mat), stringsAsFactors = FALSE),
       file.path(OUT, "EBSP_combined.rcf.data.csv"), col.names = FALSE, quote = FALSE)
cat(sprintf("[rcf] %d taxa rows x %d samples\n", length(taxa), length(keep_r)))

# === 6. Abricate: rename (replicate hits pooled) =============================
a25 <- read_csv(file.path(S25$rep, "Abricate/summary_reportClean.csv"),
                col_types = cols(.default = col_character()))
a26 <- read_csv(file.path(S26$rep, "summary_reportClean.csv"),
                col_types = cols(.default = col_character()))
stopifnot(identical(names(a25), names(a26)))
abr <- bind_rows(a25, a26)
abr <- abr[abr$`#FILE` != "#FILE", ]
old_s <- sub("_kraken2ID.*$", "", abr$`#FILE`)
new_s <- unname(to_new[old_s])
keep_a <- !is.na(new_s)
abr$`#FILE`[keep_a] <- paste0(new_s[keep_a], substring(abr$`#FILE`[keep_a], nchar(old_s[keep_a]) + 1))
write_csv(abr[keep_a, ], file.path(OUT, "summary_reportClean.csv"), na = "")
cat(sprintf("[abricate] %d hits from %d samples (dropped %d hits not in the analysis set); no hits for: %s\n",
            sum(keep_a), n_distinct(new_s[keep_a]), sum(!keep_a),
            paste(setdiff(samples, new_s[keep_a]), collapse = ", ")))

cat("=== Done. Next: Rscript run_pipeline.R projects/lung_microbiome_final/config*.yaml ===\n")
