# sensitivity_inputs.R — inputs for the EB vs SP sensitivity runs (T2a, 2026-10-14).
#
# The EB-enriched taxa of the main run (C. acnes, M. restricta, ...) are concentrated in
# batch-2 exhaled breath, and batch 2 had no upstream control subtraction. Three reruns of
# config.yaml test whether EB vs SP survives that:
#   config_sens_b1.yaml     — B1 samples only   (Metadata/combined_sample_metadata_B1.csv)
#   config_sens_b2.yaml     — B2 samples only   (Metadata/combined_sample_metadata_B2.csv)
#   config_sens_noskin.yaml — all samples, skin/environmental taxa removed from the Bracken
#                             and Re-centrifuge inputs (Reports/pipeline_inputs_noskin/)
#
# Run after combine_inputs.R:
#   Rscript projects/lung_microbiome_final/sensitivity_inputs.R
suppressPackageStartupMessages({ library(readr); library(dplyr) })

root <- "C:/Users/julio/Desktop/Lung_microbiome/Final_Analysis"
inp  <- file.path(root, "Reports/pipeline_inputs")
out  <- file.path(root, "Reports/pipeline_inputs_noskin")
dir.create(out, showWarnings = FALSE)

# ---- per-batch metadata ------------------------------------------------------
meta <- read_csv(file.path(root, "Metadata/combined_sample_metadata.csv"), show_col_types = FALSE)
for (b in c("B1", "B2")) {
  m <- filter(meta, batch == b)
  write_csv(m, file.path(root, sprintf("Metadata/combined_sample_metadata_%s.csv", b)), na = "")
  message(sprintf("%s: %d samples, %d patients (%s)", b, nrow(m), n_distinct(m$ID),
                  paste(names(table(m$sample_type)), table(m$sample_type), sep = "=", collapse = ", ")))
}

# ---- skin / environmental taxa: each node plus its whole subtree --------------
# Genera: Cutibacterium, Malassezia. Species: Staphylococcus epidermidis, Moraxella osloensis,
# Pseudomonas putida. (S. aureus, P. aeruginosa and the rest of their genera are kept.)
SKIN <- c(Cutibacterium = 1912216, Malassezia = 55193, `S. epidermidis` = 1282,
          `M. osloensis` = 34062, `P. putida` = 303)
tax <- read_tsv(file.path(inp, "kraken2_combined_taxonomy.txt"), col_types = "ccc", trim_ws = FALSE)
depth <- nchar(tax$name) - nchar(sub("^ +", "", tax$name))
drop <- character(0)
for (t in SKIN) {
  i <- match(as.character(t), tax$taxid)
  stopifnot(!is.na(i))
  j <- i + 1
  while (j <= nrow(tax) && depth[j] > depth[i]) j <- j + 1
  drop <- c(drop, tax$taxid[i:(j - 1)])
}
drop <- unique(drop)
message(sprintf("skin subtrees: %d taxids (%s)", length(drop), paste(names(SKIN), collapse = ", ")))

# Bracken: taxid is the second-to-last column. Text-level filter keeps every other byte intact.
br <- readLines(file.path(inp, "bracken_arranged.csv"))
br_taxid <- vapply(strsplit(br, ",", fixed = TRUE), function(f) f[length(f) - 1], "")
keep <- seq_along(br) == 1 | !br_taxid %in% drop
writeLines(br[keep], file.path(out, "bracken_arranged.csv"))
message(sprintf("bracken: removed %d of %d rows", sum(!keep), length(br) - 1))

# Re-centrifuge: three header rows, then one row per taxid in the first column.
rc <- readLines(file.path(inp, "EBSP_combined.rcf.data.csv"))
rc_id <- sub(",.*$", "", rc)
keep <- seq_along(rc) <= 3 | !rc_id %in% drop
writeLines(rc[keep], file.path(out, "EBSP_combined.rcf.data.csv"))
message(sprintf("re-centrifuge: removed %d of %d rows", sum(!keep), length(rc) - 3))

for (f in c("kraken2_combined_taxonomy.txt", "summary_reportClean.csv"))
  file.copy(file.path(inp, f), file.path(out, f), overwrite = TRUE)
