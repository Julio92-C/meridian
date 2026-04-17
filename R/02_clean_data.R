# 02_clean_data.R — merge Abricate + Kraken2 + Re-centrifuge; apply taxid fixes;
# subtract negative-control counts; produce the "LR-GEs contigs" table used by
# downstream modules. Mirrors cleanData.R from the reference studies but
# parameterised by config (no hardcoded sample IDs or paths).

clean_data <- function(inputs, cfg) {
  pipeline_log(cfg, "Cleaning + decontaminating")

  # --- Abricate summary: normalise sample and taxid columns ----------------
  summary_df <- as.data.frame(inputs$abricate)
  names(summary_df)[names(summary_df) == "#FILE"] <- "sample"
  summary_df$sample <- sub("_.*", "", gsub(".fasta", "", summary_df$sample))
  summary_df <- summary_df |>
    tidyr::separate(SEQUENCE, into = c("sequence", "taxid"), sep = "_") |>
    dplyr::mutate(taxid = sub(".*\\|", "", taxid)) |>
    dplyr::distinct()

  # --- Kraken2 taxid → name lookup ----------------------------------------
  taxid_name <- dplyr::select(inputs$kraken2, taxid, name)

  # --- Re-centrifuge contaminant counts -----------------------------------
  rcf <- inputs$recentrifuge
  rcf <- rcf[-c(1, 2), ]
  # Keep only columns that are read counts (every 3rd column after the
  # Samples/taxid column) — matches the original cleanData.R indexing.
  keep_cols <- c(1, 2, seq(5, ncol(rcf), by = 3))
  keep_cols <- keep_cols[keep_cols <= ncol(rcf)]
  rcf <- rcf[, keep_cols]
  names(rcf)[names(rcf) == "Samples"] <- "taxid"
  colnames(rcf) <- sapply(colnames(rcf), function(n) {
    n <- gsub(".*?/", "", n)
    gsub("\\..*", "", n)
  })

  # --- Apply taxid fixes (replaces noncontaminants_list[NN, 4] <- "...") ---
  rcf_named <- dplyr::left_join(rcf, taxid_name, by = "taxid")
  if (!is.null(inputs$taxid_fixes)) {
    fix_lookup <- setNames(inputs$taxid_fixes$name, inputs$taxid_fixes$taxid)
    rcf_named$name <- dplyr::coalesce(fix_lookup[rcf_named$taxid], rcf_named$name)
  }

  # --- Split contaminant vs non-contaminant counts -------------------------
  controls <- cfg$metadata$controls
  sample_cols <- setdiff(colnames(rcf_named),
                         c("taxid", "name", "Classifier", controls))

  pivoted <- rcf_named |>
    tidyr::pivot_longer(cols = dplyr::all_of(sample_cols),
                        names_to = "sample", values_to = "count") |>
    dplyr::mutate(
      count = as.numeric(count),
      dplyr::across(dplyr::all_of(controls), as.numeric)
    )

  # Use mean of control columns as contaminant baseline.
  pivoted$control_max <- apply(pivoted[, controls, drop = FALSE], 1, max, na.rm = TRUE)

  noncontaminants <- pivoted |>
    dplyr::filter(count > control_max) |>
    dplyr::filter(!grepl(
      "root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples",
      name
    ))

  # --- Link Abricate → taxa (LR-GEs contigs) -------------------------------
  abri_kraken2 <- dplyr::inner_join(summary_df, noncontaminants,
                                    by = c("taxid", "sample"))

  # --- Persist cleaned outputs --------------------------------------------
  out_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(noncontaminants, file.path(out_dir, "noncontaminants_list.csv"))
  readr::write_csv(abri_kraken2,    file.path(out_dir, "abri_kraken2_cleaned.csv"))

  list(noncontaminants = noncontaminants, abri_kraken2 = abri_kraken2)
}
