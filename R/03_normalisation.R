# 03_normalisation.R — TPM normalisation per (sample, GENE).
# Mirrors normData.R from the chicken_batch1 reference scripts:
#   meanLength = mean(|END - START| + 1) within (sample, GENE)
#   meanCount  = number of hits of that gene in that sample
#   RPK        = meanCount / (meanLength / 1000)
#   TPM        = RPK / (sum(RPK) / 1e6) within sample
# Output columns: sample, GENE, meanLength, meanCount, RPK, TPM,
# plus PRODUCT/RESISTANCE/DATABASE and any cfg$metadata$fixed_effects /
# random_effect columns that were attached upstream.

normalise_data <- function(cleaned, cfg) {
  pipeline_log(cfg, "TPM-normalising gene elements")

  if (is.null(cleaned$merged)) {
    stop("R/03 requires cleaned$merged (abricate + bracken + metadata). ",
         "Check that cfg$inputs$bracken_report points at a valid file and ",
         "that R/02 produced abri_kraken2Bracken_merged.csv.")
  }
  df <- cleaned$merged

  required <- c("sample", "GENE", "START", "END")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Required columns missing from cleaned$merged: ",
         paste(missing, collapse = ", "))
  }

  df_genetable <- df |>
    dplyr::group_by(sample, GENE) |>
    dplyr::mutate(
      GeneLength = abs(END - START) + 1,
      geneCount  = dplyr::n()
    ) |>
    dplyr::summarise(
      meanLength = round(mean(GeneLength), 1),
      meanCount  = round(mean(geneCount), 1),
      .groups = "drop"
    ) |>
    dplyr::mutate(RPK = meanCount / (meanLength / 1000)) |>
    dplyr::select(sample, GENE, meanLength, meanCount, RPK) |>
    dplyr::distinct()

  df_scaled <- df_genetable |>
    dplyr::group_by(sample) |>
    dplyr::mutate(TPM = RPK / (sum(RPK) / 1e6)) |>
    dplyr::ungroup()

  meta_cols <- intersect(
    c("PRODUCT", "RESISTANCE", "DATABASE",
      cfg$metadata$fixed_effects, cfg$metadata$random_effect),
    colnames(df)
  )
  meta_per_gene <- df |>
    dplyr::select(sample, GENE, dplyr::all_of(meta_cols)) |>
    dplyr::distinct(sample, GENE, .keep_all = TRUE)

  out_df <- merge(df_scaled, meta_per_gene, by = c("sample", "GENE"))

  out_path <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                        "genetable_normdata.csv")
  readr::write_csv(out_df, out_path)
  out_df
}
