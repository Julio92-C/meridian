# 03_normalisation.R — TPM normalisation for gene elements (ARG/VF/MGE).
# Per §2.6 of the Chicken batch 1 paper: "raw reads counts were normalised as
# transcripts per million (TPM) for each GEs".

tpm_normalise <- function(counts, gene_length_kb) {
  rpk <- counts / gene_length_kb
  scaling <- sum(rpk, na.rm = TRUE) / 1e6
  if (scaling == 0) return(rep(0, length(rpk)))
  rpk / scaling
}

normalise_data <- function(cleaned, cfg) {
  pipeline_log(cfg, "TPM-normalising gene elements")
  df <- cleaned$abri_kraken2

  # ABRicate reports COVERAGE as "start-end/length"; derive gene length (bp).
  if ("COVERAGE" %in% colnames(df)) {
    df$gene_len_bp <- suppressWarnings(as.numeric(
      sub(".*/", "", df$COVERAGE)
    ))
  } else {
    df$gene_len_bp <- NA_real_
  }
  df$gene_len_kb <- df$gene_len_bp / 1000

  df <- df |>
    dplyr::group_by(sample) |>
    dplyr::mutate(tpm = tpm_normalise(count, gene_len_kb)) |>
    dplyr::ungroup()

  out <- file.path(cfg$project_root, cfg$outputs$datasets_dir,
                   "abri_kraken2_tpm.csv")
  readr::write_csv(df, out)
  df
}
