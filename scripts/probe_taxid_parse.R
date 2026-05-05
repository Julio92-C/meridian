suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})
root <- "C:/Users/julio/OneDrive/Desktop/PC_JC_2024-11-29_Julio_Gallus"
abr <- read_csv(file.path(root, "Reports/Abricate/summary_reportClean.csv"),
                show_col_types = FALSE)
names(abr)[1] <- "sample"
abr$sample <- sub("_.*", "", gsub(".fasta", "", abr$sample))

chu <- abr |> dplyr::filter(sample == "D19", GENE == "chuA") |>
  dplyr::select(sample, GENE, SEQUENCE)

cat("--- raw Abricate D19/chuA hits ---\n")
print(chu)

cat("\n--- GT separate(sep=\"_\") result ---\n")
print(
  chu |>
    tidyr::separate(SEQUENCE, into = c("sequence", "taxid"),
                    sep = "_", extra = "drop") |>
    dplyr::mutate(taxid = sub(".*\\|", "", taxid))
)

cat("\n--- ours regex result ---\n")
print(
  chu |>
    dplyr::mutate(
      sequence = sub("_kraken:taxid\\|\\d+$", "", SEQUENCE),
      taxid    = sub(".*\\|", "", SEQUENCE)
    ) |>
    dplyr::select(sample, GENE, sequence, taxid)
)
