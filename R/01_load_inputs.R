# 01_load_inputs.R — read raw Kraken2/Bracken/ABRicate/Re-centrifuge reports
# and the study metadata. Returns a named list; downstream modules read from it.

load_inputs <- function(cfg) {
  pipeline_log(cfg, "Loading inputs")

  abricate <- readr::read_csv(
    file.path(cfg$project_root, cfg$inputs$abricate_summary),
    show_col_types = FALSE
  )

  kraken2 <- readr::read_delim(
    file.path(cfg$project_root, cfg$inputs$kraken2_report),
    delim = "\t", escape_double = FALSE, trim_ws = TRUE,
    show_col_types = FALSE
  )

  bracken <- tryCatch(
    readr::read_csv(
      file.path(cfg$project_root, cfg$inputs$bracken_report),
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )

  # Re-centrifuge path may be a glob.
  rcf_pattern <- file.path(cfg$project_root, cfg$inputs$recentrifuge_csv)
  rcf_files   <- Sys.glob(rcf_pattern)
  stopifnot("No Re-centrifuge csv matched the glob" = length(rcf_files) >= 1)
  recentrifuge <- readr::read_csv(rcf_files[[1]], show_col_types = FALSE)

  metadata <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file),
    show_col_types = FALSE
  )

  taxid_fixes <- NULL
  if (!is.null(cfg$taxid_fixes_file)) {
    fix_path <- file.path(cfg$project_root, cfg$taxid_fixes_file)
    if (file.exists(fix_path)) {
      taxid_fixes <- readr::read_csv(fix_path, show_col_types = FALSE)
    }
  }

  list(
    abricate     = abricate,
    kraken2      = kraken2,
    bracken      = bracken,
    recentrifuge = recentrifuge,
    metadata     = metadata,
    taxid_fixes  = taxid_fixes
  )
}
