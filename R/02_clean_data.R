# 02_clean_data.R — merge Abricate + Kraken2 + Re-centrifuge; apply taxid fixes;
# optionally subtract negative-control counts; produce the "LR-GEs contigs"
# table used by downstream modules. Mirrors cleanData.R from the reference
# studies but parameterised by config (no hardcoded sample IDs or paths).
# Studies without negative controls leave cfg$metadata$controls empty/null
# and the control-baseline subtraction is skipped. metadata.control_match_cols
# limits each sample's baseline to the controls of its own run/batch.

# Which negative controls each sample is compared against. With no match
# columns every sample gets every control. With match columns (e.g.
# c("seq_run", "batch")) a sample gets the controls that share its value in
# the FIRST column where any control does: its own run's controls, else its
# batch's. A sample with no match in any column gets none (not subtracted).
control_sets_for_samples <- function(meta, sid_col, samples, controls, match_cols,
                                     log = function(msg) pipeline_log(NULL, msg)) {
  if (length(match_cols) == 0) {
    return(setNames(rep(list(controls), length(samples)), samples))
  }
  missing_cols <- setdiff(match_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(sprintf("metadata.control_match_cols not in metadata columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  ids <- as.character(meta[[sid_col]])
  missing_ctrls <- setdiff(controls, ids)
  if (length(missing_ctrls) > 0) {
    stop(sprintf("Controls need a metadata row for control_match_cols: %s",
                 paste(missing_ctrls, collapse = ", ")))
  }
  value <- function(col, id) {
    v <- as.character(meta[[col]][match(id, ids)])
    ifelse(is.na(v) | v == "", NA_character_, v)
  }
  # Columns with no metadata row (e.g. Re-centrifuge's derived *_EXCLUSIVE_* /
  # *_SHARED_* summary columns) keep the all-controls baseline.
  sets <- lapply(samples, function(s) {
    if (!s %in% ids) return(controls)
    for (col in match_cols) {
      sv <- value(col, s)
      if (is.na(sv)) next
      cv  <- value(col, controls)
      hit <- controls[!is.na(cv) & cv == sv]
      if (length(hit) > 0) return(hit)
    }
    character(0)
  })
  names(sets) <- samples
  none <- samples[lengths(sets) == 0]
  if (length(none) > 0) {
    log(sprintf("No negative control shares %s with %d sample(s); not control-subtracted: %s",
                paste(match_cols, collapse = " or "), length(none),
                paste(none, collapse = ", ")))
  }
  sets
}

# Per-row control baseline for the long table (one row per taxon x sample,
# control counts as wide columns): the max over that sample's own controls,
# -Inf when it has none (so `count > control_max` keeps everything).
control_max_by_sample <- function(long, sets) {
  key <- vapply(sets, paste, character(1), collapse = "\r")
  out <- rep(-Inf, nrow(long))
  for (k in unique(key[nzchar(key)])) {
    rows <- long$sample %in% names(key)[key == k]
    cols <- strsplit(k, "\r", fixed = TRUE)[[1]]
    out[rows] <- apply(as.matrix(long[rows, cols, drop = FALSE]), 1, max, na.rm = TRUE)
  }
  out
}

clean_data <- function(inputs, cfg) {
  pipeline_log(cfg, "Cleaning + decontaminating")

  # --- Abricate summary: normalise sample and taxid columns ----------------
  # SEQUENCE values look like "<contig_id>_kraken:taxid|<taxid>", with an
  # optional contig-variant suffix between the UUID and "_kraken" (e.g.
  # "uuid_1_kraken:taxid|543"). Splitting on "_" loses the variant suffix
  # AND the actual taxid for those rows, so extract via regex on the stable
  # "_kraken:taxid|<digits>" tail instead.
  summary_df <- as.data.frame(inputs$abricate)
  names(summary_df)[names(summary_df) == "#FILE"] <- "sample"
  # Abricate #FILE is "<sample>_kraken2ID_PFP.fasta" across studies. Strip that
  # tail to recover the sample id. The previous sub("_.*") truncated at the
  # FIRST underscore, which silently mangled sample ids that themselves contain
  # underscores (e.g. wetlands_flyway GA1_1 -> GA1, breaking the metadata join).
  # Stripping the _kraken2ID... suffix is backward-compatible for underscore-free
  # ids (chicken D19, hospital 1SW) and preserves the rest.
  summary_df$sample <- sub("_kraken2ID.*$", "", summary_df$sample)
  summary_df$sample <- sub("\\.fasta$", "", summary_df$sample)
  summary_df <- summary_df |>
    dplyr::mutate(
      sequence = sub("_kraken:taxid\\|\\d+$", "", SEQUENCE),
      taxid    = sub(".*\\|", "", SEQUENCE)
    ) |>
    dplyr::select(-SEQUENCE) |>
    dplyr::distinct()

  # --- Kraken2 taxid → name lookup ----------------------------------------
  # A lookup concatenated from several reports (e.g. hospital Air + Surface)
  # can list a taxid more than once; keep one row so the joins below stay 1:1.
  taxid_name <- dplyr::select(inputs$kraken2, taxid, name) |>
    dplyr::distinct(taxid, .keep_all = TRUE)

  # --- Kraken2 taxid → ancestry lookup (phylum + genus) -------------------
  # Standard kraken2 combined-report rows are in DFS pre-order, with
  # `lvl_type` carrying the NCBI rank code (U / R / R1 / D / D1 / P / C /
  # O / F / G / S / S1 / ...). Walking rows top-to-bottom and tracking the
  # most recent P-rank and G-rank ancestors yields per-taxid ancestry
  # columns that propagate onto `noncontaminants` → `abri_kraken2` →
  # `merged`. Phylum feeds the C11 taxon→ARG→MGE Sankey (R/12); genus is
  # the rank-aware source for C4 `genus_heatmap` (R/05) and the
  # category-level Mantel rollup (R/12), replacing the rank-blind
  # first-word string heuristic that previously leaked phyla / classes /
  # orders / families into "genus" plots. Rows above each rank
  # (Domain / Root / Unclassified for phylum; everything above G for
  # genus) map to NA on that column.
  # Walk the full ordered table first (duplicates included, so each report's
  # tree context is intact), then de-duplicate.
  taxid_ancestry <- build_taxid_ancestry(inputs$kraken2) |>
    dplyr::distinct(taxid, .keep_all = TRUE)

  # --- Re-centrifuge contaminant counts -----------------------------------
  # NOTE: the by-3 column slicing below is the original cleanData.R logic
  # and is known to be too coarse for the real Re-centrifuge layout (mixes
  # ranks, can include trailing Rank/Name columns). Full rewrite is queued
  # in docs/pipeline_rework_scoping.md. Until then we drop the obvious
  # trailing character columns so as.numeric() doesn't warn on them.
  rcf <- inputs$recentrifuge
  rcf <- rcf[-c(1, 2), ]
  # Drop trailing Rank/Name columns. Their row-1 header is "Details" (so the
  # CSV reader auto-disambiguates them as "Details", "Details.1" via fread,
  # or "Details", "Details...2" via readr). Match both styles.
  rcf <- rcf[, !grepl("^Details($|\\.)", colnames(rcf))]
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
  rcf$taxid        <- as.character(rcf$taxid)
  taxid_name$taxid <- as.character(taxid_name$taxid)
  rcf_named <- dplyr::left_join(rcf, taxid_name, by = "taxid")
  if (!is.null(inputs$taxid_fixes)) {
    fix_lookup <- setNames(inputs$taxid_fixes$name, inputs$taxid_fixes$taxid)
    rcf_named$name <- dplyr::coalesce(fix_lookup[rcf_named$taxid], rcf_named$name)
  }
  taxid_ancestry$taxid <- as.character(taxid_ancestry$taxid)
  rcf_named <- dplyr::left_join(rcf_named, taxid_ancestry, by = "taxid")

  # --- Split contaminant vs non-contaminant counts -------------------------
  # Negative controls are optional. Three cases:
  #   1. controls declared and all present in rcf  -> validate, then subtract
  #   2. controls declared but ALL absent from rcf -> assume upstream
  #      subtraction (e.g. a project-specific combine script that strips
  #      control columns after per-batch filtering); keep `controls`
  #      populated downstream so R/06 / R/14 still know which metadata
  #      rows are controls, but skip this stage's filter.
  #   3. partial mismatch                          -> hard error (typo)
  declared <- cfg$metadata$controls
  if (is.null(declared)) declared <- character(0)
  ctrl_present  <- intersect(declared, colnames(rcf_named))
  missing_ctrls <- setdiff(declared, colnames(rcf_named))

  upstream_subtracted <- length(declared) > 0 && length(ctrl_present) == 0
  if (upstream_subtracted) {
    pipeline_log(cfg, sprintf(
      "Declared controls (%s) absent from Re-centrifuge columns — assuming upstream subtraction; R/02 baseline filter skipped",
      paste(declared, collapse = ", ")
    ))
    controls <- character(0)
  } else if (length(missing_ctrls) > 0) {
    stop(sprintf(
      "Negative controls declared in config but not found in Re-centrifuge data: %s",
      paste(missing_ctrls, collapse = ", ")
    ))
  } else {
    controls <- declared
  }
  has_controls <- length(controls) > 0

  if (!has_controls && !upstream_subtracted) {
    pipeline_log(cfg, "No negative controls declared — skipping control-baseline subtraction")
  }

  sample_cols <- setdiff(colnames(rcf_named),
                         c("taxid", "name",
                           "phylum", "class", "order", "family", "genus",
                           "tax_rank", "lineage_class", "Classifier", controls))

  pivoted <- rcf_named |>
    tidyr::pivot_longer(cols = dplyr::all_of(sample_cols),
                        names_to = "sample", values_to = "count") |>
    dplyr::mutate(count = as.numeric(count))

  if (has_controls) {
    pivoted <- dplyr::mutate(
      pivoted,
      dplyr::across(dplyr::all_of(controls), as.numeric)
    )
    # metadata.control_match_cols (optional) pairs each sample with the
    # controls that share its run/batch; unset = every control for every sample.
    sets <- control_sets_for_samples(inputs$metadata, cfg$metadata$sample_id_col,
                                     unique(pivoted$sample), controls,
                                     unlist(cfg$metadata$control_match_cols),
                                     log = function(msg) pipeline_log(cfg, msg))
    pivoted$control_max <- control_max_by_sample(pivoted, sets)
  } else {
    # No controls → no baseline to subtract; count > -Inf is always TRUE.
    pivoted$control_max <- -Inf
  }

  noncontaminants <- pivoted |>
    dplyr::filter(count > control_max) |>
    dplyr::filter(!grepl(
      "root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples",
      name
    ))

  # --- Optional: microbial taxa only ---------------------------------------
  # filters.microbial_only drops plant (Embryophyta), animal (Metazoa, incl.
  # human) and unplaced taxa by Kraken2 lineage, before anything downstream
  # (taxonomy, diversity, DA, gene linkage, network) sees them. Off by default.
  excluded_non_microbial <- NULL
  if (isTRUE(cfg$filters$microbial_only)) {
    microbial <- c("prokaryote", "virus", "microbial_eukaryote")
    cls <- noncontaminants$lineage_class
    cls[is.na(cls)] <- "unmatched"
    bad <- !(cls %in% microbial)
    excluded_non_microbial <- noncontaminants[bad, ] |>
      dplyr::mutate(lineage_class = cls[bad]) |>
      dplyr::filter(.data$count > 0) |>
      dplyr::group_by(.data$taxid, .data$name, .data$lineage_class) |>
      dplyr::summarise(total_count = sum(.data$count), n_samples = dplyr::n(),
                       .groups = "drop") |>
      dplyr::arrange(dplyr::desc(.data$total_count))
    pipeline_log(cfg, sprintf(
      "Microbial-only filter: dropped %d non-microbial taxa (%s reads)",
      nrow(excluded_non_microbial),
      format(sum(excluded_non_microbial$total_count), big.mark = ",")
    ))
    noncontaminants <- noncontaminants[!bad, ]
  }

  # --- Link Abricate → taxa (LR-GEs contigs) -------------------------------
  abri_kraken2 <- dplyr::inner_join(summary_df, noncontaminants,
                                    by = c("taxid", "sample"))

  # Drop duplicate hits (same contig, same coordinates) — mirrors the
  # `distinct(sequence, START, END, .keep_all = T)` step in normData.R.
  abri_kraken2 <- dplyr::distinct(abri_kraken2, sequence, START, END,
                                  .keep_all = TRUE)

  # --- Optional: Bracken merge + metadata join ----------------------------
  # Produces the abri_kraken2Bracken_merged table consumed by R/03 and other
  # downstream stages. Sample columns in the Bracken file are detected from
  # the metadata's sample_id column (no positional slicing); studies without
  # a Bracken report leave `merged` NULL.
  merged <- NULL
  if (!is.null(inputs$bracken)) {
    bracken    <- inputs$bracken
    sid_col    <- cfg$metadata$sample_id_col
    meta_ids   <- as.character(inputs$metadata[[sid_col]])
    brk_smp    <- intersect(meta_ids, colnames(bracken))
    if (length(brk_smp) == 0) {
      stop("No metadata sample IDs match Bracken column headers")
    }

    bracken_long <- bracken |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of(brk_smp),
        names_to  = "sample",
        values_to = "sampleCount"
      ) |>
      dplyr::mutate(
        sampleCount = as.numeric(sampleCount),
        taxid       = as.character(taxid)
      ) |>
      dplyr::filter(!grepl(
        "root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples",
        name
      ))

    rename_map <- c(TotalBCount = "tot_all", PercB = "#perc")
    rename_map <- rename_map[rename_map %in% colnames(bracken_long)]
    if (length(rename_map) > 0) {
      bracken_long <- dplyr::rename(bracken_long, !!!rename_map)
    }
    keep <- intersect(
      c("sample", "taxid", "name", "sampleCount", "PercB", "TotalBCount"),
      colnames(bracken_long)
    )
    bracken_long <- bracken_long[, keep, drop = FALSE]

    abri_kraken2$taxid <- as.character(abri_kraken2$taxid)
    merged <- dplyr::inner_join(abri_kraken2, bracken_long,
                                by = c("sample", "taxid", "name"))

    # Attach metadata columns declared in cfg (fixed_effects, random_effect).
    meta_attach <- unique(c(cfg$metadata$fixed_effects,
                            cfg$metadata$random_effect))
    meta_attach <- intersect(meta_attach, colnames(inputs$metadata))
    if (length(meta_attach) > 0) {
      md <- inputs$metadata[, c(sid_col, meta_attach), drop = FALSE]
      names(md)[names(md) == sid_col] <- "sample"
      md$sample      <- as.character(md$sample)
      merged$sample  <- as.character(merged$sample)
      merged <- dplyr::left_join(merged, md, by = "sample")
    }
  }

  # --- Persist cleaned outputs --------------------------------------------
  out_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(noncontaminants, file.path(out_dir, "noncontaminants_list.csv"))
  readr::write_csv(abri_kraken2,    file.path(out_dir, "abri_kraken2_cleaned.csv"))
  if (!is.null(merged)) {
    readr::write_csv(merged, file.path(out_dir, "abri_kraken2Bracken_merged.csv"))
  }
  if (!is.null(excluded_non_microbial)) {
    readr::write_csv(excluded_non_microbial,
                     file.path(out_dir, "excluded_non_microbial_taxa.csv"))
  }

  list(noncontaminants = noncontaminants,
       abri_kraken2    = abri_kraken2,
       merged          = merged)
}

# Build a per-taxid ancestry data frame from a kraken2 combined-report
# frame. Tracks every NCBI rank above species (phylum, class, order,
# family, genus). Relies on the report being in NCBI DFS pre-order with
# `lvl_type` carrying rank codes — both standard for kraken2
# combined-reports.
#
# For each ancestor rank we track a "current" value that is:
#   - cleared (NA) when the row's rank is *strictly above* that ancestor,
#   - set when the row's rank is exactly that ancestor (bare `P`, `C`,
#     `O`, `F`, or `G`),
#   - carried forward when the row's rank is at or below that ancestor
#     (so e.g. species rows inherit the genus, family, order, class, and
#     phylum their kraken2 lineage set).
#
# Reset uses the first letter of the rank code so sub-ranks (`P1`, `C1`,
# `F1`, ...) reset their respective ancestor lookups correctly.
build_taxid_ancestry <- function(kraken2) {
  needed <- c("taxid", "name", "lvl_type")
  miss <- setdiff(needed, colnames(kraken2))
  if (length(miss) > 0) {
    stop(sprintf(
      "build_taxid_ancestry: kraken2 input missing column(s) %s",
      paste(miss, collapse = ", ")
    ))
  }
  rk <- as.character(kraken2$lvl_type)
  nm <- trimws(as.character(kraken2$name))
  td <- as.character(kraken2$taxid)

  # Rank codes that sit strictly above each ancestor we track. Mapping
  # is hierarchical: phylum-reset is the strictest (only U/R/D/K above),
  # each subsequent rank adds the parent ranks above it.
  rank_resets <- list(
    phylum = c("U", "R", "D", "K"),
    class  = c("U", "R", "D", "K", "P"),
    order  = c("U", "R", "D", "K", "P", "C"),
    family = c("U", "R", "D", "K", "P", "C", "O"),
    genus  = c("U", "R", "D", "K", "P", "C", "O", "F")
  )
  rank_set_code <- c(phylum = "P", class = "C", order = "O",
                     family = "F", genus  = "G")

  cols <- lapply(names(rank_resets), function(r) character(length(rk)))
  names(cols) <- names(rank_resets)
  cur <- setNames(rep(NA_character_, length(rank_resets)),
                  names(rank_resets))

  # Lineage class (microbial vs not). Kraken2 rank codes nest strictly along a
  # path (D < D1 < K < P < P2 < ...), so a stack popped back to the row's rank
  # depth holds that row's full lineage. Same classes as the hospital polish
  # figures: Bacteria/Archaea -> prokaryote; Viruses -> virus; Embryophyta ->
  # plant; Metazoa -> animal; other Eukaryota -> microbial_eukaryote (fungi,
  # protists, microalgae); anything else (root, unclassified) -> other.
  rank_letters <- c("U", "R", "D", "K", "P", "C", "O", "F", "G", "S")
  rank_key <- function(code) {
    if (is.na(code) || !nzchar(code)) return(NA_real_)
    pos <- match(substr(code, 1L, 1L), rank_letters)
    sub <- suppressWarnings(as.numeric(substring(code, 2L)))
    if (is.na(pos)) NA_real_ else pos * 100 + ifelse(is.na(sub), 0, sub)
  }
  lin_class <- rep(NA_character_, length(rk))
  stack_key <- numeric(0)
  stack_nm  <- character(0)

  for (i in seq_along(rk)) {
    r <- rk[[i]]
    k <- rank_key(r)
    if (!is.na(k)) {
      keep      <- stack_key < k
      stack_key <- c(stack_key[keep], k)
      stack_nm  <- c(stack_nm[keep], nm[[i]])
      lin_class[[i]] <-
        if (any(stack_nm %in% c("Bacteria", "Archaea"))) "prokaryote"
        else if ("Viruses" %in% stack_nm) "virus"
        else if ("Embryophyta" %in% stack_nm) "plant"
        else if ("Metazoa" %in% stack_nm) "animal"
        else if ("Eukaryota" %in% stack_nm) "microbial_eukaryote"
        else "other"
    }
    first <- if (is.na(r)) NA_character_ else substr(r, 1L, 1L)
    for (rk_name in names(rank_resets)) {
      if (!is.na(first) && first %in% rank_resets[[rk_name]]) {
        cur[[rk_name]] <- NA_character_
      }
      if (!is.na(r) && r == rank_set_code[[rk_name]]) {
        cur[[rk_name]] <- nm[[i]]
      }
      cols[[rk_name]][[i]] <- cur[[rk_name]]
    }
  }
  # `tax_rank` is the row's own kraken2 rank code (S = species, S1 = strain
  # / subspecies, G, F, ...). R/05 keeps rank "S" only so strain rows do not
  # double-count their parent species.
  data.frame(taxid  = td,
             tax_rank = rk,
             lineage_class = lin_class,
             phylum = cols$phylum,
             class  = cols$class,
             order  = cols$order,
             family = cols$family,
             genus  = cols$genus,
             stringsAsFactors = FALSE)
}
