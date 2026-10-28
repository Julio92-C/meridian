# utils_taxa.R — Generic taxonomic-name helpers shared across plotting stages.
#
# clean_taxa_names(df, cfg) takes a long-form taxa table (must contain a
# character `name` column) and returns the same data frame with:
#   1. higher-level clade-group rows dropped (regex blacklist),
#   2. non-species rows dropped (binomial-only filter),
#   3. NCBI provisional-placement brackets stripped from `name`,
#   4. genus abbreviated to its first letter (`Genus species ...` -> `G. species ...`),
#   5. names exceeding a length cap truncated to their first N words.
#
# All behaviour is driven by `cfg$taxa_display` with built-in defaults so the
# helper is usable on any project without explicit config. Per-project
# overrides go in the YAML, e.g.:
#
#   taxa_display:
#     species_only: true
#     drop_patterns:
#       - "Terrabacteria group"
#       - "FCB group"
#     abbreviate_genus: true
#     max_name_len: 30
#     truncate_to_words: 2

clean_taxa_names <- function(df, cfg = NULL) {
  if (!is.data.frame(df)) {
    stop("clean_taxa_names: `df` must be a data frame")
  }
  if (!"name" %in% colnames(df)) {
    stop("clean_taxa_names: `df` must have a `name` column")
  }

  td <- cfg$taxa_display %||% list()
  species_only      <- td$species_only      %||% TRUE
  abbreviate_genus  <- td$abbreviate_genus  %||% TRUE
  max_name_len      <- td$max_name_len      %||% 30
  truncate_to_words <- td$truncate_to_words %||% 2
  drop_patterns     <- td$drop_patterns     %||% c(
    # NCBI rollup suffixes — informal labels above species rank.
    # Catches "X group", "X clade", "X complex", "X section", "X subgroup",
    # "X serogroup", "X cluster" at end of name (e.g. "Terrabacteria group",
    # "Treponema group", "Mycobacterium tuberculosis complex"). End-anchored
    # so valid species containing these words mid-name are not affected.
    " (group|clade|complex|section|subgroup|serogroup|cluster)$",
    # Non-species placeholder labels that pass the binomial word count
    # filter but are not real species (uncultured isolates, environmental
    # samples, "X bacterium" / "X archaeon" generic names).
    "\\buncultured\\b",
    "\\bunidentified\\b",
    "\\bmetagenome\\b",
    "\\benvironmental\\b",
    "\\bbacterium\\b",
    "\\barchaeon\\b",
    "\\beukaryote\\b",
    "\\bfungus\\b"
  )

  out <- df
  out$name <- as.character(out$name)

  if (length(drop_patterns) > 0) {
    pattern <- paste(drop_patterns, collapse = "|")
    # perl = TRUE so \b word boundaries in default placeholder patterns work.
    out <- out[!grepl(pattern, out$name, perl = TRUE), , drop = FALSE]
  }

  if (isTRUE(species_only)) {
    word_count <- vapply(strsplit(out$name, "\\s+"), length, integer(1))
    out <- out[word_count > 1, , drop = FALSE]
  }

  # Strip NCBI provisional-placement brackets so "[Clostridium] x" abbreviates
  # to "C. x" via the generic rule below, rather than needing per-name cases.
  out$name <- gsub("\\[|\\]", "", out$name)

  if (isTRUE(abbreviate_genus)) {
    full <- out$name
    abbr <- sub("^(\\w)\\w+\\s+(\\w+)", "\\1. \\2", full)
    # Different taxa that share an abbreviation (Acinetobacter soli and
    # Agromyces soli -> "A. soli") would be merged downstream; keep their full
    # genus so every stage still counts them as separate taxa.
    n_full <- vapply(split(full, abbr), function(x) length(unique(x)), integer(1))
    clash <- names(n_full)[n_full > 1]
    out$name <- ifelse(abbr %in% clash, full, abbr)
    substr(out$name, 1, 1) <- toupper(substr(out$name, 1, 1))
  }

  long_idx <- nchar(out$name) > max_name_len
  if (any(long_idx)) {
    trunc_pat <- sprintf("^((?:\\S+\\s+){%d}).*$", truncate_to_words)
    out$name[long_idx] <- sub(trunc_pat, "\\1", out$name[long_idx])
    out$name <- trimws(out$name)
  }

  out
}

# community_filter(df, cfg, stage) — optional community-level filters for the
# taxa stages (R/04 taxonomy, R/05 relative abundance, R/06 alpha, R/07 beta,
# R/08 taxa DA). Off unless `cfg$taxa_filter` is set; gene-based stages do
# not call it. Applied to a long-form (sample, name, count[, tax_rank]) table:
#   taxa_filter.rank             keep rows whose Kraken2 rank code is in this
#                                set (e.g. "S"), using the real `tax_rank`
#                                column from R/02 instead of the binomial
#                                word-count heuristic, so strain / subspecies
#                                rows (S1, S2) are not double-counted against
#                                their species and two-word higher ranks
#                                ("Candidatus Saccharibacteria") do not leak in.
#   taxa_filter.min_sample_reads drop samples whose summed reads at that rank,
#                                over taxa above filters.min_count_per_sample,
#                                fall below N. Near-empty samples otherwise
#                                appear as 100%-one-taxon profiles.
community_filter <- function(df, cfg, stage = "Taxa") {
  tf <- cfg$taxa_filter
  if (is.null(tf) || is.null(df) || nrow(df) == 0) return(df)

  if (!is.null(tf$rank) && "tax_rank" %in% colnames(df)) {
    n0 <- nrow(df)
    rk <- trimws(as.character(df$tax_rank))
    df <- df[!is.na(rk) & rk %in% tf$rank, , drop = FALSE]
    pipeline_log(cfg, sprintf("%s: rank filter (%s) kept %d/%d rows",
                              stage, paste(tf$rank, collapse = ","), nrow(df), n0))
  }

  min_reads <- tf$min_sample_reads
  if (!is.null(min_reads) && min_reads > 0) {
    floor_n <- cfg$filters$min_count_per_sample %||% 0
    depth <- df |>
      dplyr::filter(.data$count > floor_n) |>
      dplyr::group_by(.data$sample) |>
      dplyr::summarise(n = sum(.data$count), .groups = "drop")
    keep <- depth$sample[depth$n >= min_reads]
    drop <- setdiff(unique(df$sample), keep)
    if (length(drop) > 0) {
      pipeline_log(cfg, sprintf(
        "%s: dropped %d sample(s) with < %s reads (taxa_filter.min_sample_reads): %s",
        stage, length(drop), min_reads, paste(sort(drop), collapse = ", ")))
    }
    df <- df[df$sample %in% keep, , drop = FALSE]
  }
  df
}
