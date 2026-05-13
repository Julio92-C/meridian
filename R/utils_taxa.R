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
    out$name <- sub("^(\\w)\\w+\\s+(\\w+)", "\\1. \\2", out$name)
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
