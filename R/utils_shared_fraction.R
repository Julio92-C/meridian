# utils_shared_fraction.R — pure helpers for the taxa shared-fraction stage
# (R/04_shared_fraction.R) and the taxonomic tree (R/04_taxa_tree.R).
#
# Works for any number of groups G >= 2. Definitions (see
# docs/superpowers/specs/2026-10-30-taxa-shared-fraction-design.md):
#   n_g      samples of group g where the species is detected
#   prev_g   n_g / samples in g
#   abund_g  mean relative abundance (%) over all samples of g
#   pattern  groups with n_g > 0, joined by "+"
#   category "unique:<g>", "shared:<g1+g2..>" (subset, G > 2) or "shared:all"
#   tier     shared species only, over the groups of its pattern:
#            widespread  n_g >= min_samples and prev_g >= min_prev in every group
#            recurrent   not widespread, n_g >= 2 in every group
#            sporadic    otherwise (one sample in at least one group)

SF_TIERS <- c("widespread", "recurrent", "sporadic")
SF_TIER_COLOURS <- c(widespread = "#CB181D", recurrent = "#FCBBA1", sporadic = "#D9D9D9")

# Threshold settings from cfg$shared_fraction with template defaults.
sf_settings <- function(cfg) {
  s <- cfg$shared_fraction %||% list()
  list(
    min_samples   = s$widespread_min_samples %||% 3,
    min_prev      = s$widespread_min_prev    %||% 0.10,
    heatmap_top_n = s$heatmap_top_n          %||% 30,
    elements_min_reads = s$elements_min_reads %||% 2,
    linkage_vars  = unlist(s$linkage_vars %||% list()),
    linkage_bins  = s$linkage_bins           %||% 3,
    linkage_min_shared = s$linkage_min_shared %||% 1,
    linkage_scope = s$linkage_scope          %||% "between",
    linkage_keep_na = isTRUE(s$linkage_keep_na),
    linkage_test  = isTRUE(s$linkage_test),
    linkage_n_perm = s$linkage_n_perm        %||% 999
  )
}

# Tier of one species given its per-group sample counts and prevalences,
# restricted to the groups where it is present. NA for species in one group.
sf_tier <- function(n, prev, min_samples, min_prev) {
  present <- n > 0
  if (sum(present) < 2) return(NA_character_)
  n <- n[present]; prev <- prev[present]
  if (all(n >= min_samples & prev >= min_prev)) return("widespread")
  if (all(n >= 2)) return("recurrent")
  "sporadic"
}

# Per-species table. `long`: data frame with sample, name, count, group
# (one row per detected species x sample, already floored). `group_n`: named
# integer vector of samples per group (all analysed samples, including any
# with no detected species). Returns one row per species.
sf_species_table <- function(long, group_n, min_samples = 3, min_prev = 0.10) {
  groups <- names(group_n)
  pct <- long |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(pct = .data$count / sum(.data$count) * 100) |>
    dplyr::ungroup()
  by_g <- pct |>
    dplyr::group_by(.data$name, .data$group) |>
    dplyr::summarise(n = dplyr::n_distinct(.data$sample), pct_sum = sum(.data$pct), .groups = "drop")
  out <- tibble::tibble(name = sort(unique(by_g$name)))
  for (g in groups) {
    s <- by_g[by_g$group == g, , drop = FALSE]
    i <- match(out$name, s$name)
    n <- ifelse(is.na(i), 0L, s$n[i])
    out[[paste0("n_", g)]]     <- as.integer(n)
    out[[paste0("prev_", g)]]  <- n / group_n[[g]]
    out[[paste0("abund_", g)]] <- ifelse(is.na(i), 0, s$pct_sum[i] / group_n[[g]])
  }
  nmat <- as.matrix(out[paste0("n_", groups)])
  pmat <- as.matrix(out[paste0("prev_", groups)])
  present <- nmat > 0
  out$n_groups <- rowSums(present)
  out$pattern  <- apply(present, 1, function(p) paste(groups[p], collapse = "+"))
  out$category <- ifelse(out$n_groups == 1, paste0("unique:", out$pattern),
                  ifelse(out$n_groups == length(groups), "shared:all", paste0("shared:", out$pattern)))
  out$tier <- vapply(seq_len(nrow(out)), function(i)
    sf_tier(nmat[i, ], pmat[i, ], min_samples, min_prev), character(1))
  out$prev_all   <- rowSums(nmat) / sum(group_n)
  out$abund_all  <- as.numeric(as.matrix(out[paste0("abund_", groups)]) %*% group_n) / sum(group_n)
  out$min_prev_present <- apply(ifelse(present, pmat, NA), 1, min, na.rm = TRUE)
  lin <- intersect(c("phylum", "class", "order", "family", "genus"), colnames(long))
  if (length(lin) > 0) {
    l1 <- long[!duplicated(long$name), c("name", lin), drop = FALSE]
    out <- dplyr::left_join(out, l1, by = "name")
  }
  out
}

# Pairwise view: for every pair of groups, the species present in both,
# tiered on those two groups only.
sf_pairwise <- function(tab, groups, min_samples = 3, min_prev = 0.10) {
  if (length(groups) < 2) return(NULL)
  prs <- utils::combn(groups, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(prs, function(p) {
    ni <- tab[[paste0("n_", p[1])]]; nj <- tab[[paste0("n_", p[2])]]
    keep <- ni > 0 & nj > 0
    if (!any(keep)) return(NULL)
    t <- tab[keep, , drop = FALSE]
    tibble::tibble(
      pair = paste(p[1], "vs", p[2]), group_i = p[1], group_j = p[2], name = t$name,
      n_i = t[[paste0("n_", p[1])]], n_j = t[[paste0("n_", p[2])]],
      prev_i = t[[paste0("prev_", p[1])]], prev_j = t[[paste0("prev_", p[2])]],
      abund_i = t[[paste0("abund_", p[1])]], abund_j = t[[paste0("abund_", p[2])]],
      tier = vapply(seq_len(nrow(t)), function(k)
        sf_tier(c(t[[paste0("n_", p[1])]][k], t[[paste0("n_", p[2])]][k]),
                c(t[[paste0("prev_", p[1])]][k], t[[paste0("prev_", p[2])]][k]),
                min_samples, min_prev), character(1))
    )
  }))
}

# Element class from an ABRicate DATABASE value.
sf_element_class <- function(db) {
  db <- tolower(as.character(db))
  dplyr::case_when(
    grepl("card|resfinder|ncbi|argannot|megares", db) ~ "ARG",
    grepl("vfdb|victors", db)                         ~ "VF",
    grepl("plasmid|mge|isfinder", db)                 ~ "MGE",
    TRUE                                              ~ NA_character_
  )
}

# Genes carried on the reads of the given species. `hits`: element-carrying
# reads with sample, GENE, DATABASE, name (and optionally sequence = read id).
# Returns one row per species x gene with read and sample counts and a
# low_support flag (reads < min_reads).
sf_elements <- function(hits, species, samples, min_reads = 2) {
  if (is.null(hits) || nrow(hits) == 0) return(NULL)
  h <- hits
  h$name <- gsub("\\[|\\]", "", as.character(h$name))
  h <- h[h$name %in% species & h$sample %in% samples, , drop = FALSE]
  if (nrow(h) == 0) return(NULL)
  h$class <- sf_element_class(h$DATABASE)
  read_id <- if ("sequence" %in% colnames(h)) h$sequence else seq_len(nrow(h))
  h$read_id <- read_id
  h |>
    dplyr::filter(!is.na(.data$class)) |>
    dplyr::group_by(.data$name, .data$class, gene = .data$GENE) |>
    dplyr::summarise(reads = dplyr::n_distinct(.data$read_id),
                     samples = dplyr::n_distinct(.data$sample), .groups = "drop") |>
    dplyr::mutate(low_support = .data$reads < min_reads) |>
    dplyr::arrange(.data$name, .data$class, dplyr::desc(.data$reads))
}

# Levels of a linkage variable: categorical as is; numeric binned into
# quantiles (labels low/mid/high for 3 bins, Q1..Qk otherwise).
sf_bin_variable <- function(x, bins = 3) {
  if (!is.numeric(x)) return(as.character(x))
  ok <- !is.na(x)
  if (length(unique(x[ok])) <= bins) return(as.character(x))
  br <- unique(stats::quantile(x[ok], probs = seq(0, 1, length.out = bins + 1), names = FALSE))
  labs <- if (length(br) - 1 == 3) c("low", "mid", "high") else paste0("Q", seq_len(length(br) - 1))
  as.character(cut(x, breaks = br, labels = labs, include.lowest = TRUE))
}

# Species x node presence and node-pair edges for one linkage variable.
# `long`: sample, name, group (floored); `levels`: named character vector
# sample -> variable level (NA allowed); `species`: species to link (the
# widespread shared species). Nodes are "<group> · <level>".
sf_linkage <- function(long, levels, species, group_of, min_shared = 1,
                       scope = "between", keep_na = FALSE) {
  smp <- names(levels)
  lv <- levels
  if (!keep_na) smp <- smp[!is.na(lv[smp])]
  lv[is.na(lv)] <- "NA"
  node_of <- stats::setNames(paste(group_of[smp], lv[smp], sep = " · "), smp)
  d <- long[long$sample %in% smp & long$name %in% species, c("sample", "name"), drop = FALSE]
  d$node <- node_of[d$sample]
  nodes <- tibble::tibble(node = unique(node_of)) |>
    dplyr::mutate(group = group_of[match(.data$node, node_of)],
                  level = lv[names(node_of)[match(.data$node, node_of)]],
                  n_samples = as.integer(table(node_of)[.data$node]))
  pres <- unique(d[c("node", "name")])
  nodes$n_species <- as.integer(table(factor(pres$node, levels = nodes$node)))
  nodes <- nodes[order(nodes$group, nodes$level), , drop = FALSE]
  sets <- split(pres$name, factor(pres$node, levels = nodes$node))
  if (nrow(nodes) < 2) return(list(nodes = nodes, edges = NULL, presence = pres))
  prs <- utils::combn(nodes$node, 2, simplify = FALSE)
  edges <- dplyr::bind_rows(lapply(prs, function(p) {
    gi <- nodes$group[nodes$node == p[1]]; gj <- nodes$group[nodes$node == p[2]]
    if (scope == "between" && gi == gj) return(NULL)
    if (scope == "within" && gi != gj) return(NULL)
    k <- length(intersect(sets[[p[1]]], sets[[p[2]]]))
    if (k < min_shared) return(NULL)
    tibble::tibble(from = p[1], to = p[2], shared = k)
  }))
  list(nodes = nodes, edges = edges, presence = pres)
}

# Permutation test for one linkage variable: shuffle the variable's levels
# among samples within each group and compare each edge weight (and the total)
# with the null. Returns edges with p_perm and a global p.
sf_linkage_test <- function(long, levels, species, group_of, link, n_perm = 999,
                            seed = 20261005L, ...) {
  if (is.null(link$edges) || nrow(link$edges) == 0) return(link)
  set.seed(seed)
  key <- paste(link$edges$from, link$edges$to)
  obs <- link$edges$shared; obs_tot <- sum(obs)
  ge <- numeric(length(obs)); ge_tot <- 0
  smp <- names(levels)
  for (b in seq_len(n_perm)) {
    perm <- levels
    for (g in unique(group_of[smp])) {
      idx <- smp[group_of[smp] == g]
      perm[idx] <- sample(levels[idx])
    }
    pl <- sf_linkage(long, perm, species, group_of, min_shared = 0, ...)
    k <- if (is.null(pl$edges)) numeric(0) else stats::setNames(pl$edges$shared, paste(pl$edges$from, pl$edges$to))
    kb <- ifelse(is.na(k[key]), 0, k[key])
    ge <- ge + (kb >= obs)
    ge_tot <- ge_tot + (sum(kb) >= obs_tot)
  }
  link$edges$p_perm <- (ge + 1) / (n_perm + 1)
  link$p_global <- (ge_tot + 1) / (n_perm + 1)
  link
}

# Tips of the taxonomic tree: the `top_n` species by mean abundance among those
# detected in at least `min_samples` samples (all groups together), plus every
# widespread shared species. The sample floor keeps a species seen once, in a
# low-depth sample, from entering on its relative abundance alone.
sf_tree_tips <- function(tab, groups, top_n = 60, min_samples = 2, widespread = character(0)) {
  n_tot <- rowSums(as.matrix(tab[paste0("n_", groups)]))
  elig <- tab[n_tot >= min_samples, , drop = FALSE]
  top <- elig$name[order(-elig$abund_all)][seq_len(min(top_n, nrow(elig)))]
  unique(c(widespread, top))
}

# Genera excluded from a stage: `<key>.exclude_genera` (shared_fraction or
# taxa_tree), falling back to `differential_abundance.exclude_genera` so one
# list serves every taxa stage. An explicit empty list (`[]`) disables it.
sf_exclude_genera <- function(cfg, key) {
  own <- cfg[[key]]$exclude_genera
  ex <- if (!is.null(own)) own else cfg$differential_abundance$exclude_genera
  as.character(unlist(ex %||% character(0)))
}

# Genus of a full species name: its first word (same rule as da_taxa_prefilter).
sf_genus_of <- function(name) sub("^\\s*(\\S+).*$", "\\1", as.character(name))

# Figure output for the shared-fraction and tree stages, following the
# pipeline's conventions:
#   - PNG + sibling .rds via save_panel_ggplot(): the PNG is the manifest
#     artefact and the .rds lets R/13 re-render the figure natively inside a
#     composite panel (so any of these figures can be placed in a panel later);
#   - extra siblings from `figures.formats` (tiff, pdf, svg). When unset they
#     follow the composite switches `panels.tiff` / `panels.pdf`, so a project
#     that publishes TIFF gets TIFF here too. TIFF = LZW 300 dpi (as R/13); PDF =
#     cairo_pdf (fonts embedded, text editable in Illustrator); SVG = base cairo.
sf_set_figure_formats <- function(cfg) {
  f <- cfg$figures$formats
  fmts <- if (!is.null(f)) tolower(as.character(unlist(f))) else
    c(if (isTRUE(cfg$panels$tiff)) "tiff", if (isTRUE(cfg$panels$pdf)) "pdf")
  options(meridian.figure_formats = setdiff(fmts, "png"))
}

sf_save <- function(p, file, width, height, dpi = 300) {
  # A PNG held open by an image viewer cannot be overwritten on Windows; warn
  # and still write the other formats rather than failing the stage.
  tryCatch(save_panel_ggplot(file, p, width = width, height = height, dpi = dpi),
           error = function(e) message(sprintf("PNG export failed for %s (file open elsewhere?): %s",
                                               basename(file), conditionMessage(e))))
  stem <- sub("[.]png$", "", file)
  for (fmt in getOption("meridian.figure_formats", character(0))) {
    out <- paste0(stem, ".", fmt)
    tryCatch(switch(fmt,
      tiff = ggplot2::ggsave(out, p, width = width, height = height, dpi = dpi, bg = "white",
                             device = grDevices::tiff, compression = "lzw"),
      pdf  = ggplot2::ggsave(out, p, width = width, height = height, device = grDevices::cairo_pdf),
      svg  = ggplot2::ggsave(out, p, width = width, height = height, device = grDevices::svg),
      message(sprintf("figures.formats: unknown format '%s' ignored", fmt))),
      error = function(e) message(sprintf("%s export failed for %s: %s", fmt, basename(file), conditionMessage(e))))
  }
  invisible(file)
}

# Short display label: "Genus species" -> "G. species".
sf_short <- function(x) sub("^(\\w)\\w+\\s+", "\\1. ", x)
