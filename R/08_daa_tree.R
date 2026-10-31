# 08_daa_tree.R — circular taxonomic tree of the top differentially abundant
# taxa, one per ALDEx2 group pair, from the tables R/08 writes to
# <datasets_dir>/aldex2/taxa/<g1>_vs_<g2>.tsv. Same lineage tree as
# R/04_taxa_tree.R (NCBI lineage, not evolutionary distance).
#   tips   top `daa_tree.top_n` taxa by |effect| (default 30); with
#          daa_tree.balance (default true) top_n/2 from each direction
#   colour group in which the taxon is more abundant (higher rab.win.<group>)
#   size   |ALDEx2 effect|
#   red    significant (wi.eBH, or we.eBH where Wilcoxon was dropped, < alpha)
#   ring   direction of enrichment, opacity scaled by |effect|
#
# Output: <figures_dir>/daa_tree/daa_tree_<pair>.png (+ formats),
#         <datasets_dir>/daa_tree/daa_tree_<pair>.csv

run_daa_tree <- function(cleaned, cfg) {
  sf_set_figure_formats(cfg)
  pipeline_log(cfg, "DAA tree")
  da_dir <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "aldex2", "taxa")
  files <- list.files(da_dir, pattern = "_vs_.*[.]tsv$", full.names = TRUE)
  files <- files[!grepl("^(merged|candidates)", basename(files))]
  if (length(files) == 0) { pipeline_log(cfg, "DAA tree: no ALDEx2 taxa tables — skipped"); return(invisible(NULL)) }
  lin <- daa_tree_lineage(cleaned, cfg)
  if (is.null(lin)) { pipeline_log(cfg, "DAA tree: lineage columns missing — skipped"); return(invisible(NULL)) }
  top_n <- cfg$daa_tree$top_n %||% 30
  alpha <- cfg$differential_abundance$alpha %||% cfg$stats$alpha %||% 0.05
  groups_all <- cfg_group_levels(cfg)
  pal <- group_palette(groups_all, cfg)
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "daa_tree")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "daa_tree")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  # Widespread shared species from the shared_fraction stage (if it ran), so
  # the DAA tree shows which differentially abundant taxa are also widespread
  # in every group.
  sf_csv <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "shared_fraction", "shared_fraction_species.csv")
  ws <- if (file.exists(sf_csv)) {
    s <- readr::read_csv(sf_csv, show_col_types = FALSE); s$name[!is.na(s$tier) & s$tier == "widespread"]
  } else character(0)
  out <- list()
  for (f in files) {
    cn <- sub("[.]tsv$", "", basename(f))
    t <- daa_tree_table(readr::read_tsv(f, show_col_types = FALSE), lin, top_n, alpha,
                        balance = !isFALSE(cfg$daa_tree$balance),
                        min_effect = cfg$daa_tree$min_effect, max_tips = cfg$daa_tree$max_tips %||% 80)
    if (!is.null(t) && !is.na(attr(t, "n_pass") %||% NA) && attr(t, "n_pass") > (cfg$daa_tree$max_tips %||% 80))
      pipeline_log(cfg, sprintf("DAA tree %s: %d taxa with |effect| >= %s; capped at daa_tree.max_tips = %d",
                                cn, attr(t, "n_pass"), cfg$daa_tree$min_effect, cfg$daa_tree$max_tips %||% 80))
    if (is.null(t) || nrow(t) < 3) { pipeline_log(cfg, sprintf("DAA tree %s: fewer than 3 taxa with lineage — skipped", cn)); next }
    t$widespread_shared <- t$name %in% ws
    readr::write_csv(t[setdiff(colnames(t), c("x", "y"))], file.path(ds_dir, paste0("daa_tree_", cn, ".csv")))
    p <- daa_tree_plot(t, attr(t, "groups"), pal, alpha, cn)
    sf_save(p, file.path(fig_dir, paste0("daa_tree_", cn, ".png")), width = 9, height = 9.6)
    um <- attr(t, "unmatched_names")
    pipeline_log(cfg, sprintf("DAA tree %s: %d taxa (%d significant at %.2g), %d unmatched to a lineage%s",
                              cn, nrow(t), sum(t$significant), alpha, attr(t, "unmatched"),
                              if (length(um)) sprintf(" (ambiguous short name: %s)", paste(um, collapse = ", ")) else ""))
    out[[cn]] <- t
  }
  invisible(out)
}

# Display name (as written by R/08, via clean_taxa_names) -> full name and
# lineage. Cleaning the full species list with the same cfg reproduces the
# names R/08 used, including the collision-safe abbreviation.
daa_tree_lineage <- function(cleaned, cfg) {
  nc <- cleaned$noncontaminants
  ranks <- intersect(c("phylum", "class", "order", "family", "genus"), colnames(nc))
  if (!"phylum" %in% ranks) return(NULL)
  u <- nc[!duplicated(nc$name), c("name", ranks), drop = FALSE]
  u$full <- gsub("[][]", "", as.character(u$name))
  cl <- clean_taxa_names(u, cfg)
  cl <- cl[!duplicated(cl$name), , drop = FALSE]
  attr(cl, "ranks") <- ranks
  cl
}

# Top taxa of one ALDEx2 pair with lineage, enrichment and significance.
# balance = TRUE takes the top top_n/2 by |effect| in each direction, so both
# groups are represented even when one side has all the largest effects.
# min_effect (if set) replaces the top-N rule: every taxon with
# |effect| >= min_effect, capped at max_tips (largest effects kept).
daa_tree_table <- function(res, lin, top_n = 30, alpha = 0.05, balance = TRUE,
                           min_effect = NULL, max_tips = 80) {
  feat <- colnames(res)[1]
  rab <- grep("^rab[.]win[.]", colnames(res), value = TRUE)
  if (length(rab) != 2 || !"effect" %in% colnames(res)) return(NULL)
  groups <- sub("^rab[.]win[.]", "", rab)
  pcol <- if ("wi.eBH" %in% colnames(res)) "wi.eBH" else if ("we.eBH" %in% colnames(res)) "we.eBH" else NA
  res <- res[order(-abs(res$effect)), , drop = FALSE]
  n_pass <- NA_integer_
  if (!is.null(min_effect)) {
    res <- res[abs(res$effect) >= min_effect, , drop = FALSE]
    n_pass <- nrow(res)
    res <- res[seq_len(min(max_tips, nrow(res))), , drop = FALSE]
  } else if (balance) {
    side <- ifelse(res[[rab[1]]] >= res[[rab[2]]], 1L, 2L)
    k <- ceiling(top_n / 2)
    res <- res[unlist(lapply(1:2, function(s) utils::head(which(side == s), k))), , drop = FALSE]
    res <- res[order(-abs(res$effect)), , drop = FALSE]
  } else {
    res <- res[seq_len(min(top_n, nrow(res))), , drop = FALSE]
  }
  m <- match(res[[feat]], lin$name)
  # Fallback: R/08 abbreviates a genus when the clash is absent from its own
  # (filtered) taxon set, while the full species list keeps the full genus.
  # Match such names on the abbreviation only when exactly one full name fits.
  if (anyNA(m)) {
    abbr <- sub("^(\\w)\\w+\\s+(\\w+)", "\\1. \\2", lin$full)
    for (i in which(is.na(m))) {
      hit <- which(abbr == res[[feat]][i])
      if (length(hit) == 1) m[i] <- hit
    }
  }
  unmatched <- sum(is.na(m))
  unmatched_names <- res[[feat]][is.na(m)]
  keep <- !is.na(m)
  ranks <- attr(lin, "ranks")
  t <- tibble::tibble(display = res[[feat]][keep], name = lin$full[m[keep]])
  for (r in ranks) t[[r]] <- lin[[r]][m[keep]]
  t$effect <- res$effect[keep]
  t$abs_effect <- abs(t$effect)
  t$enriched <- ifelse(res[[rab[1]]][keep] >= res[[rab[2]]][keep], groups[1], groups[2])
  t$p_adj <- if (is.na(pcol)) NA_real_ else res[[pcol]][keep]
  t$p_col <- pcol
  t$significant <- !is.na(t$p_adj) & t$p_adj < alpha
  attr(t, "groups") <- groups
  attr(t, "unmatched") <- unmatched
  attr(t, "unmatched_names") <- unmatched_names
  attr(t, "ranks") <- ranks
  attr(t, "n_pass") <- n_pass
  attr(t, "min_effect") <- min_effect
  t
}

daa_tree_plot <- function(t, groups, pal, alpha, cn) {
  ranks <- attr(t, "ranks")
  tl <- sf_tree_layout(t, ranks, n_rings = 1)
  leaf <- tl$leaf
  ring <- data.frame(x = 1.08 * cos(leaf$ang), y = 1.08 * sin(leaf$ang), enriched = leaf$enriched,
                     a = leaf$abs_effect / max(leaf$abs_effect))
  pcol <- t$p_col[1]
  ggraph::ggraph(tl$layout) +
    ggraph::geom_edge_diagonal(colour = "grey55", width = 0.3) +
    ggplot2::geom_point(data = ring, ggplot2::aes(.data$x, .data$y, colour = .data$enriched, alpha = .data$a),
                        shape = 15, size = 2.6, show.legend = FALSE) +
    ggplot2::scale_alpha(range = c(0.25, 1), limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::geom_point(data = leaf, ggplot2::aes(.data$x, .data$y, fill = .data$enriched, size = .data$abs_effect),
                        shape = 21, colour = "grey20", stroke = 0.3) +
    ggplot2::geom_point(data = leaf[leaf$significant, ], ggplot2::aes(.data$x, .data$y), shape = 21, fill = NA,
                        colour = SF_TIER_COLOURS[["widespread"]], size = 4.2, stroke = 0.9) +
    # Black triangle between the ring and the label: widespread shared species.
    ggplot2::geom_point(data = leaf[leaf$widespread_shared, ],
                        ggplot2::aes(1.135 * cos(.data$ang), 1.135 * sin(.data$ang)), shape = 17, size = 1.6) +
    ggplot2::geom_text(data = leaf[!leaf$significant, ], ggplot2::aes(.data$lx, .data$ly, label = sf_short(.data$name),
                                                                      angle = .data$angle, hjust = .data$hjust), size = 2.3, fontface = "italic") +
    ggplot2::geom_text(data = leaf[leaf$significant, ], ggplot2::aes(.data$lx, .data$ly, label = sf_short(.data$name),
                                                                     angle = .data$angle, hjust = .data$hjust), size = 2.4,
                       fontface = "bold.italic", colour = SF_TIER_COLOURS[["widespread"]]) +
    ggplot2::scale_fill_manual(values = pal, name = "More abundant in", drop = TRUE) +
    ggplot2::scale_size_area(max_size = 5, name = "|ALDEx2 effect|") +
    ggplot2::coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.75, 1.75), clip = "off") +
    ggplot2::labs(title = sprintf("Differentially abundant taxa: %s vs %s (%s)", groups[1], groups[2],
                                  paste(sprintf("%d more abundant in %s", as.integer(table(factor(t$enriched, groups))), groups), collapse = ", ")),
                  subtitle = sprintf(paste0("NCBI lineage (not evolutionary distance). %s; outer ring: direction, opacity = |effect|.\n",
                                            "Red ring and label: %s < %.2g (n = %d).%s"),
                                     if (is.null(attr(t, "min_effect"))) "Top taxa by |ALDEx2 effect|"
                                     else sprintf("All taxa with |ALDEx2 effect| >= %s", attr(t, "min_effect")),
                                     if (is.na(pcol)) "adjusted p" else pcol, alpha, sum(leaf$significant),
                                     if (any(leaf$widespread_shared)) sprintf(" Black triangle: widespread shared species (n = %d).",
                                                                              sum(leaf$widespread_shared)) else "")) +
    sf_tree_phylum_labels(leaf) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom",
                   plot.margin = ggplot2::margin(10, 10, 10, 10))
}
