# 04_taxa_tree.R — circular taxonomic tree of the most abundant species
# (plus every widespread shared species), built from the NCBI lineage the
# pipeline carries (phylum > class > order > family > genus > species), in
# the style of GraPhlAn / MetaPhlAn cladograms. Branch lengths are not
# evolutionary distances. Tips: colour = group membership, size = mean
# abundance; outer rings = prevalence per group; black dot = reads of the
# species carry an ARG, VF or MGE.
#
# Output: <figures_dir>/taxa_tree/taxa_tree.png, <datasets_dir>/taxa_tree/taxa_tree_nodes.csv

run_taxa_tree <- function(cleaned, cfg) {
  sf_set_figure_formats(cfg)
  pipeline_log(cfg, "Taxonomic tree")
  inp <- sf_prepare_inputs(cleaned, cfg, "Taxa tree", key = "taxa_tree")
  if (is.null(inp)) { pipeline_log(cfg, "Taxa tree: no species — skipped"); return(invisible(NULL)) }
  if (!all(c("phylum", "genus") %in% colnames(inp$long))) {
    pipeline_log(cfg, "Taxa tree: lineage columns missing — skipped"); return(invisible(NULL))
  }
  st <- sf_settings(cfg)
  top_n <- cfg$taxa_tree$top_n %||% 60
  min_samples <- cfg$taxa_tree$min_samples %||% 2
  tab <- sf_species_table(inp$long, inp$group_n, st$min_samples, st$min_prev)
  ws <- tab$name[!is.na(tab$tier) & tab$tier == "widespread"]
  sel <- sf_tree_tips(tab, names(inp$group_n), top_n, min_samples, ws)
  n_rare <- sum(rowSums(as.matrix(tab[paste0("n_", names(inp$group_n))])) < min_samples)
  pipeline_log(cfg, sprintf("Taxa tree: %d species in fewer than %d samples not eligible as tips (taxa_tree.min_samples)",
                            n_rare, min_samples))
  t <- tab[tab$name %in% sel, , drop = FALSE]
  ranks <- intersect(c("phylum", "class", "order", "family", "genus"), colnames(t))
  hits <- cleaned$abri_kraken2
  carriers <- if (!is.null(hits)) unique(gsub("[][]", "", hits$name[!is.na(sf_element_class(hits$DATABASE))])) else character(0)
  t$carrier <- t$name %in% carriers
  groups <- inp$groups; pal <- group_palette(groups, cfg)
  tl <- sf_tree_layout(t, ranks, n_rings = length(groups))
  lay <- tl$layout; leaf <- tl$leaf; t <- tl$table
  cats <- unique(leaf$category); cc <- sf_category_colours(cats, pal); names(cc) <- sf_category_label(cats, groups)
  leaf$cat <- sf_category_label(leaf$category, groups)
  leaf$ws  <- leaf$name %in% ws
  rings <- dplyr::bind_rows(lapply(seq_along(groups), function(i) {
    r <- 1.08 + 0.07 * (i - 1)
    data.frame(x = r * cos(leaf$ang), y = r * sin(leaf$ang), prev = leaf[[paste0("prev_", groups[i])]], group = groups[i])
  }))

  phy_names <- unique(t$phylum)

  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_diagonal(colour = "grey55", width = 0.3) +
    ggplot2::geom_point(data = rings, ggplot2::aes(.data$x, .data$y, alpha = .data$prev, colour = .data$group),
                        shape = 15, size = 2.2, show.legend = FALSE) +
    ggplot2::scale_alpha(range = c(0.05, 1), limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::geom_point(data = leaf, ggplot2::aes(.data$x, .data$y, fill = .data$cat, size = .data$abund_all),
                        shape = 21, colour = "grey20", stroke = 0.3) +
    # Widespread shared species: fixed-size red ring (visible however small the
    # tip) and a red bold label, as in the prevalence-abundance figure.
    ggplot2::geom_point(data = leaf[leaf$ws, ], ggplot2::aes(.data$x, .data$y), shape = 21, fill = NA,
                        colour = SF_TIER_COLOURS[["widespread"]], size = 4, stroke = 0.9) +
    ggplot2::geom_point(data = leaf[leaf$carrier, ], ggplot2::aes(.data$x * 0.96, .data$y * 0.96), size = 0.7) +
    ggplot2::geom_text(data = leaf[!leaf$ws, ], ggplot2::aes(.data$lx, .data$ly, label = sf_short(.data$name), angle = .data$angle,
                                                             hjust = .data$hjust), size = 2.1, fontface = "italic") +
    ggplot2::geom_text(data = leaf[leaf$ws, ], ggplot2::aes(.data$lx, .data$ly, label = sf_short(.data$name), angle = .data$angle,
                                                            hjust = .data$hjust), size = 2.2, fontface = "bold.italic",
                       colour = SF_TIER_COLOURS[["widespread"]]) +
    ggplot2::scale_fill_manual(values = cc, name = NULL) +
    ggplot2::scale_size_area(max_size = 5, name = "Mean abundance (%)") +
    # Room for the rotated tip labels outside the rings, so they do not run into the title.
    ggplot2::coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.75, 1.75), clip = "off") +
    ggplot2::labs(title = "Taxonomic tree of the most abundant and widespread shared species",
                  subtitle = sprintf(paste0("NCBI lineage (not evolutionary distance). Outer rings, inside out: prevalence in %s.\n",
                                            "Red ring and label: widespread shared species (n = %d). Black dot: carries ARG/VF/MGE."),
                                     paste(groups, collapse = ", "), sum(leaf$ws))) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom",
                   plot.margin = ggplot2::margin(10, 10, 10, 10))
  p <- p + sf_tree_phylum_labels(leaf)

  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "taxa_tree")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "taxa_tree")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  sf_save(p, file.path(fig_dir, "taxa_tree.png"), width = 9, height = 9.6, dpi = 300)
  readr::write_csv(t[c("name", ranks, "category", "tier", "abund_all", paste0("prev_", groups), "carrier")],
                   file.path(ds_dir, "taxa_tree_nodes.csv"))
  pipeline_log(cfg, sprintf("Taxa tree: %d species (%d widespread shared), %d phyla", nrow(t), length(ws), length(phy_names)))
  invisible(t)
}
