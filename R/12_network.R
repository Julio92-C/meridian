# 12_network.R — Tripartite (sample × taxon × gene) network + chord + sankey.
# Mirrors the GT scripts tripartiteNetwork.R, chordDiagram_v1.2.1.R,
# chordDiagramNorm_v1.2.1.R, and VFsProfile_sankeyDiagram.R, but drives
# everything from cfg — no hand-tuned row-index cell edits, no hardcoded
# treatments / sample IDs / palettes. Uses clean_taxa_names() from
# R/utils_taxa.R, classify_resistance() from R/09_resistome.R, and
# extract_vf_function() from R/10_virulome.R.
#
# Outputs (under <project>/<figures_dir>/network/):
#   network.png             ggraph FR-layout tripartite network
#   chord_overall.png       Chord of sample → taxon → gene → drug class (all)
#   chord_<group>.png       One chord per group level (toggleable)
#   sankey_overall.html     4-tier sankey sample→taxon→gene→category (all)
#   sankey_<group>.html     One sankey per group level (toggleable)
#
# Datasets (under <project>/<datasets_dir>/network/):
#   gephi_edges.csv         Source / Target / Weight rows for Gephi import
#   gephi_nodes.csv         Id / node_kind / <group> / geneCategory /
#                           Function_group / distance_cluster
#                           (node_kind ∈ {sample, taxa, gene}; kept distinct
#                           from any user metadata column that might happen
#                           to be called `type`.)
#   topology.csv            node / kind / degree / betweenness / module
#   sample_clusters.csv     sample → Bray-Curtis hierarchical cluster id
#   chord_long.csv          (sample, name, GENE, RESISTANCE_class, <group>)
#   sankey_long.csv         (sample, name, GENE, category, <group>)

`%||%` <- function(a, b) if (is.null(a)) b else a

run_network <- function(cleaned, cfg) {
  pipeline_log(cfg, "Network analysis (tripartite + chord)")
  ncfg <- cfg$network %||% list()

  df_all <- .network_load_table(cleaned, cfg)
  if (is.null(df_all) || nrow(df_all) == 0) {
    pipeline_log(cfg, "Network: no usable rows — skipping")
    return(invisible(NULL))
  }

  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "network")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir,  "network")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  group <- cfg$metadata$group_cols[[1]]

  # --- Tripartite network ----------------------------------------------
  min_count_net <- ncfg$min_sample_count %||% 10
  df_net <- df_all |>
    dplyr::filter(.data$sampleCount > min_count_net) |>
    dplyr::distinct()
  if (nrow(df_net) > 0) {
    .network_build_tripartite(df_net, group, ncfg, fig_dir, ds_dir, cfg)
  } else {
    pipeline_log(cfg, sprintf(
      "Network: no rows with sampleCount > %g — tripartite skipped",
      min_count_net
    ))
  }

  # --- Chord diagrams ---------------------------------------------------
  chord_cfg <- ncfg$chord %||% list()
  if (!isFALSE(chord_cfg$enabled %||% TRUE)) {
    .network_build_chord(df_all, group, chord_cfg, fig_dir, ds_dir, cfg)
  }

  # --- Sankey diagrams --------------------------------------------------
  sankey_cfg <- ncfg$sankey %||% list()
  if (!isFALSE(sankey_cfg$enabled %||% TRUE)) {
    .network_build_sankey(df_all, group, sankey_cfg, fig_dir, ds_dir, cfg)
  }

  # --- 4-omics Mantel triangle (PIPELINE_V2_GAPS C10) -------------------
  if (!isFALSE(ncfg$mantel %||% TRUE)) {
    .network_build_mantel_triangle(cleaned, cfg, fig_dir, ds_dir)
  }

  # --- Mobile ARG fraction bar (PIPELINE_V2_GAPS C12) -------------------
  if (!isFALSE(ncfg$mobile_fraction %||% TRUE)) {
    .network_build_mobile_fraction(cleaned, cfg, group, fig_dir, ds_dir)
  }

  # --- Taxon → ARG → MGE sankey (PIPELINE_V2_GAPS C11) ------------------
  if (!isFALSE(ncfg$sankey_taxon_arg_mge %||% TRUE)) {
    .network_build_sankey_taxon_arg_mge(cleaned, cfg, fig_dir, ds_dir)
  }

  invisible(NULL)
}

# ============================================================================
# Data load
# ============================================================================

.network_load_table <- function(cleaned, cfg) {
  src <- cleaned$merged
  if (is.null(src) || nrow(src) == 0) {
    pipeline_log(cfg, "Network: cleaned$merged unavailable — skipping")
    return(NULL)
  }
  needed <- c("sample", "name", "GENE", "DATABASE", "sampleCount")
  miss <- setdiff(needed, colnames(src))
  if (length(miss)) {
    pipeline_log(cfg, sprintf(
      "Network: cleaned$merged missing column(s) %s — skipping",
      paste(miss, collapse = ", ")
    ))
    return(NULL)
  }

  # Gene-label renames are safe for every downstream view, but taxa-name
  # cleaning (species-only filter + genus abbreviation) is intentionally
  # applied *inside* the tripartite and chord builds — the sankey defaults
  # to GT's permissive view that keeps higher-rank taxa rows.
  renames <- c(
    .default_network_gene_renames(),
    as.list(cfg$network$gene_renames %||% list())
  )
  src$GENE <- .apply_substring_renames(src$GENE, renames)
  src
}

.default_network_gene_renames <- function() {
  list(
    "Bifidobacterium_adolescentis_rpoB_mutants_conferring_resistance_to_rifampicin" = "rpoB_mutants",
    "Bifidobacterium_bifidum_ileS_conferring_resistance_to_mupirocin"               = "ileS",
    "vanW_gene_in_vanB_cluster" = "vanWB",
    "vanR_gene_in_vanB_cluster" = "vanRB",
    "vanX_gene_in_vanB_cluster" = "vanXB",
    "vanY_gene_in_vanB_cluster" = "vanYB",
    "vanS_gene_in_vanB_cluster" = "vanSB",
    "vanH_gene_in_vanB_cluster" = "vanHB",
    "vanR_gene_in_vanA_cluster" = "vanRA",
    "vanH_gene_in_vanA_cluster" = "vanHA",
    "vanS_gene_in_vanA_cluster" = "vanSA",
    "Escherichia_coli_ampC_beta-lactamase"       = "ampC",
    "Escherichia_coli_mdfA"                      = "mdfA",
    "Escherichia_coli_emrE"                      = "emrE",
    "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein" = "aac(6')-Ie/aph(2'')-Ia"
  )
}

.apply_substring_renames <- function(x, renames) {
  if (length(renames) == 0) return(x)
  x <- as.character(x)
  for (k in names(renames)) {
    x <- ifelse(grepl(k, x, fixed = TRUE), renames[[k]], x)
  }
  x
}

.network_db_to_category <- function(db) {
  db <- tolower(as.character(db))
  dplyr::case_when(
    grepl("card",          db) ~ "AMR",
    grepl("vfdb",          db) ~ "VFs",
    grepl("plasmidfinder", db) ~ "MGEs",
    TRUE                       ~ db
  )
}

.network_assign_function_group <- function(taxa, groups) {
  out <- rep(NA_character_, length(taxa))
  for (g in names(groups)) {
    hits <- taxa %in% groups[[g]]
    out[hits & is.na(out)] <- g
  }
  out[is.na(out)] <- "Other/unknown"
  out
}

# ============================================================================
# Tripartite network
# ============================================================================

.network_build_tripartite <- function(df, group, ncfg, fig_dir, ds_dir, cfg) {
  df <- clean_taxa_names(df, cfg)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Network: no rows after clean_taxa_names — tripartite skipped")
    return(invisible(NULL))
  }
  # Edge list: sample → taxa (weighted by Σ sampleCount), taxa → gene (binary)
  e_st <- df |>
    dplyr::transmute(from = .data$sample, to = .data$name,
                     weight = .data$sampleCount) |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(weight = sum(.data$weight, na.rm = TRUE),
                     .groups = "drop")
  e_tg <- df |>
    dplyr::distinct(.data$name, .data$GENE) |>
    dplyr::transmute(from = .data$name, to = .data$GENE, weight = 1)
  edges <- dplyr::bind_rows(e_st, e_tg)

  # Node metadata --------------------------------------------------------
  sample_meta <- df |>
    dplyr::distinct(.data$sample, !!rlang::sym(group)) |>
    dplyr::rename(name = .data$sample)

  taxa_meta <- dplyr::distinct(df, .data$name)
  if (!is.null(ncfg$taxa_groups)) {
    taxa_meta$Function_group <- .network_assign_function_group(
      taxa_meta$name, ncfg$taxa_groups
    )
  }

  gene_meta <- df |>
    dplyr::distinct(.data$GENE, .data$DATABASE) |>
    dplyr::mutate(geneCategory = .network_db_to_category(.data$DATABASE)) |>
    dplyr::rename(name = .data$GENE) |>
    dplyr::select(.data$name, .data$geneCategory) |>
    dplyr::distinct(.data$name, .keep_all = TRUE)

  # `node_kind` (not `type`) labels each node as sample/taxa/gene. The
  # column is intentionally NOT called `type` because user metadata may
  # carry a column with that exact name (hospital_microbiome's group_col
  # is `type`) — a name collision in the left_joins below silently turns
  # the kind column into type.x/type.y and breaks every downstream filter.
  nodes <- dplyr::tibble(name = unique(c(edges$from, edges$to))) |>
    dplyr::mutate(node_kind = dplyr::case_when(
      .data$name %in% df$sample ~ "sample",
      .data$name %in% df$name   ~ "taxa",
      TRUE                       ~ "gene"
    )) |>
    dplyr::left_join(sample_meta, by = "name") |>
    dplyr::left_join(taxa_meta,   by = "name") |>
    dplyr::left_join(gene_meta,   by = "name") |>
    dplyr::distinct(.data$name, .keep_all = TRUE)

  # Sample hierarchical clusters (Bray-Curtis on relative abundance) -----
  cluster_k      <- ncfg$cluster_k      %||% 4
  cluster_method <- ncfg$cluster_method %||% "average"
  sample_clusters <- .network_sample_clusters(df, cluster_k, cluster_method)
  if (length(sample_clusters)) {
    nodes$distance_cluster <- sample_clusters[match(nodes$name,
                                                    names(sample_clusters))]
    readr::write_csv(
      data.frame(sample           = names(sample_clusters),
                 distance_cluster = unname(sample_clusters)),
      file.path(ds_dir, "sample_clusters.csv")
    )
  } else {
    nodes$distance_cluster <- NA_integer_
  }

  # Gephi exports --------------------------------------------------------
  edges_gephi <- dplyr::rename(edges, Source = from, Target = to,
                               Weight = weight)
  nodes_gephi <- dplyr::rename(nodes, Id = name)
  for (col in c(group, "geneCategory", "Function_group", "distance_cluster")) {
    if (col %in% colnames(nodes_gephi)) {
      nodes_gephi[[col]][is.na(nodes_gephi[[col]])] <- ""
    }
  }
  readr::write_csv(edges_gephi, file.path(ds_dir, "gephi_edges.csv"))
  readr::write_csv(nodes_gephi, file.path(ds_dir, "gephi_nodes.csv"))

  pipeline_log(cfg, sprintf(
    "Network: %d edges, %d nodes (samples=%d, taxa=%d, genes=%d)",
    nrow(edges), nrow(nodes),
    sum(nodes$node_kind == "sample"),
    sum(nodes$node_kind == "taxa"),
    sum(nodes$node_kind == "gene")
  ))

  # Connectivity Venns (Frontiers Fig S13 panels A/B) -------------------
  # Both are derived from the just-built tripartite network member lists
  # (NOT the raw taxonomy / abricate tables) so they reflect what's
  # actually connected in the graph after the sampleCount > min cut. The
  # taxa Venn uses the primary grouping; the gene-sets Venn uses the
  # AMR/VFs/MGEs partition coming from .network_db_to_category().
  if (group %in% colnames(df)) {
    group_levels <- sort(unique(as.character(df[[group]])))
    pal_group <- resolve_top_level_colors(group, group_levels, cfg) %||%
                 ge_named_palette(group_levels, ncfg$group_colors,
                                   ncfg$group_palette %||% "ggsci::default_nejm")
    ge_plot_category_venn(
      df, category_col = "name", group = group, pal_group = pal_group,
      file = file.path(fig_dir, "connectivity_venn_taxa.png"),
      main_title = sprintf("Network-connected taxa across %s groups", group),
      log_label = "Network connectivity taxa", cfg = cfg,
      value_col = "sampleCount"
    )
  } else {
    pipeline_log(cfg, sprintf(
      "Network connectivity taxa: group column '%s' missing — skipped", group
    ))
  }

  gene_nodes <- nodes |> dplyr::filter(.data$node_kind == "gene",
                                       !is.na(.data$geneCategory))
  if (nrow(gene_nodes) > 0 &&
      dplyr::n_distinct(gene_nodes$geneCategory) >= 2) {
    geneset_levels <- sort(unique(as.character(gene_nodes$geneCategory)))
    pal_geneset <- c(AMR = "#d62728", VFs = "#9467bd", MGEs = "#2ca02c")
    # Fill in any non-canonical categories with hcl colours.
    extras <- setdiff(geneset_levels, names(pal_geneset))
    if (length(extras) > 0) {
      pal_geneset <- c(pal_geneset,
                       setNames(grDevices::hcl.colors(length(extras),
                                                       palette = "Dark 3"),
                                extras))
    }
    gene_nodes$.count <- 1L
    ge_plot_category_venn(
      gene_nodes, category_col = "name", group = "geneCategory",
      pal_group = pal_geneset[geneset_levels],
      file = file.path(fig_dir, "connectivity_venn_genesets.png"),
      main_title = "Network-connected gene elements (ARGs / VFs / MGEs)",
      log_label = "Network connectivity gene-sets", cfg = cfg,
      value_col = ".count"
    )
  } else {
    pipeline_log(cfg,
      "Network connectivity gene-sets: <2 gene categories — skipped")
  }

  # Topology + ggraph plot ----------------------------------------------
  if (!requireNamespace("igraph", quietly = TRUE)) {
    pipeline_log(cfg,
                 "Network: igraph not available — topology/ggraph skipped")
    return(invisible(NULL))
  }
  g <- igraph::graph_from_data_frame(edges, vertices = nodes,
                                     directed = FALSE)
  cluster <- igraph::cluster_louvain(g)
  topo <- data.frame(
    node        = igraph::V(g)$name,
    kind        = igraph::V(g)$node_kind,
    degree      = igraph::degree(g),
    betweenness = igraph::betweenness(g, normalized = TRUE),
    module      = cluster$membership
  )
  readr::write_csv(topo, file.path(ds_dir, "topology.csv"))
  pipeline_log(cfg, sprintf(
    "Network topology: modularity = %.3f",
    igraph::modularity(cluster)
  ))

  # Degree distribution PNG — top-N nodes ranked by degree, coloured by
  # node kind (sample / taxa / gene). Replaces the report's in-page plotly
  # recompute so figures stay pipeline-owned.
  top_n_nodes <- ncfg$degree_distribution_top_n %||% 50
  topo_top <- topo[order(-topo$degree), , drop = FALSE]
  topo_top <- utils::head(topo_top, top_n_nodes)
  dd <- ggplot2::ggplot(topo_top,
        ggplot2::aes(x = stats::reorder(.data$node, -.data$degree),
                     y = .data$degree, fill = .data$kind)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "Degree",
                  title = sprintf("Top %d nodes by degree", nrow(topo_top)),
                  fill = "Node kind") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      text        = ggplot2::element_text(size = 12)
    )
  save_panel_ggplot(file.path(fig_dir, "degree_distribution.png"),
                  dd, width = max(8, 0.18 * nrow(topo_top) + 4),
                  height = 5.5, dpi = 200, bg = "white")

  if (requireNamespace("ggraph", quietly = TRUE)) {
    .network_plot_ggraph(g, nodes, ncfg, fig_dir)
  } else {
    pipeline_log(cfg, "Network: ggraph not available — network.png skipped")
  }
}

.network_sample_clusters <- function(df, k, method) {
  if (!requireNamespace("vegan", quietly = TRUE)) return(integer(0))
  mat <- df |>
    dplyr::select(.data$sample, .data$name, .data$sampleCount) |>
    dplyr::group_by(.data$sample, .data$name) |>
    dplyr::summarise(sampleCount = sum(.data$sampleCount, na.rm = TRUE),
                     .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$name,
                       values_from = .data$sampleCount,
                       values_fill = 0)
  if (nrow(mat) < 2) return(integer(0))
  m <- as.matrix(mat[, -1, drop = FALSE])
  rownames(m) <- mat$sample
  rs <- rowSums(m)
  rs[rs == 0] <- 1
  m <- sweep(m, 1, rs, "/")
  d <- vegan::vegdist(m, method = "bray")
  hc <- hclust(d, method = method)
  k <- max(2, min(k, nrow(m) - 1))
  cutree(hc, k = k)
}

.network_plot_ggraph <- function(g, nodes, ncfg, fig_dir) {
  layout           <- ncfg$layout            %||% "fr"
  pt_size          <- ncfg$node_size         %||% 5
  label_size       <- ncfg$label_size        %||% 3
  label_top_n      <- ncfg$label_top_n       %||% 30
  edge_width_range <- ncfg$edge_width_range  %||% c(0.3, 2.5)
  edge_alpha       <- ncfg$edge_alpha        %||% 0.2
  node_alpha       <- ncfg$node_alpha        %||% 0.75

  # Limit labels to the top-N nodes by degree. With 100+ nodes from the
  # tripartite (sample x taxon x gene) graph, labeling everything packs
  # the centre into illegible mud. Top-N by degree keeps the most
  # connected hubs visible and leaves the periphery as unlabelled dots.
  deg_vec  <- igraph::degree(g)
  n_label  <- min(as.integer(label_top_n), length(deg_vec))
  top_nodes <- if (n_label > 0L) {
    names(sort(deg_vec, decreasing = TRUE))[seq_len(n_label)]
  } else character()
  igraph::V(g)$label_show <- igraph::V(g)$name %in% top_nodes

  col_key <- ifelse(
    nodes$node_kind == "gene",   as.character(nodes$geneCategory),
    ifelse(nodes$node_kind == "sample",
           paste0("cluster_", nodes$distance_cluster),
           "taxa")
  )
  col_key[is.na(col_key) | col_key == "cluster_NA"] <- "taxa"
  igraph::V(g)$col_key <- col_key[match(igraph::V(g)$name, nodes$name)]

  default_pal <- list(
    cluster_1 = "#1f77b4", cluster_2 = "#ff7f0e",
    cluster_3 = "#e9c200", cluster_4 = "#2ca02c",
    cluster_5 = "#9467bd", cluster_6 = "#17becf",
    AMR  = "#d62728", VFs  = "#9467bd", MGEs = "#2ca02c",
    taxa = "#8c564b"
  )
  pal <- unlist(modifyList(default_pal, as.list(ncfg$node_palette %||% list())))

  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = .data$weight), alpha = edge_alpha
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = .data$node_kind, color = .data$col_key),
      size = pt_size, alpha = node_alpha
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = ifelse(.data$label_show, .data$name, "")),
      repel = TRUE, size = label_size, max.overlaps = 20,
      segment.color = "grey60", min.segment.length = 0.5
    ) +
    ggraph::scale_edge_width(range = edge_width_range) +
    ggplot2::scale_color_manual(values = pal, na.value = "grey60",
                                  name  = "Node type") +
    ggplot2::scale_shape_manual(
      values = c(sample = 16, taxa = 17, gene = 15),
      name   = "Node shape"
    ) +
    ggplot2::theme_void() +
    # theme_void() blanks axes, but R/14_panels.R adds axis.text and
    # axis.title back via .panels_composite_theme so the .rds re-render
    # ends up showing x/y coords ("0", "-2.5", ...). Force them blank
    # here so the composite stays clean regardless of theme overrides.
    ggplot2::theme(legend.position = "right",
                   axis.title      = ggplot2::element_blank(),
                   axis.text       = ggplot2::element_blank(),
                   axis.ticks      = ggplot2::element_blank())

  save_panel_ggplot(
    file.path(fig_dir, "network.png"),
    p, width = 14, height = 12, dpi = 200, bg = "white"
  )
}

# ============================================================================
# Chord diagrams
# ============================================================================

.network_build_chord <- function(df_all, group, chord_cfg, fig_dir, ds_dir,
                                  cfg) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    pipeline_log(cfg, "Network chord: circlize not available — skipping")
    return(invisible(NULL))
  }

  database <- tolower(chord_cfg$database %||% "card")
  df <- dplyr::filter(df_all, tolower(.data$DATABASE) %in% database)
  if (nrow(df) == 0 || !"RESISTANCE" %in% colnames(df)) {
    pipeline_log(cfg, sprintf(
      "Network chord: no rows in DATABASE == '%s' — skipping",
      paste(database, collapse = "/")
    ))
    return(invisible(NULL))
  }
  df <- clean_taxa_names(df, cfg)
  if (nrow(df) == 0) {
    pipeline_log(cfg, "Network chord: no rows after clean_taxa_names — skipping")
    return(invisible(NULL))
  }

  min_count  <- chord_cfg$min_sample_count %||% 50
  mls_rollup <- isTRUE(chord_cfg$mls_rollup %||% TRUE)

  df$RESISTANCE_raw <- df$RESISTANCE
  if (mls_rollup) {
    df$RESISTANCE <- stringr::str_replace_all(
      df$RESISTANCE,
      stringr::fixed(
        "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B"
      ),
      "MLS"
    )
  }
  df$RESISTANCE_class <- vapply(df$RESISTANCE, classify_resistance,
                                FUN.VALUE = character(1))
  df$RESISTANCE_class <- ifelse(
    df$RESISTANCE_class %in% c("Mls", "mls"), "MLS", df$RESISTANCE_class
  )

  df <- df |>
    dplyr::filter(!is.na(.data$RESISTANCE_class)) |>
    dplyr::filter(.data$sampleCount > min_count) |>
    dplyr::distinct(.data$sample, .data$name, .data$GENE,
                    .data$RESISTANCE_class, !!rlang::sym(group))
  if (nrow(df) == 0) {
    pipeline_log(cfg, sprintf(
      "Network chord: no rows after sampleCount > %g + class — skipping",
      min_count
    ))
    return(invisible(NULL))
  }

  palette_name <- chord_cfg$palette       %||% "ggsci::default_ucscgb"
  start_degree <- chord_cfg$start_degree  %||% 152

  pipeline_log(cfg, sprintf(
    "Network chord: %d rows, %d samples, %d taxa, %d genes, %d classes",
    nrow(df), dplyr::n_distinct(df$sample), dplyr::n_distinct(df$name),
    dplyr::n_distinct(df$GENE), dplyr::n_distinct(df$RESISTANCE_class)
  ))

  readr::write_csv(df, file.path(ds_dir, "chord_long.csv"))

  .network_render_chord(df, file.path(fig_dir, "chord_overall.png"),
                        palette_name, start_degree)

  if (!isFALSE(chord_cfg$per_group %||% TRUE) && !is.null(group) &&
      group %in% colnames(df)) {
    for (lvl in sort(unique(df[[group]]))) {
      sub <- dplyr::filter(df, .data[[group]] == lvl)
      if (nrow(sub) == 0) next
      safe <- gsub("[^A-Za-z0-9_-]+", "_", lvl)
      .network_render_chord(
        sub, file.path(fig_dir, sprintf("chord_%s.png", safe)),
        palette_name, start_degree
      )
    }
  }
}

.network_render_chord <- function(df, png_path, palette_name, start_degree) {
  samples <- sort(unique(df$sample))
  taxa    <- sort(unique(df$name))
  genes   <- sort(unique(df$GENE))
  classes <- sort(unique(df$RESISTANCE_class))
  cats    <- c(samples, taxa, genes, classes)

  mat <- matrix(0, nrow = length(cats), ncol = length(cats),
                dimnames = list(cats, cats))
  for (i in seq_len(nrow(df))) {
    s <- df$sample[i]; n <- df$name[i]
    g <- df$GENE[i];   r <- df$RESISTANCE_class[i]
    mat[s, n] <- mat[s, n] + 1
    mat[n, g] <- mat[n, g] + 1
    mat[g, r] <- mat[g, r] + 1
  }

  cols <- c(
    .chord_colors(length(samples), palette_name),
    .chord_colors(length(taxa),    palette_name),
    .chord_colors(length(genes),   palette_name),
    .chord_colors(length(classes), palette_name)
  )
  names(cols) <- cats

  grDevices::png(png_path, width = 3200, height = 3200, res = 300, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  circlize::circos.clear()
  # Drop canvas.xlim/ylim — they pin the canvas to the unit circle and
  # clip any label that extends past r=1. circle.margin alone expands
  # the canvas symmetrically (here ~45 % on each side) so even the long
  # taxon / drug-class strings fit. 2026-06-16.
  circlize::circos.par(
    track.height = 0.1, start.degree = start_degree, gap.degree = 2,
    circle.margin = c(0.30, 0.30, 0.30, 0.30),
    unit.circle.segments = 500
  )
  circlize::chordDiagram(
    mat, transparency = 0.5, annotationTrack = "grid",
    scale = FALSE, directional = 1, diffHeight = circlize::mm_h(3),
    grid.col = cols,
    preAllocateTracks = list(track.height = 0.1, unit.circle.segments = 100,
                             start.degree = 90, scale = TRUE)
  )
  circlize::circos.trackPlotRegion(
    track.index = 1, panel.fun = function(x, y) {
      sector.index <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.text(
        circlize::get.cell.meta.data("xcenter"),
        circlize::get.cell.meta.data("ylim")[1],
        sector.index, facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.8), cex = 0.7
      )
    }, bg.border = NA
  )
  circlize::circos.clear()
}

.chord_colors <- function(n, palette_name) {
  if (n == 0) return(character(0))
  base <- tryCatch(
    as.character(paletteer::paletteer_d(palette_name)),
    error = function(e) NULL
  )
  if (is.null(base) || length(base) == 0) {
    base <- grDevices::rainbow(max(8, n))
  }
  if (length(base) >= n) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

# ============================================================================
# Sankey diagrams (4 tiers: sample → taxon → gene → category)
# ============================================================================

.network_build_sankey <- function(df_all, group, scfg, fig_dir, ds_dir, cfg) {
  if (!requireNamespace("networkD3", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    pipeline_log(cfg, "Network sankey: networkD3/htmlwidgets not available — skipping")
    return(invisible(NULL))
  }

  database <- tolower(scfg$database %||% "vfdb")
  df <- dplyr::filter(df_all, tolower(.data$DATABASE) %in% database)
  if (nrow(df) == 0) {
    pipeline_log(cfg, sprintf(
      "Network sankey: no rows in DATABASE == '%s' — skipping",
      paste(database, collapse = "/")
    ))
    return(invisible(NULL))
  }

  if (isTRUE(scfg$clean_taxa %||% FALSE)) {
    df <- clean_taxa_names(df, cfg)
    if (nrow(df) == 0) {
      pipeline_log(cfg, "Network sankey: no rows after clean_taxa_names — skipping")
      return(invisible(NULL))
    }
  }

  category_src <- scfg$category %||% "auto"
  df <- .network_sankey_assign_category(df, database, category_src, scfg)
  if (is.null(df) || nrow(df) == 0) {
    pipeline_log(cfg, "Network sankey: no rows after category assignment — skipping")
    return(invisible(NULL))
  }

  min_count <- scfg$min_sample_count %||% 0
  if (min_count > 0) {
    df <- dplyr::filter(df, .data$sampleCount > min_count)
  }

  # --- Per-gene prevalence + top-N filter -------------------------------
  # Reduces sankey clutter by keeping only the most prevalent / abundant
  # genes BEFORE rendering. Whitelist is computed once on the full df and
  # then applied to per-group sub-sankeys too, so group panels stay
  # comparable (same gene set across panels).
  min_prev   <- scfg$min_gene_prevalence %||% 0
  top_n_gene <- scfg$top_genes           %||% Inf
  rank_by    <- scfg$rank_by             %||% "prevalence_x_abundance"
  if ((min_prev > 0 || is.finite(top_n_gene)) && nrow(df) > 0) {
    total_samples <- dplyr::n_distinct(df$sample)
    gene_stats <- df |>
      dplyr::group_by(.data$GENE) |>
      dplyr::summarise(
        prevalence = dplyr::n_distinct(.data$sample) /
                     max(total_samples, 1L),
        abundance  = sum(.data$sampleCount, na.rm = TRUE),
        .groups    = "drop"
      ) |>
      dplyr::mutate(
        score = dplyr::case_when(
          rank_by == "mean_abundance" ~ .data$abundance,
          rank_by == "prevalence"     ~ .data$prevalence,
          TRUE                         ~ .data$prevalence * .data$abundance
        )
      )
    n_before <- nrow(gene_stats)
    gene_stats <- dplyr::filter(gene_stats, .data$prevalence >= min_prev)
    if (is.finite(top_n_gene)) {
      gene_stats <- gene_stats |>
        dplyr::arrange(dplyr::desc(.data$score)) |>
        dplyr::slice_head(n = as.integer(top_n_gene))
    }
    df <- dplyr::filter(df, .data$GENE %in% gene_stats$GENE)
    pipeline_log(cfg, sprintf(
      "Sankey gene filter: %d -> %d genes (min_prev=%.2f, top_n=%s, rank_by=%s)",
      n_before, nrow(gene_stats), min_prev,
      if (is.finite(top_n_gene)) as.character(top_n_gene) else "Inf",
      rank_by
    ))
  }

  df <- df |>
    dplyr::filter(!is.na(.data$category), nzchar(.data$category)) |>
    dplyr::distinct(.data$sample, .data$name, .data$GENE, .data$category,
                    !!rlang::sym(group))
  if (nrow(df) == 0) {
    pipeline_log(cfg, sprintf(
      "Network sankey: no rows after filters (sampleCount > %g) — skipping",
      min_count
    ))
    return(invisible(NULL))
  }

  pipeline_log(cfg, sprintf(
    "Network sankey: %d rows, %d samples, %d taxa, %d genes, %d categories",
    nrow(df), dplyr::n_distinct(df$sample), dplyr::n_distinct(df$name),
    dplyr::n_distinct(df$GENE), dplyr::n_distinct(df$category)
  ))
  readr::write_csv(df, file.path(ds_dir, "sankey_long.csv"))

  .network_render_sankey(df, file.path(fig_dir, "sankey_overall.html"), scfg)
  # Companion static PNG for the manuscript (Frontiers Fig S9). The HTML
  # stays for the supplementary website; this PNG is a ggalluvial render
  # of the same 4-tier data so the figure can be embedded in a journal
  # PDF without requiring headless-Chrome.
  .network_render_sankey_png(
    df, file.path(fig_dir, "sankey_overall.png"), scfg, cfg
  )

  if (!isFALSE(scfg$per_group %||% TRUE) && !is.null(group) &&
      group %in% colnames(df)) {
    for (lvl in sort(unique(df[[group]]))) {
      sub <- dplyr::filter(df, .data[[group]] == lvl)
      if (nrow(sub) == 0) next
      safe <- gsub("[^A-Za-z0-9_-]+", "_", lvl)
      .network_render_sankey(
        sub, file.path(fig_dir, sprintf("sankey_%s.html", safe)), scfg
      )
    }
  }
}

# ggalluvial-based static PNG render of the 4-tier sankey. Skipped (with a
# log line) when ggalluvial isn't installed — keeps the HTML render path
# unaffected. Style is intentionally muted: strata are white boxes with
# grey borders + labels; alluvia are alpha-blended bands coloured by the
# rightmost tier (category) so flows are easy to trace. Width/height scale
# with cardinality to keep labels legible across study sizes.
.network_render_sankey_png <- function(df, png_path, scfg, cfg) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    pipeline_log(cfg, "Network sankey PNG: ggalluvial not available — skipping")
    return(invisible(NULL))
  }

  df_alluv <- df |>
    dplyr::distinct(.data$sample, .data$name, .data$GENE, .data$category) |>
    dplyr::mutate(value = 1L)
  if (nrow(df_alluv) == 0) {
    pipeline_log(cfg, "Network sankey PNG: empty alluvial frame — skipping")
    return(invisible(NULL))
  }

  # ggalluvial sizes each stratum proportionally to its flow count, so the
  # smallest sample (e.g. one passing through a single category) gets a
  # stratum of height ~1/nrow(df_alluv). Scaling the plot height by total
  # flow count (not max-axis cardinality) gives those small strata enough
  # room for their labels to clear the neighbours.
  total_flows <- nrow(df_alluv)
  ph <- scfg$png_height %||% max(10, 0.07 * total_flows + 5)
  pw <- scfg$png_width  %||% 14
  label_size <- scfg$png_label_size %||% 2.0

  p <- ggplot2::ggplot(
        df_alluv,
        ggplot2::aes(axis1 = .data$sample, axis2 = .data$name,
                     axis3 = .data$GENE,   axis4 = .data$category,
                     y     = .data$value)
      ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = .data$category), alpha = 0.5, width = 1/8
    ) +
    ggalluvial::geom_stratum(width = 1/8, fill = "white", colour = "grey40") +
    # ggalluvial registers the "stratum" stat on package load. When we
    # only `requireNamespace()` the package (not library()), bare
    # stat = "stratum" isn't in ggplot2's registry — reference the
    # StatStratum ggproto object directly so it resolves regardless.
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = label_size
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Sample", "Taxon", "Gene", "Category"),
      # Extra right margin so long category labels like
      # "Antimicrobial activity/Competitive advantage" don't clip the
      # canvas edge after the Category stratum.
      expand = ggplot2::expansion(mult = c(0.05, 0.20))
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Category") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box.just = "left",
      text            = ggplot2::element_text(size = 11)
    ) +
    # Wrap categories onto 2 rows so a long label like
    # "Antimicrobial activity/Competitive advantage" doesn't push the
    # last entry off-canvas in the composite cell.
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  save_panel_ggplot(png_path, p, width = pw, height = ph, dpi = 300, bg = "white")
  pipeline_log(cfg, sprintf("Network sankey PNG: %s", basename(png_path)))
  invisible(NULL)
}

# Resolve the GENE → category mapping for a given database. "auto" picks
# VF-function for vfdb, drug-class for card, and the GENE itself otherwise
# (so the sankey degenerates into sample → taxon → gene → gene). Explicit
# values: "vf_function", "drug_class", "gene".
.network_sankey_assign_category <- function(df, database, category_src, scfg) {
  src <- category_src
  if (src == "auto") {
    src <- dplyr::case_when(
      "vfdb" %in% database ~ "vf_function",
      "card" %in% database ~ "drug_class",
      TRUE                  ~ "gene"
    )
  }

  if (src == "vf_function") {
    if (!"PRODUCT" %in% colnames(df)) return(NULL)
    df$category <- vapply(df$PRODUCT, extract_vf_function,
                          FUN.VALUE = character(1))
    renames <- as.list(scfg$category_renames %||% list())
    if (length(renames) > 0) {
      df$category <- .apply_literal_renames(df$category, renames)
    }
  } else if (src == "drug_class") {
    if (!"RESISTANCE" %in% colnames(df)) return(NULL)
    if (isTRUE(scfg$mls_rollup %||% TRUE)) {
      df$RESISTANCE <- stringr::str_replace_all(
        df$RESISTANCE,
        stringr::fixed(
          "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B"
        ),
        "MLS"
      )
    }
    df$category <- vapply(df$RESISTANCE, classify_resistance,
                          FUN.VALUE = character(1))
    df$category <- ifelse(df$category %in% c("Mls", "mls"), "MLS", df$category)
  } else {
    df$category <- df$GENE
  }
  df
}

# Best-effort literal-substring rename (mirrors R/10's .apply_literal_renames
# but local — keep R/12 self-contained).
.apply_literal_renames <- function(x, renames) {
  if (length(renames) == 0) return(x)
  x <- as.character(x)
  for (k in names(renames)) {
    x <- ifelse(grepl(k, x, fixed = TRUE), renames[[k]], x)
  }
  x
}

.network_render_sankey <- function(df, html_path, scfg) {
  node_cols <- modifyList(
    list(sample = "#175709", taxon = "#825cdb",
         gene   = "#fc3503", class = "#b5b5b5"),
    as.list(scfg$node_colors %||% list())
  )
  link_cols <- modifyList(
    list(sample_taxon = "#98ed85", taxon_gene = "#4fc1e3",
         gene_class   = "#d1d0c24D"),
    as.list(scfg$link_colors %||% list())
  )

  # ---- Build node table -----------------------------------------------
  samples <- sort(unique(df$sample))
  taxa    <- sort(unique(df$name))
  genes   <- sort(unique(df$GENE))
  cats    <- sort(unique(df$category))

  nodes <- dplyr::bind_rows(
    data.frame(name = samples, tier = "sample", color = node_cols$sample,
               stringsAsFactors = FALSE),
    data.frame(name = taxa,    tier = "taxon",  color = node_cols$taxon,
               stringsAsFactors = FALSE),
    data.frame(name = genes,   tier = "gene",   color = node_cols$gene,
               stringsAsFactors = FALSE),
    data.frame(name = cats,    tier = "class",  color = node_cols$class,
               stringsAsFactors = FALSE)
  )
  nodes$id <- seq_len(nrow(nodes)) - 1L  # networkD3 is 0-indexed

  idx <- setNames(nodes$id, nodes$name)

  # ---- Build link table (three tiers of edges) ------------------------
  l0 <- df |>
    dplyr::distinct(.data$sample, .data$name) |>
    dplyr::transmute(source = idx[.data$sample], target = idx[.data$name],
                     color = link_cols$sample_taxon)
  l1 <- df |>
    dplyr::distinct(.data$name, .data$GENE) |>
    dplyr::transmute(source = idx[.data$name], target = idx[.data$GENE],
                     color = link_cols$taxon_gene)
  l2 <- df |>
    dplyr::distinct(.data$GENE, .data$category) |>
    dplyr::transmute(source = idx[.data$GENE], target = idx[.data$category],
                     color = link_cols$gene_class)
  links <- dplyr::bind_rows(l0, l1, l2) |>
    dplyr::distinct(.data$source, .data$target, .keep_all = TRUE)
  links$value <- 1

  # ---- Render and save -------------------------------------------------
  width     <- scfg$width      %||% 1500
  height    <- scfg$height     %||% 750
  font_size <- scfg$font_size  %||% 14
  font_fam  <- scfg$font_family %||% "arial"

  sk <- networkD3::sankeyNetwork(
    Links = as.data.frame(links), Nodes = as.data.frame(nodes),
    Source = "source", Target = "target", Value = "value", NodeID = "name",
    NodeGroup = "tier", LinkGroup = "color",
    fontSize = font_size, fontFamily = font_fam,
    width = width, height = height, sinksRight = TRUE,
    margin = list(top = 10, right = 10, bottom = 10, left = 10)
  )

  selfcontained <- isTRUE(scfg$selfcontained %||% TRUE)
  # As with R/05 / R/07: saveWidget(selfcontained=TRUE) still leaves a
  # benign <name>_files/ directory after rendering; suppressWarnings keeps
  # the run log clean on re-runs. Falls through if pandoc is missing.
  suppressWarnings(
    htmlwidgets::saveWidget(sk, normalizePath(html_path, mustWork = FALSE),
                            selfcontained = selfcontained)
  )
  invisible(NULL)
}

# ============================================================================
# 4-omics Mantel triangle (PIPELINE_V2_GAPS C10)
# ============================================================================

# Pairwise Mantel correlations between Bray-Curtis distance matrices for
# the four omics layers (taxonomy / resistome / virulome / mobilome). The
# cells display Mantel r + significance stars; colour intensity = sign+|r|.
# Domains without enough samples or features are dropped silently with a
# log line; if fewer than 2 survive, the whole figure is skipped.
.network_build_mantel_triangle <- function(cleaned, cfg, fig_dir, ds_dir) {
  if (!requireNamespace("vegan", quietly = TRUE) ||
      !requireNamespace("pheatmap", quietly = TRUE)) {
    pipeline_log(cfg, "Mantel triangle: vegan/pheatmap missing — skipping")
    return(invisible(NULL))
  }

  # Level switch (PIPELINE_V2_GAPS C10 follow-up). Default "category"
  # — Bray-Curtis at the ecologically-meaningful rollup (genus / drug
  # class / VF function / replicon family) so the 4 domains are scale-
  # matched and noise from rare features is averaged out. "feature"
  # falls back to the gene/species level original behaviour.
  level <- cfg$network$mantel$level %||% "category"
  if (!level %in% c("category", "feature")) {
    pipeline_log(cfg, sprintf(
      "Mantel triangle: unknown level '%s' — falling back to 'category'",
      level
    ))
    level <- "category"
  }

  # Helper: long-form (sample, id, val) -> Bray-Curtis dist on sample x id.
  build_dist <- function(df, id_col, val_col, label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    wide <- df |>
      dplyr::filter(!is.na(.data[[id_col]])) |>
      dplyr::group_by(.data$sample, .data[[id_col]]) |>
      dplyr::summarise(v = sum(.data[[val_col]], na.rm = TRUE),
                       .groups = "drop") |>
      tidyr::pivot_wider(names_from = dplyr::all_of(id_col),
                         values_from = "v", values_fill = 0)
    if (nrow(wide) < 3) return(NULL)
    mat <- as.matrix(wide[, -1, drop = FALSE])
    rownames(mat) <- as.character(wide$sample)
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    if (nrow(mat) < 3) {
      pipeline_log(cfg, sprintf(
        "Mantel %s: only %d non-empty samples — domain skipped",
        label, nrow(mat)
      ))
      return(NULL)
    }
    vegan::vegdist(mat, method = "bray")
  }

  dists <- list()

  # ---- Taxonomy --------------------------------------------------------
  if (level == "category") {
    # Genus rollup uses the kraken2-derived `genus` column attached in
    # R/02 (`build_taxid_ancestry`). Higher-rank rows (phylum / class /
    # order / family) and unclassified rows have genus = NA and are
    # dropped, so the Bray-Curtis distance is computed on true genera
    # only — matches the C4 heatmap recipe in R/05.
    tx <- cleaned$noncontaminants
    if (!is.null(tx) && nrow(tx) > 0 && "genus" %in% colnames(tx)) {
      tx <- dplyr::filter(tx, !is.na(.data$genus), nzchar(.data$genus))
      dists$Taxonomy <- build_dist(tx, "genus", "count", "Taxonomy")
    }
  } else {
    dists$Taxonomy <- build_dist(cleaned$noncontaminants, "name", "count",
                                  "Taxonomy")
  }

  if (level == "category") {
    # Resistome / Virulome / Mobilome from the per-sample category-TPM
    # CSVs written by R/09 / R/10 / R/11. Mantel-comparable scale: ~20
    # drug classes vs ~5-10 VF functions vs ~3-5 replicon families.
    ds_root <- file.path(cfg$project_root, cfg$outputs$datasets_dir)
    cat_sources <- list(
      list(label = "Resistome", path = "resistome/drug_class_per_sample_TPM.csv",
           id = "DRUG"),
      list(label = "Virulome",  path = "virulome/vf_function_per_sample_TPM.csv",
           id = "Functions"),
      list(label = "Mobilome",  path = "mobilome/replicon_family_per_sample_TPM.csv",
           id = "Replicon_Family")
    )
    for (src in cat_sources) {
      f <- file.path(ds_root, src$path)
      if (!file.exists(f)) {
        pipeline_log(cfg, sprintf(
          "Mantel %s: %s not found — domain skipped (run R/09-11 first)",
          src$label, src$path
        ))
        next
      }
      df <- readr::read_csv(f, show_col_types = FALSE)
      dists[[src$label]] <- build_dist(df, src$id, "TPM", src$label)
    }
  } else {
    norm <- ge_load_norm_table(cfg, "Mantel triangle")
    if (!is.null(norm)) {
      db_map <- list(Resistome = "card", Virulome = "vfdb",
                     Mobilome  = "plasmidfinder")
      for (label in names(db_map)) {
        sub <- norm[tolower(norm$DATABASE) == db_map[[label]], , drop = FALSE]
        dists[[label]] <- build_dist(sub, "GENE", "TPM", label)
      }
    }
  }

  dists <- Filter(Negate(is.null), dists)
  if (length(dists) < 2) {
    pipeline_log(cfg,
      "Mantel triangle: <2 domains with usable data — skipping")
    return(invisible(NULL))
  }

  dom_names <- names(dists)
  k     <- length(dom_names)
  r_mat <- matrix(NA_real_, k, k, dimnames = list(dom_names, dom_names))
  p_mat <- r_mat
  perms <- cfg$stats$permanova_permutations %||% 9999

  for (i in seq_len(k - 1)) {
    for (j in seq.int(i + 1, k)) {
      d1 <- dists[[i]]; d2 <- dists[[j]]
      s1 <- attr(d1, "Labels"); s2 <- attr(d2, "Labels")
      common <- intersect(s1, s2)
      if (length(common) < 3) {
        pipeline_log(cfg, sprintf(
          "Mantel %s vs %s: only %d shared samples — pair skipped",
          dom_names[i], dom_names[j], length(common)
        ))
        next
      }
      d1_sub <- stats::as.dist(as.matrix(d1)[common, common])
      d2_sub <- stats::as.dist(as.matrix(d2)[common, common])
      res <- tryCatch(
        vegan::mantel(d1_sub, d2_sub, method = "spearman",
                       permutations = perms),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        r_mat[i, j] <- res$statistic
        r_mat[j, i] <- res$statistic
        p_mat[i, j] <- res$signif
        p_mat[j, i] <- res$signif
        pipeline_log(cfg, sprintf(
          "Mantel %s vs %s: r = %.3f, p = %.4g (n=%d samples)",
          dom_names[i], dom_names[j], res$statistic, res$signif,
          length(common)
        ))
      }
    }
  }
  diag(r_mat) <- 1

  # Family-wide multiple-testing correction for the upper-triangle pairs
  # (one test per omics pair, 6 tests for k = 4 layers). Mantel p-values
  # for omics layers built from the same samples are clearly dependent —
  # cfg$stats$padjust_method drives this (default BH; BY is the
  # arbitrary-dependence-safe choice if you want to be conservative).
  pad_method <- padjust_method(cfg)
  upper_idx  <- which(upper.tri(p_mat))
  raw_p      <- p_mat[upper_idx]
  if (any(!is.na(raw_p))) {
    adj_p <- padjust_p(raw_p, cfg)
    padj_mat <- p_mat
    padj_mat[upper_idx] <- adj_p
    padj_mat[lower.tri(padj_mat)] <- t(padj_mat)[lower.tri(padj_mat)]
  } else {
    padj_mat <- p_mat
  }

  # Save both the r matrix and the (raw + adjusted) pair-wise p-values.
  out_df <- cbind(data.frame(domain = dom_names),
                  as.data.frame(round(r_mat, 4)))
  readr::write_csv(out_df,
                    file.path(ds_dir, "mantel_correlation_triangle.csv"))
  if (length(upper_idx) > 0) {
    pair_idx <- which(upper.tri(p_mat), arr.ind = TRUE)
    pairs_df <- data.frame(
      domain_a = dom_names[pair_idx[, 1]],
      domain_b = dom_names[pair_idx[, 2]],
      mantel_r = round(r_mat[upper_idx], 4),
      p_raw    = signif(p_mat[upper_idx], 4),
      p_adj    = signif(padj_mat[upper_idx], 4),
      method   = pad_method,
      stringsAsFactors = FALSE
    )
    readr::write_csv(pairs_df,
                      file.path(ds_dir, "mantel_pairwise_padj.csv"))
  }

  stars <- ifelse(is.na(padj_mat), "",
    ifelse(padj_mat < 0.001, "***",
    ifelse(padj_mat < 0.01,  "**",
    ifelse(padj_mat < 0.05,  "*", ""))))
  labels <- ifelse(is.na(r_mat), "",
                   paste0(sprintf("%.2f", r_mat), stars))
  diag(labels) <- "—"

  pheatmap::pheatmap(
    r_mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    color  = grDevices::colorRampPalette(
      c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = labels,
    number_color    = "black",
    fontsize_number = 10,
    border_color    = "grey70",
    fontsize_row = 10, fontsize_col = 10, fontsize = 10,
    main = sprintf(
      "Mantel triangle (Spearman, %s-level, %d perms, p adj %s)",
      level, perms, pad_method
    ),
    filename = file.path(fig_dir, "mantel_correlation_triangle.png"),
    width = 7, height = 6
  )
  pipeline_log(cfg, sprintf(
    "Mantel triangle: %d domains x %d, %d pairs computed, p adjusted via %s",
    k, k, sum(!is.na(r_mat[upper.tri(r_mat)])), pad_method
  ))
  invisible(NULL)
}

# ============================================================================
# Mobile ARG fraction (PIPELINE_V2_GAPS C12)
# ============================================================================

# Per-treatment stacked bar of % ARG TPM that's "mobile". A contig is
# considered mobile when it carries both a CARD hit and a PlasmidFinder
# hit in the same sample. An ARG gene's TPM is then assigned to the
# mobile bucket if any of its contigs in that sample are mobile.
# Approximation: TPM is per-(sample, GENE), not per-(sample, contig, GENE),
# so a gene with multiple contigs (some mobile, some not) contributes its
# full TPM to "mobile" if ANY of them are mobile. Defensible upper bound
# on the mobile fraction.
.network_build_mobile_fraction <- function(cleaned, cfg, group,
                                            fig_dir, ds_dir) {
  if (is.null(cleaned$abri_kraken2)) {
    pipeline_log(cfg,
      "Mobile fraction: abri_kraken2 missing — skipping")
    return(invisible(NULL))
  }
  required <- c("sample", "sequence", "DATABASE", "GENE")
  miss <- setdiff(required, colnames(cleaned$abri_kraken2))
  if (length(miss) > 0) {
    pipeline_log(cfg, sprintf(
      "Mobile fraction: abri_kraken2 missing column(s) %s — skipping",
      paste(miss, collapse = ", ")
    ))
    return(invisible(NULL))
  }

  abri <- cleaned$abri_kraken2
  contig_db <- abri |>
    dplyr::filter(!is.na(.data$DATABASE), !is.na(.data$sequence)) |>
    dplyr::mutate(db = tolower(.data$DATABASE)) |>
    dplyr::distinct(.data$sample, .data$sequence, .data$db)

  mobile_contigs <- contig_db |>
    dplyr::group_by(.data$sample, .data$sequence) |>
    dplyr::summarise(
      has_card = any(.data$db == "card"),
      has_pf   = any(.data$db == "plasmidfinder"),
      .groups  = "drop"
    ) |>
    dplyr::filter(.data$has_card & .data$has_pf) |>
    dplyr::select("sample", "sequence")

  if (nrow(mobile_contigs) == 0) {
    pipeline_log(cfg,
      "Mobile fraction: no contigs co-harbouring CARD + PlasmidFinder — skipping")
    return(invisible(NULL))
  }

  card_hits <- abri |>
    dplyr::filter(tolower(.data$DATABASE) == "card",
                  !is.na(.data$sequence), !is.na(.data$GENE)) |>
    dplyr::select("sample", "GENE", "sequence") |>
    dplyr::left_join(
      dplyr::mutate(mobile_contigs, is_mobile = TRUE),
      by = c("sample", "sequence")
    ) |>
    dplyr::mutate(is_mobile = !is.na(.data$is_mobile)) |>
    dplyr::group_by(.data$sample, .data$GENE) |>
    dplyr::summarise(is_mobile = any(.data$is_mobile), .groups = "drop")

  norm <- ge_load_norm_table(cfg, "Mobile fraction")
  if (is.null(norm)) return(invisible(NULL))
  card_tpm <- norm |>
    dplyr::filter(tolower(.data$DATABASE) == "card") |>
    dplyr::select("sample", "GENE", "TPM")

  per_sample <- dplyr::inner_join(card_tpm, card_hits,
                                   by = c("sample", "GENE")) |>
    dplyr::group_by(.data$sample) |>
    dplyr::summarise(
      mobile_TPM     = sum(.data$TPM[.data$is_mobile],  na.rm = TRUE),
      non_mobile_TPM = sum(.data$TPM[!.data$is_mobile], na.rm = TRUE),
      total_TPM      = sum(.data$TPM, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(mobile_pct = ifelse(.data$total_TPM > 0,
                                       100 * .data$mobile_TPM / .data$total_TPM,
                                       0))

  if (nrow(per_sample) == 0) {
    pipeline_log(cfg, "Mobile fraction: no per-sample TPM rows — skipping")
    return(invisible(NULL))
  }

  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE)
  per_sample <- dplyr::inner_join(per_sample,
                                   meta[, c(sid, group)],
                                   by = c("sample" = sid))
  per_sample[[group]] <- as.character(per_sample[[group]])

  readr::write_csv(per_sample,
                    file.path(ds_dir, "mobile_arg_fraction_per_sample.csv"))

  kw <- tryCatch(
    kruskal.test(stats::reformulate(group, "mobile_pct"), data = per_sample),
    error = function(e) NULL
  )
  kw_label <- if (!is.null(kw))
    sprintf("Kruskal-Wallis p = %.3g", kw$p.value)
  else "Kruskal-Wallis p = NA"

  per_group <- per_sample |>
    dplyr::group_by(.data[[group]]) |>
    dplyr::summarise(
      mobile_pct_mean     = mean(.data$mobile_pct, na.rm = TRUE),
      non_mobile_pct_mean = 100 - mean(.data$mobile_pct, na.rm = TRUE),
      .groups = "drop"
    )
  # Use `mobile_class` (not `type`) for the pivot's names column — when
  # the user's group_col is itself `type` (hospital_microbiome), a `type`
  # names_to collides with the already-present group column and trips
  # pivot_longer's duplicate-name guard.
  per_group_long <- tidyr::pivot_longer(
    per_group,
    cols      = c("mobile_pct_mean", "non_mobile_pct_mean"),
    names_to  = "mobile_class",
    values_to = "pct"
  )
  per_group_long$mobile_class <- ifelse(
    per_group_long$mobile_class == "mobile_pct_mean", "Mobile", "Non-mobile"
  )
  per_group_long$mobile_class <- factor(per_group_long$mobile_class,
                                         levels = c("Non-mobile", "Mobile"))

  p <- ggplot2::ggplot(per_group_long,
        ggplot2::aes(x = .data[[group]], y = .data$pct,
                     fill = .data$mobile_class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`Non-mobile` = "#4575b4",
                                           Mobile      = "#d73027")) +
    ggplot2::labs(
      x = group, y = "% of total ARG TPM",
      fill = NULL,
      title = sprintf("ARG mobile fraction by treatment (%s)", kw_label)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "top",
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      text        = ggplot2::element_text(size = 12)
    )
  save_panel_ggplot(
    file.path(fig_dir, "mobile_arg_fraction_bar.png"),
    p,
    width  = max(7, 1.2 * dplyr::n_distinct(per_group_long[[group]]) + 4),
    height = 6, dpi = 300
  )
  pipeline_log(cfg, sprintf(
    "Mobile ARG fraction: %d samples, %d treatment groups, %d mobile contigs",
    nrow(per_sample), dplyr::n_distinct(per_sample[[group]]),
    nrow(mobile_contigs)
  ))
  invisible(NULL)
}

# ============================================================================
# Taxon → ARG → MGE sankey (PIPELINE_V2_GAPS C11)
# ============================================================================
# Three-axis alluvial flow: bacterial phylum (left) → ARG drug class (centre)
# → MGE replicon family (right). Uses the same CARD + PlasmidFinder
# contig-co-occurrence join as `.network_build_mobile_fraction` to identify
# mobile ARGs, then groups by the phylum the contig was assigned to (via
# kraken2 taxid, carried on `cleaned$abri_kraken2` after the R/02 phylum
# join). When `ge_load_norm_table` is available the ribbon weight is CARD
# TPM summed per (phylum, drug, replicon); otherwise it falls back to the
# Bracken sampleCount. Skips with a log line if no co-occurring contigs
# exist (chicken_batch1 currently triggers this — same condition as
# mobile_fraction).
.network_build_sankey_taxon_arg_mge <- function(cleaned, cfg, fig_dir, ds_dir) {
  if (is.null(cleaned$abri_kraken2)) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: abri_kraken2 missing — skipping")
    return(invisible(NULL))
  }
  required <- c("sample", "sequence", "DATABASE", "GENE",
                "RESISTANCE", "phylum")
  miss <- setdiff(required, colnames(cleaned$abri_kraken2))
  if (length(miss) > 0) {
    pipeline_log(cfg, sprintf(
      "Taxon→ARG→MGE sankey: abri_kraken2 missing column(s) %s — skipping",
      paste(miss, collapse = ", ")
    ))
    return(invisible(NULL))
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: ggalluvial not available — skipping")
    return(invisible(NULL))
  }

  abri <- cleaned$abri_kraken2

  # 1. Mobile contigs (same definition as mobile_fraction_bar).
  contig_db <- abri |>
    dplyr::filter(!is.na(.data$DATABASE), !is.na(.data$sequence)) |>
    dplyr::mutate(db = tolower(.data$DATABASE)) |>
    dplyr::distinct(.data$sample, .data$sequence, .data$db)

  mobile_contigs <- contig_db |>
    dplyr::group_by(.data$sample, .data$sequence) |>
    dplyr::summarise(
      has_card = any(.data$db == "card"),
      has_pf   = any(.data$db == "plasmidfinder"),
      .groups  = "drop"
    ) |>
    dplyr::filter(.data$has_card & .data$has_pf) |>
    dplyr::select("sample", "sequence")

  if (nrow(mobile_contigs) == 0) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: no contigs co-harbouring CARD + PlasmidFinder — skipping")
    return(invisible(NULL))
  }

  # 2. CARD hits on mobile contigs — annotate DRUG (via classify_resistance)
  #    and inherit `phylum` from the R/02 kraken2-derived map.
  mls_rollup <- isTRUE(cfg$network$sankey$mls_rollup %||% TRUE)
  card_mobile <- abri |>
    dplyr::filter(tolower(.data$DATABASE) == "card",
                  !is.na(.data$sequence), !is.na(.data$GENE)) |>
    dplyr::inner_join(mobile_contigs, by = c("sample", "sequence"))
  if (nrow(card_mobile) == 0) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: no CARD hits on mobile contigs — skipping")
    return(invisible(NULL))
  }
  card_mobile$RESISTANCE <- as.character(card_mobile$RESISTANCE)
  if (mls_rollup) {
    card_mobile$RESISTANCE <- stringr::str_replace_all(
      card_mobile$RESISTANCE,
      "lincosamide|macrolide|streptogramin", "MLS"
    )
  }
  card_mobile$DRUG <- vapply(card_mobile$RESISTANCE, classify_resistance,
                              FUN.VALUE = character(1))

  # 3. PlasmidFinder hits on mobile contigs — assign Replicon_Family then
  #    pick the dominant family per (sample, sequence) so each contig
  #    contributes one MGE-type edge instead of fanning out.
  patterns <- cfg$mobilome$family_patterns %||% .default_replicon_patterns()
  pf_mobile <- abri |>
    dplyr::filter(tolower(.data$DATABASE) == "plasmidfinder",
                  !is.na(.data$sequence), !is.na(.data$GENE)) |>
    dplyr::inner_join(mobile_contigs, by = c("sample", "sequence")) |>
    dplyr::mutate(Replicon_Family = classify_replicon_family(.data$GENE,
                                                              patterns)) |>
    dplyr::filter(!is.na(.data$Replicon_Family))
  if (nrow(pf_mobile) == 0) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: no classified PlasmidFinder hits on mobile contigs — skipping")
    return(invisible(NULL))
  }
  pf_by_contig <- pf_mobile |>
    dplyr::group_by(.data$sample, .data$sequence, .data$Replicon_Family) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data$sample, .data$sequence) |>
    dplyr::slice_max(.data$n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("sample", "sequence", "Replicon_Family")

  # 4. Flow rows: (sample, sequence, phylum, DRUG, GENE, Replicon_Family).
  flows <- card_mobile |>
    dplyr::inner_join(pf_by_contig, by = c("sample", "sequence")) |>
    dplyr::filter(!is.na(.data$phylum), !is.na(.data$DRUG),
                  nzchar(.data$phylum), nzchar(.data$DRUG))
  if (nrow(flows) == 0) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: empty flow frame after phylum/drug filter — skipping")
    return(invisible(NULL))
  }

  # 5. Weight ribbons by CARD TPM when normalisation is available; fall
  #    back to Bracken `sampleCount` (always present on cleaned$abri_kraken2).
  norm <- ge_load_norm_table(cfg, "Taxon→ARG→MGE sankey")
  weight_source <- "sampleCount"
  if (!is.null(norm) && all(c("sample", "GENE", "TPM", "DATABASE") %in%
                              colnames(norm))) {
    card_tpm <- norm |>
      dplyr::filter(tolower(.data$DATABASE) == "card") |>
      dplyr::select("sample", "GENE", "TPM")
    flows <- dplyr::left_join(flows, card_tpm, by = c("sample", "GENE"))
    if (any(!is.na(flows$TPM))) {
      flows$weight <- ifelse(is.na(flows$TPM), 0, flows$TPM)
      weight_source <- "TPM"
    } else {
      flows$weight <- as.numeric(flows$sampleCount %||% 1L)
    }
  } else {
    flows$weight <- as.numeric(flows$sampleCount %||% 1L)
  }

  agg <- flows |>
    dplyr::group_by(.data$phylum, .data$DRUG, .data$Replicon_Family) |>
    dplyr::summarise(weight = sum(.data$weight, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(.data$weight > 0)
  if (nrow(agg) == 0) {
    pipeline_log(cfg,
      "Taxon→ARG→MGE sankey: zero-weight aggregate — skipping")
    return(invisible(NULL))
  }

  readr::write_csv(agg,
                    file.path(ds_dir, "sankey_taxon_arg_mge_long.csv"))

  p <- ggplot2::ggplot(
        agg,
        ggplot2::aes(axis1 = .data$phylum,
                     axis2 = .data$DRUG,
                     axis3 = .data$Replicon_Family,
                     y     = .data$weight)
      ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = .data$DRUG), alpha = 0.6, width = 1/8
    ) +
    ggalluvial::geom_stratum(width = 1/8, fill = "white", colour = "grey40") +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = 2.8
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Phylum", "ARG drug class", "MGE type"),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Drug class") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "bottom",
      text            = ggplot2::element_text(size = 11)
    )

  ph <- max(8, 0.35 * nrow(agg) + 5)
  save_panel_ggplot(
    file.path(fig_dir, "sankey_taxon_arg_mge.png"),
    p, width = 12, height = ph, dpi = 300, bg = "white"
  )
  pipeline_log(cfg, sprintf(
    "Taxon→ARG→MGE sankey: %d flows, %d phyla, %d drug classes, %d MGE types (weight=%s)",
    nrow(agg), dplyr::n_distinct(agg$phylum),
    dplyr::n_distinct(agg$DRUG), dplyr::n_distinct(agg$Replicon_Family),
    weight_source
  ))
  invisible(NULL)
}
