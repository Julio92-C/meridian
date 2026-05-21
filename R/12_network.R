# 12_network.R — Tripartite (sample × taxon × gene) network + chord diagrams.
# Mirrors the GT scripts tripartiteNetwork.R, chordDiagram_v1.2.1.R, and
# chordDiagramNorm_v1.2.1.R, but drives everything from cfg — no hand-tuned
# row-index cell edits, no hardcoded treatments / sample IDs / palettes.
# Uses clean_taxa_names() from R/utils_taxa.R and classify_resistance() from
# R/09_resistome.R.
#
# Outputs (under <project>/<figures_dir>/network/):
#   network.png             ggraph FR-layout tripartite network
#   chord_overall.png       Chord of sample → taxon → gene → drug class (all)
#   chord_<group>.png       One chord per group level (toggleable)
#
# Datasets (under <project>/<datasets_dir>/network/):
#   gephi_edges.csv         Source / Target / Weight rows for Gephi import
#   gephi_nodes.csv         Id / type / <group> / geneCategory /
#                           Function_group / distance_cluster
#   topology.csv            node / kind / degree / betweenness / module
#   sample_clusters.csv     sample → Bray-Curtis hierarchical cluster id

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
  src <- clean_taxa_names(src, cfg)

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

  nodes <- dplyr::tibble(name = unique(c(edges$from, edges$to))) |>
    dplyr::mutate(type = dplyr::case_when(
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
    sum(nodes$type == "sample"),
    sum(nodes$type == "taxa"),
    sum(nodes$type == "gene")
  ))

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
    kind        = igraph::V(g)$type,
    degree      = igraph::degree(g),
    betweenness = igraph::betweenness(g, normalized = TRUE),
    module      = cluster$membership
  )
  readr::write_csv(topo, file.path(ds_dir, "topology.csv"))
  pipeline_log(cfg, sprintf(
    "Network topology: modularity = %.3f",
    igraph::modularity(cluster)
  ))

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
  edge_width_range <- ncfg$edge_width_range  %||% c(0.3, 2.5)
  edge_alpha       <- ncfg$edge_alpha        %||% 0.2
  node_alpha       <- ncfg$node_alpha        %||% 0.75

  col_key <- ifelse(
    nodes$type == "gene",   as.character(nodes$geneCategory),
    ifelse(nodes$type == "sample",
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
      ggplot2::aes(shape = .data$type, color = .data$col_key),
      size = pt_size, alpha = node_alpha
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      repel = TRUE, size = label_size, max.overlaps = Inf,
      segment.color = "grey60", min.segment.length = 0.5
    ) +
    ggraph::scale_edge_width(range = edge_width_range) +
    ggplot2::scale_color_manual(values = pal, na.value = "grey60") +
    ggplot2::scale_shape_manual(
      values = c(sample = 16, taxa = 17, gene = 15)
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(
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

  grDevices::png(png_path, width = 2400, height = 2400, res = 200, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  circlize::circos.clear()
  circlize::circos.par(
    track.height = 0.1, start.degree = start_degree, gap.degree = 2,
    canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1),
    circle.margin = c(1, 1), unit.circle.segments = 500
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
        adj = c(0, 0.8), cex = 0.8
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
