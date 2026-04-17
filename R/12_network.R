# 12_network.R — tripartite network (Treatment × taxa × GE) for Gephi + an
# in-R igraph rendering. Computes degree, betweenness, modularity as per the
# paper's topology metrics.

run_network <- function(cleaned, cfg) {
  pipeline_log(cfg, "Network analysis (tripartite)")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "network")
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "network")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ds_dir,  recursive = TRUE, showWarnings = FALSE)

  meta <- readr::read_csv(
    file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE
  )
  sid   <- cfg$metadata$sample_id_col
  group <- cfg$metadata$group_cols[[1]]

  df <- dplyr::inner_join(cleaned$abri_kraken2, meta,
                          by = c("sample" = sid))

  # Edge list: treatment → taxon, taxon → GE ------------------------------
  e_tt <- df |>
    dplyr::distinct(.data[[group]], name) |>
    dplyr::rename(source = !!group, target = name) |>
    dplyr::mutate(type = "treatment-taxon")
  e_tg <- df |>
    dplyr::distinct(name, GENE) |>
    dplyr::rename(source = name, target = GENE) |>
    dplyr::mutate(type = "taxon-gene")

  edges <- dplyr::bind_rows(e_tt, e_tg)
  nodes <- dplyr::tibble(
    id = unique(c(edges$source, edges$target))
  ) |>
    dplyr::mutate(
      kind = dplyr::case_when(
        id %in% unique(df[[group]]) ~ "treatment",
        id %in% unique(df$name)      ~ "taxon",
        TRUE                         ~ "gene"
      )
    )

  readr::write_csv(edges, file.path(ds_dir, "gephi_edges.csv"))
  readr::write_csv(nodes, file.path(ds_dir, "gephi_nodes.csv"))

  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
    topo <- data.frame(
      node       = igraph::V(g)$name,
      kind       = igraph::V(g)$kind,
      degree     = igraph::degree(g),
      betweenness = igraph::betweenness(g, normalized = TRUE)
    )
    cluster <- igraph::cluster_louvain(g)
    topo$module <- cluster$membership
    readr::write_csv(topo, file.path(ds_dir, "topology.csv"))
    pipeline_log(cfg, sprintf("Network: %d nodes, %d edges, modularity = %.3f",
                              igraph::vcount(g), igraph::ecount(g),
                              igraph::modularity(cluster)))
  }

  invisible(list(edges = edges, nodes = nodes))
}
