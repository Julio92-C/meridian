# 04_shared_fraction.R — the species shared between groups, how widespread
# they are, which genes their reads carry, and how they link the levels of
# metadata variables (ward, room, clinical status ...). Any number of groups.
#
# Uses the same species table as R/04_taxonomy.R (community_filter ->
# clean_taxa_names -> taxonomy.min_count_per_sample), so the unique / shared
# counts agree with the taxonomy Venn of the same run. Names are kept
# unabbreviated internally (the tree needs the genus) and shortened for display.
#
# Outputs: <datasets_dir>/shared_fraction/*.csv, <figures_dir>/shared_fraction/*.png
# Definitions and config: utils_shared_fraction.R and the design spec
# docs/superpowers/specs/2026-10-30-taxa-shared-fraction-design.md.

# Species x sample table (floored, species rank, NA groups dropped) plus the
# analysed samples per group. Shared by the shared-fraction and tree stages.
sf_prepare_inputs <- function(cleaned, cfg, stage = "Shared fraction") {
  cfg_full <- cfg
  cfg_full$taxa_display$abbreviate_genus <- FALSE
  cfg_full$taxa_display$max_name_len <- 1000
  df <- clean_taxa_names(community_filter(cleaned$noncontaminants, cfg, stage), cfg_full)
  if (nrow(df) == 0) return(NULL)
  group <- cfg$metadata$group_cols[[1]]
  sid   <- cfg$metadata$sample_id_col
  meta  <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file), show_col_types = FALSE)
  meta[[sid]] <- as.character(meta[[sid]])
  df$sample <- as.character(df$sample)
  meta <- meta[meta[[sid]] %in% unique(df$sample) & !is.na(meta[[group]]), , drop = FALSE]
  df <- df[df$sample %in% meta[[sid]], , drop = FALSE]
  group_of <- stats::setNames(as.character(meta[[group]]), meta[[sid]])
  lv <- levels(factor(group_of))
  group_n <- vapply(lv, function(g) sum(group_of == g), integer(1))
  min_count <- cfg$taxonomy$min_count_per_sample %||% 5
  long <- df[df$count >= min_count, , drop = FALSE]
  long <- long |>
    dplyr::group_by(.data$sample, .data$name) |>
    dplyr::summarise(count = sum(.data$count),
                     dplyr::across(dplyr::any_of(c("phylum", "class", "order", "family", "genus")), dplyr::first),
                     .groups = "drop")
  long$group <- group_of[long$sample]
  list(long = long, group_of = group_of, group_n = group_n, groups = lv,
       meta = meta, sid = sid, group = group, min_count = min_count)
}

# Display colour for each category: unique -> group colour; shared -> blend
# of the colours of the groups in the pattern.
sf_category_colours <- function(categories, pal) {
  blend <- function(cols) {
    m <- grDevices::col2rgb(cols)
    grDevices::rgb(t(rowMeans(m)), maxColorValue = 255)
  }
  vapply(categories, function(ct) {
    grp <- sub("^(unique|shared):", "", ct)
    gs <- if (grp == "all") names(pal) else strsplit(grp, "+", fixed = TRUE)[[1]]
    if (startsWith(ct, "unique:")) pal[[gs]] else blend(pal[gs])
  }, character(1))
}

sf_category_label <- function(ct, groups) {
  ifelse(ct == "shared:all", if (length(groups) == 2) "Shared" else "Shared by all groups",
  ifelse(startsWith(ct, "unique:"), paste(sub("^unique:", "", ct), "only"),
         paste("Shared:", gsub("+", " + ", sub("^shared:", "", ct), fixed = TRUE))))
}

sf_theme <- function() {
  ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                   legend.key.size = ggplot2::unit(0.35, "cm"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"))
}

run_shared_fraction <- function(cleaned, cfg) {
  pipeline_log(cfg, "Shared fraction (overlap tiers, elements, linkage)")
  st <- sf_settings(cfg)
  inp <- sf_prepare_inputs(cleaned, cfg)
  if (is.null(inp) || length(inp$groups) < 2) {
    pipeline_log(cfg, "Shared fraction: fewer than 2 groups or no species — skipped")
    return(invisible(NULL))
  }
  ds_dir  <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "shared_fraction")
  fig_dir <- file.path(cfg$project_root, cfg$outputs$figures_dir, "shared_fraction")
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  groups <- inp$groups; pal <- group_palette(groups, cfg)

  tab <- sf_species_table(inp$long, inp$group_n, st$min_samples, st$min_prev)
  pw  <- sf_pairwise(tab, groups, st$min_samples, st$min_prev)
  ws  <- tab$name[!is.na(tab$tier) & tab$tier == "widespread"]
  el  <- sf_elements(cleaned$abri_kraken2, ws, names(inp$group_of), st$elements_min_reads)
  readr::write_csv(tab, file.path(ds_dir, "shared_fraction_species.csv"))
  if (!is.null(pw)) readr::write_csv(pw, file.path(ds_dir, "shared_fraction_pairwise.csv"))
  if (!is.null(el)) readr::write_csv(el, file.path(ds_dir, "shared_fraction_widespread_elements.csv"))
  cnt <- table(tab$category)
  pipeline_log(cfg, sprintf("Shared fraction: %d species; %s; tiers among shared: %s; widespread (n>=%s, prev>=%s): %d",
    nrow(tab), paste(sprintf("%s=%d", names(cnt), as.integer(cnt)), collapse = ", "),
    paste(sprintf("%s=%d", SF_TIERS, as.integer(table(factor(tab$tier, SF_TIERS)))), collapse = ", "),
    st$min_samples, st$min_prev, length(ws)))

  sf_plot_prevalence_abundance(tab, groups, pal, ws, file.path(fig_dir, "shared_fraction_prevalence_abundance.png"))
  if (!is.null(pw)) sf_plot_pairwise(pw, groups, file.path(fig_dir, "shared_fraction_pairwise.png"))
  sf_plot_upset(tab, groups, pal, file.path(fig_dir, "shared_fraction_upset.png"))
  sf_plot_tiers(tab, file.path(fig_dir, "shared_fraction_tiers.png"))
  if (length(ws) > 0) sf_plot_widespread(tab[tab$name %in% ws, ], groups, pal, el,
                                         file.path(fig_dir, "shared_fraction_widespread.png"), st$elements_min_reads)
  sf_plot_heatmap(inp, tab, st, file.path(fig_dir, "shared_fraction_heatmap.png"))

  # ---- linkage by metadata variable ----------------------------------------
  edges_all <- list()
  for (v in st$linkage_vars) {
    if (!v %in% colnames(inp$meta)) {
      pipeline_log(cfg, sprintf("Shared fraction: linkage variable '%s' not in metadata — skipped", v)); next
    }
    if (length(ws) == 0) { pipeline_log(cfg, "Shared fraction: no widespread species — linkage skipped"); break }
    lv <- stats::setNames(sf_bin_variable(inp$meta[[v]], st$linkage_bins), inp$meta[[inp$sid]])
    link <- sf_linkage(inp$long, lv, ws, inp$group_of, st$linkage_min_shared, st$linkage_scope, st$linkage_keep_na)
    if (st$linkage_test) {
      link <- sf_linkage_test(inp$long, lv, ws, inp$group_of, link, st$linkage_n_perm, cluster_seed(cfg),
                              scope = st$linkage_scope, keep_na = st$linkage_keep_na)
      pipeline_log(cfg, sprintf("Shared fraction: linkage '%s' permutation p (total shared) = %.3f", v, link$p_global %||% NA))
    }
    if (!is.null(link$edges) && nrow(link$edges) > 0) edges_all[[v]] <- dplyr::mutate(link$edges, variable = v, .before = 1)
    sf_plot_linkage(link, v, groups, pal, ws, file.path(fig_dir, sprintf("shared_fraction_linkage_%s.png", make.names(v))))
  }
  if (length(edges_all) > 0)
    readr::write_csv(dplyr::bind_rows(edges_all), file.path(ds_dir, "shared_fraction_linkage_edges.csv"))
  invisible(list(table = tab, pairwise = pw, elements = el))
}

# ---- plots -------------------------------------------------------------------

sf_plot_prevalence_abundance <- function(tab, groups, pal, ws, file) {
  cats <- unique(tab$category)
  cc <- sf_category_colours(cats, pal); names(cc) <- sf_category_label(cats, groups)
  d <- tab |> dplyr::mutate(cat = sf_category_label(.data$category, groups),
                            y = log10(pmax(.data$abund_all, 1e-4)), x = .data$prev_all * 100)
  lab <- d |> dplyr::filter(.data$name %in% ws | rank(-.data$abund_all) <= 4) |>
    dplyr::slice_max(.data$abund_all, n = 10)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$cat), size = 1, alpha = 0.6) +
    ggplot2::geom_point(data = d[d$name %in% ws, ], shape = 21, size = 2.6, stroke = 0.8,
                        colour = SF_TIER_COLOURS[["widespread"]], fill = NA) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = sf_short(.data$name)), size = 2.6,
                             fontface = "italic", max.overlaps = 30, segment.size = 0.2) +
    ggplot2::scale_colour_manual(values = cc, name = NULL) +
    ggplot2::labs(x = "Prevalence, all samples (%)", y = expression(log[10]~"mean relative abundance (%)"),
                  title = "Prevalence and abundance of all species",
                  subtitle = sprintf("Red circles: widespread shared species (n = %d)", length(ws))) +
    sf_theme() + ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file, p, width = 6.5, height = 5, dpi = 300)
}

sf_plot_pairwise <- function(pw, groups, file) {
  pw$tier <- factor(pw$tier, SF_TIERS)
  np <- length(unique(pw$pair))
  p <- ggplot2::ggplot(pw, ggplot2::aes(.data$prev_i * 100, .data$prev_j * 100)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$tier, size = .data$abund_i + .data$abund_j),
                        shape = 21, colour = "grey40", stroke = 0.2, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = SF_TIER_COLOURS, name = "Tier", drop = FALSE) +
    ggplot2::scale_size_area(max_size = 6, name = "Mean abundance (%)") +
    ggplot2::facet_wrap(~ pair, ncol = min(3, np)) +
    ggplot2::labs(x = "Prevalence in first group (%)", y = "Prevalence in second group (%)",
                  title = "Shared species: prevalence in each pair of groups") +
    sf_theme()
  ggplot2::ggsave(file, p, width = 3.6 * min(3, np) + 1.6, height = 3.4 * ceiling(np / 3) + 0.8, dpi = 300)
}

sf_plot_upset <- function(tab, groups, pal, file, top = 20) {
  # Bars of species found in one group only take that group's colour; shared
  # patterns are stacked by prevalence tier.
  only <- paste(groups, "only")
  d <- tab |> dplyr::mutate(tier = ifelse(is.na(.data$tier), paste(.data$pattern, "only"), .data$tier))
  pat <- d |> dplyr::count(.data$pattern, name = "n") |> dplyr::slice_max(.data$n, n = top, with_ties = FALSE)
  d <- d[d$pattern %in% pat$pattern, ]
  d$pattern <- factor(d$pattern, levels = pat$pattern[order(-pat$n)])
  d$tier <- factor(d$tier, c(only, SF_TIERS))
  fills <- c(stats::setNames(unname(pal[groups]), only), SF_TIER_COLOURS)
  bars <- ggplot2::ggplot(d, ggplot2::aes(.data$pattern, fill = .data$tier)) +
    ggplot2::geom_bar(width = 0.7) +
    ggplot2::scale_fill_manual(values = fills, name = NULL, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Species", title = "Species by group membership") +
    sf_theme() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                                legend.position = "top")
  mm <- expand.grid(pattern = levels(d$pattern), group = groups, stringsAsFactors = FALSE)
  mm$on <- mapply(function(p, g) g %in% strsplit(p, "+", fixed = TRUE)[[1]], mm$pattern, mm$group)
  mm$pattern <- factor(mm$pattern, levels(d$pattern)); mm$group <- factor(mm$group, rev(groups))
  dots <- ggplot2::ggplot(mm, ggplot2::aes(.data$pattern, .data$group)) +
    ggplot2::geom_line(data = mm[mm$on, ], ggplot2::aes(group = .data$pattern), colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(colour = ifelse(.data$on, as.character(.data$group), "off")), size = 3) +
    ggplot2::scale_colour_manual(values = c(pal, off = "grey88"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) + sf_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
  p <- patchwork::wrap_plots(bars, dots, ncol = 1, heights = c(3, 0.35 + 0.25 * length(groups)))
  ggplot2::ggsave(file, p, width = max(4.5, 0.35 * nlevels(d$pattern) + 2), height = 4 + 0.25 * length(groups), dpi = 300)
}

sf_plot_tiers <- function(tab, file) {
  d <- tab[!is.na(tab$tier), , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  phy <- if ("phylum" %in% colnames(d)) d$phylum else rep("Unknown", nrow(d))
  phy[is.na(phy)] <- "Unknown"
  top <- names(sort(table(phy), decreasing = TRUE))[1:min(6, length(unique(phy)))]
  d$phy <- factor(ifelse(phy %in% top, phy, "Other"), c(top, "Other"))
  d$tier <- factor(d$tier, rev(SF_TIERS))
  tol <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677", "#BBBBBB")
  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$tier, fill = .data$phy)) + ggplot2::geom_bar(width = 0.7) +
    ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count), group = .data$tier, fill = NULL),
                       hjust = -0.2, size = 3) +
    ggplot2::scale_fill_manual(values = stats::setNames(tol[seq_along(levels(d$phy))], levels(d$phy)), name = "Phylum") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = "Shared species", y = NULL, title = "Shared species by tier and phylum") + sf_theme()
  ggplot2::ggsave(file, p, width = 6.5, height = 3, dpi = 300)
}

sf_plot_widespread <- function(wt, groups, pal, el, file, min_reads = 2) {
  d <- tidyr::pivot_longer(wt[c("name", paste0("prev_", groups))], -"name", names_to = "group", values_to = "prev")
  d$group <- factor(sub("^prev_", "", d$group), groups)
  ord <- wt$name[order(wt$min_prev_present)]
  d$name <- factor(sf_short(d$name), sf_short(ord))
  lab <- tibble::tibble(name = factor(sf_short(ord), sf_short(ord)), txt = "0/0/0")
  if (!is.null(el)) {
    k <- el[!el$low_support, ] |> dplyr::count(.data$name, .data$class)
    lab$txt <- vapply(ord, function(n) paste(vapply(c("ARG", "VF", "MGE"), function(c)
      as.character(sum(k$n[k$name == n & k$class == c])), character(1)), collapse = "/"), character(1))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$prev * 100, .data$name, fill = .data$group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = 104, y = .data$name, label = .data$txt), inherit.aes = FALSE,
                       hjust = 0, size = 2.6) +
    ggplot2::annotate("text", x = 104, y = length(ord) + 0.8, label = "ARG/VF/MGE", hjust = 0, size = 2.6, fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 125), breaks = seq(0, 100, 25)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Prevalence (% of samples)", y = NULL, title = "Widespread shared species",
                  subtitle = sprintf("ARG/VF/MGE: distinct genes on the species' reads, each supported by at least %d reads", min_reads)) +
    sf_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom")
  ggplot2::ggsave(file, p, width = 6, height = 1.5 + 0.28 * length(ord) * max(1, length(groups) / 2), dpi = 300)
}

sf_plot_heatmap <- function(inp, tab, st, file) {
  sh <- tab[!is.na(tab$tier), , drop = FALSE]
  if (nrow(sh) == 0) return(invisible(NULL))
  top <- sh |> dplyr::arrange(dplyr::desc(.data$min_prev_present), dplyr::desc(.data$abund_all)) |>
    dplyr::slice_head(n = st$heatmap_top_n)
  pct <- inp$long |> dplyr::group_by(.data$sample) |>
    dplyr::mutate(pct = .data$count / sum(.data$count) * 100) |> dplyr::ungroup() |>
    dplyr::filter(.data$name %in% top$name)
  smp <- names(inp$group_of)
  v <- st$linkage_vars[st$linkage_vars %in% colnames(inp$meta)][1]
  sub_lv <- if (!is.na(v %||% NA)) stats::setNames(sf_bin_variable(inp$meta[[v]], st$linkage_bins), inp$meta[[inp$sid]]) else
    stats::setNames(rep("", length(smp)), smp)
  grid <- expand.grid(sample = smp, name = top$name, stringsAsFactors = FALSE) |>
    dplyr::left_join(pct[c("sample", "name", "pct")], by = c("sample", "name"))
  grid$group <- inp$group_of[grid$sample]; grid$sub <- sub_lv[grid$sample]
  grid$name <- factor(sf_short(grid$name), rev(sf_short(top$name)))
  ord <- smp[order(inp$group_of[smp], sub_lv[smp], smp)]
  grid$sample <- factor(grid$sample, ord)
  ycol <- ifelse(top$tier[match(levels(grid$name), sf_short(top$name))] == "widespread", "#CB181D", "black")
  p <- ggplot2::ggplot(grid, ggplot2::aes(.data$sample, .data$name, fill = log10(.data$pct))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"),
                                  na.value = "grey95", name = expression(log[10]~"%"), limits = c(-3, 2), oob = scales::squish) +
    ggplot2::facet_grid(~ group + sub, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Top %d shared species by prevalence in every group", nrow(top)),
                  subtitle = "Red labels: widespread tier") +
    sf_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
                                axis.text.y = ggplot2::element_text(face = "italic", colour = ycol),
                                panel.spacing.x = ggplot2::unit(1, "pt"),
                                strip.text = ggplot2::element_text(size = 6.5, angle = 90, hjust = 0))
  ggplot2::ggsave(file, p, width = min(16, 3 + 0.12 * length(smp)), height = 1.8 + 0.17 * nrow(top), dpi = 300)
}

sf_plot_linkage <- function(link, v, groups, pal, ws, file) {
  nd <- link$nodes
  if (nrow(nd) == 0) return(invisible(NULL))
  nd$gx <- match(nd$group, groups)
  nd <- nd |> dplyr::group_by(.data$group) |>
    dplyr::mutate(y = 1 - (dplyr::row_number() - 0.5) / dplyr::n()) |> dplyr::ungroup()
  e <- link$edges
  p <- ggplot2::ggplot()
  if (!is.null(e) && nrow(e) > 0) {
    e <- e |> dplyr::left_join(nd[c("node", "gx", "y")], by = c("from" = "node")) |>
      dplyr::left_join(nd[c("node", "gx", "y")], by = c("to" = "node"), suffix = c("", "_to"))
    p <- p + ggplot2::geom_segment(data = e, ggplot2::aes(.data$gx, .data$y, xend = .data$gx_to, yend = .data$y_to,
                                                           linewidth = .data$shared, alpha = .data$shared),
                                   colour = SF_TIER_COLOURS[["widespread"]], lineend = "round") +
      ggplot2::scale_linewidth(range = c(0.2, 3), name = "Widespread species shared",
                               breaks = scales::breaks_pretty(n = 4), labels = scales::label_number(accuracy = 1)) +
      ggplot2::scale_alpha(range = c(0.15, 0.7), guide = "none")
  }
  # Labels sit beside the node, on the outer side (left of the first group,
  # right of the others), offset so they clear the largest node.
  nd$lx <- ifelse(nd$gx == 1, nd$gx - 0.1, nd$gx + 0.1)
  p <- p + ggplot2::geom_point(data = nd, ggplot2::aes(.data$gx, .data$y, size = .data$n_species, fill = .data$group),
                               shape = 21, colour = "black") +
    ggplot2::geom_text(data = nd, ggplot2::aes(.data$lx, .data$y, label = sprintf("%s (%d; n = %d)", .data$level, .data$n_species, .data$n_samples),
                                               hjust = ifelse(.data$gx == 1, 1, 0)), size = 2.6) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_size_area(max_size = 7, name = "Species detected",
                             breaks = scales::breaks_pretty(n = 4), labels = scales::label_number(accuracy = 1)) +
    ggplot2::scale_x_continuous(breaks = seq_along(groups), labels = groups,
                                expand = ggplot2::expansion(add = 0.9)) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Widespread shared species linking %s levels", v),
                  subtitle = sprintf("%d widespread species; schematic layout, not to scale", length(ws))) +
    sf_theme() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file, p, width = 3 + 2.6 * length(groups), height = max(4, 0.35 * max(table(nd$group)) + 2), dpi = 300)
}
