#!/usr/bin/env Rscript
# polish_taxa_overview.R — hospital_microbiome Figure 1: microbial community
# overview on the same A–F grid as the genetic-element overview figures
# (resistome / virulome / mobilome composites), built from pipeline outputs:
#   A  species richness per sample (violin), Kruskal–Wallis + BH
#   B  species Shannon diversity per sample (violin), Kruskal–Wallis + BH
#   C  species Bray–Curtis NMDS (relative abundances), PERMANOVA / PERMDISP
#   D  Venn of species detected in air vs surface samples
#   E  top-10 genera | top-10 species, mean relative abundance per group (1x2)
#   F  class-level mean relative abundance per group
# Community rules as in the text: Kraken2 rank S (C for panel F), > 10 reads
# per sample, samples with < 100 species-level reads excluded. Plasmodium
# (likely false positive) is excluded from E and F as in the composition
# figures; A–D keep the community tables used for the reported statistics.
#
# Inputs : <datasets_dir>/{alpha_diversity.csv, noncontaminants_list.csv,
#          relative_abundance/genus_mean_relative_abundance_by_group.csv,
#          relative_abundance/species_top30_long.csv,
#          relative_abundance/relative_abundance_colour_key.csv}
#          (the last three from polish_relative_abundance.R — run it first)
# Outputs: <figures_dir>/taxa_overview/taxa_overview.{png,tiff}
#          <datasets_dir>/taxa_overview/taxa_overview_stats.csv
# Run from the repo root: Rscript projects/hospital_microbiome/polish_taxa_overview.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(VennDiagram)
  library(patchwork)
})
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
`%||%` <- function(a, b) if (is.null(a)) b else a

script_dir <- "projects/hospital_microbiome"
cfg     <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
root    <- cfg$project_root
ds_dir  <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "taxa_overview")
out_dir <- file.path(ds_dir, "taxa_overview")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols      <- unlist(cfg$colors$type)       # Air grey, Surface pink
floor_n   <- cfg$filters$min_count_per_sample %||% 0
min_reads <- cfg$taxa_filter$min_sample_reads %||% 0
nperm     <- cfg$stats$permanova_permutations %||% 9999
group_lab <- cfg$panels$group_label %||% "Group"
excl_gen  <- "Plasmodium"
base_sz   <- 12

meta <- read_csv(file.path(root, cfg$metadata$file), show_col_types = FALSE) |>
  filter(type %in% c("Air", "Surface")) |>
  select(sample, type)

nc <- read_csv(file.path(ds_dir, "noncontaminants_list.csv"),
               col_select = c("sample", "name", "tax_rank", "count"),
               show_col_types = FALSE) |>
  mutate(name = trimws(name), tax_rank = trimws(tax_rank)) |>
  inner_join(meta, by = "sample")

# ---- Species table + analysed samples -------------------------------------------
sp <- filter(nc, tax_rank == "S", count > floor_n)
depth <- sp |> group_by(sample) |> summarise(n = sum(count), .groups = "drop")
keep_samples <- depth$sample[depth$n >= min_reads]
sp <- filter(sp, sample %in% keep_samples)

theme_panel <- function() {
  theme_classic(base_size = base_sz) +
    theme(plot.title = element_text(face = "bold", size = 16),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
}
fp <- function(x) if (x < 0.001) "p < 0.001" else sprintf("p = %.3f", x)

# ---- A, B: alpha diversity violins ------------------------------------------------
alpha <- read_csv(file.path(ds_dir, "alpha_diversity.csv"), show_col_types = FALSE) |>
  filter(sample %in% keep_samples) |>
  select(sample, richness, shannon, simpson, chao1) |>
  inner_join(meta, by = "sample") |>
  mutate(type = factor(type, levels = c("Air", "Surface")))
kw <- sapply(c("richness", "shannon", "simpson", "chao1"), function(v) kruskal.test(alpha[[v]], alpha$type)$p.value)
kw_adj <- p.adjust(kw, method = "BH")   # BH over the four indices, as the pipeline (R/06)

violin <- function(v, ylab, title) {
  lab <- sprintf("Kruskal–Wallis p = %.3f | p_adj (BH) = %.2f", kw[[v]], kw_adj[[v]])
  ggplot(alpha, aes(type, .data[[v]], fill = type)) +
    geom_violin(trim = TRUE, alpha = 0.6, colour = "grey30", linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    geom_jitter(width = 0.08, size = 1.3, alpha = 0.8) +
    annotate("text", x = 1.5, y = Inf, label = lab, vjust = 1.2, size = 3.6) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = group_lab, y = ylab, title = title) +
    theme_panel()
}
pA <- violin("richness", "Species richness", "A")
pB <- violin("shannon", "Shannon index", "B")

# ---- C: species Bray–Curtis NMDS ------------------------------------------------
w <- sp |> group_by(sample, name) |> summarise(v = sum(count), .groups = "drop") |>
  pivot_wider(names_from = name, values_from = v, values_fill = 0)
m <- as.matrix(w[, -1]); rownames(m) <- w$sample
m <- m / rowSums(m)
grp <- factor(meta$type[match(rownames(m), meta$sample)], levels = c("Air", "Surface"))
d <- vegdist(m, method = "bray")
set.seed(12345); nm <- metaMDS(d, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE)
set.seed(12345); pm <- adonis2(d ~ grp, permutations = nperm)
set.seed(12345); pd <- permutest(betadisper(d, grp), permutations = nperm)$tab$`Pr(>F)`[1]
pts <- as.data.frame(scores(nm, display = "sites")); pts$type <- grp
labC <- sprintf("PERMANOVA R² = %.3f, %s\nPERMDISP %s\nNMDS stress = %.3f",
                pm$R2[1], fp(pm$`Pr(>F)`[1]), fp(pd), nm$stress)
pC <- ggplot(pts, aes(NMDS1, NMDS2, colour = type, fill = type)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, alpha = 0.15, show.legend = FALSE) +
  geom_point(size = 2.2, alpha = 0.9) +
  annotate("text", x = -Inf, y = Inf, label = labC, hjust = -0.04, vjust = 1.15, size = 3.6,
           lineheight = 1.05) +
  scale_colour_manual(values = cols, name = group_lab) +
  scale_fill_manual(values = cols, name = group_lab) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(title = "C") +
  theme_panel() + theme(legend.position = "right")

# ---- D: species Venn -------------------------------------------------------------
sets <- split(sp$name, sp$type) |> lapply(unique)
ord <- names(sort(lengths(sets[c("Air", "Surface")]), decreasing = TRUE))
g <- draw.pairwise.venn(
  area1 = length(sets[[ord[1]]]), area2 = length(sets[[ord[2]]]),
  cross.area = length(intersect(sets$Air, sets$Surface)),
  category = c("", ""), fill = unname(cols[ord]),
  rotation.degree = if (ord[1] == "Surface") 180 else 0,
  print.mode = c("raw", "percent"), sigdigs = 3, alpha = 0.6, lty = "blank",
  cex = 1.2, fontfamily = "sans", scaled = FALSE, ind = FALSE)
labs_v <- grid::gList(
  grid::textGrob("Air", x = 0.30, y = 0.92, gp = grid::gpar(fontsize = 14)),
  grid::textGrob("Surface", x = 0.70, y = 0.92, gp = grid::gpar(fontsize = 14)))
pD <- wrap_elements(full = grid::gTree(children = grid::gList(g, labs_v))) +
  labs(title = "D") + theme(plot.title = element_text(face = "bold", size = 16))

# ---- E: top-10 genera | top-10 species per group ---------------------------------
key <- read_csv(file.path(ds_dir, "relative_abundance", "relative_abundance_colour_key.csv"),
                show_col_types = FALSE)
gen <- read_csv(file.path(ds_dir, "relative_abundance", "genus_mean_relative_abundance_by_group.csv"),
                show_col_types = FALSE) |>
  filter(!taxon %in% excl_gen)
spl <- read_csv(file.path(ds_dir, "relative_abundance", "species_top30_long.csv"),
                show_col_types = FALSE) |>
  filter(sample %in% keep_samples, taxon != "Other", !sub(" .*", "", taxon) %in% excl_gen)
n_grp <- tibble(sample = keep_samples) |> inner_join(meta, by = "sample") |> count(type, name = "n")
spm <- spl |> group_by(type, taxon) |> summarise(s = sum(prop), .groups = "drop") |>
  inner_join(n_grp, by = "type") |> mutate(mean_prop = s / n) |> select(type, taxon, mean_prop)

top_bars <- function(df, level, title, italic_lab) {
  top <- df |> group_by(type) |> slice_max(mean_prop, n = 10, with_ties = FALSE) |> pull(taxon) |> unique()
  dd <- df |> mutate(taxon = ifelse(taxon %in% top, taxon, "Other")) |>
    group_by(type, taxon) |> summarise(p = 100 * sum(mean_prop), .groups = "drop")
  dd <- bind_rows(dd, dd |> group_by(type) |> summarise(p = 100 - sum(p)) |>
                    mutate(taxon = "Other")) |>
    group_by(type, taxon) |> summarise(p = sum(p), .groups = "drop")
  lev <- dd |> filter(taxon != "Other") |> group_by(taxon) |> summarise(t = sum(p)) |>
    arrange(desc(t)) |> pull(taxon)
  dd$taxon <- factor(dd$taxon, levels = c(lev, "Other"))
  k <- key |> filter(level == !!level)
  pal <- setNames(k$colour[match(lev, k$taxon)], lev)
  pal[is.na(pal)] <- grDevices::hcl.colors(sum(is.na(pal)), "Set 2")
  pal <- c(pal, Other = "grey85")
  ggplot(dd, aes(type, p, fill = taxon)) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = title) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
    labs(x = group_lab, y = "Mean relative abundance (%)") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    theme_panel() +
    theme(legend.text = element_text(size = 9, face = "italic"),
          legend.key.size = grid::unit(0.38, "cm"))
}
pE1 <- top_bars(gen, "genus", "Genus", TRUE) + labs(title = "E")
pE2 <- top_bars(spm, "species", "Species", TRUE) + labs(y = NULL)
pE <- wrap_elements(full = pE1 | pE2)

# ---- F: class-level composition -------------------------------------------------
cls <- nc |> filter(sample %in% keep_samples, tax_rank == "C", count > floor_n)
plas <- nc |> filter(sample %in% keep_samples, tax_rank == "G", name %in% excl_gen) |>
  group_by(sample) |> summarise(pl = sum(count), .groups = "drop")
cls <- cls |> left_join(plas, by = "sample") |>
  mutate(count = ifelse(name == "Aconoidasida", pmax(count - coalesce(pl, 0), 0), count)) |>
  filter(count > 0) |>
  group_by(sample, type) |> mutate(prop = count / sum(count)) |> ungroup()
clm <- cls |> group_by(type, name) |> summarise(s = sum(prop), .groups = "drop") |>
  inner_join(n_grp, by = "type") |> mutate(p = 100 * s / n)
top_cls <- clm |> group_by(name) |> summarise(t = sum(p)) |> arrange(desc(t)) |> slice_head(n = 15) |> pull(name)
clm <- clm |> mutate(name = ifelse(name %in% top_cls, name, "Others")) |>
  group_by(type, name) |> summarise(p = sum(p), .groups = "drop") |>
  mutate(name = factor(name, levels = c(top_cls, "Others")), x = as.numeric(factor(type)))
igv <- c("#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD",
         "#802268", "#6BD76B", "#D595A7", "#924822", "#837B8D", "#C75127", "#D58F5C", "#7A65A5")
pal_c <- c(setNames(igv[seq_along(top_cls)], top_cls), Others = "grey80")
pF <- ggplot(clm, aes(x, p, fill = name)) +
  geom_area(colour = "white", linewidth = 0.2) +
  scale_x_continuous(breaks = 1:2, labels = c("Air", "Surface"), expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = pal_c, name = "Class") +
  labs(x = group_lab, y = "Relative abundance (%)", title = "F") +
  theme_panel() + theme(legend.key.size = grid::unit(0.38, "cm"), legend.text = element_text(size = 9))

# ---- Compose + save -------------------------------------------------------------
fig <- wrap_plots(pA, pB, pC, pD, pE, pF, design = "ABC\nDEF", widths = c(1, 1.75, 1.1))
ggsave(file.path(fig_dir, "taxa_overview.png"), fig, width = 18, height = 10.5, dpi = 300)
ggsave(file.path(fig_dir, "taxa_overview.tiff"), fig, width = 18, height = 10.5, dpi = 300,
       compression = "lzw")

vc <- c(air_only = length(setdiff(sets$Air, sets$Surface)),
        surface_only = length(setdiff(sets$Surface, sets$Air)),
        shared = length(intersect(sets$Air, sets$Surface)))
med <- alpha |> group_by(type) |> summarise(rich = median(richness), shan = median(shannon))
stats <- data.frame(
  n_air = sum(grp == "Air"), n_surface = sum(grp == "Surface"),
  richness_median_air = med$rich[1], richness_median_surface = med$rich[2],
  richness_kw_p = kw[["richness"]], richness_p_adj = kw_adj[["richness"]],
  shannon_median_air = med$shan[1], shannon_median_surface = med$shan[2],
  shannon_kw_p = kw[["shannon"]], shannon_p_adj = kw_adj[["shannon"]],
  permanova_r2 = pm$R2[1], permanova_p = pm$`Pr(>F)`[1], permdisp_p = pd, nmds_stress = nm$stress,
  t(vc), shared_pct = 100 * vc[["shared"]] / sum(vc))
write_csv(stats, file.path(out_dir, "taxa_overview_stats.csv"))
write_csv(select(clm, type, class = name, mean_pct = p), file.path(out_dir, "class_mean_by_group.csv"))
cat(sprintf("A/B: n = %d air / %d surface; richness median %.0f vs %.0f (p_adj %.2f); Shannon median %.2f vs %.2f (p_adj %.2f)\n",
            stats$n_air, stats$n_surface, med$rich[1], med$rich[2], kw_adj[["richness"]],
            med$shan[1], med$shan[2], kw_adj[["shannon"]]))
cat(sprintf("C: PERMANOVA R2 = %.3f, p = %.4f, PERMDISP p = %.4f, stress = %.3f\n",
            pm$R2[1], pm$`Pr(>F)`[1], pd, nm$stress))
cat(sprintf("D: air-only %d, surface-only %d, shared %d (%.1f%%)\n", vc[1], vc[2], vc[3], stats$shared_pct))
cat("F classes:", paste(top_cls, collapse = ", "), "\n")
cat("Wrote", file.path(fig_dir, "taxa_overview.png"), "(+ .tiff)\n")
