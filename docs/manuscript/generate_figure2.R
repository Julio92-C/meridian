suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Per-stage wall time (seconds) and peak R-process heap (MB) across four
# cohorts, all run on MERIDIAN v1.2.0 (tag v1.2.0, 91fb5df) with each
# project's shipped config, single-thread, sequentially on one machine
# (2026-10-29/30). Values parsed from each run's pipeline.log.
# B1 = chicken_batch1 (n=18, 8m 08s); B2 = chicken_batch2 (n=32, 11m 46s);
# HW = hospital_wastewater_eg (n=20, 3m 23s); WL = wetlands_flyway (n=24, 6m 24s).
stages <- tibble::tribble(
  ~stage,                  ~order, ~b1_sec, ~b1_mb, ~b2_sec, ~b2_mb, ~hw_sec, ~hw_mb, ~wl_sec, ~wl_mb,
  "load_inputs",                1,    0.9,   168,    1.5,   205,    0.7,   166,    0.7,   181,
  "clean_data",                 2,    2.6,   247,    7.4,   644,    2.1,   273,    3.1,   377,
  "normalisation",              3,    0.5,   163,    0.6,   207,    0.3,   160,    0.4,   207,
  "taxonomy",                   4,    2.5,   226,    2.4,   350,    3.4,   290,    4.1,   375,
  "relative_abundance",         5,   17.6,   341,   17.2,   492,   18.4,   376,   19.4,   512,
  "alpha_diversity",            6,   28.9,   479,   31.1,   617,   14.5,   561,   56.6,   633,
  "beta_diversity",             7,    5.0,   466,    4.8,   522,    5.3,   514,    7.4,   596,
  "differential_abundance",     8,  276.0,  1078,  437.0,   953,   69.0,  2090,  136.0,  3529,
  "resistome",                  9,   12.9,  1029,   13.8,   881,    6.1,  1101,    9.4,  1181,
  "virulome",                  10,    9.3,  1022,   10.0,  1113,    7.7,  1086,    7.6,  1204,
  "mobilome",                  11,    9.0,  1019,    9.4,  1136,    5.3,  1061,    1.2,   673,
  "network",                   12,   50.2,  1181,   82.0,  1447,   21.0,  1117,   68.0,  1432,
  "report",                    13,   28.3,   612,   39.4,   648,   21.1,   616,   24.1,   650,
  "panels",                    14,   41.7,  1731,   45.6,  1566,   25.9,  1331,   43.7,  1438,
  "manifest",                  15,    2.8,   920,    3.4,  1041,    2.1,   736,    2.1,   875
)

cohort_colors <- c(B1 = "#2c7fb8", B2 = "#e6550d",
                   HW = "#31a354", WL = "#756bb1")
cohort_labels <- c(
  B1 = "B1  Chicken batch 1 (n = 18)",
  B2 = "B2  Chicken batch 2 (n = 32)",
  HW = "HW  Hospital wastewater (n = 20)",
  WL = "WL  Wetlands (n = 24)"
)
cohort_levels <- c("B1", "B2", "HW", "WL")

pretty_stage <- function(s) {
  s <- gsub("_", " ", s)
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

fmt_duration <- function(sec) {
  vapply(sec, function(s) {
    if (s < 60) {
      if (s < 10) sprintf("%.1f s", s) else sprintf("%.0f s", s)
    } else {
      m <- floor(s / 60); r <- round(s - 60 * m)
      sprintf("%dm %02ds", m, r)
    }
  }, character(1))
}

stage_order_by_b2 <- stages %>%
  arrange(desc(b2_sec)) %>%
  pull(stage)

# ----- Panel A: wall time per stage -----
runtime_long <- stages %>%
  select(stage, B1 = b1_sec, B2 = b2_sec, HW = hw_sec, WL = wl_sec) %>%
  pivot_longer(all_of(cohort_levels), names_to = "cohort", values_to = "sec") %>%
  mutate(
    stage = factor(pretty_stage(stage), levels = rev(pretty_stage(stage_order_by_b2))),
    cohort = factor(cohort, levels = cohort_levels)
  )

panel_a <- ggplot(runtime_long, aes(x = sec, y = stage, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_x_log10(
    breaks = c(0.5, 1, 5, 10, 30, 60, 120, 300),
    labels = c("0.5 s", "1 s", "5 s", "10 s", "30 s", "1 m", "2 m", "5 m"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_manual(values = cohort_colors, labels = cohort_labels) +
  labs(
    title = "A. Wall time per stage",
    x = "wall time (log scale)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

# ----- Panel B: peak memory per stage (same stage order as A) -----
mem_long <- stages %>%
  select(stage, B1 = b1_mb, B2 = b2_mb, HW = hw_mb, WL = wl_mb) %>%
  pivot_longer(all_of(cohort_levels), names_to = "cohort", values_to = "mb") %>%
  mutate(
    stage = factor(pretty_stage(stage), levels = rev(pretty_stage(stage_order_by_b2))),
    cohort = factor(cohort, levels = cohort_levels)
  )

panel_b <- ggplot(mem_long, aes(x = mb, y = stage, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000),
    labels = function(x) paste0(x, " MB"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_manual(values = cohort_colors, labels = cohort_labels) +
  labs(
    title = "B. Peak R heap per stage",
    x = "peak resident MB",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

# ----- Panel C: cumulative runtime in pipeline order -----
cum_runtime <- stages %>%
  arrange(order) %>%
  mutate(B1 = cumsum(b1_sec), B2 = cumsum(b2_sec),
         HW = cumsum(hw_sec), WL = cumsum(wl_sec)) %>%
  select(order, stage, all_of(cohort_levels)) %>%
  pivot_longer(all_of(cohort_levels), names_to = "cohort", values_to = "cum_sec") %>%
  mutate(cohort = factor(cohort, levels = cohort_levels))

end_labels_c <- cum_runtime %>%
  group_by(cohort) %>%
  filter(order == max(order)) %>%
  ungroup() %>%
  mutate(label = fmt_duration(cum_sec))

panel_c <- ggplot(cum_runtime, aes(x = order, y = cum_sec, color = cohort)) +
  geom_step(linewidth = 0.7) +
  geom_point(size = 1.2) +
  geom_text(
    data = end_labels_c,
    aes(label = label),
    hjust = 1.1, vjust = -0.6, size = 2.7, fontface = "bold", show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = stages$order,
    labels = pretty_stage(stages$stage),
    expand = expansion(add = c(0.3, 0.3))
  ) +
  scale_y_continuous(
    breaks = c(0, 60, 120, 180, 240, 300, 360, 420, 480),
    labels = function(s) ifelse(s == 0, "0", fmt_duration(s)),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_color_manual(values = cohort_colors, labels = cohort_labels) +
  labs(
    title = "C. Cumulative wall time",
    x = NULL,
    y = "cumulative time",
    color = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ----- Panel D: running peak memory in pipeline order -----
run_peak <- stages %>%
  arrange(order) %>%
  mutate(B1 = cummax(b1_mb), B2 = cummax(b2_mb),
         HW = cummax(hw_mb), WL = cummax(wl_mb)) %>%
  select(order, stage, all_of(cohort_levels)) %>%
  pivot_longer(all_of(cohort_levels), names_to = "cohort", values_to = "max_mb") %>%
  mutate(cohort = factor(cohort, levels = cohort_levels))

end_labels_d <- run_peak %>%
  group_by(cohort) %>%
  filter(order == max(order)) %>%
  ungroup() %>%
  mutate(label = sprintf("%d MB", max_mb))

panel_d <- ggplot(run_peak, aes(x = order, y = max_mb, color = cohort)) +
  geom_step(linewidth = 0.7) +
  geom_point(size = 1.2) +
  geom_text(
    data = end_labels_d,
    aes(label = label),
    hjust = 1.1, vjust = -0.6, size = 2.7, fontface = "bold", show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = stages$order,
    labels = pretty_stage(stages$stage),
    expand = expansion(add = c(0.3, 0.3))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, " MB"),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_color_manual(values = cohort_colors, labels = cohort_labels) +
  labs(
    title = "D. Running peak memory",
    x = NULL,
    y = "running max heap",
    color = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ----- Assemble -----
shared_legend <- ggplot(runtime_long, aes(x = sec, y = stage, fill = cohort)) +
  geom_col() +
  scale_fill_manual(values = cohort_colors, labels = cohort_labels) +
  labs(fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8)
  )

legend_grob <- cowplot::get_legend(shared_legend)

fig2 <- (
  (panel_a | panel_b) /
  (panel_c | panel_d)
) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

fig2_with_legend <- patchwork::wrap_elements(legend_grob) / fig2 +
  plot_layout(heights = c(0.08, 1))

out_dir <- "C:/Users/julio/Desktop/Metagenomics_pipeline_automation/docs/manuscript"
png_path  <- file.path(out_dir, "figure2_pipeline_performance.png")
tiff_path <- file.path(out_dir, "figure2_pipeline_performance.tiff")

ggsave(
  png_path, fig2_with_legend,
  width = 190, height = 205, units = "mm", dpi = 300
)
ggsave(
  tiff_path, fig2_with_legend,
  width = 190, height = 205, units = "mm", dpi = 300,
  compression = "lzw"
)

message("Wrote ", png_path)
message("Wrote ", tiff_path)
