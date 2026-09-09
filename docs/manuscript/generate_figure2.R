suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Per-stage wall time (seconds) and peak R-process heap (MB) across four
# cohorts. B1/B2 from the manuscript Table 1 run (commit eee122c, 2026-06-23).
# HW = hospital_wastewater_eg (n=20, run_hww 2026-09-06, 4m 08s).
# WL = wetlands_flyway (n=24, standard run rerun_wet2 2026-09-09, 7m 37s;
#      the complete-ABRicate standard analysis, not the covariate exploration).
stages <- tibble::tribble(
  ~stage,                  ~order, ~b1_sec, ~b1_mb, ~b2_sec, ~b2_mb, ~hw_sec, ~hw_mb, ~wl_sec, ~wl_mb,
  "load_inputs",                1,    0.7,   168,    1.3,   205,    0.7,   167,    1.3,   181,
  "clean_data",                 2,    2.1,   245,    6.4,   522,    2.2,   266,    3.4,   337,
  "normalisation",              3,    0.4,   160,    0.5,   205,    0.4,   157,    0.5,   203,
  "taxonomy",                   4,    2.2,   230,    2.3,   308,    2.5,   280,    3.9,   338,
  "relative_abundance",         5,   27.6,   362,   29.6,   574,   26.2,   391,   32.0,   466,
  "alpha_diversity",            6,   32.6,   472,   31.4,   587,   21.5,   539,   65.0,   646,
  "beta_diversity",             7,    5.9,   472,    5.2,   482,    6.4,   479,    9.1,   601,
  "differential_abundance",     8,  127.0,  1028,  198.0,  1133,   72.0,  1982,  150.0,  3590,
  "resistome",                  9,   29.3,  1028,   25.9,  1102,   16.9,  1194,   23.1,  1240,
  "virulome",                  10,   16.0,  1028,   17.3,  1132,   11.8,  1137,   13.8,  1234,
  "mobilome",                  11,   15.3,  1028,   16.1,  1100,    9.3,  1146,    1.3,   665,
  "network",                   12,   53.3,  1183,   86.0,  1418,   26.1,  1185,   76.0,  1435,
  "panels",                    13,   49.0,  1330,   49.0,  1505,   25.4,  1301,   43.0,  1451,
  "manifest",                  14,    2.8,   922,    3.1,  1040,    2.2,   734,    2.3,   873,
  "report",                    15,   36.9,   611,   27.4,   647,   24.2,   614,   32.1,   646
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
