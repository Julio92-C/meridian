# 14_panels.R — assemble publication-ready multi-panel figures + a
# supplementary tables XLSX workbook from the per-stage PNGs the pipeline
# emits plus the kind-tagged manifest. Driven by templates/frontiers_v2_slots.yaml
# (Frontiers v2 layout) by default; override via cfg$panels$slots_yaml.
#
# Outputs (under <project>/<figures_dir>/panels/):
#   main/fig01_<slug>.png    ... fig10_<slug>.png    — 10 main figures
#   supplementary/figS01_<slug>.png ... figS17_<slug>.png — 17 supplementary
#
# Outputs (under <project>/<datasets_dir>/panels/):
#   supplementary_tables.xlsx — multi-tab workbook with one tab per
#   `table_entry` kind in the manifest, plus an Index sheet.

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

run_panels <- function(cfg) {
  pipeline_log(cfg, "Publication panels + supplementary tables")
  pcfg <- cfg$panels %||% list()
  if (isFALSE(pcfg$enabled %||% TRUE)) {
    pipeline_log(cfg, "Panels stage disabled — skipping")
    return(invisible(NULL))
  }

  # Default slot layout ships in templates/; override via cfg$panels$slots_yaml.
  default_yaml <- file.path(.panels_repo_root(), "templates",
                             "frontiers_v2_slots.yaml")
  yaml_path <- pcfg$slots_yaml %||% default_yaml
  if (!file.exists(yaml_path)) {
    pipeline_log(cfg, sprintf(
      "Panels: slot YAML not found at %s — skipping", yaml_path
    ))
    return(invisible(NULL))
  }

  # Panels runs BEFORE manifest in the current pipeline ordering (so the
  # manifest can catalogue what was actually rendered). Instead of round-
  # tripping through manifest.json, ask the per-stage builders in R/13 for
  # the same in-memory dict via build_stages_index(). Empty stage_times is
  # fine — panels only needs figures[] / tables[], not status/duration_s.
  slots         <- yaml::read_yaml(yaml_path)
  stages_index  <- build_stages_index(cfg)
  manifest_like <- list(stages = stages_index)
  figures_index <- .panels_index_figures(manifest_like)
  tables_index  <- .panels_index_tables(manifest_like)

  fig_root <- file.path(cfg$project_root, cfg$outputs$figures_dir, "panels")
  dir.create(file.path(fig_root, "main"),          recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(fig_root, "supplementary"), recursive = TRUE,
             showWarnings = FALSE)

  # Wipe any stale panel PNGs from previous runs / previous YAML versions
  # so the directory only ever reflects the current slot layout. Without
  # this, renaming a slot or moving one between main/supplementary leaves
  # an orphaned file behind.
  for (sub in c("main", "supplementary")) {
    stale <- list.files(file.path(fig_root, sub),
                         pattern = "\\.(png|pdf|tiff?)$", full.names = TRUE)
    if (length(stale) > 0) file.remove(stale)
  }

  built_main <- character()
  for (slug in names(slots$main %||% list())) {
    out_path <- file.path(fig_root, "main", paste0(slug, ".png"))
    built <- .panels_build_slot(slug, slots$main[[slug]], figures_index,
                                cfg, out_path)
    if (!is.null(built)) built_main <- c(built_main, built)
  }
  built_supp <- character()
  for (slug in names(slots$supplementary %||% list())) {
    out_path <- file.path(fig_root, "supplementary", paste0(slug, ".png"))
    built <- .panels_build_slot(slug, slots$supplementary[[slug]],
                                 figures_index, cfg, out_path)
    if (!is.null(built)) built_supp <- c(built_supp, built)
  }
  pipeline_log(cfg, sprintf(
    "Panels: %d/%d main + %d/%d supplementary slots assembled",
    length(built_main), length(slots$main %||% list()),
    length(built_supp), length(slots$supplementary %||% list())
  ))

  # ---- Supplementary tables XLSX -----------------------------------------
  ds_root <- file.path(cfg$project_root, cfg$outputs$datasets_dir, "panels")
  dir.create(ds_root, recursive = TRUE, showWarnings = FALSE)
  xlsx_path <- file.path(ds_root, "supplementary_tables.xlsx")
  .panels_build_supplementary_xlsx(tables_index, cfg, xlsx_path)

  invisible(list(
    main_panels          = built_main,
    supplementary_panels = built_supp,
    supplementary_xlsx   = xlsx_path
  ))
}

# ---------------------------------------------------------------------------
# Repo-root discovery (same trick as run_pipeline.R / smoketest scripts).
# ---------------------------------------------------------------------------

.panels_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_dir <- normalizePath(dirname(sub("^--file=", "", file_arg[1])))
    # Entry can be the repo-root run_pipeline.R OR a one-level-deep script
    # in R/ or scripts/. Probe for templates/ next to the script first; if
    # not there, walk up one level.
    if (dir.exists(file.path(script_dir, "templates"))) return(script_dir)
    return(normalizePath(dirname(script_dir)))
  }
  normalizePath(getwd())
}

# ---------------------------------------------------------------------------
# Composite theme override
# ---------------------------------------------------------------------------

# Applied to every ggplot we re-render natively from .rds. Per-stage plots
# are tuned for standalone viewing (12pt base text, full-size legend); in
# a composite cell those defaults overlap with neighbours. This override
# shrinks text + legend, drops the per-plot title (the composite uses
# A/B/C tags instead), and adds a small plot margin.
#
# Knobs (all optional, sensible defaults):
#   cfg$panels$composite_base_size    (default 9)   base text pt
#   cfg$panels$composite_legend_size  (default 7)   legend text pt
#   cfg$panels$composite_margin_pt    (default 6)   plot margin pt on all sides
#   cfg$panels$composite_unified      (default FALSE) when TRUE, force axis.text =
#                                                    axis.title = base and legend.title =
#                                                    legend.text (no ±1pt offsets), so every
#                                                    glyph in a ggplot panel renders at one of
#                                                    just two sizes (base / legend).
#   cfg$panels$strip_plot_title       (default TRUE) drop per-plot title
.panels_composite_theme <- function(cfg) {
  base <- cfg$panels$composite_base_size   %||% 9
  legd <- cfg$panels$composite_legend_size %||% 7
  marg <- cfg$panels$composite_margin_pt   %||% 6
  unified <- isTRUE(cfg$panels$composite_unified %||% FALSE)
  strip_title <- isTRUE(cfg$panels$strip_plot_title %||% TRUE)
  axis_text_size  <- if (unified) base else max(6, base - 1)
  legend_title_sz <- if (unified) legd else legd + 1
  th <- ggplot2::theme(
    text           = ggplot2::element_text(size = base),
    axis.text      = ggplot2::element_text(size = axis_text_size),
    axis.title     = ggplot2::element_text(size = base),
    strip.text     = ggplot2::element_text(size = base),
    plot.tag       = ggplot2::element_text(size = base),
    legend.text    = ggplot2::element_text(size = legd),
    legend.title   = ggplot2::element_text(size = legend_title_sz),
    legend.key.size  = grid::unit(0.4, "cm"),
    legend.spacing   = grid::unit(0.15, "cm"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin    = ggplot2::margin(marg, marg, marg, marg, "pt")
  )
  if (strip_title) {
    th <- th + ggplot2::theme(
      plot.title    = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  }
  th
}

# ---------------------------------------------------------------------------
# Per-slot / per-panel composite overrides (from the slots YAML)
# ---------------------------------------------------------------------------

# Recognised keys in a slot's `composite:` block (or a panel's `composite:`
# block). Panel-level beats slot-level when both are set.
#
#   legend_position : "right" | "left" | "top" | "bottom" | "none"
#   legend_ncol     : integer  (forces guide_legend ncol)
#   legend_nrow     : integer  (forces guide_legend nrow)
#   axis_text_size  : numeric  pt size override for axis.text
#   strip_axis_text : "x" | "y" | "both" | "none"   (drop labels entirely)
#
# Anything not present is left untouched so the global theme override and
# the .rds's own theme apply.

.panels_merge_composite <- function(slot_block, panel_block) {
  base <- if (is.null(slot_block)) list() else slot_block
  over <- if (is.null(panel_block)) list() else panel_block
  modifyList(base, over)
}

.panels_apply_overrides <- function(plot_obj, overrides) {
  if (length(overrides) == 0) return(plot_obj)

  th <- list()
  if (!is.null(overrides$legend_position)) {
    th$legend.position <- overrides$legend_position
  }
  if (!is.null(overrides$legend_text_size)) {
    th$legend.text  <- ggplot2::element_text(
      size = as.numeric(overrides$legend_text_size)
    )
    th$legend.title <- ggplot2::element_text(
      size = as.numeric(overrides$legend_text_size) + 1
    )
  }
  if (!is.null(overrides$legend_key_size_cm)) {
    th$legend.key.size <- grid::unit(
      as.numeric(overrides$legend_key_size_cm), "cm"
    )
  }
  if (!is.null(overrides$axis_text_size)) {
    th$axis.text <- ggplot2::element_text(
      size = as.numeric(overrides$axis_text_size)
    )
  }
  if (!is.null(overrides$strip_axis_text)) {
    s <- tolower(overrides$strip_axis_text)
    if (s %in% c("x", "both", "all")) {
      th$axis.text.x <- ggplot2::element_blank()
    }
    if (s %in% c("y", "both", "all")) {
      th$axis.text.y <- ggplot2::element_blank()
    }
    # "all" also blanks axis titles + ticks — needed for ggraph-based
    # panels (fig09 network) where theme_void() in the source gets
    # overwritten by the composite axis.title default and the layout
    # x/y coords end up bleeding into the panel.
    if (s == "all") {
      th$axis.title  <- ggplot2::element_blank()
      th$axis.ticks  <- ggplot2::element_blank()
    }
  }
  if (length(th) > 0) {
    plot_obj <- plot_obj + do.call(ggplot2::theme, th)
  }

  if (!is.null(overrides$legend_ncol) || !is.null(overrides$legend_nrow)) {
    plot_obj <- plot_obj + ggplot2::guides(
      fill   = ggplot2::guide_legend(ncol = overrides$legend_ncol,
                                       nrow = overrides$legend_nrow),
      colour = ggplot2::guide_legend(ncol = overrides$legend_ncol,
                                       nrow = overrides$legend_nrow)
    )
  }
  plot_obj
}

# ---------------------------------------------------------------------------
# `per_group: true` slot expansion
# ---------------------------------------------------------------------------

# Expand any panel spec carrying `per_group: true` into one entry per level
# of cfg$metadata$group_cols[[1]] (controls excluded — sourced from
# cfg_group_levels() in R/13_manifest.R). Each expanded entry inherits the
# original kind / domain / composite block and gets a `group:` filter equal
# to the sanitised level name. Sanitisation must match the convention used
# when writing the underlying PNG — R/12_network.R::.network_build_chord
# does `gsub("[^A-Za-z0-9_-]+", "_", lvl)`, so we apply the same here so
# the filter matches the figures_index `group` field (which is derived from
# the filename by R/13).
#
# If no levels are resolvable (missing metadata, empty group column), the
# `per_group` entry is dropped silently — the slot's `required` semantics
# then drive whether the whole slot skips or partially renders.
.panels_expand_per_group <- function(panels, cfg, slug, sec_id) {
  out <- list()
  expanded <- 0L
  for (p in panels) {
    if (isTRUE(p$per_group %||% FALSE)) {
      levels <- tryCatch(cfg_group_levels(cfg), error = function(e) character(0))
      if (length(levels) == 0) {
        pipeline_log(cfg, sprintf(
          "Panels[%s/%s]: per_group spec has no resolvable levels — dropped",
          slug, sec_id
        ))
        next
      }
      p$per_group <- NULL
      p$tag       <- NULL  # let the resolver auto-assign LETTERS in order
      for (lvl in levels) {
        slot_p       <- p
        slot_p$group <- gsub("[^A-Za-z0-9_-]+", "_", lvl)
        out[[length(out) + 1]] <- slot_p
      }
      expanded <- expanded + length(levels)
    } else {
      out[[length(out) + 1]] <- p
    }
  }
  if (expanded > 0) {
    pipeline_log(cfg, sprintf(
      "Panels[%s/%s]: per_group expanded to %d panel(s)",
      slug, sec_id, expanded
    ))
  }
  out
}

# ---------------------------------------------------------------------------
# Manifest indexing
# ---------------------------------------------------------------------------

# Flatten every stage.figures entry into a single list, preserving the
# `kind`, `path`, and any filterable metadata (domain/metric/pair/...) the
# pipeline tagged it with.
.panels_index_figures <- function(manifest) {
  out <- list()
  stages <- manifest$stages %||% list()
  for (stage_name in names(stages)) {
    figs <- stages[[stage_name]]$figures %||% list()
    for (f in figs) {
      f$stage <- stage_name
      out[[length(out) + 1]] <- f
    }
  }
  out
}

.panels_index_tables <- function(manifest) {
  out <- list()
  stages <- manifest$stages %||% list()
  for (stage_name in names(stages)) {
    tbls <- stages[[stage_name]]$tables %||% list()
    for (t in tbls) {
      t$stage <- stage_name
      out[[length(out) + 1]] <- t
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Slot composition
# ---------------------------------------------------------------------------

.panels_build_slot <- function(slug, spec, figures_index, cfg, out_path) {
  if (!requireNamespace("cowplot", quietly = TRUE) ||
      !requireNamespace("magick",  quietly = TRUE)) {
    pipeline_log(cfg, sprintf(
      "Panels[%s]: cowplot/magick not available — skipping", slug
    ))
    return(NULL)
  }

  use_labels <- !identical(spec$labels, FALSE)

  # Two layout modes:
  # 1. Flat (legacy)  — spec$panels + spec$rows + spec$cols. Single grid.
  # 2. Sections (new) — spec$sections is a list of {panels, rows, cols,
  #    rel_height} groups stacked vertically. Used when one panel needs
  #    to span the full row width above a sub-grid (e.g. the diet-effects
  #    supplementaries: A on top, B-E in 2x2 below).
  if (!is.null(spec$sections) && length(spec$sections) > 0) {
    section_blocks <- list()
    section_rel    <- numeric()
    section_rows   <- integer()
    total_cols     <- 1L
    total_panels   <- 0L
    for (i in seq_along(spec$sections)) {
      sec    <- spec$sections[[i]]
      block  <- .panels_build_section_grid(slug, sec, figures_index, cfg,
                                            use_labels, sprintf("sec%d", i))
      if (is.null(block)) return(NULL)
      section_blocks[[length(section_blocks) + 1]] <- block$grid
      section_rel <- c(section_rel,
                       as.numeric(sec$rel_height %||% block$rows))
      section_rows <- c(section_rows, block$rows)
      total_cols   <- max(total_cols, as.integer(block$cols))
      total_panels <- total_panels + block$n_panels
    }
    composed <- cowplot::plot_grid(
      plotlist     = section_blocks,
      ncol         = 1,
      rel_heights  = section_rel
    )
    # Sum the *actual* row counts across sections so the output canvas
    # gets enough vertical room (rel_height is a relative split within
    # the canvas, NOT a row-count substitute). For 1+2x2 layout this
    # gives rows_eq = 3, matching the legacy 3x2 layout footprint.
    rows_eq <- sum(section_rows)
    cols    <- total_cols
    n_panels_for_log <- total_panels
  } else {
    block <- .panels_build_section_grid(slug, spec, figures_index, cfg,
                                         use_labels, "main")
    if (is.null(block)) return(NULL)
    composed         <- block$grid
    rows_eq          <- block$rows
    cols             <- block$cols
    n_panels_for_log <- block$n_panels
  }

  # Slot title is opt-in (default off). Journal figures get a caption,
  # not a baked-in italic title; including one steals canvas height and
  # clashes with the panel labels (A/B/C/...). Re-enable via
  # cfg$panels$show_title = TRUE.
  if (isTRUE(cfg$panels$show_title %||% FALSE)) {
    title <- spec$title %||% slug
    composed_final <- cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(title, fontface = "italic",
                                                size = 13),
      composed,
      ncol = 1,
      rel_heights = c(0.05, 0.95)
    )
  } else {
    composed_final <- composed
  }

  dpi    <- cfg$panels$dpi    %||% 300
  width  <- cfg$panels$width  %||% (5.5 * cols)
  height <- cfg$panels$height %||% (4.5 * rows_eq + 0.4)

  # Output format toggles. PNG is the default container (and the path returned
  # to the caller / logged); TIFF + PDF are optional siblings. Disable PNG via
  # cfg$panels$png = FALSE when you only need a publication TIFF.
  emit_png  <- isTRUE(cfg$panels$png  %||% TRUE)
  emit_tiff <- isTRUE(cfg$panels$tiff %||% FALSE)
  emit_pdf  <- isTRUE(cfg$panels$pdf  %||% FALSE)
  if (!emit_png && !emit_tiff && !emit_pdf) emit_png <- TRUE  # safety net

  if (emit_png) {
    ggplot2::ggsave(out_path, composed_final,
                     width = width, height = height,
                     dpi = dpi, bg = "white")
  }

  if (emit_tiff) {
    tiff_path <- sub("\\.png$", ".tiff", out_path)
    ggplot2::ggsave(tiff_path, composed_final,
                     width = width, height = height,
                     dpi = dpi, bg = "white",
                     device = grDevices::tiff,
                     compression = "lzw")
  }

  if (emit_pdf) {
    pdf_path <- sub("\\.png$", ".pdf", out_path)
    ggplot2::ggsave(pdf_path, composed_final,
                     width = width, height = height,
                     device = grDevices::cairo_pdf)
  }

  pipeline_log(cfg, sprintf(
    "Panels[%s]: %d panel(s) assembled -> %s",
    slug, n_panels_for_log, basename(out_path)
  ))
  out_path
}

# Build a single panel grid from a `panels` spec block. Used by both the
# legacy flat layout and the new section-based layout. Returns a list of
# (grid, rows, cols, n_panels) on success, NULL on failure / empty.
.panels_build_section_grid <- function(slug, sec_spec, figures_index, cfg,
                                        use_labels, sec_id) {
  panels <- sec_spec$panels %||% list()
  if (length(panels) == 0) return(NULL)

  # Expand any `per_group: true` spec into one panel per primary-grouping
  # level. Lets a slot like figS_chord_by_treatment adapt from
  # chicken_batch1's 3 treatments (Dulce / Reference_diet / Soyabean_meal)
  # to lung_microbiome's 2 (Exhale / Sputum) without hardcoded group names.
  panels <- .panels_expand_per_group(panels, cfg, slug, sec_id)

  resolved <- vector("list", length(panels))
  missing_required <- character()
  for (i in seq_along(panels)) {
    p <- panels[[i]]
    hit <- .panels_find_match(p, figures_index, cfg)
    if (is.null(hit)) {
      if (isFALSE(p$required %||% TRUE)) {
        resolved[[i]] <- NULL
      } else {
        missing_required <- c(
          missing_required,
          sprintf("%s (kind=%s)", p$tag %||% LETTERS[i], p$kind)
        )
      }
    } else {
      resolved[[i]] <- list(panel = p, entry = hit)
    }
  }
  if (length(missing_required) > 0) {
    pipeline_log(cfg, sprintf(
      "Panels[%s/%s]: missing required panels %s — slot skipped",
      slug, sec_id, paste(missing_required, collapse = ", ")
    ))
    return(NULL)
  }
  resolved <- Filter(Negate(is.null), resolved)
  if (length(resolved) == 0) {
    pipeline_log(cfg, sprintf(
      "Panels[%s/%s]: no usable panels — section skipped", slug, sec_id
    ))
    return(NULL)
  }

  plots <- list()
  tags  <- character()
  n_native <- 0L
  n_raster <- 0L
  for (item in resolved) {
    abs_path <- .panels_resolve_path(cfg, item$entry$path)
    if (!file.exists(abs_path)) {
      pipeline_log(cfg, sprintf(
        "Panels[%s/%s]: PNG missing on disk (%s) — slot skipped",
        slug, sec_id, abs_path
      ))
      return(NULL)
    }
    # Prefer a sibling .rds containing the original ggplot object so we can
    # re-render at the composite canvas resolution. Fall back to reading the
    # PNG via magick when no .rds is available (non-ggplot outputs like
    # pheatmap / circlize / UpSet are still raster).
    rds_path <- sub("\\.(png|pdf)$", ".rds", abs_path,
                    ignore.case = TRUE)
    plot_obj <- NULL
    if (rds_path != abs_path && file.exists(rds_path)) {
      plot_obj <- tryCatch(readRDS(rds_path), error = function(e) NULL)
      if (is.null(plot_obj) || !inherits(plot_obj, "ggplot")) {
        plot_obj <- NULL  # fall through to raster
      }
    }
    if (is.null(plot_obj)) {
      img <- tryCatch(magick::image_read(abs_path), error = function(e) NULL)
      if (is.null(img)) {
        pipeline_log(cfg, sprintf(
          "Panels[%s/%s]: failed to read %s — slot skipped",
          slug, sec_id, abs_path
        ))
        return(NULL)
      }
      plot_obj <- cowplot::ggdraw() + cowplot::draw_image(img)
      n_raster <- n_raster + 1L
    } else {
      # Apply the composite theme override so per-stage plots tuned for
      # standalone viewing (12pt text, big legend) shrink to fit cleanly
      # at composite scale.
      plot_obj <- plot_obj + .panels_composite_theme(cfg)
      # Per-slot or per-panel composite overrides from the slots YAML
      # (e.g. legend_position: right, legend_ncol: 2). Slot-level applies
      # to every panel in the slot; panel-level wins when both are set.
      overrides <- .panels_merge_composite(sec_spec$composite,
                                            item$panel$composite)
      plot_obj <- .panels_apply_overrides(plot_obj, overrides)
      n_native <- n_native + 1L
    }
    plots[[length(plots) + 1]] <- plot_obj
    tags <- c(tags, item$panel$tag %||% LETTERS[length(plots)])
  }
  if (n_native + n_raster > 0) {
    pipeline_log(cfg, sprintf(
      "Panels[%s/%s]: %d native (.rds) + %d raster (.png)",
      slug, sec_id, n_native, n_raster
    ))
  }

  rows <- max(1L, as.integer(sec_spec$rows %||% 1L))
  cols <- max(1L, as.integer(sec_spec$cols %||% length(plots)))
  if (rows * cols < length(plots)) {
    rows <- ceiling(length(plots) / cols)
  }

  grid <- cowplot::plot_grid(
    plotlist       = plots,
    labels         = if (use_labels) tags else NULL,
    label_size     = 16,
    label_fontface = "bold",
    nrow           = rows,
    ncol           = cols
  )
  list(grid = grid, rows = rows, cols = cols, n_panels = length(plots))
}

# ---------------------------------------------------------------------------
# Slot panel spec → manifest entry
# ---------------------------------------------------------------------------

# Filter the figures index by every field present on the panel spec (kind,
# domain, metric, pair, organism, rank). First match wins.
.panels_find_match <- function(panel_spec, figures_index, cfg) {
  filters <- panel_spec[setdiff(names(panel_spec),
                                 c("tag", "required", "composite"))]
  if (is.null(filters$kind)) return(NULL)
  for (entry in figures_index) {
    match <- TRUE
    for (f in names(filters)) {
      want <- filters[[f]]
      have <- entry[[f]]
      if (is.null(have)) { match <- FALSE; break }
      if (f == "pair") {
        # pair is a 2-element list; compare as set so order doesn't matter
        if (length(want) != length(have) ||
            !all(sort(unlist(want)) == sort(unlist(have)))) {
          match <- FALSE; break
        }
      } else if (length(want) > 1) {
        if (length(have) != length(want) ||
            !all(unlist(want) == unlist(have))) {
          match <- FALSE; break
        }
      } else {
        if (as.character(want) != as.character(have)) {
          match <- FALSE; break
        }
      }
    }
    if (match) return(entry)
  }
  NULL
}

# Resolve a manifest-relative path back to an absolute path on disk. The
# pipeline writes paths as `<datasets_dir>/...` or `<figures_dir>/...`
# under project_root; manifest stores them as relative-from-project-root.
.panels_resolve_path <- function(cfg, rel_path) {
  if (file.exists(rel_path)) return(normalizePath(rel_path, mustWork = FALSE))
  candidate <- file.path(cfg$project_root, rel_path)
  if (file.exists(candidate)) return(normalizePath(candidate, mustWork = FALSE))
  candidate
}

# ---------------------------------------------------------------------------
# Supplementary tables XLSX
# ---------------------------------------------------------------------------

# Reads every (table_entry kind, csv path) pair the manifest exposes,
# builds a multi-tab workbook with one tab per kind. Sheet names are
# truncated + sanitised (Excel cap = 31 chars, no /\?*[]: characters).
# Adds an Index sheet listing (kind, sheet, description, n_rows, source).
.panels_build_supplementary_xlsx <- function(tables_index, cfg, out_path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    pipeline_log(cfg,
      "Panels XLSX: openxlsx not available — skipping supplementary_tables.xlsx")
    return(invisible(NULL))
  }
  if (length(tables_index) == 0) {
    pipeline_log(cfg, "Panels XLSX: no table entries in manifest — skipping")
    return(invisible(NULL))
  }

  wb <- openxlsx::createWorkbook()
  index_rows <- list()
  used_sheet_names <- character()

  for (t in tables_index) {
    abs_csv <- .panels_resolve_path(cfg, t$path %||% "")
    if (!file.exists(abs_csv)) next
    df <- tryCatch(
      readr::read_csv(abs_csv, show_col_types = FALSE, progress = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) next

    sheet <- .panels_safe_sheet_name(t$kind %||% basename(abs_csv),
                                       used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet)

    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df, headerStyle = openxlsx::createStyle(
      textDecoration = "bold", fgFill = "#E8EEF7"
    ))
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet,
                            cols = seq_len(ncol(df)), widths = "auto")

    index_rows[[length(index_rows) + 1]] <- data.frame(
      kind         = t$kind %||% NA_character_,
      sheet        = sheet,
      description  = t$description %||% NA_character_,
      n_rows       = nrow(df),
      source_path  = t$path %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (length(index_rows) == 0) {
    pipeline_log(cfg, "Panels XLSX: no readable tables found — skipping")
    return(invisible(NULL))
  }

  index_df <- do.call(rbind, index_rows)
  openxlsx::addWorksheet(wb, "Index")
  openxlsx::worksheetOrder(wb) <- c(length(used_sheet_names) + 1L,
                                     seq_len(length(used_sheet_names)))
  openxlsx::writeData(wb, "Index", index_df,
                      headerStyle = openxlsx::createStyle(
                        textDecoration = "bold", fgFill = "#D9E2F3"
                      ))
  openxlsx::freezePane(wb, "Index", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Index", cols = seq_len(ncol(index_df)),
                          widths = "auto")

  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  pipeline_log(cfg, sprintf(
    "Panels XLSX: %d tables -> %s", nrow(index_df), basename(out_path)
  ))
  invisible(out_path)
}

.panels_safe_sheet_name <- function(raw, used) {
  s <- gsub("[\\\\/\\?\\*\\[\\]:]", "_", as.character(raw))
  if (nchar(s) > 31) s <- substr(s, 1, 31)
  # Deduplicate against already-used names (Excel sheet names are unique).
  base <- s
  i <- 2L
  while (s %in% used) {
    suffix <- sprintf("_%d", i)
    s <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1L
  }
  s
}
