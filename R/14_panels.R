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

  manifest_path <- file.path(cfg$project_root, dirname(cfg$outputs$log_file),
                              "manifest.json")
  if (!file.exists(manifest_path)) {
    pipeline_log(cfg, sprintf(
      "Panels: manifest %s missing — run the manifest stage first; skipping",
      manifest_path
    ))
    return(invisible(NULL))
  }

  slots <- yaml::read_yaml(yaml_path)
  manifest <- jsonlite::read_json(manifest_path)
  figures_index <- .panels_index_figures(manifest)
  tables_index  <- .panels_index_tables(manifest)

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
                         pattern = "\\.(png|pdf)$", full.names = TRUE)
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
    return(normalizePath(dirname(dirname(sub("^--file=", "", file_arg[1])))))
  }
  # Sourced interactively — best effort: current wd.
  normalizePath(getwd())
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
  panels <- spec$panels %||% list()
  if (length(panels) == 0) return(NULL)

  resolved <- vector("list", length(panels))
  missing_required <- character()
  for (i in seq_along(panels)) {
    p <- panels[[i]]
    hit <- .panels_find_match(p, figures_index, cfg)
    if (is.null(hit)) {
      if (isFALSE(p$required %||% TRUE)) {
        resolved[[i]] <- NULL  # optional panel just absent
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
      "Panels[%s]: missing required panels %s — slot skipped",
      slug, paste(missing_required, collapse = ", ")
    ))
    return(NULL)
  }
  resolved <- Filter(Negate(is.null), resolved)
  if (length(resolved) == 0) {
    pipeline_log(cfg, sprintf(
      "Panels[%s]: no usable panels — slot skipped", slug
    ))
    return(NULL)
  }

  plots <- list()
  tags  <- character()
  for (item in resolved) {
    abs_path <- .panels_resolve_path(cfg, item$entry$path)
    if (!file.exists(abs_path)) {
      pipeline_log(cfg, sprintf(
        "Panels[%s]: PNG missing on disk (%s) — slot skipped",
        slug, abs_path
      ))
      return(NULL)
    }
    img <- tryCatch(magick::image_read(abs_path), error = function(e) NULL)
    if (is.null(img)) {
      pipeline_log(cfg, sprintf(
        "Panels[%s]: failed to read %s — slot skipped", slug, abs_path
      ))
      return(NULL)
    }
    plots[[length(plots) + 1]] <- cowplot::ggdraw() +
      cowplot::draw_image(img)
    tags <- c(tags, item$panel$tag %||% LETTERS[length(plots)])
  }

  rows <- max(1L, as.integer(spec$rows %||% 1L))
  cols <- max(1L, as.integer(spec$cols %||% length(plots)))
  # cowplot tiles plots in row-major order. We allow rows*cols < length(plots)
  # by overflowing to extra rows (rare; only when slot YAML is misconfigured).
  if (rows * cols < length(plots)) {
    rows <- ceiling(length(plots) / cols)
  }

  # Standalone "no-label" slots (`labels: false` in YAML) get a single
  # full-width composition with no A/B/C tag overlay — the slot title at
  # the top is the only annotation. Used for relative-abundance stacked
  # bars, heatmaps, chord, sankey, network.
  use_labels <- !identical(spec$labels, FALSE)
  composed <- cowplot::plot_grid(
    plotlist = plots,
    labels   = if (use_labels) tags else NULL,
    label_size = 16,
    label_fontface = "bold",
    nrow = rows, ncol = cols
  )

  title <- spec$title %||% slug
  composed_with_title <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(title, fontface = "italic",
                                              size = 13),
    composed,
    ncol = 1,
    rel_heights = c(0.05, 0.95)
  )

  dpi    <- cfg$panels$dpi    %||% 300
  width  <- cfg$panels$width  %||% (5.5 * cols)
  height <- cfg$panels$height %||% (4.5 * rows + 0.4)

  ggplot2::ggsave(out_path, composed_with_title,
                   width = width, height = height,
                   dpi = dpi, bg = "white")

  if (isTRUE(cfg$panels$pdf %||% FALSE)) {
    pdf_path <- sub("\\.png$", ".pdf", out_path)
    ggplot2::ggsave(pdf_path, composed_with_title,
                     width = width, height = height,
                     device = grDevices::cairo_pdf)
  }

  pipeline_log(cfg, sprintf(
    "Panels[%s]: %d panel(s) assembled -> %s",
    slug, length(plots), basename(out_path)
  ))
  out_path
}

# ---------------------------------------------------------------------------
# Slot panel spec → manifest entry
# ---------------------------------------------------------------------------

# Filter the figures index by every field present on the panel spec (kind,
# domain, metric, pair, organism, rank). First match wins.
.panels_find_match <- function(panel_spec, figures_index, cfg) {
  filters <- panel_spec[setdiff(names(panel_spec),
                                 c("tag", "required"))]
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
