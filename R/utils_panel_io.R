# utils_panel_io.R — Save a ggplot to PNG AND a sibling .rds so the
# publication panels stage (R/13_panels.R) can re-render the figure
# natively at the composite canvas resolution instead of rasterising the
# saved PNG via magick.
#
# Used as a drop-in replacement for ggplot2::ggsave() at every site whose
# output may end up inside a panel slot. The PNG is still emitted (and
# remains the canonical artefact recorded in the manifest); the .rds is
# additive and ignored by anything that doesn't know to look for it.

# Axis / legend / title text for a metadata column: "sample_type" -> "Sample type".
# Display only; the column itself keeps its name.
pretty_var <- function(x) {
  s <- gsub("_", " ", x, fixed = TRUE)
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# y-scale expansion that leaves room above the data for an in-panel stats
# label drawn at y = Inf with `n_lines` lines of text, so the label never
# covers points, whiskers or ellipses.
stats_headroom <- function(n_lines = 1) {
  ggplot2::expansion(mult = c(0.05, 0.06 + 0.1 * n_lines))
}

# Number of text lines in a label.
n_text_lines <- function(label) length(strsplit(label, "\n", fixed = TRUE)[[1]])

save_panel_ggplot <- function(filename, plot, width, height,
                              dpi = 300, bg = "white", ...,
                              emit_rds = TRUE) {
  # figures.set (utils_figures.R): a figure outside the run's set is not written.
  if (exists("figure_wanted", mode = "function") && !figure_wanted(filename)) {
    return(invisible(NULL))
  }
  ggplot2::ggsave(filename, plot,
                   width = width, height = height,
                   dpi = dpi, bg = bg, limitsize = FALSE, ...)
  if (isTRUE(emit_rds)) {
    rds_path <- sub("\\.(png|pdf|jpg|jpeg|svg|tif|tiff)$", ".rds",
                    filename, ignore.case = TRUE)
    if (rds_path != filename) {
      tryCatch(saveRDS(plot, rds_path),
               error = function(e) {
                 message(sprintf(
                   "save_panel_ggplot: failed to write %s (%s)",
                   rds_path, conditionMessage(e)
                 ))
               })
    }
  }
  invisible(filename)
}

# Save an htmlwidget, tolerating hosts without pandoc. selfcontained = TRUE
# needs pandoc (via rmarkdown); when it is missing saveWidget() throws, which
# used to abort the whole stage. Retry with selfcontained = FALSE (assets in a
# sibling <name>_files/ dir) and log a line; if that fails too, log and skip
# only this widget. Returns TRUE/FALSE invisibly.
save_widget_safe <- function(widget, file, selfcontained = TRUE, cfg = NULL) {
  say <- function(msg) {
    if (exists("pipeline_log", mode = "function")) pipeline_log(cfg, msg)
    else message(msg)
  }
  try_save <- function(sc) {
    tryCatch({
      suppressWarnings(htmlwidgets::saveWidget(widget, file, selfcontained = sc))
      TRUE
    }, error = function(e) conditionMessage(e))
  }
  res <- try_save(selfcontained)
  if (isTRUE(res)) return(invisible(TRUE))
  if (isTRUE(selfcontained)) {
    say(sprintf("saveWidget(selfcontained = TRUE) failed for %s (%s) — retrying with selfcontained = FALSE",
                basename(file), res))
    res <- try_save(FALSE)
    if (isTRUE(res)) return(invisible(TRUE))
  }
  say(sprintf("saveWidget failed for %s (%s) — skipping html", basename(file), res))
  invisible(FALSE)
}
