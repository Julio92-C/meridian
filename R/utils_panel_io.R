# utils_panel_io.R — Save a ggplot to PNG AND a sibling .rds so the
# publication panels stage (R/13_panels.R) can re-render the figure
# natively at the composite canvas resolution instead of rasterising the
# saved PNG via magick.
#
# Used as a drop-in replacement for ggplot2::ggsave() at every site whose
# output may end up inside a panel slot. The PNG is still emitted (and
# remains the canonical artefact recorded in the manifest); the .rds is
# additive and ignored by anything that doesn't know to look for it.

save_panel_ggplot <- function(filename, plot, width, height,
                              dpi = 300, bg = "white", ...,
                              emit_rds = TRUE) {
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
