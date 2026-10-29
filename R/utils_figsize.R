# utils_figsize.R — one place to decide how big a figure may be and how it is written to disk.
#
# Why this exists. Heatmap call sites sized their canvas from the data and nothing else, for example
# `width = max(8, 0.55 * n_cols + 4)` with `fontsize_row = 10` fixed. With 60 columns that is a 37-inch
# canvas; placed on a 6.5-inch page it is scaled to about 18%, so 10 pt text renders under 2 pt. Some
# figures were also written through a bare png() with no `res`, which Word imports at 96 dpi, making the
# figure enormous on the page and then shrinking it further. Both faults are visible in the lung cohort:
# taxa_heatmap.png was 19290 x 3900 px with no DPI recorded.
#
# The rules here:
#   * a canvas never exceeds `fig_max_width_in` (default 2x the printable width), so nothing is scaled
#     into illegibility;
#   * font size follows the cell size, so text stays readable at whatever size the canvas ends up;
#   * `res` is always set, so the DPI is recorded in the file;
#   * when the data do not fit, `fig_split_plan()` says how to split them rather than squeezing.
#
# Config (any project, all optional; defaults suit an A4 portrait thesis page):
#   figures:
#     page_width_in:     6.5    # printable text width of the target document
#     max_width_factor:  2      # a figure may be this many page widths before it must split
#     dpi:               300
#     min_fontsize:      6      # never render text smaller than this
#     max_fontsize:      12
#     heatmap:
#       max_samples:       40   # above this many columns, split (mirrors relative_abundance.per_sample)
#       samples_per_panel: 30
#       split_by:          null # metadata column; default is the group column
#       max_features:      null # optional cap on rows; the remainder goes to a supplementary table

.fig_cfg <- function(cfg) {
  f <- cfg$figures %||% list()
  list(page_width_in    = f$page_width_in    %||% 6.5,
       max_width_factor = f$max_width_factor %||% 2,
       dpi              = f$dpi              %||% 300,
       min_fontsize     = f$min_fontsize     %||% 6,
       max_fontsize     = f$max_fontsize     %||% 12,
       heatmap          = f$heatmap          %||% list())
}

# Largest canvas we allow, in inches.
fig_max_width <- function(cfg) {
  f <- .fig_cfg(cfg)
  f$page_width_in * f$max_width_factor
}

# Canvas size and matching font size for a matrix of n_rows x n_cols cells.
#
# `per_cell_*` are the ideal inches per cell (the old call sites' coefficients); `pad_*` the margin for
# labels, legend and dendrogram. The ideal size is computed first, then capped: if the cap bites, the
# cells shrink and the font shrinks with them, down to `min_fontsize`. Returns width, height, fontsize
# and `fits`, which is FALSE when the cap bit — the caller should then split (see fig_split_plan).
fig_heatmap_size <- function(n_rows, n_cols, cfg,
                             per_cell_w = 0.45, per_cell_h = 0.22,
                             pad_w = 4, pad_h = 2.5,
                             min_w = 8, min_h = 6) {
  f <- .fig_cfg(cfg)
  max_w <- fig_max_width(cfg)
  max_h <- max_w * 1.4                       # keep the aspect within something a page can take
  ideal_w <- max(min_w, per_cell_w * n_cols + pad_w)
  ideal_h <- max(min_h, per_cell_h * n_rows + pad_h)
  w <- min(ideal_w, max_w)
  h <- min(ideal_h, max_h)
  # how much the cells were squeezed; font follows the tighter of the two axes
  shrink <- min(1, (w - pad_w) / max(1e-6, ideal_w - pad_w),
                  (h - pad_h) / max(1e-6, ideal_h - pad_h))
  base <- f$max_fontsize * shrink
  list(width = w, height = h,
       fontsize = max(f$min_fontsize, min(f$max_fontsize, base)),
       dpi = f$dpi,
       fits = (ideal_w <= max_w && ideal_h <= max_h))
}

# How to split a heatmap whose columns do not fit on one canvas.
#
# Mirrors the contract relative_abundance.per_sample already uses (max_samples / layout), so the two
# behave the same way: above `max_samples` columns the figure is drawn once per level of `split_by`,
# and within a level in panels of `samples_per_panel`.
# Returns list(split = TRUE/FALSE, by = <column or NULL>, per_panel = <int>, reason = <text>).
fig_split_plan <- function(n_cols, cfg, group_col = NULL) {
  h <- .fig_cfg(cfg)$heatmap
  max_samples <- h$max_samples %||% 40
  if (n_cols <= max_samples) {
    return(list(split = FALSE, by = NULL, per_panel = n_cols,
                reason = sprintf("%d columns (<= max_samples = %d)", n_cols, max_samples)))
  }
  list(split = TRUE,
       by = h$split_by %||% group_col,
       per_panel = h$samples_per_panel %||% 30,
       reason = sprintf("%d columns (> max_samples = %d)", n_cols, max_samples))
}

# Rows kept when a cap is configured; the caller writes the remainder to a supplementary table.
fig_cap_features <- function(x, cfg, by = NULL) {
  cap <- .fig_cfg(cfg)$heatmap$max_features
  if (is.null(cap) || nrow(x) <= cap) return(list(keep = x, dropped = NULL))
  ord <- if (is.null(by)) order(-rowSums(abs(as.matrix(x)), na.rm = TRUE)) else order(-by)
  list(keep = x[utils::head(ord, cap), , drop = FALSE],
       dropped = x[utils::tail(ord, nrow(x) - cap), , drop = FALSE])
}

# A png() device that always records its resolution. Use instead of a bare png().
fig_png <- function(file, width_in, height_in, cfg) {
  grDevices::png(file, width = width_in, height = height_in, units = "in",
                 res = .fig_cfg(cfg)$dpi, bg = "white")
}

# Draw a pheatmap at a size the page can take, splitting when the columns do not fit.
#
# `mat` is the matrix; `...` is passed to pheatmap unchanged, so every call site keeps its own colours,
# breaks, clustering and annotations. fontsize_row / fontsize_col / fontsize are supplied from the
# computed size unless the caller sets them explicitly.
# `groups` (named by column) drives the split; when absent, the columns are split into equal panels.
# Returns the files written.
fig_save_heatmap <- function(mat, file, cfg, groups = NULL,
                             per_cell_w = 0.45, per_cell_h = 0.22,
                             pad_w = 4, pad_h = 2.5, log_label = "heatmap", ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) return(character(0))
  plan <- fig_split_plan(ncol(mat), cfg, group_col = if (is.null(groups)) NULL else "group")
  panels <- list()
  if (!plan$split) {
    panels[[basename(file)]] <- seq_len(ncol(mat))
  } else {
    idx <- if (!is.null(groups)) split(seq_len(ncol(mat)), groups[colnames(mat)]) else
      list(all = seq_len(ncol(mat)))
    for (lv in names(idx)) {
      cols <- idx[[lv]]
      chunks <- split(cols, ceiling(seq_along(cols) / plan$per_panel))
      for (k in seq_along(chunks)) {
        tag <- if (length(idx) > 1 && length(chunks) > 1) sprintf("_%s_%d", .fig_slug(lv), k)
               else if (length(idx) > 1) sprintf("_%s", .fig_slug(lv))
               else sprintf("_%d", k)
        panels[[sub("\\.png$", paste0(tag, ".png"), basename(file))]] <- chunks[[k]]
      }
    }
  }
  out <- character(0)
  for (nm in names(panels)) {
    sub <- mat[, panels[[nm]], drop = FALSE]
    sz <- fig_heatmap_size(nrow(sub), ncol(sub), cfg, per_cell_w, per_cell_h, pad_w, pad_h)
    args <- list(sub, filename = file.path(dirname(file), nm),
                 width = sz$width, height = sz$height,
                 fontsize = sz$fontsize, fontsize_row = sz$fontsize, fontsize_col = sz$fontsize)
    extra <- list(...)
    args[names(extra)] <- extra                      # caller's settings win
    do.call(pheatmap::pheatmap, args)
    out <- c(out, file.path(dirname(file), nm))
  }
  if (!is.null(cfg$.log_fn) || exists("pipeline_log")) {
    msg <- sprintf("%s: %d file(s); %s", log_label, length(out), plan$reason)
    try(pipeline_log(cfg, msg), silent = TRUE)
  }
  out
}

.fig_slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
