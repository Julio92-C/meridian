# utils_palette.R — Centralised group-level color resolver.
#
# Every stage that colors a categorical metadata variable (typically
# group_cols[[1]]) routes through resolve_top_level_colors() so colors
# stay consistent across alpha-violin, beta-PCoA, DA effect plots, and
# resistome/virulome/mobilome bars.
#
# Lookup order for a stage's per-group palette helper:
#   1. cfg$colors[[variable]]          — explicit per-variable override (map or
#                                         paletteer palette name)
#   2. cfg$<stage>$group_colors        — legacy per-stage override
#   3. cfg$taxonomy$treatment_colors   — legacy (chicken_batch1)
#   4. paletteer palette name from cfg$<stage>$palette (or stage-specific
#      default, e.g. ggsci::default_nejm).

`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve a named palette for `levels` of a metadata `variable` using the
# top-level cfg$colors block. Returns a named character vector of hex colors,
# or NULL if the top-level block has no usable entry for this variable
# (callers then fall through to legacy per-stage lookup).
#
# cfg$colors entry conventions:
#   colors:
#     sample_type:                  # named map level -> color
#       Exhale: "#1F77B4"
#       Sputum: "#FF7F0E"
#     ID: "ggsci::default_jama"     # OR a paletteer_d palette name (string)
resolve_top_level_colors <- function(variable, levels, cfg) {
  if (is.null(variable) || is.null(cfg$colors)) return(NULL)
  entry <- cfg$colors[[variable]]
  if (is.null(entry)) return(NULL)

  # Case A: single string -> treat as paletteer palette name
  if (is.character(entry) && length(entry) == 1L && is.null(names(entry))) {
    pal <- tryCatch(
      as.character(paletteer::paletteer_d(entry)),
      error = function(e) NULL
    )
    if (is.null(pal) || length(pal) == 0L) {
      pal <- grDevices::hcl.colors(length(levels), palette = "Dark 3")
    }
    pal <- rep_len(pal, length(levels))
    names(pal) <- levels
    return(pal)
  }

  # Case B: named map level -> hex. Accepts list or named character vector.
  user_map <- unlist(entry)
  if (is.null(names(user_map))) return(NULL)
  pal <- user_map[levels]
  if (length(pal) == length(levels) && all(!is.na(pal))) {
    names(pal) <- levels
    return(pal)
  }
  # Incomplete map -> let callers handle the fallback (so partial overrides
  # don't silently break a stage). Logging is left to the caller, which
  # has the stage context.
  NULL
}
