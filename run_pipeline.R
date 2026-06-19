#!/usr/bin/env Rscript
# run_pipeline.R — master runner.
# Usage:
#   Rscript run_pipeline.R projects/chicken_batch1/config.yaml

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript run_pipeline.R <config.yaml>")
cfg_path <- args[[1]]

ofile    <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
script   <- ofile %||% (if (length(file_arg)) file_arg[[1]] else NULL)
repo_root <- if (!is.null(script)) dirname(normalizePath(script)) else getwd()

# Source every R/NN_*.R stage module in order, plus any R/utils_*.R helpers.
# Stage files sort numerically; helpers are sourced afterwards so they can
# rely on anything defined in 00_setup.R.
mods <- sort(list.files(file.path(repo_root, "R"),
                        pattern = "^\\d+_.*\\.R$",
                        full.names = TRUE))
utils <- sort(list.files(file.path(repo_root, "R"),
                         pattern = "^utils_.*\\.R$",
                         full.names = TRUE))
invisible(lapply(c(mods, utils), source))

cfg <- load_config(cfg_path)

# Format a duration in seconds as a human-readable string.
# < 60s  -> "12.3s"
# < 1h   -> "3m 12s"
# >= 1h  -> "1h 23m 12s"
fmt_duration <- function(secs) {
  if (is.na(secs) || secs < 0) return(NA_character_)
  if (secs < 60) return(sprintf("%.1fs", secs))
  h <- as.integer(secs %/% 3600)
  m <- as.integer((secs %% 3600) %/% 60)
  s <- as.integer(round(secs %% 60))
  if (h > 0) sprintf("%dh %dm %ds", h, m, s) else sprintf("%dm %ds", m, s)
}

pipeline_start <- Sys.time()
pipeline_log(cfg, sprintf("Starting pipeline for study '%s'", cfg$study$id))

# Cross-module stats accumulators: clear at startup so re-runs don't
# append to stale rows. Each surface is filled by ge_alpha_kw / similar
# helpers and finalised by R/13 when building the manifest.
padj_summary_reset(cfg, "ge_alpha_kw")

# Build the ordered list of stages that will actually run, so we can show
# "[i/n  pct%] stage_name" progress in the log as each one starts/finishes.
enabled_names <- c(
  "load_inputs",
  if (isTRUE(cfg$stages$clean_data))             "clean_data",
  if (isTRUE(cfg$stages$normalisation))          "normalisation",
  if (isTRUE(cfg$stages$taxonomy))               "taxonomy",
  if (isTRUE(cfg$stages$relative_abundance))     "relative_abundance",
  if (isTRUE(cfg$stages$alpha_diversity))        "alpha_diversity",
  if (isTRUE(cfg$stages$beta_diversity))         "beta_diversity",
  if (isTRUE(cfg$stages$differential_abundance)) "differential_abundance",
  if (isTRUE(cfg$stages$resistome))              "resistome",
  if (isTRUE(cfg$stages$virulome))               "virulome",
  if (isTRUE(cfg$stages$mobilome))               "mobilome",
  if (isTRUE(cfg$stages$network))                "network",
  if (isTRUE(cfg$stages$report))                 "report",
  # Panels runs BEFORE manifest so the manifest catalogues the rendered
  # composites (the build_stages_index() helper in R/13 scans the output
  # dirs panels just wrote to).
  if (isTRUE(cfg$stages$panels   %||% TRUE))     "panels",
  if (isTRUE(cfg$stages$manifest %||% TRUE))     "manifest"
)
total_stages <- length(enabled_names)
stage_idx    <- 0L
# Accumulator for the per-stage breakdown printed at the end of the run.
stage_times  <- list()
pipeline_log(cfg, sprintf("Plan: %d stage(s) — %s",
                          total_stages, paste(enabled_names, collapse = ", ")))

run_stage <- function(name, fn) {
  if (!(name %in% enabled_names)) return(invisible(NULL))
  stage_idx <<- stage_idx + 1L
  pct <- round(100 * (stage_idx - 1) / total_stages)
  pipeline_log(cfg, sprintf("[%d/%d %3d%%] %s — start",
                            stage_idx, total_stages, pct, name))
  t0 <- Sys.time()
  result <- fn()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  stage_times[[name]] <<- elapsed
  pipeline_log(cfg, sprintf("[%d/%d %3d%%] %s — done in %s",
                            stage_idx, total_stages, pct, name,
                            fmt_duration(elapsed)))
  result
}

inputs  <- run_stage("load_inputs", function() load_inputs(cfg))
cleaned <- run_stage("clean_data",  function() clean_data(inputs, cfg))

if (!is.null(cleaned)) {
  run_stage("normalisation",          function() normalise_data(cleaned, cfg))
  run_stage("taxonomy",               function() run_taxonomy(cleaned, cfg))
  run_stage("relative_abundance",     function() run_relative_abundance(cleaned, cfg))
  run_stage("alpha_diversity",        function() run_alpha_diversity(cleaned, cfg))
  run_stage("beta_diversity",         function() run_beta_diversity(cleaned, cfg))
  run_stage("differential_abundance", function() run_aldex(cleaned, cfg))
  run_stage("resistome",              function() run_resistome(cleaned, cfg))
  run_stage("virulome",               function() run_virulome(cleaned, cfg))
  run_stage("mobilome",               function() run_mobilome(cleaned, cfg))
  run_stage("network",                function() run_network(cleaned, cfg))
}

run_stage("report", function() {
  report_src <- file.path(repo_root, "templates", "report.qmd")
  report_out <- file.path(cfg$project_root, cfg$outputs$report_html)
  if (file.exists(report_src) && nzchar(Sys.which("quarto"))) {
    pipeline_log(cfg, "Rendering Quarto HTML report")
    cfg_abs <- normalizePath(cfg_path, winslash = "/", mustWork = TRUE)
    dir.create(dirname(report_out), recursive = TRUE, showWarnings = FALSE)
    # Render in place — passing --output-dir / --output on Windows breaks
    # Quarto's lookup of its own report_files/ support dir during embed.
    # Format (`dashboard`) is set in the .qmd YAML; don't pass --to here.
    system2("quarto", c("render", report_src,
                        "-P", paste0("config=", cfg_abs)))
    stem      <- tools::file_path_sans_ext(basename(report_src))
    default_h <- file.path(dirname(report_src), paste0(stem, ".html"))
    files_dir <- file.path(dirname(report_src), paste0(stem, "_files"))
    if (file.exists(default_h)) {
      file.copy(default_h, report_out, overwrite = TRUE)
      file.remove(default_h)
    }
    # embed-resources inlines everything; the staging _files/ dir is not
    # needed once the HTML is copied to its destination.
    if (dir.exists(files_dir)) unlink(files_dir, recursive = TRUE, force = TRUE)
  } else {
    pipeline_log(cfg, "Quarto not found or report.qmd missing — skipping HTML")
  }
})

# Publication panels + supplementary tables XLSX. Runs BEFORE manifest so
# the manifest can catalogue the rendered composites (R/14 builds its
# figures_index in-process via build_stages_index() — no manifest.json
# round-trip). Gated by cfg$stages$panels (default TRUE).
run_stage("panels", function() run_panels(cfg))

pipeline_end <- Sys.time()

# Emit manifest.json after panels so its build_stage_panels() entry lists
# the actual TIFF/PNG/PDF files that were just written. Gated by
# cfg$stages$manifest (default TRUE); consults stage_times to mark each
# stage complete / skipped / failed.
run_stage("manifest", function() {
  write_manifest(
    cfg,
    stage_times        = stage_times,
    pipeline_start     = pipeline_start,
    pipeline_end       = pipeline_end,
    pipeline_repo_root = repo_root
  )
})

total_elapsed <- as.numeric(difftime(Sys.time(), pipeline_start, units = "secs"))
pipeline_log(cfg, sprintf("Pipeline finished — %d/%d stages complete in %s",
                          stage_idx, total_stages, fmt_duration(total_elapsed)))

# Per-stage breakdown, sorted slowest -> fastest, with a Total footer.
# Helps spot bottlenecks without scrolling the run log line by line.
if (length(stage_times) > 0) {
  st_secs <- unlist(stage_times)
  ord     <- order(st_secs, decreasing = TRUE)
  pct     <- 100 * st_secs / total_elapsed
  width   <- max(nchar(c(names(st_secs), "Total")))
  lines   <- sprintf("  %-*s  %8s  %5.1f%%",
                     width, names(st_secs)[ord],
                     vapply(st_secs[ord], fmt_duration, character(1)),
                     pct[ord])
  divider <- paste0("  ", strrep("-", width + 2 + 8 + 2 + 6))
  total   <- sprintf("  %-*s  %8s  %5.1f%%",
                     width, "Total",
                     fmt_duration(total_elapsed), 100.0)
  pipeline_log(cfg, paste0("Stage breakdown (slowest first):\n",
                           paste(c(lines, divider, total), collapse = "\n")))
}
