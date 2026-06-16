# Smoke-test the native-render panels path end-to-end against chicken_batch1.
#
# Strategy: copy projects/chicken_batch1/config.yaml to a sibling override
# that disables the report stage (needs pandoc) and the differential_abundance
# stage (~2 min ALDEx2 burn that isn't on the rds critical path), then invoke
# run_pipeline.R with the override config.
#
# Verifies after the run:
#   1. .rds siblings exist next to .png in the figure subdirs we expect.
#   2. The panels stage logged at least one "N native (.rds) + M raster (.png)"
#      line with N > 0 (proving the new path was taken).

src_cfg  <- "projects/chicken_batch1/config.yaml"
test_cfg <- "projects/chicken_batch1/config_smoke_panels.yaml"

stopifnot(file.exists(src_cfg))
lines <- readLines(src_cfg)
lines <- sub("^(\\s*differential_abundance:)\\s*true",
             "\\1 false", lines)
lines <- sub("^(\\s*report:)\\s*true", "\\1 false", lines)
writeLines(lines, test_cfg)

cat("Smoke run with overrides: differential_abundance=false, report=false\n")
cat(sprintf("Config: %s\n", test_cfg))

# Invoke run_pipeline.R via Rscript so it sees a clean commandArgs() and the
# repo-root sniffing in run_pipeline.R picks up the right --file= path.
status <- system2(file.path(R.home("bin"), "Rscript.exe"),
                   c("run_pipeline.R", test_cfg))
cat(sprintf("\nrun_pipeline exit status: %d\n", status))

# ---- Verify .rds emission ------------------------------------------------
project_root <- "C:/Users/julio/Desktop/PC_JC_2024-11-29_Julio_Gallus"
fig_root     <- file.path(project_root, "test_run", "Figures")

probe_dirs <- c(
  file.path(fig_root, "alpha_diversity"),
  file.path(fig_root, "beta_diversity"),
  file.path(fig_root, "relative_abundance"),
  file.path(fig_root, "resistome"),
  file.path(fig_root, "virulome"),
  file.path(fig_root, "mobilome")
)
for (d in probe_dirs) {
  if (!dir.exists(d)) {
    cat(sprintf("  [missing] %s\n", d)); next
  }
  pngs <- list.files(d, pattern = "\\.png$", full.names = FALSE)
  rdss <- list.files(d, pattern = "\\.rds$", full.names = FALSE)
  cat(sprintf("  %s : %d png, %d rds\n", basename(d), length(pngs), length(rdss)))
}

# ---- Grep pipeline.log for the native/raster summary lines ---------------
log_path <- file.path(project_root, "test_run", "Reports", "pipeline.log")
if (file.exists(log_path)) {
  hits <- grep("native \\(\\.rds\\)", readLines(log_path), value = TRUE)
  cat(sprintf("\nNative/raster panel lines in %s:\n", basename(log_path)))
  for (h in tail(hits, 30)) cat("  ", h, "\n", sep = "")
}
