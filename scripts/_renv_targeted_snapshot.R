options(renv.config.startup.quiet = TRUE)
source("renv/activate.R")

lf_before <- renv::lockfile_read()
cat("Before: lockfile has", length(lf_before$Packages), "packages\n")

# renv::record() is additive — appends entries for the named packages to the
# existing lockfile rather than recomputing a minimal one (which is what
# renv::snapshot(packages=) does and which would destroy unrelated entries).
renv::record(c("magick", "openxlsx", "zip"))

lf_after <- renv::lockfile_read()
cat("After:  lockfile has", length(lf_after$Packages), "packages\n")

new_pkgs <- setdiff(names(lf_after$Packages), names(lf_before$Packages))
cat("Newly recorded packages:", paste(new_pkgs, collapse = ", "), "\n")

for (p in c("magick", "openxlsx", "zip")) {
  entry <- lf_after$Packages[[p]]
  if (is.null(entry)) {
    cat(sprintf("  %s: MISSING\n", p))
  } else {
    cat(sprintf("  %s: version=%s source=%s\n", p,
                entry$Version %||% "?", entry$Source %||% "?"))
  }
}

cat("\n--- status after ---\n")
renv::status()
