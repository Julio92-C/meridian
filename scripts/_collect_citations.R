options(renv.config.startup.quiet = TRUE)
source("renv/activate.R")

pkgs <- c("vegan", "ggraph", "igraph", "networkD3", "quarto")
for (p in pkgs) {
  cat("===", p, "===\n")
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  cat("version:", v, "\n")
  cit <- tryCatch(
    paste(format(citation(p), style = "text"), collapse = "\n\n"),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  cat(cit, "\n\n")
}

cat("=== R ===\n")
cat(paste(format(citation(), style = "text"), collapse = "\n\n"), "\n")
cat("R version:", R.version.string, "\n")
