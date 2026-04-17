# 10_virulome.R — VF profile analysis. Uses VFDB rows from ABRicate output.

run_virulome <- function(cleaned, cfg) {
  pipeline_log(cfg, "Virulome (VFDB)")
  .run_ge_block(cleaned, cfg, db = "vfdb", label = "VF",
                out_subdir = "virulome")
}
