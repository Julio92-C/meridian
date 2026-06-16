# Lightweight parse check: walks R/ and reports any file that fails to
# parse. Used as a smoke test after sweeps that touch many files.
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
results <- vapply(files, function(f) {
  tryCatch({
    parse(file = f)
    TRUE
  }, error = function(e) {
    cat(sprintf("PARSE FAIL: %s :: %s\n", f, conditionMessage(e)))
    FALSE
  })
}, logical(1))
cat(sprintf("parse_ok=%s  files=%d  ok=%d  fail=%d\n",
            all(results), length(results), sum(results), sum(!results)))
if (!all(results)) quit(status = 1)
