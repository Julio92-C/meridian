# Functional smoke test for save_panel_ggplot() and the R/14 .rds-prefer
# path. Builds a trivial ggplot, saves via the helper, verifies the .rds
# sibling exists and round-trips back to an identical-looking plot.
source("R/utils_panel_io.R")
suppressPackageStartupMessages(library(ggplot2))

tmp <- tempfile(fileext = ".png")
p <- ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) +
  geom_point(size = 2.5) +
  labs(title = "Smoke test") +
  theme_classic()

save_panel_ggplot(tmp, p, width = 6, height = 4, dpi = 300)

rds <- sub("\\.png$", ".rds", tmp)
stopifnot(file.exists(tmp))
stopifnot(file.exists(rds))

p2 <- readRDS(rds)
stopifnot(inherits(p2, "ggplot"))

cat(sprintf("PNG bytes: %d\n", file.size(tmp)))
cat(sprintf("RDS bytes: %d\n", file.size(rds)))
cat("OK: save_panel_ggplot wrote PNG + sibling RDS; RDS round-trips to ggplot.\n")
