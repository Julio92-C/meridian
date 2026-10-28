#!/usr/bin/env Rscript
# validate_cluster_perm.R — unit + simulation checks for R/utils_cluster_perm.R.
# Run from the repo root:  Rscript scripts/validate_cluster_perm.R
# Exits non-zero if any check fails.

suppressPackageStartupMessages(library(vegan))
`%||%` <- function(a, b) if (is.null(a)) b else a
source("R/utils_cluster_perm.R")

fails <- 0L
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
  if (!isTRUE(ok)) fails <<- fails + 1L
}

set.seed(1)

# --- 1. Pseudo-F / R2 identical to adonis2 ---------------------------------
X  <- matrix(rpois(40 * 25, 5), 40)
g  <- rep(c("A", "B", "C", "D"), each = 10)
d  <- vegdist(decostand(X, "hellinger"), "bray")
ad <- adonis2(d ~ g, data = data.frame(g = g), permutations = 0)
cp <- cluster_permanova(d, g, clusters = as.character(seq_along(g)), nperm = 99)
check("pseudo-F matches adonis2", abs(cp$statistic - ad$F[1]) < 1e-8,
      sprintf("(%.6f vs %.6f)", cp$statistic, ad$F[1]))
check("R2 matches adonis2", abs(cp$R2 - ad$R2[1]) < 1e-8,
      sprintf("(%.6f vs %.6f)", cp$R2, ad$R2[1]))

# --- 2. No clustering (singletons) ~ classical tests ------------------------
x  <- c(rnorm(15), rnorm(15, 0.8))
g2 <- rep(c("A", "B"), each = 15)
kw_cl <- cluster_kw(x, g2, clusters = as.character(seq_along(x)), nperm = 9999)
kw_as <- kruskal.test(x, factor(g2))
check("singleton KW H matches kruskal.test (no ties)",
      abs(kw_cl$statistic - unname(kw_as$statistic)) < 1e-8)
check("singleton KW perm p ~ asymptotic p",
      abs(kw_cl$p.value - kw_as$p.value) < 0.01,
      sprintf("(%.4f vs %.4f)", kw_cl$p.value, kw_as$p.value))
check("singleton design is 'between'", kw_cl$design == "between")

ad_p <- adonis2(d ~ g, data = data.frame(g = g), permutations = 4999)$`Pr(>F)`[1]
cp_p <- cluster_permanova(d, g, as.character(seq_along(g)), nperm = 4999)$p.value
check("singleton PERMANOVA p ~ adonis2 p", abs(ad_p - cp_p) < 0.03,
      sprintf("(%.4f vs %.4f)", cp_p, ad_p))

# --- 3. Permutation structure ----------------------------------------------
cl  <- rep(paste0("P", 1:6), times = c(1, 2, 3, 4, 5, 6))       # unbalanced
gw  <- unlist(lapply(c(1, 2, 3, 4, 5, 6), function(n) rep(c("EB", "SP"), length.out = n)))
plw <- cluster_perm_labels(gw, cl, nperm = 200)
check("within design detected", plw$design == "within")
same_counts <- all(apply(plw$perms, 1, function(p)
  all(tapply(p, cl, function(v) paste(sort(v), collapse = "")) ==
      tapply(plw$observed, cl, function(v) paste(sort(v), collapse = "")))))
check("within: per-cluster label multiset preserved", same_counts)

gb  <- c("CVID", "SID", "CVID", "XLA", "SID", "HG")[match(cl, paste0("P", 1:6))]
plb <- cluster_perm_labels(gb, cl, nperm = 200)
check("between design detected", plb$design == "between")
const <- all(apply(plb$perms, 1, function(p) all(tapply(p, cl, function(v) length(unique(v))) == 1)))
check("between: label constant within each cluster in every permutation", const)
check("between: unequal cluster sizes accepted", plb$n_clusters == 6)

# units: rows nested in samples
units <- rep(paste0("S", 1:6), each = 3)
clu   <- rep(c("P1", "P1", "P2", "P2", "P3", "P3"), each = 3)
gu    <- rep(c("EB", "SP", "EB", "SP", "EB", "SP"), each = 3)
plu   <- cluster_perm_labels(gu, clu, units = units, nperm = 100)
check("units: rows of one sample always share a label",
      all(apply(plu$perms, 1, function(p) all(tapply(p, units, function(v) length(unique(v))) == 1))))

# --- 4. RNG of the caller is untouched -------------------------------------
set.seed(99); a <- runif(1); set.seed(99)
invisible(cluster_perm_labels(gb, cl, nperm = 10)); b <- runif(1)
check("caller RNG stream unchanged", identical(a, b))

# --- 5. Type-I error under clustering (between-subject, no true effect) ----
sim_between <- function(nsim = 300, nperm = 199) {
  rej_naive <- rej_cl <- 0
  for (s in seq_len(nsim)) {
    np   <- 20
    size <- sample(1:8, np, replace = TRUE)
    pid  <- rep(paste0("P", seq_len(np)), size)
    peff <- rnorm(np, 0, 2)[match(pid, paste0("P", seq_len(np)))]   # strong patient effect
    y    <- peff + rnorm(length(pid), 0, 0.5)
    grpP <- sample(rep(c("CVID", "SID", "XLA", "HG"), length.out = np))
    grp  <- grpP[match(pid, paste0("P", seq_len(np)))]
    rej_naive <- rej_naive + (kruskal.test(y, factor(grp))$p.value < 0.05)
    rej_cl    <- rej_cl + (cluster_kw(y, grp, pid, nperm = nperm, seed = s)$p.value < 0.05)
  }
  c(naive = rej_naive / nsim, cluster = rej_cl / nsim)
}
set.seed(2024)
te <- sim_between()
cat(sprintf("      type-I (between, 300 sims): naive KW = %.3f, cluster KW = %.3f\n",
            te["naive"], te["cluster"]))
check("between: naive KW is anti-conservative under clustering", te["naive"] > 0.15)
check("between: cluster KW keeps nominal 5% (<= 0.085)", te["cluster"] <= 0.085)

sim_within <- function(nsim = 300, nperm = 199) {
  rej_cl <- 0
  for (s in seq_len(nsim)) {
    np   <- 20
    size <- sample(1:8, np, replace = TRUE)
    pid  <- rep(paste0("P", seq_len(np)), size)
    peff <- rnorm(np, 0, 2)[match(pid, paste0("P", seq_len(np)))]
    y    <- peff + rnorm(length(pid), 0, 0.5)
    grp  <- ifelse(runif(length(pid)) < 0.3, "SP", "EB")         # varies within patient
    res  <- cluster_kw(y, grp, pid, nperm = nperm, seed = s)
    if (!is.null(res)) rej_cl <- rej_cl + (res$p.value < 0.05)
  }
  rej_cl / nsim
}
set.seed(7)
tw <- sim_within()
cat(sprintf("      type-I (within, 300 sims): cluster KW = %.3f\n", tw))
check("within: cluster KW keeps nominal 5% (<= 0.085)", tw <= 0.085)

# Power sanity: a real within-patient shift is still detected.
set.seed(11)
pid <- rep(paste0("P", 1:20), each = 4)
y   <- rnorm(20, 0, 2)[match(pid, paste0("P", 1:20))] +
       rep(c(0, 0, 1.5, 1.5), 20) + rnorm(80, 0, 0.5)
grp <- rep(c("EB", "EB", "SP", "SP"), 20)
check("within: real paired shift detected (p < 0.01)",
      cluster_kw(y, grp, pid, nperm = 999)$p.value < 0.01)

# --- 6. PERMDISP + pooling ---------------------------------------------------
pd <- cluster_permdisp(d, g, as.character(seq_along(g)), nperm = 999)
pt <- permutest(betadisper(d, factor(g)), permutations = 999)
check("singleton PERMDISP F matches permutest F",
      abs(pd$statistic - pt$tab$F[1]) < 1e-8,
      sprintf("(%.5f vs %.5f)", pd$statistic, pt$tab$F[1]))

mat  <- matrix(1:12, 2, 6, dimnames = list(c("t1", "t2"), paste0("S", 1:6)))
meta <- data.frame(sample = paste0("S", 1:6), grp = c("EB", "SP", "EB", "EB", "SP", "SP"),
                   ID = c("P1", "P1", "P1", "P2", "P2", "P3"))
pool <- cluster_pool_counts(mat, meta, "sample", "grp", "ID")
check("pool: within design -> patient x group units", pool$design == "within" &&
      setequal(colnames(pool$mat), c("P1|EB", "P1|SP", "P2|EB", "P2|SP", "P3|SP")))
check("pool: counts summed", pool$mat["t1", "P1|EB"] == mat["t1", "S1"] + mat["t1", "S3"])
check("pool: total counts preserved", sum(pool$mat) == sum(mat))

meta$dx <- c("CVID", "CVID", "CVID", "SID", "SID", "XLA")
pool_b <- cluster_pool_counts(mat, meta, "sample", "dx", "ID")
check("pool: between design -> one unit per patient",
      pool_b$design == "between" && setequal(colnames(pool_b$mat), c("P1", "P2", "P3")))

cat(sprintf("\n%d check(s) failed\n", fails))
quit(status = if (fails > 0) 1L else 0L)
