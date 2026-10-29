# utils_cluster_perm.R — patient-aware (cluster-respecting) permutation tests.
#
# Studies where one subject contributes several samples (repeated visits,
# paired exhale/sputum, pens/blocks) violate the independence assumption of
# Kruskal-Wallis / Wilcoxon / PERMANOVA. These helpers keep the familiar
# test statistics but build their null distribution by permuting group
# labels in a way that respects the clustering column
# (cfg$metadata$random_effect):
#
#   * "within"  design — the group varies within at least one cluster
#     (e.g. sample_type: a patient gave both EB and SP). Labels are shuffled
#     only among each cluster's own units (restricted permutation). Clusters
#     carrying a single level contribute no exchangeable information, which
#     is the correct within-subject null.
#   * "between" design — the group is constant within every cluster
#     (e.g. immunodeficiency diagnosis: one label per patient). Whole
#     clusters' labels are shuffled among clusters, so the effective sample
#     size is the number of clusters, not samples. Unequal cluster sizes are
#     fine (permute::how(plots=) would require balanced plots, which is why
#     this is hand-rolled).
#
# Opt-in: active only when cfg$stats$random_effect_tests is TRUE AND
# cfg$metadata$random_effect names a column. Otherwise every call site keeps
# its original (independence-assuming) test, so other projects' results are
# unchanged.
#
# `units`: when rows are sub-sample observations (e.g. one row per gene hit),
# pass the sample id per row — labels are permuted at the unit (sample)
# level and expanded back to rows, so hits from one sample always move
# together.
#
# p-values are (1 + #{null >= observed}) / (nperm + 1).

cluster_tests_enabled <- function(cfg) {
  re <- cfg$metadata$random_effect
  isTRUE(cfg$stats$random_effect_tests) && !is.null(re) && nzchar(re)
}

cluster_nperm <- function(cfg) {
  as.integer(cfg$stats$cluster_permutations %||%
             cfg$stats$permanova_permutations %||% 9999)
}

cluster_seed <- function(cfg) cfg$stats$cluster_seed %||% 20261005L

# Cluster id per sample (aligned to `samples`), read from the metadata file.
# Samples with a missing cluster id become their own singleton cluster.
cluster_vector <- function(cfg, samples) {
  re   <- cfg$metadata$random_effect
  sid  <- cfg$metadata$sample_id_col
  meta <- readr::read_csv(file.path(cfg$project_root, cfg$metadata$file),
                          show_col_types = FALSE, progress = FALSE)
  if (!re %in% colnames(meta)) {
    stop(sprintf("random_effect column '%s' not in metadata", re))
  }
  lk <- setNames(as.character(meta[[re]]), as.character(meta[[sid]]))
  cl <- unname(lk[as.character(samples)])
  miss <- is.na(cl)
  cl[miss] <- paste0(".singleton_", as.character(samples)[miss])
  cl
}

# Evaluate `expr` under a fixed seed without disturbing the caller's RNG.
.cp_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}

cluster_design <- function(groups, clusters) {
  varies <- tapply(as.character(groups), clusters,
                   function(g) length(unique(g)) > 1)
  if (any(varies)) "within" else "between"
}

# Integer-coded permuted labels: nperm x length(groups) matrix + metadata.
cluster_perm_labels <- function(groups, clusters, units = NULL,
                                nperm = 9999L, seed = 20261005L) {
  groups   <- as.character(groups)
  clusters <- as.character(clusters)
  n <- length(groups)
  stopifnot(length(clusters) == n)
  if (is.null(units)) units <- seq_len(n)
  units <- as.character(units)

  lev    <- sort(unique(groups))
  u      <- unique(units)
  row_u  <- match(units, u)
  first  <- match(u, units)
  u_grp  <- match(groups[first], lev)
  u_cl   <- clusters[first]
  # A unit (sample) must carry one label; a sample split across clusters is a data error.
  chk <- tapply(groups, units, function(g) length(unique(g)))
  if (any(chk > 1)) stop("cluster_perm_labels: a unit carries more than one group label")
  chk_cl <- tapply(clusters, units, function(g) length(unique(g)))
  if (any(chk_cl > 1)) stop("cluster_perm_labels: a unit belongs to more than one cluster")

  design <- cluster_design(u_grp, u_cl)
  perms  <- matrix(0L, nrow = nperm, ncol = length(u))
  .cp_with_seed(seed, {
    if (design == "within") {
      idx_by_cl <- split(seq_along(u), u_cl)
      idx_by_cl <- idx_by_cl[lengths(idx_by_cl) > 1]
      for (b in seq_len(nperm)) {
        lab <- u_grp
        for (ix in idx_by_cl) lab[ix] <- lab[ix][sample.int(length(ix))]
        perms[b, ] <- lab
      }
    } else {
      cl_ids <- unique(u_cl)
      cl_grp <- u_grp[match(cl_ids, u_cl)]
      pos    <- match(u_cl, cl_ids)
      for (b in seq_len(nperm)) {
        perms[b, ] <- cl_grp[sample.int(length(cl_ids))][pos]
      }
    }
  })
  list(observed   = match(groups, lev),
       perms      = perms[, row_u, drop = FALSE],
       levels     = lev,
       design     = design,
       n_clusters = length(unique(u_cl)),
       n_units    = length(u),
       nperm      = nperm)
}

.cp_pvalue <- function(obs, null) {
  (1 + sum(null >= obs - 1e-10 * max(1, abs(obs)), na.rm = TRUE)) /
    (sum(!is.na(null)) + 1)
}

# --- Statistics on integer-coded labels ------------------------------------
.cp_kw_H <- function(r, g, N) {
  n  <- tabulate(g)
  Rs <- vapply(seq_along(n), function(j) sum(r[g == j]), numeric(1))
  ok <- n > 0
  12 / (N * (N + 1)) * sum(Rs[ok]^2 / n[ok]) - 3 * (N + 1)
}

.cp_anova_F <- function(z, g) {
  n  <- tabulate(g); ok <- n > 0; k <- sum(ok); N <- length(z)
  if (k < 2 || N <= k) return(NA_real_)
  m   <- vapply(seq_along(n), function(j) if (n[j]) mean(z[g == j]) else 0, numeric(1))
  ssb <- sum(n[ok] * (m[ok] - mean(z))^2)
  ssw <- sum((z - m[g])^2)
  (ssb / (k - 1)) / (ssw / (N - k))
}

# Pseudo-F (McArdle & Anderson 2001) for a one-factor design; D2 = squared
# distance matrix. Identical to vegan::adonis2's F and R2 for `d ~ group`.
.cp_pseudo_F <- function(D2, g) {
  N  <- nrow(D2); n <- tabulate(g); ok <- n > 0; k <- sum(ok)
  sst <- sum(D2) / (2 * N)
  ssw <- sum(vapply(which(ok), function(j) {
    ix <- which(g == j); sum(D2[ix, ix]) / (2 * n[j])
  }, numeric(1)))
  list(F  = ((sst - ssw) / (k - 1)) / (ssw / (N - k)),
       R2 = (sst - ssw) / sst)
}

.cp_result <- function(stat, p, pl, method, extra = list()) {
  c(list(statistic = stat, p.value = p, method = method,
         design = pl$design, n_clusters = pl$n_clusters,
         n_units = pl$n_units, nperm = pl$nperm), extra)
}

# Short provenance tag for logs / plot annotations.
cluster_tag <- function(res) {
  if (is.null(res$design)) return("")
  sprintf("cluster perm: %s, %d clusters",
          if (res$design == "within") "labels shuffled within clusters"
          else "whole-cluster labels shuffled",
          res$n_clusters)
}

# Plot/log label for a KW result: flags the patient-aware variant so a
# figure never shows a cluster-permutation p under a plain "KW p" label.
kw_name <- function(kw, long = FALSE) {
  base <- if (long) "Kruskal-Wallis p" else "KW p"
  if (!is.null(kw$design)) paste(base, "(patient perm)") else base
}

# Kruskal-Wallis dispatcher used by every per-sample KW call site: the
# cluster-permutation test when cfg enables it, else stats::kruskal.test.
# `samples` gives the sample id per row (rows may be sub-sample
# observations; they are then permuted together at the sample level).
kw_test <- function(values, groups, samples, cfg) {
  if (cluster_tests_enabled(cfg)) {
    cluster_kw(values, groups, cluster_vector(cfg, samples), units = samples,
               nperm = cluster_nperm(cfg), seed = cluster_seed(cfg))
  } else {
    stats::kruskal.test(values, factor(groups))
  }
}

# --- Public tests ------------------------------------------------------------
# Kruskal-Wallis H with a cluster-respecting permutation null. Drop-in for
# kruskal.test(): returns $statistic and $p.value.
cluster_kw <- function(values, groups, clusters, units = NULL,
                       nperm = 9999L, seed = 20261005L) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- groups[ok]; clusters <- clusters[ok]
  if (!is.null(units)) units <- units[ok]
  if (length(unique(groups)) < 2) return(NULL)
  pl <- cluster_perm_labels(groups, clusters, units, nperm, seed)
  r  <- rank(values); N <- length(values)
  obs  <- .cp_kw_H(r, pl$observed, N)
  null <- apply(pl$perms, 1, function(g) .cp_kw_H(r, g, N))
  .cp_result(obs, .cp_pvalue(obs, null), pl,
             "Kruskal-Wallis H, cluster permutation")
}

# Pairwise two-sided Wilcoxon rank-sum tests (two-group KW H) with cluster permutation,
# one row per pair; p adjusted across pairs with `p.adjust.method`.
# Columns mirror ggpubr::compare_means (group1, group2, p, p.adj, p.signif).
cluster_wilcox_pairs <- function(values, groups, clusters, units = NULL,
                                 nperm = 9999L, seed = 20261005L,
                                 p.adjust.method = "BH") {
  groups <- as.character(groups)
  lev    <- sort(unique(groups[!is.na(groups)]))
  if (length(lev) < 2) return(NULL)
  prs <- utils::combn(lev, 2, simplify = FALSE)
  out <- lapply(prs, function(pr) {
    keep <- groups %in% pr & !is.na(values)
    x  <- values[keep]; g <- groups[keep]; cl <- clusters[keep]
    un <- if (is.null(units)) NULL else units[keep]
    if (min(table(factor(g, levels = pr))) < 1) return(NULL)
    pl <- cluster_perm_labels(g, cl, un, nperm, seed)
    r  <- rank(x)
    # Two-group Kruskal-Wallis H = the standardised rank-sum statistic. The
    # raw |U - n1*n2/2| was used before, but cluster/unit permutations change
    # the group sizes, so its scale changed from permutation to permutation
    # and the pairwise p could disagree with the overall KW for two groups.
    wstat <- function(code) {
      n1 <- sum(code == 1L)
      if (n1 == 0 || n1 == length(code)) return(NA_real_)
      .cp_kw_H(r, code, length(code))
    }
    obs  <- wstat(pl$observed)
    null <- apply(pl$perms, 1, wstat)
    data.frame(group1 = pr[1], group2 = pr[2],
               p = .cp_pvalue(obs, null), design = pl$design,
               n_clusters = pl$n_clusters, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$p.adj    <- stats::p.adjust(out$p, method = p.adjust.method)
  out$p.signif <- as.character(cut(out$p.adj,
                                   c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                                   c("****", "***", "**", "*", "ns")))
  out
}

# One-factor PERMANOVA (pseudo-F, R2) with a cluster-respecting null.
cluster_permanova <- function(d, groups, clusters, nperm = 9999L,
                              seed = 20261005L) {
  D2 <- as.matrix(d)^2
  pl <- cluster_perm_labels(groups, clusters, NULL, nperm, seed)
  obs  <- .cp_pseudo_F(D2, pl$observed)
  null <- apply(pl$perms, 1, function(g) .cp_pseudo_F(D2, g)$F)
  .cp_result(obs$F, .cp_pvalue(obs$F, null), pl,
             "PERMANOVA pseudo-F, cluster permutation",
             list(R2 = obs$R2))
}

# PERMDISP: distances to group spatial medians from vegan::betadisper (true
# groups, as vegan does), then the one-way ANOVA F on those distances is
# tested against a cluster-respecting permutation of the labels.
cluster_permdisp <- function(d, groups, clusters, nperm = 9999L,
                             seed = 20261005L) {
  bd <- vegan::betadisper(d, factor(groups))
  z  <- as.numeric(bd$distances)
  pl <- cluster_perm_labels(groups, clusters, NULL, nperm, seed)
  obs  <- .cp_anova_F(z, pl$observed)
  null <- apply(pl$perms, 1, function(g) .cp_anova_F(z, g))
  .cp_result(obs, .cp_pvalue(obs, null), pl,
             "PERMDISP ANOVA F on betadisper distances, cluster permutation",
             list(betadisper = bd))
}

# Plain-text record of a cluster-aware beta test (replaces the adonis2 /
# permutest printout in the *.txt files the manifest links to).
cluster_write_beta_txt <- function(res, file, group, cluster_name, label) {
  lines <- c(
    sprintf("%s — %s", label, res$method),
    sprintf("Formula: distance ~ %s; clusters: %s (%s design, %d clusters, %d samples)",
            group, cluster_name, res$design, res$n_clusters, res$n_units),
    sprintf("Permutations: %d", res$nperm),
    if (!is.null(res$R2)) sprintf("R2 = %.4f", res$R2),
    sprintf("F = %.4f", res$statistic),
    sprintf("p (permutation) = %.4g", res$p.value)
  )
  writeLines(lines, file)
}

# ALDEx2 pooling: sum counts within patient (between design) or within
# patient x group (within design) so each pooled unit is one replicate.
# Returns the pooled matrix, a metadata frame keyed by unit id in `sid`
# (with `group` and the cluster column), and the sample -> unit map.
cluster_pool_counts <- function(mat, meta_use, sid, group, cluster_name) {
  samples <- colnames(mat)
  m   <- meta_use[match(samples, as.character(meta_use[[sid]])), , drop = FALSE]
  grp <- as.character(m[[group]])
  cl  <- as.character(m[[cluster_name]])
  cl[is.na(cl)] <- paste0(".singleton_", samples[is.na(cl)])
  design <- cluster_design(grp, cl)
  unit <- if (design == "within") paste(cl, grp, sep = "|") else cl
  units <- unique(unit)
  pooled <- vapply(units, function(u) rowSums(mat[, unit == u, drop = FALSE]),
                   numeric(nrow(mat)))
  if (is.null(dim(pooled))) pooled <- matrix(pooled, nrow = nrow(mat))
  dimnames(pooled) <- list(rownames(mat), units)
  pm <- data.frame(units, grp[match(units, unit)], cl[match(units, unit)],
                   stringsAsFactors = FALSE)
  names(pm) <- c(sid, group, cluster_name)
  list(mat = pooled, meta = pm,
       unit_of = setNames(unit, samples), design = design)
}
