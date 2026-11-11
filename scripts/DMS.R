#!/usr/bin/env Rscript
# =============================================================================
# Differential RNA modification analysis (dRNA-seq, nanopore)
#
# Design : 3 basal + 3 luminal cell lines, 4 libraries per line (23 libraries)
# Unit   : the CELL LINE is the experimental unit -> effective n = 3 per group
# Model  : cbind(mod, unmod) ~ subtype + (1 | cell_line), family = betabinomial
#
# Variance levels explicitly modelled:
#   reads within library  -> binomial
#   library within line   -> betabinomial dispersion
#   line within subtype   -> random intercept  <-- the subtype effect is tested here
# =============================================================================
setwd("/projects/CGS_shared/vfama/BRIGHT_PROJECT/")
suppressPackageStartupMessages({
  library("glmmTMB")
  library("future.apply")
  library("data.table")
  library("parallel")
})
# =============================================================================
# 0.  SIMULATED INPUT  (delete this block once you plug in real data)
# =============================================================================
# Produces exactly the three objects the pipeline expects:
#   meta        : 23 rows, columns sample_id / cell_line / subtype
#   counts_mod  : sites x 23 integer matrix
#   counts_unmod: sites x 23 integer matrix
luminal <- c("BT483", "MCF7", "T47D")
basal   <- c("SUM159", "MDAMB231", "BT549")
bedMethylFiles <- list.files(path="/projects/CGS_shared/vfama/BRIGHT_PROJECT/",pattern="filtCov.m6A.bed",recursive=TRUE)
bedMethylFiles <- bedMethylFiles[grep("DMSO",bedMethylFiles)]
metadata <- data.frame("sample_id" = sapply(strsplit(bedMethylFiles, "/"), `[`, 2),
                       "cell_line" = factor(sapply(strsplit(bedMethylFiles, "/"), `[`, 1)))
metadata$subtype <- factor(ifelse(metadata$cell_line %in% luminal, "luminal",
                    ifelse(metadata$cell_line %in% basal,   "basal", NA)))

read_one <- function(f) {
  d <- fread(file.path("/projects/CGS_shared/vfama/BRIGHT_PROJECT/", f), header = TRUE, sep = "\t",
             select = c(1, 2, 10, 11, 12),
             col.names = c("chrom", "start", "Nvalid_cov","percent_modified","mod"))
  d[, site := paste(chrom, start, sep = ":")]
  d[, .(site,percent_modified,Nvalid_cov, mod)]
} #Read the bedMethyl files to get the columns we need

lst <- lapply(bedMethylFiles, read_one)
names(lst) <- metadata$sample_id
sites <- Reduce(union, lapply(lst, `[[`, "site")) #collect all the unique sites

get_mat <- function(col) {
  m <- sapply(lst, function(d) d[[col]][match(sites, d$site)])
  rownames(m) <- sites
  m[is.na(m)] <- 0          #site absent in that rep = 0 read
  m
} #Extract mod and unmod counts for each site

counts_mod   <- get_mat("mod")
counts_unmod <- get_mat("Nvalid_cov") - get_mat("mod")
stoich <- get_mat("percent_modified")

stopifnot(identical(colnames(counts_mod), metadata$sample_id),
          identical(dim(counts_mod), dim(counts_unmod)))

message("sites in: ", nrow(counts_mod), " | libraries: ", nrow(metadata))

# -----------------------------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------------------------
n_libs <- table(metadata$cell_line)# lib per line
PAR <- list(
  min_reads_library   = 20,        # drop a library at this site if below (unbalanced is fine)
  min_libs_per_line   = n_libs,    # require this many usable libraries (replicate) per line
  min_cell_per_line   = 20,        # pooled over its libraries, each line needs
                                   #>= this many modified AND unmodified reads
  min_native_fraction = 0.2,       #threshold on stoichiometry
  loo_sign_stable     = TRUE,      #leave-one-CELL-LINE-out must preserve sign
  loo_min_fraction    = 0.5,       # ... and >= this fraction of |beta|
  fdr                 = 0.05,
  n_workers           = 4
)

# =============================================================================
# 1.  LIBRARY-LEVEL QC   (before any per-site work)
# =============================================================================
# One bad library out of 24 can distort a whole cell line. Check this first:
# the per-library global modification rate should cluster by cell line, not
# scatter at random. An outlier library here is a basecalling or prep problem.

lib_qc <- data.table(
  sample_id = metadata$sample_id,
  cell_line = metadata$cell_line,
  subtype   = metadata$subtype,
  depth     = colSums(counts_mod + counts_unmod),
  global_mod_rate = colSums(counts_mod) /
                    colSums(counts_mod + counts_unmod)
)
print(lib_qc[, .(n_libs = .N,
                 depth_median = median(depth),
                 mod_rate_mean = round(mean(global_mod_rate), 4),
                 mod_rate_sd   = round(sd(global_mod_rate), 4)),
             by = .(cell_line, subtype)])

# Flag libraries >3 SD from their own cell line's mean
lib_qc[, z_within_line := (global_mod_rate - mean(global_mod_rate)) /
                           max(sd(global_mod_rate), 1e-9), by = cell_line]
bad_libs <- lib_qc[abs(z_within_line) > 3, sample_id]
if (length(bad_libs)) message("outlier libraries: ", paste(bad_libs, collapse = ", "))

# =============================================================================
# 2.  VECTORISED SITE FILTERS   (boolean, ANDed, all BEFORE any model fit)
# =============================================================================
# Keep sites that have coverage everywhere

tot <- counts_mod + counts_unmod
line_levels <- levels(metadata$cell_line)
line_of_col <- as.integer(metadata$cell_line)

# --- 2a. coverage: replicati usabili per linea ---
keep_libs <- rowSums(tot >= PAR$min_reads_library) == ncol(tot) #5.7M

# --- 2b. minimum cell counts, PER CELL LINE ----------------------------------
# From the Woolf formula SE(b1) = sqrt(sum 1/cell), precision is set by the
# SMALLEST cell. Enforcing this per line guarantees every biological replicate
# contributes a usable estimate -- strictly better than a stoichiometry cutoff.
# mod_line <- sapply(seq_along(line_levels), function(k)
#   rowSums(counts_mod[, line_of_col == k, drop = FALSE]))
# unm_line <- sapply(seq_along(line_levels), function(k)
#   rowSums(counts_unmod[, line_of_col == k, drop = FALSE]))
# colnames(mod_line) <- colnames(unm_line) <- line_levels

# keep_cells <- rowSums(mod_line >= PAR$min_cell_per_line &
#                       unm_line >= PAR$min_cell_per_line) == length(line_levels)

# --- 2b. stechiometria >= 20 in almeno un sottotipo ---
# luminal_cols <- c("BT483", "MCF7", "T47D")
# basal_cols   <- c("SUM159", "MDAMB231", "BT549")
# n_rep    <- as.vector(table(metadata$cell_line))[match(line_levels, names(table(metadata$cell_line)))]
# min_libs <- pmax(n_rep - 1, 2)
# idx <- which(keep_libs)
# high <- stoich[idx, , drop = FALSE] >= 20

# high_per_line <- t(rowsum(t(high) * 1, line_of_col))
# colnames(high_per_line) <- line_levels
# ok_high <- sweep(high_per_line, 2, min_libs, FUN = ">=")

# mod_lum <- rowSums(ok_high[, luminal_cols, drop = FALSE]) >= 2
# mod_bas <- rowSums(ok_high[, basal_cols,   drop = FALSE]) >= 2

# keep <- rep(FALSE, nrow(tot))
# keep[idx] <- mod_lum | mod_bas  
# attrition <- data.frame(
#   step = c("input", "coverage in tutte le librerie", "+ stoich >= 20 (>= 1 sottotipo)"),
#   n = c(nrow(tot), sum(keep_libs), sum(keep))
# )
# sites <- which(keep)
luminal_cols <- c("BT483", "MCF7", "T47D")
basal_cols   <- c("SUM159", "MDAMB231", "BT549")
n_rep    <- as.vector(table(metadata$cell_line))[match(line_levels, names(table(metadata$cell_line)))]
min_libs <- pmax(n_rep - 1, 2)

# coverage: invariata
idx <- which(keep_libs)

# stechiometria per linea, solo sui siti con coverage ok
high <- stoich[idx, , drop = FALSE] >= 20
high[is.na(high)] <- FALSE

high_per_line <- t(rowsum(t(high) * 1, line_of_col))
colnames(high_per_line) <- line_levels

# la linea "ha il sito" se stoich >= 20 in almeno n-1 replicati
has_site <- sweep(high_per_line, 2, min_libs, FUN = ">=")

# massimo una linea senza il sito
n_missing <- length(line_levels) - rowSums(has_site)

keep <- rep(FALSE, nrow(tot))
keep[idx] <- n_missing <= 1

attrition <- data.frame(
  step = c("input", "coverage in tutte le librerie", "+ stoich >= 20 in >= 5 linee su 6"),
  n = c(nrow(tot), sum(keep_libs), sum(keep))
)
print(attrition)
sites <- which(keep)
# =============================================================================
# 3.  THE MODEL
# =============================================================================
# ~ subtype + (1 | cell_line), betabinomial
#
#   subtype          fixed, the effect of interest
#   (1 | cell_line)  6 groups x 4 observations -> a REAL blocking factor.
#                    This is what makes the subtype effect be tested against
#                    between-line variance instead of between-library variance.
#                    Without it you are claiming n = 24 when you have n = 3+3.
#   betabinomial     residual overdispersion (library-to-library within a line)
#
# glmmTMB rather than glmer: lme4 has no betabinomial family, so it cannot fit
# both variance levels at once.

fit_site <- function(i, mod_mat, unm_mat, md, par) {
  d <- data.frame(
    mod       = as.numeric(mod_mat[i, ]),
    unmod     = as.numeric(unm_mat[i, ]),
    subtype   = md$subtype,
    cell_line = md$cell_line
  )
  d <- d[d$mod + d$unmod >= par$min_reads_library, , drop = FALSE]
  d$cell_line <- droplevels(d$cell_line)

  if (length(unique(d$subtype)) < 2L || nlevels(d$cell_line) < 4L) return(NULL)

  f <- try(glmmTMB(cbind(mod, unmod) ~ subtype + (1 | cell_line),
                   family = betabinomial, data = d), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)

  converged <- isTRUE(f$fit$convergence == 0) && isTRUE(f$sdr$pdHess)
  cf <- summary(f)$coefficients$cond
  if (nrow(cf) < 2L || !is.finite(cf[2, "Std. Error"])) return(NULL)

  vc <- glmmTMB::VarCorr(f)$cond
  sd_line <- if (!is.null(vc$cell_line)) sqrt(as.numeric(vc$cell_line[1, 1])) else NA_real_

  pr <- residuals(f, type = "pearson")
  list(
    site        = rownames(mod_mat)[i],
    intercept   = cf[1, "Estimate"],
    estimate    = cf[2, "Estimate"],
    std_error   = cf[2, "Std. Error"],
    z_value     = cf[2, "z value"],
    p_wald      = cf[2, "Pr(>|z|)"],
    sd_cellline = sd_line,
    bb_sigma    = tryCatch(sigma(f), error = function(e) NA_real_),
    dispersion  = sum(pr^2, na.rm = TRUE) / df.residual(f),
    n_obs       = nrow(d),
    n_lines     = nlevels(d$cell_line),
    converged   = converged
  )
}
# =============================================================================
# 4.  RUN  (parallel; chunk this across array jobs on the cluster)
# =============================================================================
res_list <- mclapply(sites, fit_site,
                     mod_mat = counts_mod, unm_mat = counts_unmod,
                     md = metadata, par = PAR,
                     mc.cores = PAR$n_workers)

res <- rbindlist(Filter(Negate(is.null), res_list))
message("fitted: ", nrow(res), " / ", length(sites),
        "  (converged: ", sum(res$converged), ")")
res <- res[converged == TRUE]
# =============================================================================
# 5.  DIAGNOSTICS   -- read these before trusting any p-value
# =============================================================================
# (a) Is the random effect actually being estimated?
#     sd_cellline == 0 means a singular fit: the model collapsed and the
#     subtype effect is being tested against between-LIBRARY variance, i.e.
#     pseudo-replication is back. If this is a large fraction of sites, switch
#     to the collapsed model in section 8.
frac_singular <- mean(res$sd_cellline < 1e-4, na.rm = TRUE)
message(sprintf("singular random effect: %.1f%% of sites", 100 * frac_singular)) #14% dei siti ha una sd quasi nulla
if (frac_singular > 0.40)
  warning("random effect rarely identifiable -- prefer the collapsed model (section 8)")

# (b) Did betabinomial absorb the residual overdispersion?
#     Pearson dispersion follows chi2(df)/df under correctness, whose MEDIAN is
#     NOT 1. Normalise by the expected median before reading it.
res[, df_res := n_obs - 2L]
res[, disp_norm := dispersion / (qchisq(0.5, df_res) / df_res)]
message(sprintf("median normalised dispersion: %.2f  (target ~1)",
                median(res$disp_norm, na.rm = TRUE)))              # 1.35 -> piccolo eccesso di variabilità residua

# (c) How big is the between-line variance? This is what limits your power.
message(sprintf("median SD(cell_line) on logit scale: %.2f",
                median(res$sd_cellline, na.rm = TRUE)))            #0.23 -> le linee sono abbastanza simili (ci piace?)

# =============================================================================
# 6.  GLOBAL SHIFT: estimate it on beta (logit), remove it, recalibrate
# =============================================================================
# Assumption: most sites do NOT change between subtypes, so the median of beta
# is the global offset (technical or biological) and not a site-level effect.
df_effective <- length(line_levels) - 2L          # 4

b0 <- median(res$estimate, na.rm = TRUE)
res[, est_c := estimate - b0]
res[, z_c   := est_c / std_error]

sigma0 <- mad(res$z_c, na.rm = TRUE)              # width of the null body
s0     <- max(sigma0, 1)                          # never shrink the null
res[, z_f := z_c / s0]

# raw (uncorrected) p-values, kept for comparison
res[, p_z  := p_wald]
res[, q_z  := p.adjust(p_z, "BH")]

# corrected p-values: normal reference and t(df) reference
res[, p_zc := 2 * pnorm(-abs(z_f))]
res[, p_tc := 2 * pt(-abs(z_f), df = df_effective)]
res[, q_zc := p.adjust(p_zc, "BH")]
res[, q_tc := p.adjust(p_tc, "BH")]

med_int <- median(res$intercept, na.rm = TRUE)
b0_pp   <- 100 * (plogis(med_int + b0) - plogis(med_int))
message(sprintf("shift b0 = %.3f logit (~%.1f pp at the median site) | MAD(z_c) = %.2f | scale used = %.2f",
                b0, b0_pp, sigma0, s0))
message(sprintf("hits FDR %.0f%%:  raw z %d | corrected z %d | corrected t(%d) %d",
                100 * PAR$fdr, sum(res$q_z < PAR$fdr),
                sum(res$q_zc < PAR$fdr), df_effective, sum(res$q_tc < PAR$fdr)))

# --- checks on the correction (look at these) --------------------------------
par(mfrow = c(1, 3))
hist(res$p_zc, breaks = 50, main = "corrected p-values", xlab = "p")   # flat + peak near 0
qqnorm(res$z_f, pch = ".", main = "QQ corrected z"); abline(0, 1, col = 2)
hist(res$z_f[abs(res$z_f) < 6], breaks = 100, main = "corrected z", xlab = "z")
par(mfrow = c(1, 1))

# old vs new hit lists and their sign balance
old <- res[q_z  < PAR$fdr, site]
new <- res[q_zc < PAR$fdr, site]
print(c(old = length(old), new = length(new), shared = length(intersect(old, new))))
print(table(sign(res[q_z  < PAR$fdr, estimate])))
print(table(sign(res[q_zc < PAR$fdr, est_c])))

# =============================================================================
# 7.  CANDIDATES = corrected z-reference at FDR 5%  (m = nrow(res))
# =============================================================================
hits <- res[q_zc < PAR$fdr]

# =============================================================================
# 8.  SANITY CHECK: collapsed model (6 rows, one per line)
# =============================================================================
# Cannot pseudo-replicate. Compare BOTH the effect and the standard error:
# pseudo-replication shows up as a much smaller SE in the nested model.
fit_collapsed <- function(site_id) {
  i <- match(site_id, rownames(counts_mod))
  keep_lib <- (counts_mod[i, ] + counts_unmod[i, ]) >= PAR$min_reads_library
  m <- tapply(as.numeric(counts_mod[i, keep_lib]),   metadata$cell_line[keep_lib], sum)
  u <- tapply(as.numeric(counts_unmod[i, keep_lib]), metadata$cell_line[keep_lib], sum)
  d <- data.frame(mod = as.numeric(m), unmod = as.numeric(u),
                  cell_line = names(m))
  d$subtype <- factor(meta_lines$subtype[match(d$cell_line, meta_lines$cell_line)],
                      levels = levels(metadata$subtype))
  d <- d[!is.na(d$mod), ]
  if (length(unique(d$subtype)) < 2L) return(NULL)
  f <- tryCatch(glmmTMB(cbind(mod, unmod) ~ subtype, family = betabinomial, data = d),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  cf <- summary(f)$coefficients$cond
  if (nrow(cf) < 2L) return(NULL)
  list(site = site_id, estimate_collapsed = cf[2, "Estimate"],
       se_collapsed = cf[2, "Std. Error"])
}
meta_lines <- unique(as.data.frame(metadata[, c("cell_line", "subtype")]))

set.seed(1)
check_sites <- sample(res$site, min(500L, nrow(res)))        # random, NOT only hits
coll <- rbindlist(Filter(Negate(is.null),
                         mclapply(check_sites, fit_collapsed, mc.cores = PAR$n_workers)))
cmp <- merge(res[, .(site, estimate, std_error)], coll, by = "site")
message(sprintf("nested vs collapsed: r(beta) = %.3f | slope = %.2f | median SE ratio nested/collapsed = %.2f (n = %d)",
                cor(cmp$estimate, cmp$estimate_collapsed, use = "complete.obs"),
                coef(lm(estimate ~ estimate_collapsed, cmp))[2],
                median(cmp$std_error / cmp$se_collapsed, na.rm = TRUE),
                nrow(cmp)))

# =============================================================================
# 9.  LEAVE-ONE-CELL-LINE-OUT on the hits, on CENTERED estimates
# =============================================================================
# Drop a whole line, refit, subtract the SAME b0, and require the centered effect
# to keep its sign and >= loo_min_fraction of the full centered effect.
# Simplification: b0 is not re-estimated for each refit (that would need all
# sites refitted without that line). Declare it in the methods.
loo_one_site <- function(site_id) {
  i <- match(site_id, rownames(counts_mod))
  vapply(line_levels, function(drop_line) {
    cols <- metadata$cell_line != drop_line
    d <- data.frame(
      mod       = as.numeric(counts_mod[i, cols]),
      unmod     = as.numeric(counts_unmod[i, cols]),
      subtype   = droplevels(metadata$subtype[cols]),
      cell_line = droplevels(metadata$cell_line[cols])
    )
    d <- d[d$mod + d$unmod >= PAR$min_reads_library, , drop = FALSE]
    d$cell_line <- droplevels(d$cell_line)
    if (length(unique(d$subtype)) < 2L) return(NA_real_)
    f <- tryCatch(glmmTMB(cbind(mod, unmod) ~ subtype + (1 | cell_line),
                          family = betabinomial, data = d),
                  error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    cf <- summary(f)$coefficients$cond
    if (nrow(cf) < 2L) NA_real_ else cf[2, "Estimate"]
  }, numeric(1))
}

loo_list <- mclapply(hits$site, function(s)
  tryCatch(loo_one_site(s), error = function(e) rep(NA_real_, length(line_levels))),
  mc.cores = PAR$n_workers)
loo_mat <- do.call(rbind, loo_list)                 # hits x 6 lines
colnames(loo_mat) <- line_levels
loo_c <- loo_mat - b0                               # centered LOO estimates

message(sprintf("LOO refits failed or not estimable: %d of %d",
                sum(!is.finite(loo_mat)), length(loo_mat)))

full_c <- hits$est_c
hits[, loo_ok := vapply(seq_len(nrow(hits)), function(k) {
  b <- loo_c[k, ]
  all(is.finite(b)) &&
    (!PAR$loo_sign_stable || all(sign(b) == sign(full_c[k]))) &&
    all(abs(b) >= PAR$loo_min_fraction * abs(full_c[k]))
}, logical(1))]
hits[, loo_worst_ratio := vapply(seq_len(nrow(hits)), function(k) {
  b <- loo_c[k, ]
  if (!all(is.finite(b))) NA_real_ else min(abs(b)) / abs(full_c[k])
}, numeric(1))]

message(sprintf("corrected hits: %d | LOO robust: %d (%.1f%% dropped)",
                nrow(hits), sum(hits$loo_ok), 100 * mean(!hits$loo_ok)))

# =============================================================================
# 10. REPORTING (percentage points from CENTERED estimates)
# =============================================================================
hits[, p_basal    := plogis(intercept)]
hits[, p_luminal  := plogis(intercept + est_c)]
hits[, delta_pp   := 100 * (p_luminal - p_basal)]
hits[, odds_ratio := exp(est_c)]

final <- hits[loo_ok == TRUE][order(p_zc)]

#fwrite(final, "dms_final_centered.tsv.gz", sep = "\t")
#fwrite(res,   "dms_all_sites.tsv.gz",      sep = "\t")

cat("\n================ SUMMARY ================\n")
cat(sprintf("sites tested (m for BH)             : %d\n", nrow(res)))
cat(sprintf("median normalised dispersion        : %.2f\n", median(res$disp_norm, na.rm = TRUE)))
cat(sprintf("median SD(cell_line), logit         : %.2f\n", median(res$sd_cellline, na.rm = TRUE)))
cat(sprintf("singular random effect              : %.1f%%\n", 100 * frac_singular))
cat(sprintf("global shift b0 (logit | pp)        : %.3f | %.1f\n", b0, b0_pp))
cat(sprintf("MAD of centered z                   : %.2f\n", sigma0))
cat(sprintf("hits raw z                          : %d\n", sum(res$q_z  < PAR$fdr)))
cat(sprintf("hits corrected z                    : %d\n", sum(res$q_zc < PAR$fdr)))
cat(sprintf("hits corrected t(%d)                 : %d\n", df_effective, sum(res$q_tc < PAR$fdr)))
cat(sprintf("after leave-one-line-out            : %d\n", nrow(final)))
cat(sprintf("median |delta| (percentage points)  : %.1f\n", median(abs(final$delta_pp))))
cat("=========================================\n")

print(head(final[, .(site, delta_pp = round(delta_pp, 1),
                     odds_ratio = round(odds_ratio, 2),
                     p_zc = signif(p_zc, 3), q_zc = signif(q_zc, 3),
                     sd_cellline = round(sd_cellline, 2),
                     loo_worst_ratio = round(loo_worst_ratio, 2))], 15))

