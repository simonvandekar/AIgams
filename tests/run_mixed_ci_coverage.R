# Coverage of the mixed-model RESI confidence intervals (RESI:::.resi_ci_mixed)
# in the egsingle hybrid plasmode (Model A). Each replicate computes intervals
# with the model-based and CR2 (school-clustered) covariances and the normal,
# QF and CF methods; coverage is scored against the Lyu-Welsh limit S_Lyu
# (g, m -> Inf) and the fixed-m target S_m (expectedDesignResi). Writes
# resiLongCalibrationSim/mixed_ci_coverage.rds. Run from the AIgams root:
#   Rscript tests/run_mixed_ci_coverage.R

devtools::load_all("RESI")

egsingle <- local({
  e <- new.env()
  utils::data("egsingle", package = "mlmRev", envir = e)
  e$egsingle
})
egsingle$size_100 <- egsingle$size / 100
egsingle$lowinc_10 <- egsingle$lowinc / 10

form_A <- math ~ year + female + black + hispanic + size_100 + lowinc_10 +
  (1 | schoolid/childid)
methods <- c("normal", "qf", "cf")

ci_rep <- function(fit, design) tryCatch(ci_rep_raw(fit, design), error = function(e) conditionMessage(e))
ci_rep_raw <- function(fit, design) {
  pe <- resi_pe(fit, anova = TRUE, vcovfunc = stats::vcov, level = 2)
  n_coef <- setNames(pe$coefficients$n, rownames(pe$coefficients))
  n_anova <- setNames(pe$anova$n, rownames(pe$anova))
  one <- function(robust) {
    Vhat <- if (robust) RESI:::.mixedCiCR2(fit, "schoolid") else as.matrix(stats::vcov(fit))
    est <- if (robust) resi_pe(fit, anova = TRUE, level = 2,
                               vcovfunc = function(object, ...) Vhat) else pe
    ci <- RESI:::.resi_ci_mixed(fit, Vhat, robust, n_coef, n_anova, methods,
                                cluster = "schoolid")
    p <- length(n_coef)
    se <- list(coefficients = ci$se[seq_len(p)], anova = ci$se[-seq_len(p)])
    do.call(rbind, lapply(methods, function(method) {
      do.call(rbind, lapply(c("coefficients", "anova"), function(tab) {
        lim <- ci[[method]][[tab]]
        data.frame(covariance = if (robust) "CR2" else "model", method = method,
                   table = tab, term = rownames(lim),
                   est = est[[tab]][rownames(lim), "RESI"], se = unname(se[[tab]]),
                   lower = lim[, 1L], upper = lim[, 2L], row.names = NULL)
      }))
    }))
  }
  rbind(one(FALSE), one(TRUE))
}

n_schools <- getOption("resi.n_schools", c(15, 100))
n_per_site <- getOption("resi.n_per_site", c(25, 200))
out <- hybridPlasmodeSim(
  form_A, egsingle, school = "schoolid", child = "childid",
  n_schools = n_schools, n_per_site = n_per_site,
  nsim = getOption("resi.nsim", 1000L),
  mc.cores = getOption("resi.mc.cores", 24L),
  seed = 2027, keep_reps = TRUE, rep_fun = ci_rep)

truth_table <- function(resi_obj) {
  rbind(data.frame(table = "coefficients", term = rownames(resi_obj$coefficients),
                   truth = resi_obj$coefficients[, "RESI"]),
        data.frame(table = "anova", term = rownames(resi_obj$anova),
                   truth = resi_obj$anova[, "RESI"]))
}
s_lyu <- truth_table(out$truth_resi)
between_terms <- c("(Intercept)", "size_100", "lowinc_10")

coverage <- do.call(rbind, lapply(seq_len(nrow(out$grid)), function(cell) {
  g <- out$grid$n_schools[cell]
  m <- out$grid$n_per_site[cell]
  s_m <- truth_table(expectedDesignResi(out$truth_fit, m))
  failed <- vapply(out$reps[[cell]], function(r) is.character(r$extra), logical(1L))
  if (any(failed)) message("Cell ", cell, ": CI computation failed in ", sum(failed), " replicates: ",
                           paste(unique(vapply(out$reps[[cell]][failed], `[[`, "", "extra")), collapse = "; "))
  reps <- out$reps[[cell]][!failed]
  long <- do.call(rbind, lapply(seq_along(reps), function(i) cbind(rep = i, reps[[i]]$extra)))
  long$S_Lyu <- s_lyu$truth[match(paste(long$table, long$term), paste(s_lyu$table, s_lyu$term))]
  long$S_m <- s_m$truth[match(paste(long$table, long$term), paste(s_m$table, s_m$term))]
  keys <- unique(long[, c("covariance", "method", "table", "term")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    d <- merge(keys[k, ], long)
    ok <- is.finite(d$lower) & is.finite(d$upper)
    cov_lyu <- mean(d$lower[ok] <= d$S_Lyu[ok] & d$S_Lyu[ok] <= d$upper[ok])
    cov_m <- mean(d$lower[ok] <= d$S_m[ok] & d$S_m[ok] <= d$upper[ok])
    data.frame(n_schools = g, n_per_site = m, keys[k, ], row.names = NULL,
               level = if (d$term[1] %in% between_terms) "between" else "within",
               n_rep = sum(ok), n_failed = sum(failed), S_Lyu = d$S_Lyu[1], S_m = d$S_m[1],
               mean_est = mean(d$est), mc_sd_est = stats::sd(d$est),
               mean_se = mean(d$se), median_width = stats::median(d$upper[ok] - d$lower[ok]),
               coverage_S_Lyu = cov_lyu, coverage_S_m = cov_m,
               mc_se_coverage = sqrt(0.95 * 0.05 / sum(ok)))
  }))
}))
out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
saveRDS(coverage, file.path(out_dir, "mixed_ci_coverage.rds"))

cat("\n=== Coverage of S_m (fixed-m target), nominal 0.95 ===\n")
print(xtabs(coverage_S_m ~ term + interaction(covariance, method) + interaction(n_schools, n_per_site),
            subset(coverage, table == "coefficients")), digits = 3)
cat("\n=== Coverage of S_Lyu, nominal 0.95 ===\n")
print(xtabs(coverage_S_Lyu ~ term + interaction(covariance, method) + interaction(n_schools, n_per_site),
            subset(coverage, table == "coefficients")), digits = 3)
cat("\n=== Mean estimated SE / MC SD of RESI ===\n")
print(with(subset(coverage, method == "normal" & table == "coefficients"),
           xtabs(mean_se / mc_sd_est ~ term + interaction(covariance, n_schools, n_per_site))), digits = 3)
