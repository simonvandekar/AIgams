# Shared driver for the coverage of the mixed-model RESI confidence intervals
# (RESI:::.resi_ci_mixed) in the hybrid plasmode (hybridPlasmodeSim). Nothing
# here is specific to a dataset: terms, their levels and the targets come from
# the model fit. Used by tests/run_mixed_ci_coverage.R (egsingle) and
# tests/run_star_validation.R (mlmRev::star).

ci_methods <- c("normal", "qf", "cf")

# Intervals for one fitted replicate with the model-based and CR2 covariances.
ci_table <- function(fit, school, engine) {
  pe <- resi_pe(fit, anova = TRUE, vcovfunc = stats::vcov, level = 2)
  n_coef <- setNames(pe$coefficients$n, rownames(pe$coefficients))
  n_anova <- setNames(pe$anova$n, rownames(pe$anova))
  n_level <- c(setNames(pe$coefficients$`n level`, paste("coefficients", rownames(pe$coefficients))),
               setNames(pe$anova$`n level`, paste("anova", rownames(pe$anova))))
  one <- function(robust) {
    Vhat <- if (robust) RESI:::.mixedCiCR2(fit, school) else as.matrix(stats::vcov(fit))
    est <- if (robust) resi_pe(fit, anova = TRUE, level = 2,
                               vcovfunc = function(object, ...) Vhat) else pe
    ci <- RESI:::.resi_ci_mixed(fit, Vhat, robust, n_coef, n_anova, ci_methods,
                                cluster = school)
    p <- length(n_coef)
    se <- list(coefficients = ci$se[seq_len(p)], anova = ci$se[-seq_len(p)])
    do.call(rbind, lapply(ci_methods, function(method) {
      do.call(rbind, lapply(c("coefficients", "anova"), function(tab) {
        lim <- ci[[method]][[tab]]
        data.frame(engine = engine, covariance = if (robust) "CR2" else "model",
                   method = method, table = tab, term = rownames(lim),
                   n_level = unname(n_level[paste(tab, rownames(lim))]),
                   est = est[[tab]][rownames(lim), "RESI"], se = unname(se[[tab]]),
                   lower = lim[, 1L], upper = lim[, 2L], row.names = NULL)
      }))
    }))
  }
  rbind(one(FALSE), one(TRUE))
}

# rep_fun for hybridPlasmodeSim(): intervals from the replicate's lmer fit
# and, if `lme = TRUE`, from an nlme::lme refit of the same model. Errors are
# returned as strings so that a failure is recorded rather than fatal.
make_ci_rep <- function(formula, school, child, lme = FALSE) {
  fixed <- reformulas::nobars(formula)
  random <- stats::as.formula(paste("~ 1 |", school, "/", child))
  function(fit, design) {
    lmer_tab <- tryCatch(ci_table(fit, school, "lmer"), error = function(e) conditionMessage(e))
    if (!lme) return(lmer_tab)
    lme_tab <- tryCatch({
      # Embed the data in the call so nlme::getData() can recover it later.
      fit_lme <- do.call(nlme::lme, list(fixed = fixed, random = random, data = design,
                                         method = "REML"))
      ci_table(fit_lme, school, "lme")
    }, error = function(e) paste("lme:", conditionMessage(e)))
    if (is.character(lmer_tab)) return(lmer_tab)
    if (is.character(lme_tab)) {
      attr(lmer_tab, "lme_error") <- lme_tab
      return(lmer_tab)
    }
    rbind(lmer_tab, lme_tab)
  }
}

truth_table <- function(resi_obj) {
  rbind(data.frame(table = "coefficients", term = rownames(resi_obj$coefficients),
                   truth = resi_obj$coefficients[, "RESI"]),
        data.frame(table = "anova", term = rownames(resi_obj$anova),
                   truth = resi_obj$anova[, "RESI"]))
}

# Coverage of S_Lyu (g, m -> Inf limit) and S_m (fixed-m expected-design RESI)
# per cell, engine, covariance, method and term. `school`/`child` name the
# lme4 grouping factors (e.g. "schoolid", "childid:schoolid").
summarize_coverage <- function(out, school, child) {
  s_lyu <- truth_table(out$truth_resi)
  key <- function(d) paste(d$table, d$term)
  do.call(rbind, lapply(seq_len(nrow(out$grid)), function(cell) {
    g <- out$grid$n_schools[cell]
    m <- out$grid$n_per_site[cell]
    s_m <- truth_table(expectedDesignResi(out$truth_fit, m, school, child))
    extras <- lapply(out$reps[[cell]], `[[`, "extra")
    failed <- vapply(extras, is.character, logical(1L))
    if (any(failed)) message("Cell ", cell, ": CI computation failed in ", sum(failed),
                             " replicates: ", paste(unique(unlist(extras[failed])), collapse = "; "))
    lme_err <- unlist(lapply(extras[!failed], attr, "lme_error"))
    if (length(lme_err)) message("Cell ", cell, ": lme failed in ", length(lme_err),
                                 " replicates: ", paste(unique(lme_err), collapse = "; "))
    long <- do.call(rbind, lapply(which(!failed), function(i) cbind(rep = i, extras[[i]])))
    long$S_Lyu <- s_lyu$truth[match(key(long), key(s_lyu))]
    long$S_m <- s_m$truth[match(key(long), key(s_m))]
    keys <- unique(long[, c("engine", "covariance", "method", "table", "term")])
    do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
      d <- merge(keys[k, ], long)
      ok <- is.finite(d$lower) & is.finite(d$upper)
      data.frame(n_schools = g, n_per_site = m, keys[k, ], row.names = NULL,
                 level = if (d$n_level[1] == school) "between" else "within",
                 n_rep = sum(ok), n_failed = sum(failed), S_Lyu = d$S_Lyu[1], S_m = d$S_m[1],
                 mean_est = mean(d$est), mc_sd_est = stats::sd(d$est),
                 mean_se = mean(d$se), median_width = stats::median(d$upper[ok] - d$lower[ok]),
                 coverage_S_Lyu = mean(d$lower[ok] <= d$S_Lyu[ok] & d$S_Lyu[ok] <= d$upper[ok]),
                 coverage_S_m = mean(d$lower[ok] <= d$S_m[ok] & d$S_m[ok] <= d$upper[ok]),
                 mc_se_coverage = sqrt(0.95 * 0.05 / sum(ok)))
    }))
  }))
}

# Per-cell maximum absolute lme - lmer difference in the RESI estimate and
# interval endpoints over replicates where both engines succeeded.
summarize_engine_agreement <- function(out) {
  do.call(rbind, lapply(seq_len(nrow(out$grid)), function(cell) {
    tabs <- Filter(function(x) is.data.frame(x) && all(c("lmer", "lme") %in% x$engine),
                   lapply(out$reps[[cell]], `[[`, "extra"))
    if (!length(tabs)) return(NULL)
    d <- do.call(rbind, lapply(tabs, function(x) {
      a <- x[x$engine == "lmer", ]; b <- x[x$engine == "lme", ]
      b <- b[match(paste(a$covariance, a$method, a$table, a$term),
                   paste(b$covariance, b$method, b$table, b$term)), ]
      data.frame(covariance = a$covariance, method = a$method,
                 d_est = abs(a$est - b$est), d_lower = abs(a$lower - b$lower),
                 d_upper = abs(a$upper - b$upper))
    }))
    agg <- aggregate(cbind(d_est, d_lower, d_upper) ~ covariance + method, d, max)
    cbind(n_schools = out$grid$n_schools[cell], n_per_site = out$grid$n_per_site[cell],
          n_rep = length(tabs), agg)
  }))
}

print_coverage <- function(coverage) {
  for (eng in unique(coverage$engine)) {
    cov_e <- subset(coverage, engine == eng & table == "coefficients")
    cat("\n=== [", eng, "] Coverage of S_m (fixed-m target), nominal 0.95 ===\n")
    print(xtabs(coverage_S_m ~ term + interaction(covariance, method) +
                  interaction(n_schools, n_per_site), cov_e), digits = 3)
    cat("\n=== [", eng, "] Mean estimated SE / MC SD of RESI ===\n")
    print(with(subset(cov_e, method == "normal"),
               xtabs(mean_se / mc_sd_est ~ term + interaction(covariance, n_schools, n_per_site))),
          digits = 3)
  }
}
