# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-04]
# [Prompt summary: Prototype Goal 5a calibration check -- does cluster
#  (subject-level) plasmode resampling from egsingle recover the full-data
#  population coefficients, variance components, and their standard errors,
#  for Model A (random intercepts only) and Model B (random slope for year)?]
#
# This is an ad-hoc, exploratory script (not yet a package function). If the
# approach looks right, it should be formalized into RESI/R/simulations.R as
# something like `longitudinalCalibrationSim()`.

suppressPackageStartupMessages({
  library(lme4)
  library(mlmRev)
})

data(egsingle, package = "mlmRev")
egsingle$childid  <- factor(egsingle$childid)
egsingle$schoolid <- factor(egsingle$schoolid)

fixed_form <- math ~ year + female + black + hispanic + size + lowinc
form_A <- update(fixed_form, . ~ . + (1 | schoolid/childid))
form_B <- update(fixed_form, . ~ . + (1 + year | childid) + (1 | schoolid))

## ---- helpers ---------------------------------------------------------------

# Flatten fixef() + residual variance + VarCorr() into one named vector so
# Monte Carlo summaries can be computed generically across model A/B.
.extract_params <- function(fit) {
  fe <- lme4::fixef(fit)
  names(fe) <- paste0("beta_", names(fe))

  vc <- as.data.frame(lme4::VarCorr(fit))
  # var1/var2 label the (co)variance parameter within each grouping factor
  vc_names <- ifelse(is.na(vc$var2),
                      paste0("var_", vc$grp, "_", vc$var1),
                      paste0("cov_", vc$grp, "_", vc$var1, "_", vc$var2))
  vc_vals <- vc$vcov
  names(vc_vals) <- vc_names

  c(fe, vc_vals)
}

# Model-based SE for the fixed effects only (variance-component SEs are not
# returned by summary(); left NA here and can be added later if needed, e.g.
# via bootMer or numerical Hessian).
.extract_fe_se <- function(fit) {
  se <- sqrt(diag(as.matrix(vcov(fit))))
  names(se) <- paste0("beta_", names(se))
  se
}

# Cluster (subject-level) plasmode resample: draw n children with replacement,
# keep each child's full set of rows (and its school membership) intact. This
# preserves the multilevel/repeated-measures structure by construction.
.resample_egsingle <- function(data, n) {
  ids <- levels(data$childid)
  draw <- sample(ids, n, replace = TRUE)
  # re-key resampled duplicates so lme4 treats repeated draws of the same
  # child as distinct clusters (else refit would only ever see <=1721 groups)
  out <- do.call(rbind, lapply(seq_along(draw), function(i) {
    rows <- data[data$childid == draw[i], ]
    rows$childid <- factor(paste0(rows$childid, "_rep", i))
    rows
  }))
  out$schoolid <- factor(out$schoolid)
  out
}

## ---- population ("truth") values from the full data -------------------------

message("Fitting full-data population models (Model A, Model B)...")
full_A <- lmer(form_A, data = egsingle, REML = FALSE)
full_B <- lmer(form_B, data = egsingle, REML = FALSE)

truth_A <- .extract_params(full_A)
truth_B <- .extract_params(full_B)
se_A    <- .extract_fe_se(full_A)
se_B    <- .extract_fe_se(full_B)

## ---- calibration simulation --------------------------------------------------

run_calibration <- function(form, data, truth, se_truth, n.vec, nsim, label) {
  results <- list()
  for (n in n.vec) {
    reps <- vector("list", nsim)
    for (i in seq_len(nsim)) {
      samp <- .resample_egsingle(data, n)
      fit <- tryCatch(
        suppressWarnings(suppressMessages(lmer(form, data = samp, REML = FALSE))),
        error = function(e) NULL)
      reps[[i]] <- if (is.null(fit)) NULL else list(
        params = .extract_params(fit),
        fe_se  = .extract_fe_se(fit)
      )
    }
    reps <- Filter(Negate(is.null), reps)
    n_success <- length(reps)
    if (n_success == 0L) { warning(label, " n=", n, ": all replicates failed"); next }

    param_names <- names(truth)
    param_mat <- t(vapply(reps, function(r) r$params[param_names], numeric(length(param_names))))
    colnames(param_mat) <- param_names
    mc_mean <- colMeans(param_mat, na.rm = TRUE)
    mc_sd   <- apply(param_mat, 2, sd, na.rm = TRUE)

    fe_names <- names(se_truth)
    fe_se_mat <- t(vapply(reps, function(r) r$fe_se[fe_names], numeric(length(fe_names))))
    colnames(fe_se_mat) <- fe_names
    mean_model_se <- colMeans(fe_se_mat, na.rm = TRUE)

    results[[as.character(n)]] <- data.frame(
      n = n, n_success = n_success, param = param_names,
      truth = truth[param_names], mc_mean = mc_mean[param_names],
      bias = truth[param_names] - mc_mean[param_names],
      mc_sd = mc_sd[param_names],
      # SE calibration ratio only defined for fixed effects (fe_se available)
      mean_model_se = ifelse(param_names %in% fe_names, mean_model_se[param_names], NA),
      se_ratio = ifelse(param_names %in% fe_names,
                         mean_model_se[param_names] / mc_sd[param_names], NA),
      row.names = NULL
    )
  }
  out <- do.call(rbind, results)
  out$model <- label
  out
}

nsim <- 20L
n.vec <- c(200, 500, 1000, 1721)

set.seed(20260904)
calib_A <- run_calibration(form_A, egsingle, truth_A, se_A, n.vec, nsim, "A_intercept_only")
calib_B <- run_calibration(form_B, egsingle, truth_B, se_B, n.vec, nsim, "B_random_slope")

calib <- rbind(calib_A, calib_B)
print(calib, digits = 3)

saveRDS(calib, file.path(tempdir(), "egsingle_calibration_prototype.rds"))
message("Prototype calibration table saved to ", file.path(tempdir(), "egsingle_calibration_prototype.rds"))
