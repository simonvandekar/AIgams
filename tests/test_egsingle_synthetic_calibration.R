# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-05]
# [Prompt summary: test hypothesis that real cross-school slope heterogeneity
#  (not modeled by Model A's random-intercept-only structure) is causing the
#  fixed-effect SE mis-calibration seen in longitudinalCalibrationSim(). Build
#  a synthetic version of egsingle where the outcome is simulated from Model
#  A's own fitted parameters (parametric bootstrap via simulate.merMod), so
#  every coefficient except the intercept is exactly the population value by
#  construction, then rerun the same two-stage hierarchical resampling +
#  refit pipeline on that synthetic data and compare se_ratio/bias.]

devtools::load_all("/media/alsobig/AIgams/RESI")

egsingle <- local({
  e <- new.env()
  utils::data("egsingle", package = "mlmRev", envir = e)
  e$egsingle
})
egsingle$schoolid  <- factor(egsingle$schoolid)
egsingle$childid   <- factor(egsingle$childid)
egsingle$size_100  <- egsingle$size / 100
egsingle$lowinc_10 <- egsingle$lowinc / 10

form_A <- math ~ year + female + black + hispanic + size_100 + lowinc_10 +
  (1 | schoolid/childid)
ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# --- Fit Model A on the REAL data: this is the generative model whose
# fixed effects + variance components become the *known* truth below. ---
truth_fit <- lme4::lmer(form_A, data = egsingle, REML = FALSE, control = ctrl)
cat("=== Real-data Model A fit (used as simulation DGP) ===\n")
print(lme4::fixef(truth_fit))
print(lme4::VarCorr(truth_fit))

# --- Simulate a new outcome from Model A's own fitted parameters:
# use.u = FALSE draws FRESH random intercepts (school + child) from the
# estimated variance components, and fresh iid residual noise -- so by
# construction every covariate's slope is exactly the single population
# value everywhere (no real cross-school/child slope heterogeneity can leak
# in), while intercept variation is preserved (drawn from the estimated
# between-school/child intercept variance). ---
set.seed(2026)
synthetic <- egsingle
synthetic$math <- simulate(truth_fit, nsim = 1, use.u = FALSE)[[1]]

# Known-truth params = exactly what simulate() used (not a re-fit)
true_params <- RESI:::.longSimExtractParams(truth_fit)
true_fe_se  <- RESI:::.longSimExtractFeSe(truth_fit)

school_child_index <- RESI:::.longSimBuildIndex(synthetic, "schoolid", "childid")

n_schools.vec <- c(15, 30, 60, 90, 120)
nsim <- 1000L
mc.cores.reps <- 24L

message("Running synthetic-DGP calibration check for Model A ...")
rows <- lapply(n_schools.vec, function(n_schools) {
  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- RESI:::.longSimResampleHierarchical(school_child_index, n_schools,
                                                 "schoolid", "childid")
    fit <- tryCatch(
      suppressWarnings(suppressMessages(
        lme4::lmer(form_A, data = samp, REML = FALSE, control = ctrl))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    list(params = RESI:::.longSimExtractParams(fit),
         fe_se  = RESI:::.longSimExtractFeSe(fit))
  }, mc.cores = mc.cores.reps)

  reps <- Filter(Negate(is.null), reps_raw)
  n_success <- length(reps)

  param_names <- names(true_params)
  param_mat <- t(vapply(reps, function(r) r$params[param_names], numeric(length(param_names))))
  colnames(param_mat) <- param_names
  mc_mean <- colMeans(param_mat, na.rm = TRUE)
  mc_sd   <- apply(param_mat, 2, stats::sd, na.rm = TRUE)

  fe_names  <- names(true_fe_se)
  fe_se_mat <- t(vapply(reps, function(r) r$fe_se[fe_names], numeric(length(fe_names))))
  colnames(fe_se_mat) <- fe_names
  mean_model_se <- colMeans(fe_se_mat, na.rm = TRUE)

  is_fe <- param_names %in% fe_names
  data.frame(
    n_schools     = n_schools,
    n_success     = n_success,
    param         = param_names,
    truth         = true_params[param_names],
    mc_mean       = mc_mean[param_names],
    bias          = true_params[param_names] - mc_mean[param_names],
    mc_sd         = mc_sd[param_names],
    mean_model_se = ifelse(is_fe, mean_model_se[param_names], NA_real_),
    se_ratio      = ifelse(is_fe, mean_model_se[param_names] / mc_sd[param_names], NA_real_),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
})

summary_table <- do.call(rbind, rows)
saveRDS(summary_table, file = "tests/synthetic_calibration_summary.rds")

cat("\n=== Fixed effects: se_ratio by n_schools (synthetic DGP, Model A) ===\n")
fe <- summary_table[grepl("^beta_", summary_table$param), ]
wide <- reshape(fe[, c("param", "n_schools", "se_ratio")],
                 idvar = "param", timevar = "n_schools", direction = "wide")
print(wide, digits = 3)

cat("\n=== Variance components: relative bias (%) by n_schools ===\n")
vc <- summary_table[!grepl("^beta_", summary_table$param), ]
vc$rel_bias_pct <- 100 * vc$bias / vc$truth
widevc <- reshape(vc[, c("param", "n_schools", "rel_bias_pct")],
                    idvar = "param", timevar = "n_schools", direction = "wide")
print(widevc, digits = 3)
