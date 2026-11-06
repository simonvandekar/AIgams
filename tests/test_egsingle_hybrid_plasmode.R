# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-24]
# [Prompt summary: hybrid plasmode simulation -- same fixed-design direct-DGP
#  setup as test_egsingle_synthetic_calibration.R (fixed effects taken from
#  the fitted Model A, deterministic covariate design, fresh residual Normal
#  errors), but the school and child random intercepts are drawn by
#  resampling WITH REPLACEMENT from the model's BLUPs (ranef() estimated
#  means) instead of from a Normal distribution. The BLUP pools are centered
#  to exactly mean zero and rescaled so their variance equals the
#  corresponding fitted variance component (BLUPs are shrunk toward zero
#  relative to the true random-effect variance, so this corrects for that
#  shrinkage before resampling).]

devtools::load_all("/media/alsobig/AIgams/RESI")

.hybridResiCalibrationSummary <- function(reps, truth_resi, settings) {
  sources <- list(
    coefficients = list(rep_name = "resi_coef", truth = truth_resi$coefficients),
    anova = list(rep_name = "resi_anova", truth = truth_resi$anova)
  )
  do.call(rbind, unlist(lapply(names(sources), function(table_name) {
    source <- sources[[table_name]]
    terms <- rownames(source$truth)
    lapply(terms, function(term) {
      values <- vapply(reps, function(rep) rep[[source$rep_name]][term, "RESI"], numeric(1L))
      first <- reps[[1L]][[source$rep_name]][term, , drop = FALSE]
      data.frame(
        as.list(settings), table = table_name, term = term,
        n = first$n, variation_level = first$`Variation level`,
        n_level = first$`n level`,
        truth_resi = source$truth[term, "RESI"],
        mc_mean_resi = mean(values, na.rm = TRUE),
        bias_resi = source$truth[term, "RESI"] - mean(values, na.rm = TRUE),
        mc_sd_resi = stats::sd(values, na.rm = TRUE),
        mc_se_resi = stats::sd(values, na.rm = TRUE) / sqrt(sum(is.finite(values))),
        n_success = sum(is.finite(values)),
        check.names = FALSE, row.names = NULL
      )
    })
  }), recursive = FALSE))
}

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

# Known-truth params are the fitted Model A population parameters.
true_params <- RESI:::.longSimExtractParams(truth_fit)
true_fe_se  <- RESI:::.longSimExtractFeSe(truth_fit)
# The target is the model-based RESI of the fitted DGP. A robust covariance
# from the real data would include real-data slope heterogeneity removed by
# this hybrid generator; replicate fits still use the separately chosen
# Use the same ML/model-based covariance estimator and student-level common
# normalizer for the full-data target and every simulated replicate.
truth_resi  <- RESI::resi_pe(truth_fit, anova = TRUE,
                              vcovfunc = stats::vcov, level = 2)

vc <- as.data.frame(lme4::VarCorr(truth_fit))
school_var <- vc$vcov[vc$grp == "schoolid" & is.na(vc$var2)]
child_var  <- vc$vcov[vc$grp == "childid:schoolid" & is.na(vc$var2)]
resid_sd   <- sigma(truth_fit)

# Center each BLUP pool to exactly mean zero, then rescale so its variance
# equals the fitted variance component (BLUPs are shrunk toward zero, so
# their raw empirical variance is smaller than the true random-effect
# variance).
rescale_blup <- function(blup, target_var) {
  centered <- blup - mean(blup)
  centered * sqrt(target_var / stats::var(centered))
}
re <- lme4::ranef(truth_fit)
school_pool <- rescale_blup(re$schoolid[[1]], school_var)
child_pool  <- rescale_blup(re[["childid:schoolid"]][[1]], child_var)
cat("school BLUP pool: n =", length(school_pool),
    " mean =", mean(school_pool), " var =", stats::var(school_pool), "\n")
cat("child BLUP pool: n =", length(child_pool),
    " mean =", mean(child_pool), " var =", stats::var(child_pool), "\n")

# Use the observed covariate design deterministically. Each replicate resamples
# fresh school/child effects from the BLUP pools (with replacement) and draws
# fresh residual errors; no rows of the real data are resampled.
fixed_design <- function(n_schools) {
  source_schools <- levels(egsingle$schoolid)
  out <- lapply(seq_len(n_schools), function(s) {
    source_school <- source_schools[((s - 1L) %% length(source_schools)) + 1L]
    rows <- egsingle[egsingle$schoolid == source_school, , drop = FALSE]
    rows$schoolid <- factor(s)
    rows$childid <- factor(paste0(s, "_", rows$childid))
    rows
  })
  do.call(rbind, out)
}

generate_outcome <- function(design) {
  school_effect <- sample(school_pool, nlevels(design$schoolid), replace = TRUE)
  child_effect  <- sample(child_pool, nlevels(design$childid), replace = TRUE)
  fixed_part <- model.matrix(~ year + female + black + hispanic + size_100 + lowinc_10,
                             data = design) %*% lme4::fixef(truth_fit)
  fixed_part + school_effect[design$schoolid] + child_effect[design$childid] +
    rnorm(nrow(design), sd = resid_sd)
}

n_schools.vec <- getOption("resi.n_schools", c(15, 30, 60, 90, 120))
nsim <- getOption("resi.nsim", 1000L)
mc.cores.reps <- getOption("resi.mc.cores.reps", 24L)

set.seed(2026)
message("Running hybrid plasmode (BLUP-resampled random effects) calibration check for Model A ...")
resi_rows <- vector("list", length(n_schools.vec))
names(resi_rows) <- as.character(n_schools.vec)
rows <- lapply(n_schools.vec, function(n_schools) {
  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- fixed_design(n_schools)
    samp$math <- generate_outcome(samp)
    fit <- tryCatch(
      suppressWarnings(suppressMessages(
        lme4::lmer(form_A, data = samp, REML = FALSE, control = ctrl))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
        point <- RESI::resi_pe(fit, anova = TRUE,
                   vcovfunc = stats::vcov, level = 2)
    list(params = RESI:::.longSimExtractParams(fit),
          fe_se  = RESI:::.longSimExtractFeSe(fit),
          resi_coef = point$coefficients,
          resi_anova = point$anova)
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

  resi_rows[[as.character(n_schools)]] <<- .hybridResiCalibrationSummary(
    reps, truth_resi, list(n_schools = n_schools))

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
saveRDS(summary_table, file = getOption("resi.output.summary", "tests/hybrid_plasmode_summary.rds"))
resi_summary <- do.call(rbind, resi_rows)
saveRDS(resi_summary, file = getOption("resi.output.resi_summary", "tests/hybrid_plasmode_resi_summary.rds"))

cat("\n=== Fixed effects: se_ratio by n_schools (hybrid plasmode, Model A) ===\n")
fe <- summary_table[grepl("^beta_", summary_table$param), ]
wide <- reshape(fe[, c("param", "n_schools", "se_ratio")],
                 idvar = "param", timevar = "n_schools", direction = "wide")
print(wide, digits = 3)

cat("\n=== Variance components: relative bias (%) by n_schools ===\n")
vc2 <- summary_table[!grepl("^beta_", summary_table$param), ]
vc2$rel_bias_pct <- 100 * vc2$bias / vc2$truth
widevc <- reshape(vc2[, c("param", "n_schools", "rel_bias_pct")],
                    idvar = "param", timevar = "n_schools", direction = "wide")
print(widevc, digits = 3)

cat("\n=== RESI calibration: full-data reference vs Monte Carlo ===\n")
print(resi_summary[, c("n_schools", "table", "term", "variation_level", "n_level", "n",
           "truth_resi", "mc_mean_resi", "bias_resi", "mc_sd_resi")],
  row.names = FALSE, digits = 3)
