# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-05]
# [Prompt summary: test the Model A fixed-effect SE calibration by generating
#  a fresh outcome on every replicate from the fitted fixed effects, school and
#  child random-intercept Normal distributions, and residual Normal
#  distribution. The covariate design is fixed and no sites or subjects are
#  resampled; compare se_ratio/bias after refitting Model A.]

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

# Known-truth params are the fitted Model A population parameters.
true_params <- RESI:::.longSimExtractParams(truth_fit)
true_fe_se  <- RESI:::.longSimExtractFeSe(truth_fit)

vc <- as.data.frame(lme4::VarCorr(truth_fit))
school_sd <- sqrt(vc$vcov[vc$grp == "schoolid" & is.na(vc$var2)])
child_sd  <- sqrt(vc$vcov[vc$grp == "childid:schoolid" & is.na(vc$var2)])
resid_sd  <- sigma(truth_fit)

# Use the observed covariate design deterministically. Each replicate draws
# fresh school effects, child effects, and residual errors; no rows are sampled.
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
  school_effect <- rnorm(nlevels(design$schoolid), sd = school_sd)
  child_effect <- rnorm(nlevels(design$childid), sd = child_sd)
  fixed_part <- model.matrix(~ year + female + black + hispanic + size_100 + lowinc_10,
                             data = design) %*% lme4::fixef(truth_fit)
  fixed_part + school_effect[design$schoolid] + child_effect[design$childid] +
    rnorm(nrow(design), sd = resid_sd)
}

n_schools.vec <- c(15, 30, 60, 90, 120)
nsim <- 1000L
mc.cores.reps <- 24L

set.seed(2026)
message("Running fresh-DGP synthetic calibration check for Model A ...")
rows <- lapply(n_schools.vec, function(n_schools) {
  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- fixed_design(n_schools)
    samp$math <- generate_outcome(samp)
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
