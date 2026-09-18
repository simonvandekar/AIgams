# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-18]
# [Prompt summary: user hypothesized that "the synthetic simulations sampling
#  should not have this problem" -- i.e. since the fixed-n_per_site design's
#  se_ratio mis-calibration was explained as a with-replacement bootstrap
#  artifact (not a DGP/model-misspecification issue), it should show up
#  identically for the REAL (plasmode) data too, not just the synthetic DGP.
#  This script reruns test_egsingle_synthetic_fixed_n.R's exact grid and
#  pipeline on the real egsingle data directly (no simulate() step) to confirm
#  the artifact is a property of the resampling design, not the synthetic DGP.]

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

# --- Fit Model A on the REAL data: this is both the truth reference AND the
# data resampled from (no synthetic simulate() step -- true plasmode). ---
truth_fit <- lme4::lmer(form_A, data = egsingle, REML = FALSE, control = ctrl)
cat("=== Real-data Model A fit (truth reference) ===\n")
print(lme4::fixef(truth_fit))
print(lme4::VarCorr(truth_fit))

true_params <- RESI:::.longSimExtractParams(truth_fit)
true_fe_se  <- RESI:::.longSimExtractFeSe(truth_fit)

school_child_index <- RESI:::.longSimBuildIndex(egsingle, "schoolid", "childid")

n_schools.vec  <- c(15, 30, 60)
n_per_site.vec <- c(25, 50, 100, 200)
grid <- expand.grid(n_schools = n_schools.vec, n_per_site = n_per_site.vec)
nsim <- 1000L
mc.cores.reps <- 24L

message("Running plasmode (real-data) fixed-n_per_site calibration check for Model A ...")
rows <- lapply(seq_len(nrow(grid)), function(g) {
  n_schools  <- grid$n_schools[g]
  n_per_site <- grid$n_per_site[g]

  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- RESI:::.longSimResampleFixedN(school_child_index, n_schools, n_per_site,
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
    n_per_site    = n_per_site,
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
saveRDS(summary_table, file = "tests/plasmode_fixed_n_summary.rds")

cat("\n=== Fixed effects: se_ratio by (n_schools, n_per_site), plasmode (real data) Model A ===\n")
fe <- summary_table[grepl("^beta_", summary_table$param), ]
for (p in unique(fe$param)) {
  cat("\n--", p, "--\n")
  sub <- fe[fe$param == p, ]
  wide <- reshape(sub[, c("n_schools", "n_per_site", "se_ratio")],
                   idvar = "n_schools", timevar = "n_per_site", direction = "wide")
  print(wide, digits = 3)
}

cat("\n=== Variance components: relative bias (%) by (n_schools, n_per_site) ===\n")
vc <- summary_table[!grepl("^beta_", summary_table$param), ]
vc$rel_bias_pct <- 100 * vc$bias / vc$truth
for (p in unique(vc$param)) {
  cat("\n--", p, "--\n")
  sub <- vc[vc$param == p, ]
  widevc <- reshape(sub[, c("n_schools", "n_per_site", "rel_bias_pct")],
                     idvar = "n_schools", timevar = "n_per_site", direction = "wide")
  print(widevc, digits = 3)
}
