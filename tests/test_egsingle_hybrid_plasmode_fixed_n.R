# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-24]
# [Prompt summary: hybrid plasmode follow-up to test_egsingle_hybrid_plasmode.R
#  -- same fixed-n_per_site deterministic design as
#  test_egsingle_synthetic_fixed_n.R, but the school and child random
#  intercepts are drawn by resampling WITH REPLACEMENT from the model's
#  rescaled BLUPs (mean zero, variance equal to the fitted variance
#  component) instead of a Normal distribution.]

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
school_var <- vc$vcov[vc$grp == "schoolid" & is.na(vc$var2)]
child_var  <- vc$vcov[vc$grp == "childid:schoolid" & is.na(vc$var2)]
resid_sd   <- sigma(truth_fit)

# Center each BLUP pool to exactly mean zero, then rescale so its variance
# equals the fitted variance component (correcting for BLUP shrinkage).
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

# Fixed covariate design with exactly n_per_site child profiles per site; only
# the random effects (resampled from BLUPs) and residuals vary by replicate.
fixed_design <- function(n_schools, n_per_site) {
  source_schools <- levels(egsingle$schoolid)[seq_len(n_schools)]
  out <- lapply(seq_along(source_schools), function(s) {
    source_rows <- egsingle[egsingle$schoolid == source_schools[s], , drop = FALSE]
    child_ids <- unique(source_rows$childid)
    selected_children <- child_ids[((seq_len(n_per_site) - 1L) %% length(child_ids)) + 1L]
    child_rows <- lapply(selected_children, function(child_id) {
      source_rows[source_rows$childid == child_id, , drop = FALSE]
    })
    rows <- do.call(rbind, child_rows)
    rows$schoolid <- factor(s)
    rows$childid <- factor(rep(paste0(s, "_", seq_len(n_per_site)),
                   vapply(child_rows, nrow, integer(1L))),
                           levels = paste0(s, "_", seq_len(n_per_site)))
    rows
  })
  out <- do.call(rbind, out)
  out$schoolid <- factor(out$schoolid)
  out$childid <- factor(out$childid)
  out
}

generate_outcome <- function(design) {
  school_effect <- sample(school_pool, nlevels(design$schoolid), replace = TRUE)
  child_effect  <- sample(child_pool, nlevels(design$childid), replace = TRUE)
  fixed_part <- model.matrix(~ year + female + black + hispanic + size_100 + lowinc_10,
                             data = design) %*% lme4::fixef(truth_fit)
  fixed_part + school_effect[design$schoolid] + child_effect[design$childid] +
    rnorm(nrow(design), sd = resid_sd)
}

n_schools.vec  <- c(15, 30, 60)
n_per_site.vec <- c(25, 50, 100, 200)
grid <- expand.grid(n_schools = n_schools.vec, n_per_site = n_per_site.vec)
nsim <- 1000L
mc.cores.reps <- 24L

set.seed(2026)
message("Running hybrid plasmode fixed-n_per_site calibration check for Model A ...")
rows <- lapply(seq_len(nrow(grid)), function(g) {
  n_schools  <- grid$n_schools[g]
  n_per_site <- grid$n_per_site[g]

  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- fixed_design(n_schools, n_per_site)
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
saveRDS(summary_table, file = "tests/hybrid_plasmode_fixed_n_summary.rds")

cat("\n=== Fixed effects: se_ratio by (n_schools, n_per_site), hybrid plasmode Model A ===\n")
fe <- summary_table[grepl("^beta_", summary_table$param), ]
for (p in unique(fe$param)) {
  cat("\n--", p, "--\n")
  sub <- fe[fe$param == p, ]
  wide <- reshape(sub[, c("n_schools", "n_per_site", "se_ratio")],
                   idvar = "n_schools", timevar = "n_per_site", direction = "wide")
  print(wide, digits = 3)
}

cat("\n=== Variance components: relative bias (%) by (n_schools, n_per_site) ===\n")
vc2 <- summary_table[!grepl("^beta_", summary_table$param), ]
vc2$rel_bias_pct <- 100 * vc2$bias / vc2$truth
for (p in unique(vc2$param)) {
  cat("\n--", p, "--\n")
  sub <- vc2[vc2$param == p, ]
  widevc <- reshape(sub[, c("n_schools", "n_per_site", "rel_bias_pct")],
                     idvar = "n_schools", timevar = "n_per_site", direction = "wide")
  print(widevc, digits = 3)
}
