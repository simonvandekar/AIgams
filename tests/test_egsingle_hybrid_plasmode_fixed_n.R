# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-24]
# [Prompt summary: hierarchical with-replacement covariate resampling of
#  schools and students, fresh BLUP-pool random effects and residuals, and
#  REML refits for the joint n_schools by n_per_site calibration.]

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
      replicate_n <- vapply(reps, function(rep) {
        rep[[source$rep_name]][term, "n"]
      }, numeric(1L))
      data.frame(
        as.list(settings), table = table_name, term = term,
        n = mean(replicate_n), n_sd = stats::sd(replicate_n),
        variation_level = first$`Variation level`,
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
# The full-sample ML fit defines the model-based RESI reference. Replicate
# fits use REML, as requested, with the same model-based covariance method.
truth_resi  <- RESI::resi_pe(truth_fit, anova = TRUE,
                              vcovfunc = stats::vcov, level = 2)

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

# Resample schools and then students within each selected school. The number
# sampled per school is proportional to its observed size and scaled so the
# expected average is n_per_site.
school_child_counts <- vapply(levels(egsingle$schoolid), function(school) {
  length(unique(egsingle$childid[egsingle$schoolid == school]))
}, integer(1L))
mean_site_size <- mean(school_child_counts)
fixed_design <- function(n_schools, n_per_site) {
  school_profiles <- levels(egsingle$schoolid)
  source_schools <- sample(school_profiles, n_schools, replace = TRUE)
  out <- lapply(seq_len(n_schools), function(s) {
    source_rows <- egsingle[egsingle$schoolid == source_schools[s], , drop = FALSE]
    child_ids <- unique(source_rows$childid)
    n_sampled_children <- max(1L, as.integer(round(
      n_per_site * length(child_ids) / mean_site_size)))
    selected_children <- sample(child_ids, n_sampled_children, replace = TRUE)
    child_rows <- lapply(selected_children, function(child_id) {
      source_rows[source_rows$childid == child_id, , drop = FALSE]
    })
    rows <- do.call(rbind, child_rows)
    rows$schoolid <- factor(s)
    rows$childid <- factor(rep(paste0(s, "_", seq_len(n_sampled_children)),
                   vapply(child_rows, nrow, integer(1L))),
                           levels = paste0(s, "_", seq_len(n_sampled_children)))
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

n_schools.vec  <- getOption("resi.n_schools", c(15, 30, 60, 100))
n_per_site.vec <- getOption("resi.n_per_site", c(25, 100, 200))
grid <- getOption("resi.grid", expand.grid(
  n_schools = n_schools.vec, n_per_site = n_per_site.vec
))
nsim <- getOption("resi.nsim", 1000L)
mc.cores.reps <- getOption("resi.mc.cores.reps", 24L)

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
message("Running hierarchical covariate-bootstrap REML calibration for Model A ...")
resi_rows <- vector("list", nrow(grid))
rows <- lapply(seq_len(nrow(grid)), function(g) {
  n_schools  <- grid$n_schools[g]
  n_per_site <- grid$n_per_site[g]

  reps_raw <- parallel::mclapply(seq_len(nsim), function(i) {
    samp <- fixed_design(n_schools, n_per_site)
    samp$math <- generate_outcome(samp)
    fit <- tryCatch(
      suppressWarnings(suppressMessages(
        lme4::lmer(form_A, data = samp, REML = TRUE, control = ctrl))),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
      point <- RESI::resi_pe(fit, anova = TRUE,
                   vcovfunc = stats::vcov, level = 2)
    list(params = RESI:::.longSimExtractParams(fit),
          fe_se  = RESI:::.longSimExtractFeSe(fit),
          n_schools = nlevels(samp$schoolid),
          n_students = nlevels(samp$childid),
          resi_coef = point$coefficients,
          resi_anova = point$anova)
  }, mc.cores = mc.cores.reps, mc.set.seed = TRUE)

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

  school_level_terms <- c("beta_(Intercept)", "beta_size_100", "beta_lowinc_10")
  n_reference <- setNames(ifelse(
    fe_names %in% school_level_terms,
    nlevels(egsingle$schoolid), length(unique(egsingle$childid))), fe_names)
  se_level <- setNames(ifelse(fe_names %in% school_level_terms,
                              "schools", "students"), fe_names)
  n_scaling <- setNames(ifelse(fe_names %in% school_level_terms,
                               n_schools, n_schools * n_per_site), fe_names)
  truth_se_scaled <- true_fe_se[fe_names] * sqrt(n_reference[fe_names])
  mc_error_sd_scaled <- mc_sd[fe_names] * sqrt(n_scaling[fe_names])
  scaled_mc_error_sd_ratio <- mc_error_sd_scaled / truth_se_scaled

  resi_rows[[g]] <<- .hybridResiCalibrationSummary(
    reps, truth_resi, list(n_schools = n_schools, n_per_site = n_per_site))

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
    se_n_level    = unname(se_level[param_names]),
    n_reference   = unname(n_reference[param_names]),
    n_scaling     = unname(n_scaling[param_names]),
    truth_se_scaled = unname(truth_se_scaled[param_names]),
    mc_error_sd_scaled = unname(mc_error_sd_scaled[param_names]),
    scaled_mc_error_sd_ratio = unname(scaled_mc_error_sd_ratio[param_names]),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
})

summary_table <- do.call(rbind, rows)
saveRDS(summary_table, file = getOption("resi.output.summary", "tests/hybrid_plasmode_resampled_covariates_summary.rds"))
resi_summary <- do.call(rbind, resi_rows)
saveRDS(resi_summary, file = getOption("resi.output.resi_summary", "tests/hybrid_plasmode_resampled_covariates_resi_summary.rds"))

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

cat("\n=== RESI calibration: full-data reference vs Monte Carlo ===\n")
print(resi_summary[, c("n_schools", "n_per_site", "table", "term",
                       "variation_level", "n_level", "n", "truth_resi", "mc_mean_resi",
                       "bias_resi", "mc_sd_resi")], row.names = FALSE, digits = 3)
