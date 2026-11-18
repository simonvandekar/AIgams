# [AI-generated: GitHub Copilot, 2026-11-15]
# [Prompt summary: coverage simulation of the mixed-model RESI CIs on egsingle
#  (lmer and lme).]

# Coverage of the mixed-model RESI confidence intervals in the egsingle hybrid
# plasmode, Model A or B (tests/egsingle_models.R). Each replicate computes
# intervals from the lmer fit and, for Model A, from an nlme::lme refit of the
# same model, with the model-based and CR2
# (school-clustered) covariances and the normal, QF and CF methods; coverage
# is scored against the Lyu-Welsh limit S_Lyu (g, m -> Inf) and the fixed-m
# target S_m (expectedDesignResi). Writes resiLongCalibrationSim/
# mixed_ci_coverage.rds and mixed_ci_engine_agreement.rds. Run from the
# AIgams root:
#   Rscript tests/run_mixed_ci_coverage.R [A|B]
# Model B writes mixed_ci_coverage_B.rds. Options: resi.lme (lme refits;
# default for Model A only), resi.cr2_t (also CR2 with the t(g - 1) reference;
# default TRUE), resi.coverage.tag (appended to file names), and
# resi.n_schools / resi.n_per_site / resi.nsim / resi.mc.cores.

devtools::load_all("RESI")
source("tests/mixed_ci_coverage_lib.R")
source("tests/egsingle_models.R")
model <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(model)) model <- "A"
suffix <- paste0(if (model == "A") "" else paste0("_", model),
                 getOption("resi.coverage.tag", ""))
use_lme <- getOption("resi.lme", model == "A")
form <- egsingle_models[[model]]

out <- hybridPlasmodeSim(
  form, egsingle, school = "schoolid", child = "childid",
  n_schools = getOption("resi.n_schools", c(15, 100)),
  n_per_site = getOption("resi.n_per_site", c(25, 200)),
  nsim = getOption("resi.nsim", 1000L),
  mc.cores = getOption("resi.mc.cores", 24L),
  seed = 2027, keep_reps = TRUE,
  rep_fun = make_ci_rep(form, "schoolid", "childid", lme = use_lme,
                        cr2_t = getOption("resi.cr2_t", TRUE)))

out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
long <- coverage_long(out, "schoolid", "childid:schoolid")
coverage <- summarize_coverage(long)
if (getOption("resi.save.reps", FALSE)) saveRDS(do.call(rbind, long), file.path(out_dir, paste0("mixed_ci_coverage", suffix, "_reps.rds")))
saveRDS(coverage, file.path(out_dir, paste0("mixed_ci_coverage", suffix, ".rds")))
print_coverage(coverage)
if (use_lme) {
  agreement <- summarize_engine_agreement(out)
  saveRDS(agreement, file.path(out_dir, paste0("mixed_ci_engine_agreement", suffix, ".rds")))
  cat("\n=== Max |lme - lmer| over replicates ===\n")
  print(agreement, digits = 3)
}
