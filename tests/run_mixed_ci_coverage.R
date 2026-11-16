# Coverage of the mixed-model RESI confidence intervals in the egsingle hybrid
# plasmode (Model A). Each replicate computes intervals from the lmer fit and
# from an nlme::lme refit of the same model, with the model-based and CR2
# (school-clustered) covariances and the normal, QF and CF methods; coverage
# is scored against the Lyu-Welsh limit S_Lyu (g, m -> Inf) and the fixed-m
# target S_m (expectedDesignResi). Writes resiLongCalibrationSim/
# mixed_ci_coverage.rds and mixed_ci_engine_agreement.rds. Run from the
# AIgams root:
#   Rscript tests/run_mixed_ci_coverage.R

devtools::load_all("RESI")
source("tests/mixed_ci_coverage_lib.R")

egsingle <- local({
  e <- new.env()
  utils::data("egsingle", package = "mlmRev", envir = e)
  e$egsingle
})
egsingle$size_100 <- egsingle$size / 100
egsingle$lowinc_10 <- egsingle$lowinc / 10

form_A <- math ~ year + female + black + hispanic + size_100 + lowinc_10 +
  (1 | schoolid/childid)

out <- hybridPlasmodeSim(
  form_A, egsingle, school = "schoolid", child = "childid",
  n_schools = getOption("resi.n_schools", c(15, 100)),
  n_per_site = getOption("resi.n_per_site", c(25, 200)),
  nsim = getOption("resi.nsim", 1000L),
  mc.cores = getOption("resi.mc.cores", 24L),
  seed = 2027, keep_reps = TRUE,
  rep_fun = make_ci_rep(form_A, "schoolid", "childid", lme = TRUE))

coverage <- summarize_coverage(out, "schoolid", "childid:schoolid")
agreement <- summarize_engine_agreement(out)
out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
saveRDS(coverage, file.path(out_dir, "mixed_ci_coverage.rds"))
saveRDS(agreement, file.path(out_dir, "mixed_ci_engine_agreement.rds"))

print_coverage(coverage)
cat("\n=== Max |lme - lmer| over replicates ===\n")
print(agreement, digits = 3)
