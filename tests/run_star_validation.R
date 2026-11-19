# [AI-generated: GitHub Copilot, 2026-11-16]
# [Prompt summary: out-of-sample calibration and CI coverage of mixed-model
#  RESI on mlmRev::star.]

# Out-of-sample validation of mixed-model RESI on mlmRev::star (Tennessee
# STAR), with the same hybrid plasmode framework as egsingle:
#   1. point-estimate calibration on the 4 x 3 (n_schools x n_per_site) grid
#      (as tests/run_hybrid_plasmode.R), seed 2026;
#   2. CI coverage on the 2 x 2 corner grid (as tests/run_mixed_ci_coverage.R,
#      lmer only), seed 2027.
# Writes resiLongCalibrationSim/star/*.rds, read by
# resiLongCalibrationSim/star/star_validation_report.Rmd. Run from the AIgams
# root:  Rscript tests/run_star_validation.R  [calibration|coverage]

devtools::load_all("RESI")
source("tests/mixed_ci_coverage_lib.R")
source("tests/star_data.R")

star <- load_star()
out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim/star")
steps <- commandArgs(trailingOnly = TRUE)
if (!length(steps)) steps <- c("calibration", "coverage")
mc_cores <- getOption("resi.mc.cores", 24L)
nsim <- getOption("resi.nsim", 1000L)

if ("calibration" %in% steps) {
  cal <- hybridPlasmodeSim(
    form_star, star, school = "sch", child = "sid",
    n_schools = getOption("resi.n_schools", c(15, 30, 60, 100)),
    n_per_site = getOption("resi.n_per_site", c(25, 100, 200)),
    nsim = nsim, mc.cores = mc_cores, seed = 2026)
  saveRDS(cal$summary, file.path(out_dir, "star_plasmode_summary.rds"))
  saveRDS(cal$resi_summary, file.path(out_dir, "star_plasmode_resi_summary.rds"))
  cat("\n=== Fixed effects: mean model SE / MC SD ===\n")
  fe <- cal$summary[grepl("^beta_", cal$summary$param), ]
  print(xtabs(se_ratio ~ param + n_schools + n_per_site, fe), digits = 3)
}

if ("coverage" %in% steps) {
  cov_out <- hybridPlasmodeSim(
    form_star, star, school = "sch", child = "sid",
    n_schools = getOption("resi.cov.n_schools", c(15, 100)),
    n_per_site = getOption("resi.cov.n_per_site", c(25, 200)),
    nsim = nsim, mc.cores = mc_cores, seed = 2027, keep_reps = TRUE,
    rep_fun = make_ci_rep(form_star, "sch", "sid"))
  long <- coverage_long(cov_out, "sch", "sid:sch")
  coverage <- summarize_coverage(long)
  if (getOption("resi.save.reps", FALSE)) saveRDS(do.call(rbind, long), file.path(out_dir, "star_ci_coverage_reps.rds"))
  saveRDS(coverage, file.path(out_dir, "star_ci_coverage.rds"))
  print_coverage(coverage)
}
