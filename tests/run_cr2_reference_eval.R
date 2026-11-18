# [AI-generated: GitHub Copilot, 2026-11-17]
# [Prompt summary: compare normal/chi-square vs t/F(g - 1) references for the
#  CR2 mixed-model RESI CIs.]

# Reference distribution for the CR2 (cluster-robust) mixed-model RESI
# intervals: normal / chi-square vs t(g - 1) / F(., g - 1), with g the number of
# clusters. Reruns the g = 15 coverage cells (lmer only) on egsingle (Model A)
# and mlmRev::star, keeping every replicate's intervals. Writes
# resiLongCalibrationSim/cr2_reference_eval_{egsingle,star}.rds (per-replicate
# rows from coverage_long()). Run from the AIgams root:
#   Rscript tests/run_cr2_reference_eval.R

devtools::load_all("RESI")
source("tests/mixed_ci_coverage_lib.R")
source("tests/star_data.R")

out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
nsim <- getOption("resi.nsim", 1000L)
mc_cores <- getOption("resi.mc.cores", 24L)
n_schools <- getOption("resi.n_schools", 15)
n_per_site <- getOption("resi.n_per_site", c(25, 200))

source("tests/egsingle_models.R")
form_A <- egsingle_models$A
star <- load_star()

runs <- list(
  egsingle = list(formula = form_A, data = egsingle, school = "schoolid", child = "childid"),
  star = list(formula = form_star, data = star, school = "sch", child = "sid"))

for (nm in names(runs)) {
  r <- runs[[nm]]
  out <- hybridPlasmodeSim(
    r$formula, r$data, school = r$school, child = r$child,
    n_schools = n_schools, n_per_site = n_per_site, nsim = nsim,
    mc.cores = mc_cores, seed = 2027, keep_reps = TRUE,
    rep_fun = make_ci_rep(r$formula, r$school, r$child, cr2_t = TRUE))
  long <- do.call(rbind, coverage_long(out, r$school, paste0(r$child, ":", r$school)))
  saveRDS(long, file.path(out_dir, paste0("cr2_reference_eval_", nm, ".rds")))
  cov <- summarize_coverage(split(long, long$n_per_site))
  cat("\n=====", nm, ": coverage of S_m, nominal 0.95 =====\n")
  cov <- subset(cov, covariance != "model")
  print(tapply(cov$coverage_S_m, list(paste(cov$table, cov$term),
                                      paste(cov$covariance, cov$method, cov$n_per_site)), mean),
        digits = 3)
}
