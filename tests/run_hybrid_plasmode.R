# [AI-generated: GitHub Copilot, 2026-11-15]
# [Prompt summary: run the egsingle hybrid plasmode calibration of mixed-model
#  RESI.]

# Hybrid plasmode calibration of mixed-model RESI on mlmRev::egsingle,
# Model A (random intercepts) or Model B (adds a student-level random slope
# for year). Writes the summaries read by
# resiLongCalibrationSim/calibration_report.Rmd. Run from the AIgams root:
#   Rscript tests/run_hybrid_plasmode.R [A|B]
# Full grid (12 cells x 1000 replicates) takes ~1.5 h on 24 cores.

devtools::load_all("RESI")
source("tests/egsingle_models.R")
model <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(model)) model <- "A"
suffix <- if (model == "A") "" else paste0("_", model)

out <- hybridPlasmodeSim(
  egsingle_models[[model]], egsingle, school = "schoolid", child = "childid",
  n_schools = getOption("resi.n_schools", c(15, 30, 60, 100)),
  n_per_site = getOption("resi.n_per_site", c(25, 100, 200)),
  nsim = getOption("resi.nsim", 1000L),
  mc.cores = getOption("resi.mc.cores", 24L),
  seed = 2026)

out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
saveRDS(out$summary, file.path(out_dir, paste0("hybrid_plasmode", suffix, "_summary.rds")))
saveRDS(out$resi_summary, file.path(out_dir, paste0("hybrid_plasmode", suffix, "_resi_summary.rds")))

cat("\n=== Fixed effects: mean model SE / MC SD ===\n")
fe <- out$summary[grepl("^beta_", out$summary$param), ]
print(xtabs(se_ratio ~ param + n_schools + n_per_site, fe), digits = 3)
cat("\n=== Variance components: relative bias (%) ===\n")
vc <- out$summary[!grepl("^beta_", out$summary$param), ]
print(xtabs(100 * bias / truth ~ param + n_schools + n_per_site, vc), digits = 3)
cat("\n=== RESI: Lyu limit vs Monte Carlo mean ===\n")
print(out$resi_summary[, c("n_schools", "n_per_site", "table", "term", "n",
                           "truth_resi", "mc_mean_resi", "mc_sd_resi")],
      row.names = FALSE, digits = 3)
