# Hybrid plasmode calibration of mixed-model RESI on mlmRev::egsingle
# (Model A). Writes the summaries read by
# resiLongCalibrationSim/calibration_report.Rmd. Run from the AIgams root:
#   Rscript tests/run_hybrid_plasmode.R
# Full grid (12 cells x 1000 replicates) takes ~1.5 h on 24 cores.

devtools::load_all("RESI")

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
  n_schools = getOption("resi.n_schools", c(15, 30, 60, 100)),
  n_per_site = getOption("resi.n_per_site", c(25, 100, 200)),
  nsim = getOption("resi.nsim", 1000L),
  mc.cores = getOption("resi.mc.cores", 24L),
  seed = 2026)

out_dir <- getOption("resi.output.dir", "resiLongCalibrationSim")
saveRDS(out$summary, file.path(out_dir, "hybrid_plasmode_summary.rds"))
saveRDS(out$resi_summary, file.path(out_dir, "hybrid_plasmode_resi_summary.rds"))

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
