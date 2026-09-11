# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-11]
# [Prompt summary: figure for the n_schools x n_per_site diagnostic
#  (test_egsingle_synthetic_fixed_n.R), for embedding in calibration_report.Rmd,
#  following the same base-graphics/PDF conventions as
#  RESI::longitudinalCalibrationFigures().]

summary_table <- readRDS("tests/synthetic_fixed_n_summary.rds")

out_dir <- "resiLongCalibrationSim/synthetic_diagnostic/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_path <- file.path(out_dir, "longitudinal_calibration_fixedN.pdf")

.sim_pal <- c(
  "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
  "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F"
)

fe_terms <- unique(summary_table$param[grepl("^beta_", summary_table$param)])
vc_terms <- unique(summary_table$param[!grepl("^beta_", summary_table$param)])
fe_cols  <- setNames(.sim_pal[seq_along(fe_terms)], fe_terms)

n_per_site.vec <- sort(unique(summary_table$n_per_site))
n_schools.vec  <- sort(unique(summary_table$n_schools))
np_cols <- setNames(.sim_pal[seq_along(n_per_site.vec)], n_per_site.vec)

grDevices::pdf(fig_path, width = 15, height = 8)
graphics::layout(
  matrix(c(1, 2, 3, 4,
           5, 6, 7, 8), nrow = 2L, byrow = TRUE),
  widths  = c(3.4, 3.4, 3.4, 1.9),
  heights = c(1, 1)
)

# --- Row 1: fixed-effect SE calibration vs n_per_site, one panel per n_schools ---
for (ns in n_schools.vec) {
  sub <- summary_table[summary_table$n_schools == ns & summary_table$param %in% fe_terms, ]
  graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
  graphics::plot(NULL, xlim = range(n_per_site.vec), ylim = range(c(sub$se_ratio, 1), na.rm = TRUE),
                 xlab = "n_per_site", ylab = "mean_model_se / mc_sd",
                 main = paste0("n_schools = ", ns), log = "x", xaxt = "n", bty = "l")
  graphics::axis(1, at = n_per_site.vec, labels = n_per_site.vec, mgp = c(1.7, 0.35, 0))
  graphics::abline(h = 1, lty = 2L, col = "gray40")
  for (term in fe_terms) {
    td <- sub[sub$param == term, ]
    td <- td[order(td$n_per_site), ]
    graphics::lines(td$n_per_site, td$se_ratio, col = fe_cols[[term]], lwd = 2L)
    graphics::points(td$n_per_site, td$se_ratio, col = fe_cols[[term]], pch = 16L, cex = 0.8)
  }
}

graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = fe_terms, col = fe_cols[fe_terms],
                  lwd = 2L, pch = 16L, bty = "n", cex = 0.85, title = "Fixed effect")

# --- Row 2: variance-component relative bias vs n_schools, one panel per component ---
for (term in vc_terms) {
  sub <- summary_table[summary_table$param == term, ]
  sub$rel_bias_pct <- 100 * sub$bias / sub$truth
  graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
  graphics::plot(NULL, xlim = range(n_schools.vec), ylim = range(c(sub$rel_bias_pct, 0), na.rm = TRUE),
                 xlab = "n_schools", ylab = "Relative bias (%)",
                 main = term, log = "x", xaxt = "n", bty = "l")
  graphics::axis(1, at = n_schools.vec, labels = n_schools.vec, mgp = c(1.7, 0.35, 0))
  graphics::abline(h = 0, lty = 2L, col = "gray40")
  for (np in n_per_site.vec) {
    td <- sub[sub$n_per_site == np, ]
    td <- td[order(td$n_schools), ]
    graphics::lines(td$n_schools, td$rel_bias_pct, col = np_cols[[as.character(np)]], lwd = 2L)
    graphics::points(td$n_schools, td$rel_bias_pct, col = np_cols[[as.character(np)]], pch = 16L, cex = 0.8)
  }
}

graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = n_per_site.vec, col = np_cols[as.character(n_per_site.vec)],
                  lwd = 2L, pch = 16L, bty = "n", cex = 0.85, title = "n_per_site")

grDevices::dev.off()
message("Saved: ", fig_path)
