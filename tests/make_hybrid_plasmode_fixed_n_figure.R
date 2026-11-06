# [AI-generated: Claude Sonnet 5 via GitHub Copilot, 2026-09-24]
# [Prompt summary: figure for the hybrid plasmode (BLUP-resampled random
#  effects) n_schools x n_per_site diagnostic
#  (test_egsingle_hybrid_plasmode_fixed_n.R), for embedding in
#  calibration_report.Rmd, following the same base-graphics/PDF conventions as
#  RESI::longitudinalCalibrationFigures() / make_fixed_n_figure.R.]

summary_table <- readRDS("tests/hybrid_plasmode_resampled_covariates_summary.rds")
resi_summary <- readRDS("tests/hybrid_plasmode_resampled_covariates_resi_summary.rds")

out_dir <- "resiLongCalibrationSim/hybrid_plasmode/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_path <- file.path(out_dir, "longitudinal_calibration_fixedN_hybrid.pdf")
png_path <- file.path(out_dir, "longitudinal_calibration_fixedN_hybrid.png")

.sim_pal <- c(
  "#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
  "#8C564B", "#E377C2", "#17BECF", "#BCBD22", "#7F7F7F"
)

fe_terms <- unique(summary_table$param[grepl("^beta_", summary_table$param)])
vc_terms <- unique(summary_table$param[!grepl("^beta_", summary_table$param)])
resi_coef_terms <- unique(resi_summary$term[resi_summary$table == "coefficients"])
resi_anova_terms <- unique(resi_summary$term[resi_summary$table == "anova"])
fe_cols <- setNames(.sim_pal[seq_along(fe_terms)], fe_terms)
vc_cols <- setNames(.sim_pal[seq_along(vc_terms)], vc_terms)
resi_coef_cols <- setNames(.sim_pal[seq_along(resi_coef_terms)], resi_coef_terms)
resi_anova_cols <- setNames(.sim_pal[seq_along(resi_anova_terms)], resi_anova_terms)

n_per_site.vec <- sort(unique(summary_table$n_per_site))
n_schools.vec  <- sort(unique(summary_table$n_schools))
row_limits <- function(values, reference = 0) {
  values <- c(values, reference)
  limits <- range(values, na.rm = TRUE)
  padding <- max(diff(limits) * 0.08, 0.05)
  limits + c(-padding, padding)
}

fe_limits <- row_limits(summary_table$se_ratio, 1)
vc_summary <- summary_table[summary_table$param %in% vc_terms, ]
vc_bias <- 100 * vc_summary$bias / vc_summary$truth
vc_limits <- row_limits(vc_bias, 0)
se_ratio_values <- summary_table$scaled_mc_error_sd_ratio[
  summary_table$param %in% fe_terms]
se_limits <- row_limits(se_ratio_values, 1)
resi_coef_bias <- resi_summary$bias_resi[resi_summary$table == "coefficients"]
resi_anova_bias <- resi_summary$bias_resi[resi_summary$table == "anova"]
resi_coef_limits <- row_limits(resi_coef_bias, 0)
resi_anova_limits <- row_limits(resi_anova_bias, 0)

grDevices::pdf(fig_path, width = 4 * length(n_schools.vec) + 3, height = 13)
graphics::layout(
  matrix(seq_len(5L * (length(n_schools.vec) + 1L)), nrow = 5L, byrow = TRUE),
  widths  = c(rep(3.4, length(n_schools.vec)), 2.0),
  heights = rep(1, 5)
)

# --- Row 1: fixed-effect SE calibration vs n_per_site, one panel per n_schools ---
for (ns in n_schools.vec) {
  sub <- summary_table[summary_table$n_schools == ns & summary_table$param %in% fe_terms, ]
  graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
  graphics::plot(NULL, xlim = range(n_per_site.vec), ylim = fe_limits,
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

# --- Row 2: all variance components, one panel per n_schools ---
for (ns in n_schools.vec) {
  sub <- summary_table[summary_table$n_schools == ns & summary_table$param %in% vc_terms, ]
  graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
  graphics::plot(NULL, xlim = range(n_per_site.vec), ylim = vc_limits,
                 xlab = "n_per_site", ylab = "Relative bias (%)",
                 main = paste0("n_schools = ", ns), log = "x", xaxt = "n", bty = "l")
  graphics::axis(1, at = n_per_site.vec, labels = n_per_site.vec, mgp = c(1.7, 0.35, 0))
  graphics::abline(h = 0, lty = 2L, col = "gray40")
  for (term in vc_terms) {
    td <- sub[sub$param == term, ]
    td$rel_bias_pct <- 100 * td$bias / td$truth
    td <- td[order(td$n_per_site), ]
    graphics::lines(td$n_per_site, td$rel_bias_pct, col = vc_cols[[term]], lwd = 2L)
    graphics::points(td$n_per_site, td$rel_bias_pct, col = vc_cols[[term]], pch = 16L, cex = 0.8)
  }
}

graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = vc_terms, col = vc_cols[vc_terms],
                  lwd = 2L, pch = 16L, bty = "n", cex = 0.8, title = "Variance component")

# --- Row 3: normalized MC error relative to the fixed full-sample ML SE ---
for (ns in n_schools.vec) {
  sub <- summary_table[summary_table$n_schools == ns & summary_table$param %in% fe_terms, ]
  graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
  graphics::plot(NULL, xlim = range(n_per_site.vec), ylim = se_limits,
                 xlab = "n_per_site", ylab = "Ratio to scaled full-sample SE",
                 main = paste0("n_schools = ", ns), log = "x", xaxt = "n", bty = "l")
  graphics::axis(1, at = n_per_site.vec, labels = n_per_site.vec, mgp = c(1.7, 0.35, 0))
  graphics::abline(h = 1, lty = 1L, lwd = 2L, col = "gray35")
  for (term in fe_terms) {
    td <- sub[sub$param == term, ]
    td <- td[order(td$n_per_site), ]
    graphics::lines(td$n_per_site, td$scaled_mc_error_sd_ratio,
                    col = fe_cols[[term]], lwd = 1.5, lty = 2L)
    graphics::points(td$n_per_site, td$scaled_mc_error_sd_ratio,
                     col = fe_cols[[term]], pch = 1L, cex = 0.75)
  }
}

graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = fe_terms, col = fe_cols[fe_terms],
                  lwd = 1.5, lty = 2L, pch = 1L, bty = "n", cex = 0.75,
                  title = "Coefficient error ratio")
graphics::legend("bottom", legend = c("Full-sample ML SE reference = 1",
                                       "Scaled MC error SD / truth SE"),
                  col = "gray35", lwd = 2L, lty = c(1L, 2L),
                  pch = c(NA, 1L), bty = "n", cex = 0.68)

plot_resi_row <- function(table_name, terms, colors, limits, ylab) {
  for (ns in n_schools.vec) {
    sub <- resi_summary[resi_summary$table == table_name &
                          resi_summary$n_schools == ns &
                          resi_summary$term %in% terms, ]
    graphics::par(mar = c(2.8, 3.2, 2.0, 0.4), mgp = c(1.9, 0.45, 0))
    graphics::plot(NULL, xlim = range(n_per_site.vec), ylim = limits,
                   xlab = "n_per_site", ylab = ylab,
                   main = paste0("n_schools = ", ns), log = "x", xaxt = "n", bty = "l")
    graphics::axis(1, at = n_per_site.vec, labels = n_per_site.vec, mgp = c(1.7, 0.35, 0))
    graphics::abline(h = 0, lty = 2L, col = "gray40")
    for (term in terms) {
      td <- sub[sub$term == term, ][order(sub[sub$term == term, ]$n_per_site), ]
      graphics::lines(td$n_per_site, td$bias_resi, col = colors[[term]], lwd = 2L)
      graphics::points(td$n_per_site, td$bias_resi, col = colors[[term]], pch = 16L, cex = 0.8)
    }
  }
}

plot_resi_row("coefficients", resi_coef_terms, resi_coef_cols, resi_coef_limits, "RESI bias")
graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = resi_coef_terms, col = resi_coef_cols[resi_coef_terms],
                  lwd = 2L, pch = 16L, bty = "n", cex = 0.8, title = "Coefficient")

plot_resi_row("anova", resi_anova_terms, resi_anova_cols, resi_anova_limits, "RESI bias")
graphics::par(mar = c(0.5, 0.3, 0.5, 0.3))
graphics::plot.new()
graphics::legend("center", legend = resi_anova_terms, col = resi_anova_cols[resi_anova_terms],
                  lwd = 2L, pch = 16L, bty = "n", cex = 0.8, title = "ANOVA term")

grDevices::dev.off()
if (nzchar(Sys.which("pdftoppm"))) {
  system2("pdftoppm", c("-png", "-r", "120", "-singlefile", fig_path,
                         sub("\\.png$", "", png_path)))
}
message("Saved: ", fig_path)
