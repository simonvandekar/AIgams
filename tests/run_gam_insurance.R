# [AI-generated: GitHub Copilot, 2026-11-18]
# [Prompt summary: resampling (pairs bootstrap) evaluation of resi_pe() for
#  mgcv gam fits on RESI::insurance, including smoker x nonlinear age (and
#  bmi) interactions; follows insurancePlasmodeSim() for lm/glm.]

# The population is the insurance data (N = 1338). Each replicate draws n rows
# with replacement and adds a pure-noise covariate z ~ U(0, 1) (a null smooth).
# The model, with cubic regression splines whose knots are fixed from the full
# data (so the basis is identical in every replicate), is
#   log10(charges) ~ smoker + sex + region + children + s(age) + s(bmi)
#                    + s(age, by = smoker) + s(bmi, by = smoker) + s(z)
# where `smoker` is an ordered factor, so the by-smooths are the smoker minus
# non-smoker difference curves (the nonlinear interactions).
#
# Target: the full-data projection onto the same basis (unpenalized fit on
# N = 1338), i.e. the regression-spline target T2(k) of notes/resi_gam.qmd.
# For a smooth or factor term with coefficients b_j,
#   S_j^2 = b_j' [(X'X)^{-1}]_jj^{-1} b_j / (N sigma^2),  sigma^2 = mean(resid^2),
# and for a coefficient S = b / sqrt(N sigma^2 [(X'X)^{-1}]_jj) (the lm
# convention of insurancePlasmodeSim()). Under pairs resampling the
# population's basis-free regression function is the empirical one (nearly
# saturated in bmi), so T1 is not a useful target; k in {10, 20} shows how
# much the target depends on the basis. S(z) = 0 exactly.
#
# Estimators: resi_pe() on gam(method = "REML") fits, penalized (mgcv) and
# unpenalized (fx = TRUE). Writes gamCalibrationSim/gam_insurance.rds.
# Run from the AIgams root:  Rscript tests/run_gam_insurance.R
# Options: resi.gam.n, resi.gam.k, resi.gam.nsim, resi.mc.cores.

suppressMessages({
  devtools::load_all("RESI", quiet = TRUE)
  library(mgcv)
  library(parallel)
})

pop <- RESI::insurance
pop$smoker <- ordered(pop$smoker, levels = c(0, 1), labels = c("no", "yes"))
pop$sex <- factor(pop$sex)
pop$region <- factor(pop$region)
N <- nrow(pop)

gam_formula <- function(k, fx = FALSE, null_term = TRUE) {
  sm <- function(v, by = NULL) sprintf("s(%s%s, bs = \"cr\", k = %d, fx = %s)", v,
                                       if (is.null(by)) "" else paste0(", by = ", by), k, fx)
  rhs <- c("smoker + sex + region + children", sm("age"), sm("bmi"),
           sm("age", "smoker"), sm("bmi", "smoker"), if (null_term) sm("z"))
  as.formula(paste("log10(charges) ~", paste(rhs, collapse = " + ")))
}

fixed_knots <- function(k) {
  list(age = place.knots(unique(pop$age), k), bmi = place.knots(unique(pop$bmi), k),
       z = seq(0, 1, length.out = k))
}

# Full-data targets for basis size k.
truths <- function(k) {
  b <- gam(gam_formula(k, fx = TRUE, null_term = FALSE), data = pop, knots = fixed_knots(k))
  X <- model.matrix(b)
  beta <- coef(b)
  sigma2 <- mean(residuals(b)^2)
  XtXi <- solve(crossprod(X))
  quad <- function(j) {
    sum(beta[j] * solve(XtXi[j, j, drop = FALSE], beta[j])) / (N * sigma2)
  }
  sm <- vapply(b$smooth, function(s) quad(s$first.para:s$last.para), numeric(1))
  names(sm) <- vapply(b$smooth, `[[`, character(1), "label")
  pterms <- attr(b$pterms, "term.labels")
  pa <- vapply(seq_along(pterms), function(i) quad(which(b$assign == i)), numeric(1))
  names(pa) <- pterms
  np <- length(b$assign)
  coefs <- beta[seq_len(np)] / sqrt(N * sigma2 * diag(XtXi)[seq_len(np)])
  z_lab <- sprintf("s(z)")
  rbind(data.frame(table = "anova", term = c(names(pa), names(sm), z_lab), k = k,
                   S_true = sqrt(c(pa, sm, 0))),
        data.frame(table = "coefficients", term = names(coefs), k = k, S_true = coefs))
}

fit_stats <- function(d, k, fx) {
  b <- gam(gam_formula(k, fx), data = d, method = "REML", knots = fixed_knots(k))
  r <- resi_pe(b)
  an <- r$anova
  ct <- r$coefficients
  rbind(data.frame(table = "anova", term = rownames(an), edf = an$edf, df = an$Df,
                   stat = an[, 3], p = an[, 4], S = an$RESI, est = NA_real_, se = NA_real_),
        data.frame(table = "coefficients", term = rownames(ct), edf = NA_real_, df = 1,
                   stat = ct[, 3], p = ct[, 4], S = ct$RESI, est = ct[, 1], se = ct[, 2]))
}

one_rep <- function(r, n, ks, seed) {
  set.seed(seed + 1e5 * match(n, n_grid) + r)
  d <- pop[sample(N, n, replace = TRUE), ]
  d$z <- runif(n)
  out <- list()
  for (k in ks) for (fx in c(FALSE, TRUE)) {
    res <- tryCatch(fit_stats(d, k, fx), error = function(e) NULL)
    if (is.null(res)) next
    out[[length(out) + 1]] <- cbind(data.frame(n = n, rep = r, k = k,
                                               estimator = if (fx) "unpenalized" else "mgcv"), res)
  }
  do.call(rbind, out)
}

n_grid <- getOption("resi.gam.n", c(200, 500, 1000, 2000, 5000))
k_grid <- getOption("resi.gam.k", c(10, 20))
nsim <- getOption("resi.gam.nsim", 500L)
cores <- getOption("resi.mc.cores", 4L)

truth <- do.call(rbind, lapply(k_grid, truths))
RNGkind("L'Ecuyer-CMRG")
tab <- do.call(rbind, lapply(n_grid, function(n) {
  t0 <- Sys.time()
  x <- mclapply(seq_len(nsim), one_rep, n = n, ks = k_grid, seed = 2029, mc.cores = cores)
  bad <- vapply(x, function(v) inherits(v, "try-error") || is.null(v), logical(1))
  message(sprintf("n = %d: %d replicates (%d failed), %.1f min", n, sum(!bad), sum(bad),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  do.call(rbind, x[!bad])
}))
tab <- merge(tab, truth, by = c("table", "term", "k"), all.x = TRUE)

out_dir <- "gamCalibrationSim"
dir.create(out_dir, showWarnings = FALSE)
saveRDS(list(tab = tab, truth = truth,
             settings = list(n = n_grid, k = k_grid, nsim = nsim, seed = 2029)),
        file.path(out_dir, "gam_insurance.rds"))
