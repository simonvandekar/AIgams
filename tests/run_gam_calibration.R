# [AI-generated: GitHub Copilot, 2026-11-17]
# [Prompt summary: Goal 3 null check + point-estimate calibration for RESI of
#  mgcv gam/bam smooth and parametric terms (notes/resi_gam.qmd).]

# Gaussian simulation for RESI point estimates from mgcv fits. One data set per
# replicate and n is fitted with gam (REML) and bam (default fREML) at basis
# sizes k, penalized (mgcv default) and unpenalized (fx = TRUE, reference
# estimator). Terms:
#   s(x1): sin(2 pi x1)            non-linear
#   s(x2): 0.5 x2                  linear effect fitted with s(); x1, x2
#                                  correlated (Gaussian copula, rho) for concurvity
#   s(x3): 0                       null smooth
#   w: 0.2 w, w0: 0                parametric effect / null
# Targets (notes/resi_gam.qmd): T1, the basis-free S_j^2 =
# E[(f_j - projection on the other terms)^2] / sigma^2, and T2, the same for
# the fixed-k regression-spline limit. Smooth RESI = f2S(F, Ref.df, rdf, n);
# coefficient RESI = t2S(t, rdf, n). Writes gamCalibrationSim/
# gam_calibration.rds. Run from the AIgams root:
#   Rscript tests/run_gam_calibration.R
# Options: resi.gam.n, resi.gam.k, resi.gam.nsim, resi.mc.cores.

suppressMessages({
  devtools::load_all("RESI", quiet = TRUE)
  library(mgcv)
  library(parallel)
})

dgp <- list(rho = 0.5, sigma = 1, beta_w = 0.2,
            f = list(x1 = function(x) sin(2 * pi * x), x2 = function(x) 0.5 * x,
                     x3 = function(x) 0 * x))
smooths <- c("x1", "x2", "x3")
params <- c("w", "w0")

gen_data <- function(n) {
  z1 <- rnorm(n)
  z2 <- dgp$rho * z1 + sqrt(1 - dgp$rho^2) * rnorm(n)
  d <- data.frame(x1 = pnorm(z1), x2 = pnorm(z2), x3 = runif(n), w = rnorm(n), w0 = rnorm(n))
  d$mu <- dgp$f$x1(d$x1) + dgp$f$x2(d$x2) + dgp$f$x3(d$x3) + dgp$beta_w * d$w
  d$y <- d$mu + dgp$sigma * rnorm(n)
  d
}

gam_formula <- function(k, fx = FALSE) {
  as.formula(sprintf("y ~ s(x1, k = %d, fx = %s) + s(x2, k = %d, fx = %s) + s(x3, k = %d, fx = %s) + w + w0",
                     k, fx, k, fx, k, fx))
}

# Unpenalized basis (constraint absorbed) of one smooth evaluated on `d`.
basis <- function(v, k, d, bs = "tp") {
  smoothCon(s(x, k = k, bs = bs), data = data.frame(x = d[[v]]), absorb.cons = TRUE)[[1]]$X
}

# Population targets from a large noiseless sample. T1 projects the true f_j
# on a rich (k = 40, cubic regression) basis of the other terms; T2 uses the
# noiseless fixed-k (tp) regression-spline fit and its own bases.
truths <- function(k_fit, N = 2e5, k_big = 40) {
  set.seed(1)
  d <- gen_data(N)
  others <- function(v, k, bs) do.call(cbind, c(list(1), lapply(setdiff(smooths, v), basis, k = k, d = d, bs = bs),
                                               list(d$w, d$w0)))
  resid_ss <- function(y, X) mean(stats::lm.fit(X, y)$residuals^2) / dgp$sigma^2
  T1 <- sapply(smooths, function(v) resid_ss(dgp$f[[v]](d[[v]]), others(v, k_big, "cr")))
  Xs <- do.call(cbind, c(list(1), lapply(smooths, basis, k = k_fit, d = d), list(d$w, d$w0)))
  b <- stats::lm.fit(Xs, d$mu)$coefficients
  cols <- split(seq_len(ncol(Xs))[-1], rep(c(smooths, params), c(rep(k_fit - 1, 3), 1, 1)))
  T2 <- sapply(smooths, function(v) resid_ss(Xs[, cols[[v]]] %*% b[cols[[v]]], Xs[, -cols[[v]]]))
  Xbig <- others(NA, k_big, "cr")
  sd_w <- function(p) sqrt(mean(stats::lm.fit(Xbig[, colnames(Xbig) != p, drop = FALSE] , d[[p]])$residuals^2))
  colnames(Xbig) <- c(rep("", ncol(Xbig) - 2), params)
  S_param <- c(w = dgp$beta_w, w0 = 0) * sapply(params, sd_w) / dgp$sigma
  rbind(data.frame(term = smooths, k = k_fit, S_T1 = sqrt(T1), S_T2 = sqrt(T2)),
        data.frame(term = params, k = k_fit, S_T1 = S_param, S_T2 = S_param))
}

grid_x <- seq(0.02, 0.98, length.out = 20)

# One fit: RESI estimates, statistics and pointwise term fits on grid_x.
fit_stats <- function(d, k, engine, fx) {
  f <- gam_formula(k, fx)
  b <- if (engine == "gam") gam(f, data = d, method = "REML") else bam(f, data = d)
  s <- summary(b)
  n <- nrow(d)
  rdf <- b$df.residual
  st <- s$s.table
  sm <- data.frame(term = smooths, stat = s$chi.sq, F = st[, "F"], ref_df = st[, "Ref.df"],
                   edf = st[, "edf"], p = st[, "p-value"],
                   S = f2S(st[, "F"], st[, "Ref.df"], rdf, n), est = NA_real_, se = NA_real_)
  pt <- s$p.table[params, , drop = FALSE]
  pa <- data.frame(term = params, stat = pt[, "t value"], F = NA_real_, ref_df = 1, edf = 1,
                   p = pt[, "Pr(>|t|)"], S = t2S(pt[, "t value"], rdf, n),
                   est = pt[, "Estimate"], se = pt[, "Std. Error"])
  nd <- data.frame(x1 = grid_x, x2 = grid_x, x3 = grid_x, w = 0, w0 = 0)
  pr <- predict(b, nd, type = "terms", se.fit = TRUE)
  curves <- do.call(rbind, lapply(smooths, function(v) {
    j <- sprintf("s(%s)", v)
    data.frame(term = v, x = grid_x, fit = pr$fit[, j], se = pr$se.fit[, j])
  }))
  list(tab = cbind(rbind(sm, pa), row.names = NULL), curves = curves)
}

one_rep <- function(r, n, ks, seed) {
  set.seed(seed + 1e5 * match(n, n_grid) + r)
  d <- gen_data(n)
  out <- list()
  for (k in ks) for (engine in c("gam", "bam")) for (fx in c(FALSE, TRUE)) {
    res <- tryCatch(fit_stats(d, k, engine, fx), error = function(e) NULL)
    if (is.null(res)) next
    lab <- data.frame(n = n, rep = r, k = k, engine = engine,
                      estimator = if (fx) "unpenalized" else "mgcv")
    out$tab[[length(out$tab) + 1]] <- cbind(lab, res$tab)
    out$curves[[length(out$curves) + 1]] <- cbind(lab, res$curves)
  }
  list(tab = do.call(rbind, out$tab), curves = do.call(rbind, out$curves))
}

n_grid <- getOption("resi.gam.n", c(200, 1000, 5000, 20000))
k_grid <- getOption("resi.gam.k", c(10, 20))
nsim <- getOption("resi.gam.nsim", 500L)
cores <- getOption("resi.mc.cores", 4L)

truth <- do.call(rbind, lapply(k_grid, truths))
RNGkind("L'Ecuyer-CMRG")
reps <- lapply(n_grid, function(n) {
  t0 <- Sys.time()
  x <- mclapply(seq_len(nsim), one_rep, n = n, ks = k_grid, seed = 2028, mc.cores = cores)
  bad <- vapply(x, inherits, logical(1), "try-error")
  message(sprintf("n = %d: %d replicates (%d failed), %.1f min", n, sum(!bad), sum(bad),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  x <- x[!bad]
  list(tab = do.call(rbind, lapply(x, `[[`, "tab")), curves = do.call(rbind, lapply(x, `[[`, "curves")))
})
tab <- merge(do.call(rbind, lapply(reps, `[[`, "tab")), truth, by = c("term", "k"))
curves <- do.call(rbind, lapply(reps, `[[`, "curves"))
# Truth of each smooth on grid_x, centred at its population mean (mgcv's
# sum-to-zero constraint; x1 and x2 have uniform margins).
f_mean <- sapply(smooths, function(v) integrate(dgp$f[[v]], 0, 1)$value)
curves$truth <- mapply(function(v, x) dgp$f[[v]](x) - f_mean[[v]], curves$term, curves$x)

out_dir <- "gamCalibrationSim"
dir.create(out_dir, showWarnings = FALSE)
saveRDS(list(tab = tab, curves = curves, truth = truth, dgp = dgp,
             settings = list(n = n_grid, k = k_grid, nsim = nsim, seed = 2028)),
        file.path(out_dir, "gam_calibration.rds"))
