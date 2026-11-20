# [AI-generated: GitHub Copilot, 2026-11-17; revised 2026-11-20]
# [Prompt summary: Goal 3 null check + point-estimate calibration for RESI of
#  mgcv gam/bam smooth and parametric terms against the fixed-basis target
#  T2(k), with parametric (Vp) and sandwich covariances, in realistic
#  scenarios (notes/resi_gam.qmd).]

# Gaussian simulation for RESI point estimates from mgcv fits. Each replicate
# draws one data set per scenario and n and fits
#   y ~ s(x1, k) + s(x2, k) + s(x3, k) + w + w0
# with gam (REML) and, in the base scenario, bam (fREML), penalized (mgcv
# default) and unpenalized (fx = TRUE). Every fit is scored twice: with mgcv's
# parametric covariance Vp (summary.gam default) and with mgcv's sandwich
# covariance (vcov.gam(sandwich = TRUE)). Scenarios (x1, x2 correlated
# through a Gaussian copula with rho = 0.5; x3, w, w0 independent):
#   base        sin(2 pi x1) + 0.5 x2 + 0.2 w, sigma = 1
#   step        0.8 * 1(x1 > 0.5) replaces sin: truncation error, T2 depends on k
#   weak_curve  0.25 (2 x2 - 1)^2 replaces 0.5 x2: weak non-linear effect,
#               REML smoothing parameter often >> n at small n
#   hetero      base mean, sd = exp(0.4 w + 0.6 (x1 - 0.5)): parametric and
#               robust targets differ
# s(x3) and w0 are null in every scenario.
# Target T2(k): the Wald noncentrality per observation of the fixed-k
# regression-spline (projection) fit, S_j^2 = b_j' [V_jj]^{-1} b_j, with
# V = sigma_k^2 A^{-1} (parametric, sigma_k^2 = E[working residual^2]) or
# V = A^{-1} E[x x' e^2] A^{-1} (robust), A = E[x x'].
# Writes gamCalibrationSim/gam_calibration.rds. Run from the AIgams root:
#   Rscript tests/run_gam_calibration.R
# Options: resi.gam.n, resi.gam.k, resi.gam.nsim, resi.mc.cores.

suppressMessages({
  devtools::load_all("RESI", quiet = TRUE)
  library(mgcv)
  library(parallel)
})

base <- list(rho = 0.5, beta_w = 0.2, engines = c("gam", "bam"),
             f = list(x1 = function(x) sin(2 * pi * x), x2 = function(x) 0.5 * x,
                      x3 = function(x) 0 * x),
             sd = function(d) rep(1, nrow(d)))
scenarios <- list(
  base = base,
  step = modifyList(base, list(engines = "gam", f = list(x1 = function(x) 0.8 * (x > 0.5)))),
  weak_curve = modifyList(base, list(engines = "gam", f = list(x2 = function(x) 0.25 * (2 * x - 1)^2))),
  hetero = modifyList(base, list(engines = "gam", sd = function(d) exp(0.4 * d$w + 0.6 * (d$x1 - 0.5))))
)
smooths <- c("x1", "x2", "x3")
params <- c("w", "w0")

gen_data <- function(n, dgp) {
  z1 <- rnorm(n)
  z2 <- dgp$rho * z1 + sqrt(1 - dgp$rho^2) * rnorm(n)
  d <- data.frame(x1 = pnorm(z1), x2 = pnorm(z2), x3 = runif(n), w = rnorm(n), w0 = rnorm(n))
  d$mu <- dgp$f$x1(d$x1) + dgp$f$x2(d$x2) + dgp$f$x3(d$x3) + dgp$beta_w * d$w
  d$y <- d$mu + dgp$sd(d) * rnorm(n)
  d
}

gam_formula <- function(k, fx = FALSE) {
  as.formula(sprintf("y ~ s(x1, k = %d, fx = %s) + s(x2, k = %d, fx = %s) + s(x3, k = %d, fx = %s) + w + w0",
                     k, fx, k, fx, k, fx))
}

# T2(k) from a large noiseless sample: projection of the mean on the default
# (tp) basis of size k, with E[e^2 | x] = sd(x)^2 + (mean - projection)^2.
truths <- function(sc, k, N = 2e5) {
  set.seed(1)
  dgp <- scenarios[[sc]]
  d <- gen_data(N, dgp)
  b0 <- gam(gam_formula(k, fx = TRUE), data = d, fit = FALSE)
  X <- b0$X
  fit <- stats::lm.fit(X, d$mu)
  e2 <- dgp$sd(d)^2 + fit$residuals^2
  Ai <- solve(crossprod(X) / N)
  V <- list(par = mean(e2) * Ai, rob = Ai %*% (crossprod(X * sqrt(e2)) / N) %*% Ai)
  cols <- c(lapply(b0$smooth, function(s) s$first.para:s$last.para), as.list(match(params, b0$term.names)))
  S <- sapply(V, function(v) sapply(cols, function(j) {
    bj <- fit$coefficients[j]
    sqrt(drop(bj %*% solve(v[j, j, drop = FALSE], bj))) * if (length(j) == 1) sign(bj) else 1
  }))
  data.frame(scenario = sc, term = c(smooths, params), k = k, S_T2 = S[, "par"], S_T2r = S[, "rob"],
             sigma_k = sqrt(mean(e2)), row.names = NULL)
}

grid_x <- seq(0.02, 0.98, length.out = 20)

# RESI estimates and statistics from summary(), for the covariance in b$Vp.
term_stats <- function(b, n) {
  s <- summary(b)
  rdf <- b$df.residual
  st <- s$s.table
  sm <- data.frame(term = smooths, F = st[, "F"], ref_df = st[, "Ref.df"], edf = st[, "edf"],
                   p = st[, "p-value"], S = f2S(st[, "F"], st[, "Ref.df"], rdf, n),
                   est = NA_real_, se = NA_real_)
  pt <- s$p.table[params, , drop = FALSE]
  pa <- data.frame(term = params, F = NA_real_, ref_df = 1, edf = 1, p = pt[, "Pr(>|t|)"],
                   S = t2S(pt[, "t value"], rdf, n), est = pt[, "Estimate"], se = pt[, "Std. Error"])
  rbind(sm, pa)
}

# One fit: RESI with Vp and with the sandwich, and pointwise term fits (Vp).
fit_stats <- function(d, k, engine, fx) {
  f <- gam_formula(k, fx)
  b <- if (engine == "gam") gam(f, data = d, method = "REML") else bam(f, data = d)
  n <- nrow(d)
  bs <- b
  bs$Vp <- vcov(b, sandwich = TRUE)
  tab <- rbind(cbind(cov = "Vp", term_stats(b, n)), cbind(cov = "sandwich", term_stats(bs, n)))
  nd <- data.frame(x1 = grid_x, x2 = grid_x, x3 = grid_x, w = 0, w0 = 0)
  pr <- predict(b, nd, type = "terms", se.fit = TRUE)
  curves <- do.call(rbind, lapply(smooths, function(v) {
    j <- sprintf("s(%s)", v)
    data.frame(term = v, x = grid_x, fit = pr$fit[, j], se = pr$se.fit[, j])
  }))
  list(tab = cbind(tab, row.names = NULL), curves = curves)
}

one_rep <- function(r, n, ks, seed) {
  out <- list()
  for (sc in names(scenarios)) {
    set.seed(seed + 1e5 * match(n, n_grid) + r)
    d <- gen_data(n, scenarios[[sc]])
    for (k in ks) for (engine in scenarios[[sc]]$engines) for (fx in c(FALSE, TRUE)) {
      res <- tryCatch(fit_stats(d, k, engine, fx), error = function(e) NULL)
      if (is.null(res)) next
      lab <- data.frame(scenario = sc, n = n, rep = r, k = k, engine = engine,
                        estimator = if (fx) "unpenalized" else "mgcv")
      out$tab[[length(out$tab) + 1]] <- cbind(lab, res$tab)
      out$curves[[length(out$curves) + 1]] <- cbind(lab, res$curves)
    }
  }
  list(tab = do.call(rbind, out$tab), curves = do.call(rbind, out$curves))
}

n_grid <- getOption("resi.gam.n", c(200, 1000, 5000, 20000))
k_grid <- getOption("resi.gam.k", c(5, 10, 20))
nsim <- getOption("resi.gam.nsim", 500L)
cores <- getOption("resi.mc.cores", 4L)

truth <- do.call(rbind, lapply(names(scenarios), function(sc) do.call(rbind, lapply(k_grid, truths, sc = sc))))
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
tab <- merge(do.call(rbind, lapply(reps, `[[`, "tab")), truth, by = c("scenario", "term", "k"))
# Pointwise curves are only summarized at k = 10 (keeps the .rds small).
curves <- do.call(rbind, lapply(reps, `[[`, "curves"))
curves <- curves[curves$k == 10, ]
# True smooth on grid_x, centred at its population mean (mgcv's sum-to-zero
# constraint; x1, x2 have uniform margins).
curves$truth <- mapply(function(sc, v, x) {
  f <- scenarios[[sc]]$f[[v]]
  f(x) - integrate(f, 0, 1, subdivisions = 1000)$value
}, curves$scenario, curves$term, curves$x)

out_dir <- "gamCalibrationSim"
dir.create(out_dir, showWarnings = FALSE)
saveRDS(list(tab = tab, curves = curves, truth = truth,
             settings = list(n = n_grid, k = k_grid, nsim = nsim, seed = 2028,
                             scenarios = names(scenarios))),
        file.path(out_dir, "gam_calibration.rds"), compress = "xz")
