# Project Instructions: RESI Development — Mixed & Smooth Model Support

> **How to use this file:** Fill in `[FILL IN]` sections before starting a
> new session. `[FIXED]` sections describe standing project context; update
> only if something major shifts.

---

## Project Context [FIXED]

**Package:** RESI (R package), version 1.5.1.9001 (dev, based on CRAN 1.5.1)
**Repository:** `RESI/` subdirectory of this workspace, cloned from
`https://github.com/statimagcoll/RESI` (branch `master`,
commit `b11560a`, 2026-08-26).
**Key authors:** Megan Jones, Kaidi Kang, Simon Vandekar, Gina Yu, Xinyu Zhang
**Reference papers:**
- Vandekar, Tao, & Blume (2020) <doi:10.1007/s11336-020-09698-2> — original RESI definition
- Jones et al. (2025) JSS software paper <doi:10.18637/jss.v112.i03>
- Kang et al. (2024) — CS-RESI/L-RESI for GEE models (supplement); template for the
  mixed-model extension here
- Zhang et al. (asymptotic RESI variance, under review) — extended influence-function
  derivation used by `resi_asymptotic.R` (`.resi_precompute_ext()` / `.resi_contrast_ext()`)
- Vandekar, Zhang, & Jones, *"A Reproducible Human-AI Workflow for Statistical Software
  Development"* — <https://github.com/simonvandekar/AIpackageUpgrade/blob/main/paper/paper.Rmd>
  — this project follows the validation structure (plasmode simulation +
  `simCalibrationSim` + method-comparison figures) described there.

**AI model used:** Claude Sonnet 4.6 via GitHub Copilot
**Log file:** `log/claudeLog.md` — update after every session

---

## Documentation

### Human-Readable Report
**Format:** This will be a paper following the same format as the one in https://github.com/simonvandekar/AIpackageUpgrade/blob/main/paper/paper.Rmd. I'll abandon that other one besides the framework information, but the actual statistical content will be about the update, here.
**Target audience:** Biostatisticians and Statisticians
**Output file:** A new `paper/paper.Rmd`

### AI-Reproducible Documentation
**Location:** `log/claudeLog.md` (prompt log) + `instructions/instructions.md`

**Code annotation convention:** Each AI-generated function begins with:
`# [AI-generated: <Model>, YYYY-MM-DD]`
`# [Prompt summary: one-line description of what was asked]`

---

## Goals

### Overarching Goal
Extend the RESI package's point-estimator and confidence-interval (CI) machinery —
currently implemented for `lm`/`glm` (bootstrap + asymptotic normal/QF/CF via
`resi_asymptotic.R`) and for GEE (`geeglm`) and mixed models (`lmerMod`, `lme`,
`glmmTMB`, point estimates only) — to cover:

1. `mgcv::gam` and `mgcv::bam` (generalized additive models, including smooth terms)
2. `mgcv::gamm` (GAMM wrapper around `nlme::lme` / `lme4`)
3. `lme4::lmer` — add asymptotic CIs (point estimates already exist via
   `resi_pe.lmerMod` / `resi.lmerMod`, but `resi.lmerMod()` currently only warns
   and returns point estimates)
4. `nlme::lme` — same gap as `lmerMod` (`resi.lme()` also warns "CI procedure not
   developed")

We will follow the same validation structure as the paper: calibration diagnostics (`simCalibrationSim`-style) comparing Monte Carlo
expectations to population/full-data targets, a plasmode simulation
pipeline, and comparison of any new CI against
a trusted reference (bootstrap, or an independent implementation) before
considering it validated.

### Goal 1: Longitudinal/clustered plasmode simulation dataset
**Status:** Done (2026-09-04).
`mlmRev::egsingle` is the primary continuous-outcome plasmode dataset:
60 schools, 1721 children, 7230 rows; 3-level nesting (occasion in child in
school; 5-89 children per school, mean 28.7). `geepack::muscatine` is the
binary-outcome companion; `lme4::InstEval` (crossed) and
`mlmRev::Contraception` (few large sites) are later stress tests. The dataset
surveys are archived in `log/claudeLog.md`.

| Variable | Level | Type | Description |
|---|---|---|---|
| `schoolid` | school (60) | factor | school id |
| `childid` | child (1721), nested in school | factor | child id |
| `year` | occasion | continuous | centred time (-2.5 to 2.5) |
| `math` | occasion | continuous | **outcome** |
| `female`, `black`, `hispanic` | child | binary factor | time invariant |
| `size`, `lowinc`, `mobility` | school | continuous | school enrollment, % low income, % mobility |

Covariates are rescaled by fixed divisors (not z-scores):
`size_100 = size / 100`, `lowinc_10 = lowinc / 10`; fits use
`lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))`.
**Model A** (the only model currently used):
`math ~ year + female + black + hispanic + size_100 + lowinc_10 + (1 | schoolid/childid)`.

### Goal 2: Extend `resi_pe`/`resi_asymptotic` machinery to `lmer`/`lme`
**Status:** Point estimator implemented and calibrated; asymptotic CIs
implemented for `lmerMod` provisionally (derivation note
`notes/resi_mixed_ci.qmd`, not yet reviewed by Simon); `lme` CIs not started.

**Point estimator (`RESI/R/resi_pe_mixed.R`).** `resi_pe.lmerMod` and
`resi_pe.lme` return one `RESI` per coefficient (`z2S`, with `unbiased`) and
per ANOVA term (`chisq2S`), with a **term-specific normalizing count**
(Lyu & Welsh, 2022):
- terms constant within the outermost grouping factor (intercept and
  school-level covariates such as `size_100`, `lowinc_10`) use `n = g`,
  the number of outermost groups (schools);
- all other terms use the number of groups at `level` (outermost = 1; for
  `egsingle`, `level = 2` gives `n = N` students). If `level` is omitted,
  level 1 is used with a message.
- A scalar `n` overrides the count for every term.
Extra columns (`include_level = TRUE`): `Term`, `n`, `Variation level` (the
lowest design level at which the term varies), `n level` (the grouping factor
whose count is `n`), and `Var prop <group>` (share of the design column's
sum of squares between units of each grouping factor, outermost first; the
remainder is within the innermost group). Covariance: `vcovfunc`
(default `clubSandwich::vcovCR`, `type = "CR3"` unless given; the calibration
uses `stats::vcov`). `cluster` sets the robust-variance cluster independently
of `n`. Tests: `RESI/tests/testthat/test-resi-mixed-scale.R`.
**Do not change the scaling factor without asking Simon first.**

**Asymptotic regime and targets.** With g schools and average m students per
school, the RESI depends on the design. Calibration targets the Lyu-Welsh
g, m -> Inf limit evaluated at the full-data ML estimates
(`lyuLimitResi()`); the exact finite-m value of the expected design is
`expectedDesignResi()` (both in `RESI/R/mixed_simulations.R`). Between-school
terms have a finite-m gap (school-level precision is
`g / (tau2_school + sigma2_w / m)`), partly clustered student covariates
(`black`, `hispanic`) a smaller one, and purely within-school terms (`year`,
`female`) essentially none. See `notes/resi_mixed_model_levels.qmd`.

**CIs (`RESI/R/resi_mixed_ci.R`, provisional).** `resi.lmerMod(...,
ci.method = c("normal", "qf", "cf", "none"), ci.design = TRUE)` adds
`2.5%`/`97.5%` columns; `resi.lme()` still returns point estimates only.
Var(R_hat) is estimated from outermost-cluster (school) influence values:
outcome noise through beta_hat, sampling of the per-school information I_s
(random design, `ci.design`), the unit count n_k when it counts units below
the school, the meat U_s U_s' (robust only), and, for the model-based
covariance, the variance-component (psi) uncertainty via a numerical (RE)ML
Hessian. These feed the existing normal/QF/CF algorithms in
`resi_asymptotic.R`. `vcovfunc = stats::vcov` gives the model-based version;
`clubSandwich::vcovCR` the robust one. CR2 for unweighted `lmerMod` fits with
a cluster uses the internal `.mixedCiCR2()` (identical to clubSandwich,
much faster). The four review questions at the end of the note were
answered with the proposed defaults because Simon was unavailable; revisit.

**Files:** `RESI/R/resi_pe_mixed.R`, `RESI/R/resi.R` (`resi.lmerMod`,
`resi.lme`), `RESI/R/resi_asymptotic.R`, `RESI/R/mixed_simulations.R`.

### Goal 3: `resi_pe`/`resi` methods for `gam` and `bam` (mgcv)
**Status:** Not started
**Description:** Add `resi_pe.gam`/`resi.gam` (and `bam` inherits from class
`c("bam","gam","glm","lm")`, so a shared method may be possible) supporting both
a coefficients-style table (parametric terms) and an ANOVA-style table
(parametric + smooth terms, using `mgcv`'s own Wald test for smooths via
`summary.gam`/`anova.gam`).
**Key design decisions:**
- How to define/report a RESI for a smooth term with $\mathrm{edf} > 1$
  degrees of freedom — treat it like an existing multi-df ANOVA term (reuses
  `chisq2S`/QF machinery) using `mgcv`'s approximate smooth-term Wald statistic.
- Which covariance to use: `mgcv`'s Bayesian posterior covariance
  (`vcov(gam_obj, unconditional = FALSE)`, the default), the "sandwich"/frequentist
  correction (`unconditional = TRUE`, i.e. `Vc`), or a true HC/cluster-robust
  sandwich (does `sandwich::vcovHC` support `gam`/`bam` objects? needs checking)? Let's do the default Bayesian approach first.
- Is `bam`'s fast/`discrete = TRUE` fitting path in scope, or only the
  standard `bam` fit? Let's start with default bam fit and revisit once we know that works.
**Expected behavior / acceptance criteria:** [FILL IN]
**Files:**
- Functions: `RESI/R/resi_pe.R`, `RESI/R/resi.R`, `RESI/R/resi_asymptotic.R`
- Output: Similar functionality to other methods in the package. There shouldn't be a CS-RESI here -- even if gam is using lmer to fit, it is not truly a mixed model and L-RESI should just be called RESI if we use lme4/nlme machinery to get the effect sizes and CIs.

### Goal 4: `resi_pe`/`resi` methods for `gamm` (mgcv, PQL via nlme/lme4)
**Status:** Not started
**Description:** `mgcv::gamm()` returns a list with a `$gam` component (for
prediction/summary of smooth terms) and an `$lme` (or, for `gamm4`, an `lmer`)
component (for the mixed-model machinery / random effects). Need a dispatch
strategy that reuses Goal 2 (lme/lmer CI) and Goal 3 (gam smooth terms).
**Key design decisions:**
- [FILL IN] Should `resi()` dispatch on the whole `gamm` list object, or should
  users pass `model.full$gam` / `model.full$lme` explicitly? What's the
  RESI target — the marginal/fixed-effect + smooth structure (`$gam`) or the
  mixed-model longitudinal structure (`$lme`), or both (CS-RESI/L-RESI as in
  Goal 2)? SNV: In this context, we probably want CS-RESI and L-RESI reported and the things we learn in the above goals will help to nail down the details, here.
- Is `gamm4::gamm4` (lme4-based) in scope in addition to `mgcv::gamm`
  (nlme-based)? Yes
**Expected behavior / acceptance criteria:** [PROPOSED — confirm or edit] By
"acceptance criteria" I meant a concrete, checkable definition of "this goal is
done." Proposed for `gamm`:
  1. `resi()`/`resi_pe()` dispatch correctly on `gamm`/`gamm4` objects and return
     both CS-RESI and L-RESI (per your note above) for parametric terms, plus a
     RESI for each smooth term (reusing Goal 3's smooth-term approach on the
     `$gam` component).
  2. Point estimates match a cross-check: e.g. CS-RESI from `gamm`'s `$gam`
     agrees with the plain `gam` fit (ignoring the random effect) on the same
     data within tolerance, and L-RESI agrees with fitting the equivalent
     `lmer`/`lme` directly on the parametric part.
  3. CI coverage for both CS-RESI and L-RESI hits nominal levels in the
     calibration/plasmode simulation from Goal 5, for both `mgcv::gamm` and
     `gamm4::gamm4` fits.
**Files:** [FILL IN]

### Goal 5: Validation (plasmode + calibration, per paper structure)
**5a. Hybrid plasmode calibration (done for the point estimator).**
`RESI::hybridPlasmodeSim()` (`RESI/R/mixed_simulations.R`), run by
`tests/run_hybrid_plasmode.R`, writes
`resiLongCalibrationSim/hybrid_plasmode_{summary,resi_summary}.rds`; the
report is `resiLongCalibrationSim/calibration_report.Rmd`.
Design per replicate:
1. Draw `n_schools` school profiles with replacement from the 60 schools;
   within each, draw students with replacement, `round(n_per_site * m_i / mean(m_i))`
   of them (proportional to observed size m_i). Each student keeps all rows
   and every draw gets a new id.
2. Generate `math = X beta + b_school + b_child + e` from the full-data ML
   Model A fit; random effects are drawn with replacement from the BLUP pools,
   each centred and rescaled to the fitted variance component (1/G
   variance); `e ~ N(0, sigma2)`.
3. Refit Model A by REML; compute `resi_pe(fit, vcovfunc = stats::vcov, level = 2)`.
Grid: `n_schools in {15,30,60,100}` x `n_per_site in {25,100,200}`,
1000 replicates per cell, `RNGkind("L'Ecuyer-CMRG")`, `set.seed(2026)`,
24 cores (~1.5 h).
Reported per cell: successful fits; coefficient and variance-component bias;
`mean_model_se / mc_sd` (target 1); MC SD scaled by `sqrt(n_k)` relative to
the Lyu-Welsh per-unit SE `sqrt(n_k V_kk)`, with design-averaged per-unit SE
ticks (`designAveragedSe()`, the finite-(g, m) reference; cached in
`resiLongCalibrationSim/design_averaged_se.rds`); RESI means relative to the
Lyu limit with fixed-m bands S_m. Figure: between-school terms vs
within-school terms in separate columns.
Why not resample real outcomes: a full two-stage bootstrap of the real data
gave inconsistent SEs (Simon's hypothesis: unmodelled slope heterogeneity);
history archived in `log/claudeLog.md`.

**5b. Trusted reference.** Coverage is evaluated against the known
simulation truth (Lyu limit and S_m), so no bootstrap reference is required
for the mixed-model CIs.

**5c. CI coverage** for the mixed-model CIs (Goal 2):
`tests/run_mixed_ci_coverage.R` (reduced grid `n_schools in {15,100}` x
`n_per_site in {25,200}`, 1000 replicates, seed 2027, ~3 h on 24 cores) writes
`resiLongCalibrationSim/mixed_ci_coverage.rds`; results are in the
"Confidence interval coverage" section of the calibration report. The full
grid has not been run.

<!--
SNV: I commented this because we won't focus on bug fixes for this project.
## Known Issues / GitHub Issues
Open issues on `statimagcoll/RESI` as of 2026-09-04 (none currently target
gam/bam/gamm/lmer/lme CIs directly — this project opens that scope):

| # | Title | Status | Notes |
|---|-------|--------|-------|
| 42 | Add support for rms models | Open | Out of scope for this project unless requested |
| 25 | resi_pe.geeglm CS-RESI replaces weights | Open | Bug; out of scope but worth being aware of since CS-RESI logic is shared/mirrored for `lmerMod` |
| 14 | Add LR test option (?) | Open | Out of scope unless requested |
| 12 | add support for robustbase models | Open | Out of scope unless requested |

--- -->

## Design Decisions Log

Record settled decisions here to avoid relitigating them in future sessions.

| Decision | Rationale | Date |
|----------|-----------|------|
| Cloned `RESI` from `https://github.com/statimagcoll/RESI` master branch (commit `b11560a`) | User-specified source of truth for the package repo | 2026-09-04 |
| Surveyed `lme4`, `nlme`, `geepack` datasets for a >=1000-subject longitudinal plasmode base; `geepack::muscatine` (4856 subjects) is the only candidate meeting the size requirement, but has a binary outcome only | Documented in Goal 1 pending user confirmation | 2026-09-04 |
| Installed `mlmRev`, `JM`, `joineR`, `merTools` locally (not added to `DESCRIPTION` yet) to search for a >=1000-subject continuous-outcome longitudinal dataset | Extended search requested by user | 2026-09-04 |
| Recommend `mlmRev::egsingle` (1721 children, 60 schools, 3-level growth-curve data, continuous `math` outcome) as the primary plasmode dataset; keep `geepack::muscatine` as the binary-outcome companion; keep `lme4::InstEval` (crossed) and `mlmRev::Contraception` (few large sites, cross-sectional) as secondary/later test cases | Resolves the "continuous outcome + >=1000 subjects" gap from the first search; pending final user confirmation and model formulas | 2026-09-04 |
| Goal 5 split into 5a (calibration-consistency of coefficients/variance components — no RESI variance yet) / 5b (validate against trusted reference) / 5c (new CI coverage); 5a is the agreed starting point for implementation | User: "I would start by writing and running calibration simulations... we don't know the RESI variance for these models yet" | 2026-09-04 |
| Confirmed `mlmRev::egsingle` + `geepack::muscatine` as the datasets; will fit both a random-intercept-only model (A) and a model with a random slope for `year` (B) on `egsingle` for the Goal 5a calibration check | User confirmed via question round | 2026-09-04 |
| Fixed Model A/B `lme4` convergence warning by standardizing `size`/`lowinc` (full-data mean/SD, fixed constants) and using `lmerControl(optimizer = "bobyqa")`; both models now fit with 0 convergence warnings | User: "consider scaling covariates first" | 2026-09-04 |
| Confirmed nested (not crossed) 3-level design and 60 total schools (5-89 children/school, mean 28.7); implemented two-stage hierarchical plasmode resampling (schools, then children within each drawn school) in `.longSimResampleHierarchical()`, replacing the earlier single-stage (pool-all-children) prototype | User: "it is nested not crossed... you might have to resample school, then resample student" | 2026-09-04 |
| Added `longitudinalCalibrationSim()` to `RESI/R/simulations.R` (Goal 5a), following `insurancePlasmodeSim()`'s conventions (`sim_raw/`, `summary_table.rds`), with `parallel::mclapply` over replicates defaulting to `mc.cores.reps = 24`; added `mlmRev` to `DESCRIPTION` Suggests; created `main.R` at the project root with copy-paste run lines; skipped random-effect variance-component SEs entirely | User: "Add it into simulations.R... 24 cores being the default... Write a main.R script... Don't care about SEs for RFX" | 2026-09-04 |
| Implemented Lyu & Welsh (2022) asymptotic scaling correction in `RESI/R/resi_pe_mixed.R` (lines 144-156); between-cluster parameters now use n=n_schools, within-cluster use n=n_schools×n_per_site. Verified with 60% reduction in systematic bias slopes; re-ran 1000-replicate calibration simulations and regenerated calibration report | Systematic linear bias growth with n_per_site violated asymptotic theory; investigation revealed single normalization factor applied to all terms instead of term-specific factors per Lyu & Welsh | 2026-10-?? |

| Removed the full-bootstrap pipeline (`longitudinalCalibrationSim()`, `main.R`, synthetic/plasmode/`make_*` scripts and their outputs); moved the hybrid plasmode and the Lyu-limit helpers into `RESI/R/mixed_simulations.R` (`hybridPlasmodeSim()`); history moved to `log/claudeLog.md` | User: delete superseded code (git keeps it), simulation into the package | 2026-11-15 |
| Correction to the 2026-10-?? row above: the term-specific `n` (g for between-school, N for within-school terms) is the Lyu-Welsh normalization; the "60% reduction in bias slopes" claim was not independently verified. RESI/SE targets are the Lyu-Welsh g, m -> Inf limit at the full-data ML estimates | Settled with Simon (checkpoints 001-002) | 2026-11-15 |
| Mixed-model CIs: model-based and robust CR2 (school-clustered) covariance; Normal/QF/CF; coverage of both the Lyu limit and S_m; reduced grid first; wire into `resi.lmerMod`/`resi.lme` | User decisions | 2026-11-15 |
| Mixed CIs implemented with the note's proposed defaults (design term included, `ci.design = FALSE` to drop; numerical-Hessian psi term for model-based; robust drops psi and uses unadjusted U_s in the meat term; normal/chi-square quantiles; `lmerMod` only) | Simon unavailable when asked; awaiting review | 2026-11-15 |
| Fixed a sign error in the large-noncentrality branch of `.resi_ci_qf` (`resi_asymptotic.R`; affects lm/glm QF CIs whenever n S^2 / Sigma_R or the observed statistic exceeds 1e4) | Found by the mixed CI coverage simulation (`year` ANOVA QF coverage 0.26) | 2026-11-15 |
