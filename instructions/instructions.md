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
`glmmTMB`, point estimates only, CI "not developed" warning) — to cover:

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
**Status:** Needs decision — see "Dataset investigation" below
**Description:** Identify or construct a large longitudinal dataset (~1000+
independent subjects/clusters) to serve as the plasmode base for simulations,
analogous to how `insurancePlasmodeSim()` uses the `insurance` dataset (N=1338).
**Key design decisions:**
- [FILL IN] Which dataset(s) to use (see options below) — need confirmation.
- [FILL IN] Whether to keep real covariates/outcome (pure plasmode) or keep real
  covariates but simulate a known outcome model on top (to get ground-truth mixed
  model / smooth-term parameters for calibration, since none of the candidate
  real datasets have both >=1000 subjects AND a continuous outcome AND enough
  covariate richness for smooth terms).

**Dataset investigation (done):**
Checked `RESI/tests/testthat/test-resi.R` — existing tests use
`nlme::Orthodont` (27 subjects) and `lme4::sleepstudy` (18 subjects) for
`lmerMod`/`lme` tests. Both are far too small for a ~1000-subject simulation.

Surveyed built-in datasets in `lme4`, `nlme`, and `geepack` (already a `Suggests`
dependency of RESI, used in existing tests) for longitudinal structure (repeated
measures per independent subject) and subject count:

| Dataset | Package | N subjects | Obs/subject | Outcome | Notes |
|---|---|---|---|---|---|
| `sleepstudy` | lme4 | 18 | 10 | continuous (Reaction) | used in existing tests; too small |
| `Orthodont` | nlme | 27 | 4 | continuous (distance) | used in existing tests; too small |
| `Quinidine` | nlme | 136 | ~11 | continuous (conc) | too small |
| `InstEval` | lme4 | 2972 students | crossed w/ 1128 instructors | ordinal rating | crossed, not a simple longitudinal/time structure |
| `ohio` | geepack | 537 | 4 | binary (wheeze) | classic GEE example; still <1000 |
| `dietox` | geepack | 72 pigs | ~12 | continuous (weight) | too small |
| **`muscatine`** | **geepack** | **4856** | **3 (ages 6,8,10 relative to `base_age`)** | **binary (`obese`)** | **largest candidate; age spans 6-18 across the sample, `gender` covariate; balanced 3-obs/subject design** |

**Recommendation (pending your confirmation):** `geepack::muscatine` is the only
built-in dataset across `lme4`/`nlme`/`geepack` with >=1000 independent subjects
and a genuine repeated-measures/longitudinal structure. It has a *binary* outcome,
which works for `glmer`/`bam(family=binomial)`/GEE-style validation but not for a
plain continuous-outcome `lmer`/`lme`/Gaussian-`gam` plasmode the way `insurance`
does for `lm`/`glm`. Options to resolve this, in order of how closely they match
the paper's plasmode design:
  (a) Use `muscatine` as-is for the binary/GLMM and `gam(family=binomial)` /
      `bam` settings (true plasmode, matches paper philosophy).
  (b) For continuous-outcome settings (`lmer`, `lme`, Gaussian `gam`/`gamm`),
      resample `muscatine`'s subject/covariate structure (age, gender, cluster
      sizes) but simulate a continuous outcome from a known mixed/smooth model —
      this is no longer a pure plasmode design (the "truth" is simulated, not the
      full-data MLE), so calibration targets change accordingly.
  (c) Find/construct a different >=1000-subject dataset with a continuous
      outcome — searched further per your request below; **resolved, see
      "Extended search results"**.

**Extended search results (done):** Installed and surveyed `mlmRev`, `JM`,
`joineR`, and `merTools` (none of these were previously `Suggests` dependencies
of RESI; installed locally for investigation only, not yet added to
`DESCRIPTION`).

| Dataset | Package | N subjects | Cluster/nesting | Occasions/subject | Outcome | Notes |
|---|---|---|---|---|---|---|
| **`egsingle`** | **mlmRev** | **1721 children** | **3-level: occasion in child in school (60 schools)** | **2-6 (unbalanced), continuous time `year`** | **continuous (`math`, growth-curve score)** | **Best match: only real dataset found with >=1000 subjects, a continuous repeated-measures outcome, and a continuous time variable suitable for a smooth term. See variable table below.** |
| `star` | mlmRev | 11598 students | 2-level: student in school (80 schools), classroom/treatment also nested | 1-4 (grades K-3) | continuous (`read`, `math`) | Larger than egsingle; TN STAR class-size experiment; classroom `cltype` (treatment) adds a 3rd nesting level (class in school); more complex, keep as a backup/stress-test dataset |
| `Chem97` | mlmRev | 31022 students | 3-level: student in school (2410) in LEA (131) | 1 (cross-sectional, not longitudinal) | continuous (`score`) | Good nested/nested (not crossed) nuisance structure but no repeated measures per student |
| **`InstEval`** | **lme4** (already installed) | 2972 students crossed with 1128 instructors | **crossed, not nested** | n/a (each student rates multiple, non-nested instructors) | ordinal (`y`, 1-5) | As you suggested: interesting *crossed* random-effects case; keep as a distinct, later test case for effect-size definition under crossed designs, not part of the primary longitudinal goal |
| **`Contraception`** | **mlmRev** | 1934 women | **60 districts, cluster sizes 2-118 (median 26)** | 1 (cross-sectional, not longitudinal) | binary (`use`) | Matches your "small number of sites, many observations per site" request; not longitudinal (one row/woman), so useful as a *clustered* (not repeated-measures) test case, e.g. for `gam`/`glmer` with a random intercept for district |
| `pbc2` / `aids` | JM | 312 / 467 | subject only | up to 16 / 5 | continuous (`serBilir`) / CD4 count | too small |
| `epileptic` / `heart.valve` / `liver` / `mental` | joineR | 605 / 256 / 488 / 150 | subject only | varies | continuous | all too small |
| `hsb` | merTools | (same underlying data as `Hsb82`/`MathAchieve`, ~7185 students / 160 schools) | 2-level | 1 (cross-sectional) | continuous | no repeated measures |

**Recommendation:** Use `mlmRev::egsingle` as the primary plasmode dataset for
the continuous-outcome, longitudinal goals (`lmer`, `lme`, Gaussian `gam`/`gamm`)
— it resolves the gap identified in the original search (>=1000 subjects *and*
continuous outcome *and* a continuous time covariate for smooth terms). Keep
`geepack::muscatine` as the binary-outcome companion dataset (for
`glmer`/`bam(family=binomial)`). Treat `lme4::InstEval` (crossed) and
`mlmRev::Contraception` (few-large-sites, non-longitudinal cluster) as
secondary/stress-test datasets for later, per your comments, not part of the
initial plasmode.

**Variable / nesting-structure table for `egsingle` (for specifying model
formulas):**

| Variable | Level | Varies within | Type | Description |
|---|---|---|---|---|
| `schoolid` | 3 (school, n=60) | — | factor | school identifier |
| `childid` | 2 (child, n=1721) | school | factor | child identifier, nested in school |
| `year` | 1 (occasion, n=7230 rows) | child | continuous | time, centered (range -2.5 to 2.5, integer steps) |
| `grade` | 1 (occasion) | child | continuous (0-4) | grade level, collinear with `year` |
| `math` | 1 (occasion) | child | continuous | **outcome**: math achievement score |
| `retained` | 1 (occasion) | child (mostly time-invariant: 1366/1721 children never change; time-varying for 355 children) | binary factor | held back a grade |
| `female` | 2 (child) | — (time-invariant) | binary factor | child sex |
| `black`, `hispanic` | 2 (child) | — (time-invariant) | binary factor | child race/ethnicity indicators |
| `size` | 3 (school) | — (time-invariant) | continuous | school enrollment size |
| `lowinc` | 3 (school) | — (time-invariant) | continuous | % low income at school |
| `mobility` | 3 (school) | — (time-invariant) | continuous | % student mobility at school |

**Confirmed (2026-09-04):** Use `mlmRev::egsingle` as the primary continuous-
outcome dataset and `geepack::muscatine` as the binary-outcome companion. Fit
formulas **with and without a random slope for `year`** (per your answer) —
see Goal 5a for the two working model variants.

**Files:** [FILL IN once formulas are set, e.g. `RESI/R/simulations.R`
additions, `longitudinalPlasmodeSim()`]

### Goal 2: Extend `resi_pe`/`resi_asymptotic` machinery to `lmer`/`lme`
**Status:** Not started
**Description:** Add asymptotic CI methods (QF and CF — mirroring
`resi_asymptotic.R`'s `.resi_precompute_ext()`/`.resi_contrast_ext()` framework)
for the existing `resi_pe.lmerMod` (L-RESI/CS-RESI) point estimates, and the
analogous `resi_pe.lme` (nlme) if/when that exists. Currently `resi.lmerMod()`
and `resi.lme()` explicitly warn that no CI procedure exists and return point
estimates only.
**Key design decisions:**
- Does the extended influence-function derivation in the paper (which
  treats $A_\theta$, $B_\theta$ as free parameters via the Delta method) carry
  over directly when the "meat" is a **cluster-robust** (`clubSandwich::vcovCR`)
  covariance instead of an HC (heteroscedasticity-consistent) covariance? Need a
  math derivation note in `notes/` before implementing.
- What is the effective sample size $n$ for the L-RESI CI (number of
  clusters, as used for the point estimate) vs. CS-RESI CI (number of
  observations, per the existing GEE CS-RESI convention)?
- For the L-RESI the unit really depends on the study design. For small number of sites with large number of observations at each site the unit is really the individuals at each site. For many individuals with repeated measures in individual, the unit is the individual as well. Because the individual occurs at different levels in the mixed model, there really needs to be an option for the user to specify the unit and for the code to try to figure it out if they don't and report clearly what the unit is. This affects the scaling factor "n" for computing the RESI effect size.
- `glmer` should (generalized linear mixed models) be in scope for this
  goal, or only Gaussian `lmer`/`lme`. But let's hold off for now until we've validate mixed models and gams.
**Expected behavior / acceptance criteria:** CI coverage
validated via `simCalibrationSim`-style study against a bootstrap reference,
similar to how QF/CF were validated against the existing bootstrap for `lm`/`glm`. Before that though, the calibration style validation simulations.
**Files:**
- Functions: `RESI/R/resi_asymptotic.R`, `RESI/R/resi.R` (`resi.lmerMod`,
  `resi.lme`), `RESI/R/resi_pe.R` (`resi_pe.lmerMod`, new `resi_pe.lme`?)
- Output: Similar to the existing function outputs with estimates and CIs.

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
**Status:** **5a is the immediate priority / starting point for this project**
(per your note); 5b/5c depend on Goal 1 dataset+formula confirmation.
**Description:** For each new method (Goals 2-4), follow the paper's three-step
verification.

**5a. Simulation-internal-consistency / calibration check (START HERE).**
Before any new RESI-variance formula is derived (we don't yet have one for
`lmer`/`lme`/`gam`/`gamm`), confirm the plasmode resampling machinery itself is
sound: draw resamples from the full `egsingle` (and/or `muscatine`) data, refit
the model in each replicate, and check that as replicate sample size grows,
  (i) the Monte Carlo average of the fixed-effect coefficients converges to the
      full-data (population) coefficient estimates,
  (ii) the Monte Carlo average of the variance components (residual variance,
       random-effect variances/covariances, or smoothing parameters for `gam`)
       converges to the full-data values,
  (iii) the Monte Carlo *standard error* (empirical SD across replicates) of
       each coefficient/variance-component estimator converges to the
       full-data model-based (or sandwich/cluster-robust) standard error.
This is exactly analogous to `simCalibrationSim`'s existing calibration checks
for `lm`/`glm`, but does **not** yet include a RESI-variance check, since that
requires Goal 2/3's (not-yet-derived) asymptotic variance of the RESI for these
model types.

**Working model formulas for `egsingle` (confirmed 2026-09-04 — fit with and
without a random slope for `year`):**
- Fixed effects (shared): `math ~ year + female + black + hispanic + size_100 + lowinc_10`,
  where `size_100 = size / 100` (school enrollment in hundreds of students) and
  `lowinc_10 = lowinc / 10` (percent low-income in units of 10 percentage
  points) are fixed-divisor rescalings (not z-scores), applied identically to
  every replicate, so coefficients stay directly interpretable. Raw `size` is on
  a very different scale than the other covariates (range 113-1486, SD 312.6,
  vs. e.g. `year` SD 1.39) and caused an `lme4` convergence warning
  (`max|grad| = 0.0023`); after rescaling `size`/`lowinc` and switching to
  `lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))`, both Model A
  and Model B fit with 0 convergence warnings on the full data.
- Model A (random intercepts only): `lmer(... + (1 | schoolid/childid))`;
  `lme(..., random = ~1 | schoolid/childid)`; `gamm(math ~ s(year) + female + ...,
  random = list(schoolid = ~1, childid = ~1))`
- Model B (random slope for `year` at the child level, random intercept at the
  school level): `lmer(... + (1 + year | childid) + (1 | schoolid))`;
  `lme(..., random = list(schoolid = ~1, childid = ~1 + year))`; `gamm(math ~
  s(year) + female + ..., random = list(schoolid = ~1, childid = ~1 + year))`
- [FILL IN if wrong] A fully-crossed double-nested random slope
  (`(1+year|schoolid/childid)`) was avoided as likely over-parameterized/hard to
  converge with only 60 schools; flag if you wanted that instead of the
  intercept-only school + slope-at-child form above.
- `egsingle` has 60 schools total (children/school: min 5, median 24, mean 28.7,
  max 89). The nested 3-level structure (measurement occasion in child in
  school) means plasmode resampling must be **two-stage**: draw `n_schools`
  schools with replacement, then within each drawn school draw (with
  replacement) as many children as that source school originally had, keeping
  each child's full row set (all occasions) intact. Implemented as
  `.longSimBuildIndex()`/`.longSimResampleHierarchical()` in `RESI/R/simulations.R`;
  a single-stage (pool-all-children) version was tried first and is
  methodologically wrong for this nesting — replaced.
- Random-effect (variance component) standard errors are **not** computed/tracked
  (per your "Don't care about SEs for RFX"); only fixed-effect SEs are used for
  the Monte-Carlo-SD-vs-model-SE calibration ratio. Variance-component point
  estimates are still tracked (bias vs. full-data truth).

**5b. Validate against a trusted reference** (e.g., bootstrap RESI for
`lmerMod`/`lme`, or the existing `geeglm` CS-RESI/L-RESI as a cross-check for
the `lmer` case, mirroring `test-resi.R`'s `"geeglm (exchangeable, positive
rho) L-RESI matches lmerMod (CR0) L-RESI"` test).

**5c. Evaluate new CI coverage/width/bias** against the validated pipeline,
once Goal 2/3/4's asymptotic variance derivations exist.
**Key design decisions:**
- Directory naming: follow `<MethodName>Sim/` convention (e.g.
  `resiLmerAsympSim/`, `resiGamAsympSim/`). Sounds good.
- Validation output: create a **new** article (separate from
  `RESI/vignettes/articles/ci_validation.Rmd`) for mixed models and gams.
**Files:** `RESI/R/simulations.R` additions; `<MethodName>Sim/` output directories
(`sim_raw/`, `summary_table.rds`, `figures/`)

---

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
