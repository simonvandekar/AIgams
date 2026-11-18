# [AI-generated: GitHub Copilot, 2026-11-17]
# [Prompt summary: shared egsingle data preparation and Model A / Model B
#  formulas for the hybrid plasmode scripts.]

# mlmRev::egsingle prepared for the hybrid plasmode simulations, with
# school-level covariates rescaled (size in hundreds of students, lowinc in
# tens of percentage points).
egsingle <- local({
  e <- new.env()
  utils::data("egsingle", package = "mlmRev", envir = e)
  e$egsingle
})
egsingle$size_100 <- egsingle$size / 100
egsingle$lowinc_10 <- egsingle$lowinc / 10

# Model A: school and student random intercepts. Model B: adds a
# student-level random slope for year (school level intercept only).
egsingle_models <- list(
  A = math ~ year + female + black + hispanic + size_100 + lowinc_10 +
    (1 | schoolid/childid),
  B = math ~ year + female + black + hispanic + size_100 + lowinc_10 +
    (1 | schoolid) + (year | childid:schoolid))
