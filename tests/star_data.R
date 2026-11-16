# mlmRev::star (Tennessee STAR) prepared for the mixed-model RESI validation
# simulation: school (sch) > student > yearly measurement, the same nesting as
# egsingle. Complete cases on the model variables. Rare categories are merged
# so that resampled designs rarely lose a factor level: ethnicity becomes a
# `black` indicator (only ~65 students are neither white nor black) and
# school type becomes city (inner-city + urban; 7 schools are urban) /
# suburb / rural. The student id is made unique within school (446 students
# appear in more than one school); fl_10 is the school's percentage of
# free-lunch measurements in units of 10 percentage points (a continuous
# between-school covariate).
load_star <- function() {
  star <- local({
    e <- new.env()
    utils::data("star", package = "mlmRev", envir = e)
    e$star
  })
  vars <- c("math", "yrs", "sx", "eth", "ses", "cltype", "schtype", "sch", "id")
  star <- stats::na.omit(star[vars])
  star$black <- factor(as.integer(star$eth == "B"))
  star$schtype3 <- factor(ifelse(star$schtype %in% c("inner", "urban"), "city",
                                 as.character(star$schtype)),
                          levels = c("city", "suburb", "rural"))
  star$ses <- stats::relevel(factor(star$ses), ref = "N")
  star$cltype <- stats::relevel(factor(star$cltype), ref = "reg")
  star$sch <- droplevels(star$sch)
  star$sid <- factor(paste(star$sch, star$id, sep = "_"))
  star$fl_10 <- 10 * stats::ave(star$ses == "F", star$sch)
  rownames(star) <- NULL
  star
}

form_star <- math ~ yrs + sx + black + ses + cltype + schtype3 + fl_10 + (1 | sch/sid)
