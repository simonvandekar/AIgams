# RESI Sample Size Scaling Analysis

## Summary of Finding

The RESI bias in the first row of the calibration figure increases **linearly** with `n_per_site`, which indicates the scaling factor `n` is **incorrect** for school-level (between-cluster) covariates.

## Data Analysis (n_schools=100)

### Observed Bias Slopes vs n_per_site

**School-level terms** (should vary only with n_schools):
- **(Intercept)**: slope = +0.000134, ratio_to_effect = 0.00332
- **size_100**: slope = -0.0000605, ratio_to_effect = -0.00329  
- **lowinc_10**: slope = -0.000363, ratio_to_effect = -0.00332

These three terms show consistent behavior with slope/effect ≈ ±0.0033 per unit n_per_site.

**Student/Child-level terms** (should vary with n_schools × n_per_site):
- **year**: slope = +0.00000733, ratio_to_effect ≈ 0
- **female**: slope ≈ 0, ratio_to_effect = 0.00106
- **black**: slope = -0.0000709, ratio_to_effect = -0.000457
- **hispanic**: slope = -0.0000102, ratio_to_effect = -0.000130

## Root Cause

The current RESI implementation uses **a single normalization factor** `n` for ALL terms in a mixed model:

```R
# From resi_pe_mixed.R
term_ns <- setNames(rep(selected_n, length(term_levels)), names(term_levels))
```

With `level = 2` (used in the simulation), all terms receive `n = n_students ≈ n_schools × n_per_site`.

This is **incorrect** according to the **Lyu & Welsh (2022)** paper, which shows that:

## Lyu & Welsh (2022) Theory: Different Normalizations Needed

From "Increasing cluster size asymptotics for nested error regression models":

> "We need to treat within cluster parameters (i.e. coefficients of unit-level covariates that vary within clusters and the within cluster variance σₑ²) **differently from between cluster parameters** (i.e. coefficients of cluster-level covariates that are constant within clusters and the between cluster variance σₐ²) because they require **different normalizations**."

**Key insight**: 
- **Within-cluster parameters** need normalization by √(g × m) where g = number of clusters, m = cluster size
- **Between-cluster parameters** need normalization by √g only (where g = number of clusters)

## Correct Scaling Factor Formula

For each coefficient or ANOVA term:

| Parameter Type | Varies | Correct n | Current (Wrong) n |
|---|---|---|---|
| **Intercept** | Between-cluster only | n_schools | n_schools × n_per_site |
| **School-level covariates** (size_100, lowinc_10) | Between-cluster only | n_schools | n_schools × n_per_site |
| **Observation-level covariates** (female, black, hispanic, year) | Within-cluster | n_schools × n_per_site | n_schools × n_per_site |

## Why School-Level Terms Show Linear Bias with n_per_site

When RESI computes S = z / √n, the Z-statistic scales with √(n_schools), but:

- Currently using n = n_schools × n_per_site
- Should use n = n_schools
- The ratio error: √(n_schools × n_per_site) / √(n_schools) = √(n_per_site)

This means for school-level terms:
- S is divided by √(n_per_site) too much
- As n_per_site increases, RESI estimates become artificially smaller
- Bias grows approximately as: bias ≈ -effect_size × (1 - 1/√(n_per_site)) ≈ -effect_size / (2√(n_per_site))
- At first order, this is approximately linear in n_per_site

## Implementation Fix Required

The `resiMixedScaleMap()` function in `RESI/R/resi_pe_mixed.R` needs to:

1. **Identify variation level** for each term (already stored in the design)
2. **Assign different n values** based on term type:
   - Between-cluster terms → use n_clusters (the selected level at hierarchy level 1)
   - Within-cluster terms → use current logic (selected level count)

### Pseudocode Fix:

```R
.resiMixedScaleMap <- function(design, n = NULL, level = NULL) {
  # ... existing code ...
  
  # Instead of:
  # term_ns <- setNames(rep(selected_n, length(term_levels)), names(term_levels))
  
  # Do:
  term_ns <- setNames(
    ifelse(
      is_between_cluster_term,  # Based on variation_level
      hierarchy_counts[[1L]],    # Use top-level cluster count
      selected_n                 # Use current selection for within-cluster
    ),
    names(term_levels)
  )
}
```

## Expected Improvement

With this fix, RESI bias for school-level terms should:
- No longer vary with n_per_site (flat line)
- Show bias that's consistent across the n_schools × n_per_site grid
- Be properly calibrated as sample size increases

## References

- Lyu, Z., & Welsh, A.H. (2022). "Increasing cluster size asymptotics for nested error regression models." *Journal of Statistical Planning and Inference*, 217, 52-68.
- Current calibration simulation: `tests/test_egsingle_hybrid_plasmode_fixed_n.R`
- RESI implementation: `RESI/R/resi_pe_mixed.R`
