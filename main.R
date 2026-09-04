# main.R
# ------------------------------------------------------------------
# Copy-paste-ready lines to run the Goal 5a longitudinal calibration
# simulation (RESI/R/simulations.R: longitudinalCalibrationSim()).
# See instructions/instructions.md (Goal 5a) and log/claudeLog.md.
#
# Run this with the working directory set to the AIgams project root
# (the folder that contains this file and the RESI/ subdirectory),
# e.g. `Rscript main.R` from this directory, or open this file in
# RStudio/VS Code with the project root as the working directory.
# ------------------------------------------------------------------

devtools::load_all("RESI")

# --- Quick smoke test (~10-30s): confirms the pipeline runs end to end ---
# out <- longitudinalCalibrationSim(
#   nsim              = 10,
#   n_schools.vec     = c(30, 60),
#   model             = "both",
#   output.dir        = "resiLongCalibrationSim_test",
#   mc.cores.settings = 1,
#   mc.cores.reps     = 4
# )
# print(out)

# --- Full run ---
# Uses 24 cores by default (mc.cores.reps). Reduce mc.cores.reps if your
# machine/server has fewer cores available.
out <- longitudinalCalibrationSim(
  nsim              = 1000,
  n_schools.vec     = c(15, 30, 60, 90, 120),
  model             = "both",
  output.dir        = "resiLongCalibrationSim",
  mc.cores.settings = 1,
  mc.cores.reps     = 24
)

print(out)
