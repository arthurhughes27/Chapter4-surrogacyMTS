# Master simulation script: runs all simulation scripts
# Main manuscript
source(fs::path("analysis", "simulation", "main", "simulation_parametric_calibration_main.R"))
source(fs::path("analysis", "simulation", "main", "simulation_parametric_power_main.R"))

# Supporting information
source(fs::path("analysis", "simulation", "supplementary", "simulation_parametric_power_trueMean.R"))
source(fs::path("analysis", "simulation", "supplementary", "simulation_parametric_calibration_trueMean.R"))
source(fs::path("analysis", "simulation", "supplementary", "simulation_parametric_calibration_metaAnalysisModel.R"))
source(fs::path("analysis", "simulation", "supplementary", "simulation_nonparametric_nstudies.R"))
source(fs::path("analysis", "simulation", "supplementary", "simulation_nonparametric_nsamples.R"))