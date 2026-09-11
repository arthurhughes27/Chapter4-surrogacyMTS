# Master application script: runs all application scripts
# Main manuscript
source(fs::path("analysis", "application", "main", "application_lowDim.R"))
source(fs::path("analysis", "application", "main", "application_highDim.R"))

# Supporting information
source(fs::path("analysis", "application", "supplementary", "application_supplementary_FE.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_geneset.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_influenza_timepoints.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_influenza_timepoints_pi.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_ebola.R"))
