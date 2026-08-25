# Master application script: runs all application scripts
# Main manuscript
source(fs::path("analysis", "application", "main", "application_lowDim.R"))
source(fs::path("analysis", "application", "main", "application_highDim.R"))
source(fs::path("analysis", "application", "main", "application_ebola.R"))

# Supporting information
source(fs::path("analysis", "application", "supplementary", "application_supplementary_FE.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_geneset.R"))
source(fs::path("analysis", "application", "supplementary", "application_supplementary_day7.R"))
