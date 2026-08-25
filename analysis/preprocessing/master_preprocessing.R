# Master script to run all preprocessing files in order

# Preprocessing the original datasets
source(fs::path("analysis", "preprocessing", "preprocessing_hamburg.R"))
source(fs::path("analysis", "preprocessing", "preprocessing_ebovac2.R"))
source(fs::path("analysis", "preprocessing", "preprocessing_prevac.R"))
source(fs::path("analysis", "preprocessing", "preprocessing_is2_master.R"))

# Harmonising and merging the data across studies
source(fs::path("analysis", "preprocessing", "preprocessing_clinical_harmonisation.R"))
source(fs::path("analysis", "preprocessing", "preprocessing_merging_harmonisation.R"))

# Preprocessing geneset data
source(fs::path("analysis", "preprocessing", "preprocessing_BTM.R"))
source(fs::path("analysis", "preprocessing", "preprocessing_BG3M.R"))

