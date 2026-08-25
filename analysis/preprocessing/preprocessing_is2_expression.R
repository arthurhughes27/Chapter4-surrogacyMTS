# File for preprocessing gene expression data from IS2 for downstream use
# The end result should be a dataframe with columns as variables and rows as samples
# The first two columns will give the participant identifier and time of sample, allowing for unique identification

# Load necessary packages
library(fs)
library(dplyr)
library(stringr)

# Specify folder within folder root where the raw data lives
raw_data_folder = fs::path("data-raw", "is2")

# Specify folder within folder root where the processed data lives
processed_data_folder = fs::path("data", "is2")

# Use fs::path() to specify the data paths robustly
p_load_allnoNorm <- fs::path(raw_data_folder, "all_noNorm_eset.rds")

p_load_clinical <- fs::path(processed_data_folder, "is2_clinical.rds")

# Read in the rds file
all_noNorm_eset <- readRDS(p_load_allnoNorm)

is2_clinical = readRDS(p_load_clinical)

# Load the expression data
is2_expr = all_noNorm_eset@assayData[["exprs"]] %>%
  t() %>%
  as.data.frame()

# Make the column names lowercase
colnames(is2_expr) = is2_expr %>%
  colnames() %>%
  tolower()

# In the non-noNormalised data : remove genes with NA in some studies
# This reduces the number of genes to the same amount as the noNormalised data
# young_nonoNorm_expr = young_nonoNorm_expr[, colSums(is.na(young_nonoNorm_expr)) == 0]

# Now extract information for the first two columns (participant_id and study_time_collected)
sample_info_all_noNorm = rownames(is2_expr)


# Use `stringr` and regular expressions to extract the participant_id and study_time_collected from the unique identifiers
matches_all_noNorm <- str_match(sample_info_all_noNorm , "^(SUB[0-9.]+)_(-?[0-9.]+)_Days")


# Participant ids
participant_id_all_noNorm <- matches_all_noNorm[, 2]

# Study times (numeric, rounded)
study_time_collected_all_noNorm <- matches_all_noNorm[, 3] %>%
  as.numeric() %>%
  round(2)

# Insert the identifying information as the first two columns
pids = is2_clinical$participant_id %>%
  unique()

is2_expr <- is2_expr %>%
  mutate(participant_id = participant_id_all_noNorm, study_time_collected = study_time_collected_all_noNorm) %>%
  dplyr::select(participant_id, study_time_collected, everything()) %>%
  filter(participant_id %in% pids)



# Save processed dataframes

# Use fs::path() to specify the data path robustly
p_save_all_noNorm <- fs::path(processed_data_folder, "is2_expr.rds")

# Save dataframe
saveRDS(is2_expr, file = p_save_all_noNorm)

rm(list = ls())
