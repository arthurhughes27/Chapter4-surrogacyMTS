# Libraries
library(tidyverse)
library(stringr)

# Load data
processed_data_path <- fs::path("data", "is2")

path_is2_clinical = fs::path(processed_data_path, "is2_clinical.rds")
path_is2_immResp = fs::path(processed_data_path, "is2_immResp.rds")
path_is2_expr = fs::path(processed_data_path, "is2_expr.rds")

is2_clinical = readRDS(path_is2_clinical) %>% 
  dplyr::select(-group_long) %>% 
  distinct()
is2_immResp = readRDS(path_is2_immResp)
is2_expr = readRDS(path_is2_expr)

# Merge clinical and immune response
is2_clinical_immResp = is2_clinical %>%
  left_join(is2_immResp, by = c("participant_id","study_accession"))

is2_merged = is2_clinical_immResp %>% 
  right_join(is2_expr, by = c("participant_id", "study_time_collected"))

saveRDS(is2_merged,
        file = fs::path(processed_data_path, "is2_merged.rds"))

saveRDS(is2_clinical_immResp,
        file = fs::path(processed_data_path, "is2_clinical_immResp.rds"))

rm(list = ls())
