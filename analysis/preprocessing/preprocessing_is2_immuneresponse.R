# File for preprocessing the raw immune response data files

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)
library(tidyr)

# Specify folder within folder root where the raw data lives
raw_data_folder       <- fs::path("data-raw", "is2")
processed_data_folder <- fs::path("data", "is2")

p_load_nAb      <- fs::path(raw_data_folder, "neut_ab_titer_2025-01-10_01-13-22.xlsx")
p_load_hai      <- fs::path(raw_data_folder, "hai_2025-01-10_01-13-41.xlsx")
p_load_clinical <- fs::path(processed_data_folder, "is2_clinical.rds")

# ---- Load raw data ----
response_hai <- read_excel(p_load_hai) %>% clean_names()
response_nAb <- read_excel(p_load_nAb) %>% clean_names()
is2_clinical <- readRDS(p_load_clinical)

# ---- Restrict to participants with gene expression data ----
participants <- is2_clinical %>% pull(participant_id) %>% unique()

response_hai <- response_hai %>% filter(participant_id %in% participants)
response_nAb <- response_nAb %>% filter(participant_id %in% participants)

# ---- Standardise columns and tag assay ----
response_hai <- response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai") %>%
  dplyr::select(participant_id, gender, race, cohort, study_time_collected,
                study_time_collected_unit, response_strain_analyte, value_preferred, assay)

response_nAb <- response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb") %>%
  dplyr::select(participant_id, gender, race, cohort, study_time_collected,
                study_time_collected_unit, response_strain_analyte, value_preferred, assay)

# ---- Merge assays ----
response_raw_merged <- bind_rows(response_nAb, response_hai) %>%
  arrange(participant_id) %>%
  distinct()

# ---- Bring in study/vaccine info ----
is2_studies <- is2_clinical %>%
  dplyr::select(participant_id, study_accession, group_long) %>%
  distinct()

response_raw_merged_studies <- response_raw_merged %>%
  inner_join(is2_studies, by = "participant_id")  # equivalent to merge(..., all = F), but keeps dplyr pipe + avoids row-order surprises

# ---- Assign each raw collection day to a nominal immune response timepoint ----
# Timepoints are defined to match the harmonised `ab_p_DAY` columns used
# downstream (see preprocessing_clinical_harmonisation.R). A raw collection
# day is assigned to the nearest nominal day if it falls within a +-7 day
# (inclusive) window of it; otherwise it is dropped.
nominal_days <- c(0, 7, 14, 28, 56, 63, 84, 180, 365)
day_window   <- 7

nearest_nominal_day <- function(day, targets, window) {
  diffs <- abs(day - targets)
  best  <- which.min(diffs)
  if (diffs[best] <= window) targets[best] else NA_real_
}

response_raw_merged_studies <- response_raw_merged_studies %>%
  mutate(
    study_time_collected_raw = study_time_collected,
    study_time_collected = vapply(
      study_time_collected_raw,
      nearest_nominal_day,
      numeric(1),
      targets = nominal_days,
      window  = day_window
    )
  ) %>%
  filter(!is.na(study_time_collected))

# If more than one raw collection day falls within the window of the same
# nominal day for a given participant/assay, keep only the raw day closest
# to the nominal day (consistently across all strains for that participant).
response_raw_merged_studies <- response_raw_merged_studies %>%
  group_by(participant_id, assay, study_time_collected) %>%
  mutate(day_dist = abs(study_time_collected_raw - study_time_collected)) %>%
  filter(study_time_collected_raw == study_time_collected_raw[which.min(day_dist)]) %>%
  ungroup() %>%
  dplyr::select(-day_dist, -study_time_collected_raw)

# ---- Clean strain names ONCE on distinct values, then join back ----
# (janitor::make_clean_names() de-duplicates repeated inputs, which corrupts
#  a data column with many repeated strain names if applied directly to it)
strain_lookup <- response_raw_merged_studies %>%
  distinct(response_strain_analyte) %>%
  mutate(strain_clean = janitor::make_clean_names(response_strain_analyte))

response_raw_merged_studies <- response_raw_merged_studies %>%
  left_join(strain_lookup, by = "response_strain_analyte")

# ---- Summarise across strains within participant/timepoint/assay ----
response_summary <- response_raw_merged_studies %>%
  group_by(participant_id, study_time_collected, assay, study_accession, group_long) %>%
  summarise(
    response_mean = mean(value_preferred, na.rm = TRUE),
    n_analytes    = n_distinct(response_strain_analyte),
    .groups = "drop"
  )

# ---- Choose one assay per study, preferring nAb over hai ----
assay_priority <- c("nAb", "hai")

study_assay_choice <- response_summary %>%
  filter(!is.na(response_mean)) %>%
  distinct(study_accession, assay) %>%
  mutate(priority = match(assay, assay_priority)) %>%
  filter(!is.na(priority)) %>%
  group_by(study_accession) %>%
  summarise(chosen_assay = assay[which.min(priority)], .groups = "drop")

# ---- Filter both summary and strain-level data to the chosen assay per study ----
# (do this once, reuse for both pivots, to avoid recomputing the join twice)
response_chosen_assay <- response_raw_merged_studies %>%
  left_join(study_assay_choice, by = "study_accession") %>%
  filter(assay == chosen_assay)

# ---- Pivot 1: mean response per participant x timepoint ----
response_wide <- response_chosen_assay %>%
  distinct(participant_id, study_accession, group_long, study_time_collected, assay) %>%
  left_join(response_summary,
            by = c("participant_id", "study_accession", "group_long",
                   "study_time_collected", "assay")) %>%
  dplyr::select(participant_id, study_accession, group_long, study_time_collected, response_mean) %>%
  mutate(study_time_collected = as.character(study_time_collected)) %>%
  pivot_wider(
    names_from  = study_time_collected,
    values_from = response_mean,
    names_prefix = "ab_p_"
  )

# ---- Pivot 2: individual strain-level response per participant x timepoint x strain ----
response_strain_wide <- response_chosen_assay %>%
  dplyr::select(participant_id, study_accession, group_long,
                study_time_collected, strain_clean, value_preferred) %>%
  mutate(study_time_collected = as.character(study_time_collected)) %>%
  pivot_wider(
    id_cols     = c(participant_id, study_accession, group_long),
    names_from  = c(study_time_collected, strain_clean),
    values_from = value_preferred,
    names_glue  = "ab_p_{study_time_collected}_{strain_clean}"
    # no values_fn: if duplicate rows exist per participant/time/strain this
    # will now correctly error out rather than silently averaging them
  )

# ---- Order columns: means by time, strain-level by time then strain ----
order_ab_cols <- function(cols) {
  time_vals <- as.numeric(str_extract(cols, "(?<=ab_p_)[0-9.]+"))
  cols[order(time_vals, cols)]
}

ab_cols_sorted        <- order_ab_cols(grep("^ab_p_", names(response_wide), value = TRUE))
strain_cols_sorted    <- order_ab_cols(grep("^ab_p_", names(response_strain_wide), value = TRUE))

response_wide        <- response_wide %>%
  dplyr::select(participant_id, study_accession, group_long, all_of(ab_cols_sorted))

response_strain_wide <- response_strain_wide %>%
  dplyr::select(participant_id, study_accession, group_long, all_of(strain_cols_sorted))

# ---- Combine summary means + strain-level columns into final dataframe ----
is2_immResp <- response_wide %>%
  left_join(response_strain_wide, by = c("participant_id", "study_accession", "group_long"))

# ---- Save ----
p_save <- fs::path(processed_data_folder, "is2_immResp.rds")
saveRDS(is2_immResp, file = p_save)

rm(list = ls())
