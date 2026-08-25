# Script to preprocess clinical data from prevac, ebovac2, and hamburg to harmonise naming conventions, data types etc.

# Libraries
library(tidyverse)
library(stringr)

folder_ebovac2 = fs::path("data", "ebovac2")
folder_hamburg = fs::path("data", "hamburg")
folder_prevac = fs::path("data", "prevac")
folder_is2 = fs::path("data", "is2")

# Paths
path_ebovac2 = fs::path(folder_ebovac2, "ebovac2_clinical.rds")
path_hamburg = fs::path(folder_hamburg, "hamburg_clinical.rds")
path_prevac = fs::path(folder_prevac, "prevac_clinical.RData")
path_is2 = fs::path(folder_is2, "is2_clinical_immResp.rds")

# Load the data
df_clinical_ebovac2 = readRDS(path_ebovac2)
df_clinical_hamburg = readRDS(path_hamburg)
df_clinical_prevac = readRDS(path_prevac)
df_clinical_is2 = readRDS(path_is2)


# Harmonise the column names and types of the major data of interest:
# Participant id,
# GE timepoint,
# vaccine group,
# age, sex,
# antibody measurements across time,
df_clinical_ebovac2_mutated = df_clinical_ebovac2 %>%
  mutate(
    participant_id = as.character(pid),
    study_accession = "ebovac2",
    sex = ifelse(sex == "F", "Female", "Male"),
    race = race,
    ethnicity = ethnicity,
    group = ifelse(arm == "active", "Ad26MVA", "placebo"),
    group_long = actarm,
    time = recode(time, "P+0" = "P+0D"),
    ab_p_0 = ab_p_0d,
    ab_p_7 = NA_real_,
    ab_p_14 = NA_real_,
    ab_p_28 = NA_real_,
    ab_p_56 = NA_real_,
    ab_p_63 = NA_real_,
    ab_p_84 = NA_real_,
    ab_p_180 = NA_real_,
    ab_p_365 = ab_p_364d,
    ab_b_0 = ab_b_0d,
    ab_b_21 = ab_b_21d#,
    # ab_p_0_a_brisbane_10_2007 = NA_real_,
    # ab_p_0_a_brisbane_59_2007 = NA_real_,
    # ab_p_0_b_florida_4_2006 = NA_real_,
    # ab_p_14_a_brisbane_10_2007 = NA_real_,
    # ab_p_14_a_brisbane_59_2007 = NA_real_,
    # ab_p_14_b_florida_4_2006 = NA_real_,
    # ab_p_28_a_brisbane_10_2007 = NA_real_,
    # ab_p_28_a_brisbane_59_2007 = NA_real_,
    # ab_p_28_b_florida_4_2006 = NA_real_
  )

df_clinical_hamburg_mutated = df_clinical_hamburg %>%
  mutate(
    participant_id = as.character(pid),
    study_accession = "hamburg",
    sex = ifelse(sex == "female", "Female", "Male"),
    race = ifelse(race == "Caucasian","White", "Asian"),
    ethnicity = "Unknown",
    age = age_baseline,
    group = "rVSV",
    group_long = "rVSV",
    time = recode(
      time,
      "d1" = "P+1D",
      "d3" = "P+3D",
      "d7" = "P+7D",
      "d0" = "P+0D"
    ),
    ab_p_0 = gp_specific_ab_gueckedou_0,
    ab_p_7 = gp_specific_ab_gueckedou_7,
    ab_p_14 = gp_specific_ab_gueckedou_14,
    ab_p_28 = gp_specific_ab_gueckedou_28,
    ab_p_56 = gp_specific_ab_gueckedou_56,
    ab_p_63 = NA_real_,
    ab_p_84 = gp_specific_ab_gueckedou_84,
    ab_p_180 = gp_specific_ab_gueckedou_180,
    ab_p_365 = NA_real_,
    ab_b_0 = NA_real_,
    ab_b_21 = NA_real_#,
    # ab_p_0_a_brisbane_10_2007 = NA_real_,
    # ab_p_0_a_brisbane_59_2007 = NA_real_,
    # ab_p_0_b_florida_4_2006 = NA_real_,
    # ab_p_14_a_brisbane_10_2007 = NA_real_,
    # ab_p_14_a_brisbane_59_2007 = NA_real_,
    # ab_p_14_b_florida_4_2006 = NA_real_,
    # ab_p_28_a_brisbane_10_2007 = NA_real_,
    # ab_p_28_a_brisbane_59_2007 = NA_real_,
    # ab_p_28_b_florida_4_2006 = NA_real_
  )

df_clinical_prevac_mutated = df_clinical_prevac %>%
  mutate(
    participant_id = as.character(pid),
    study_accession = "prevac",
    sex = sex,
    race = "Unknown",
    ethnicity = "Unknown",
    age = age_at_enrollement,
    group = arm_short,
    group_long = arm,
    time = recode(
      timepoint,
      "J0"    = "P+0D",
      "J03h"  = "P+3H",
      "J7"    = "P+7D",
      "J56"   = "B+0D",
      "J563h" = "B+3H",
      "J63"   = "B+7D"
    ),
    ab_p_0 = as.numeric(ab_d0),
    ab_p_7 = as.numeric(ab_d7),
    ab_p_14 = as.numeric(ab_d14),
    ab_p_28 = as.numeric(ab_d28),
    ab_p_56 = as.numeric(ab_d56),
    ab_p_63 = as.numeric(ab_d63),
    ab_p_84 = as.numeric(ab_m3),
    ab_p_180 = as.numeric(ab_m6),
    ab_p_365 = as.numeric(ab_m12),
    ab_b_0 = as.numeric(ab_d56),
    ab_b_21 = NA_real_#,
    # ab_p_0_a_brisbane_10_2007 = NA_real_,
    # ab_p_0_a_brisbane_59_2007 = NA_real_,
    # ab_p_0_b_florida_4_2006 = NA_real_,
    # ab_p_14_a_brisbane_10_2007 = NA_real_,
    # ab_p_14_a_brisbane_59_2007 = NA_real_,
    # ab_p_14_b_florida_4_2006 = NA_real_,
    # ab_p_28_a_brisbane_10_2007 = NA_real_,
    # ab_p_28_a_brisbane_59_2007 = NA_real_,
    # ab_p_28_b_florida_4_2006 = NA_real_
  )


df_clinical_is2_mutated = df_clinical_is2 %>%
  mutate(
    sex = gender,
    age = age_imputed,
    group = vaccine_name,
    group_long = vaccine_name,
    time = ifelse(
      study_time_collected %% 1 == 0,
      paste0("P+", as.integer(study_time_collected), "D"),
      paste0("P+", round(study_time_collected * 24), "H")
    ),
    ab_p_7 = NA_real_,
    ab_p_56 = NA_real_,
    ab_p_365 = NA_real_,
    ab_p_63 = NA_real_,
    ab_p_84 = NA_real_,
    ab_p_180 = NA_real_,
    ab_b_0 = NA_real_,
    ab_b_21 = NA_real_
  )

cols_to_select = c(
  "participant_id",
  "study_accession",
  "group",
  "group_long",
  "age",
  "sex",
  "race",
  "ethnicity",
  "time",
  "ab_p_0",
  "ab_p_7",
  "ab_p_14",
  "ab_p_28",
  "ab_p_56",
  "ab_p_63",
  "ab_p_84",
  "ab_p_180",
  "ab_p_365",
  "ab_b_0",
  "ab_b_21"#,
  # "ab_p_0_a_brisbane_10_2007",
  # "ab_p_0_a_brisbane_59_2007",
  # "ab_p_0_b_florida_4_2006",
  # "ab_p_14_a_brisbane_10_2007",
  # "ab_p_14_a_brisbane_59_2007",
  # "ab_p_14_b_florida_4_2006",
  # "ab_p_28_a_brisbane_10_2007",
  # "ab_p_28_a_brisbane_59_2007",
  # "ab_p_28_b_florida_4_2006"
)

df_clinical_ebovac2_mutated = df_clinical_ebovac2_mutated %>%
  dplyr::select(all_of(cols_to_select))

df_clinical_hamburg_mutated = df_clinical_hamburg_mutated %>%
  dplyr::select(all_of(cols_to_select))

df_clinical_prevac_mutated = df_clinical_prevac_mutated %>%
  dplyr::select(all_of(cols_to_select))

df_clinical_is2_mutated = df_clinical_is2_mutated %>%
  dplyr::select(all_of(cols_to_select))

df_clinical_all = bind_rows(bind_rows(
  bind_rows(df_clinical_ebovac2_mutated, df_clinical_hamburg_mutated),
  df_clinical_prevac_mutated
),
df_clinical_is2_mutated)  %>%
  mutate(study_vaccine = paste0(study_accession, "-", group)) %>%
  relocate(study_vaccine, .after = group)

# Define unified pre and post-vaccination responses
## This is necessary as not all studies have post-vaccination responses measured at the same timepoints
## In particular, since ebovac2 (Ad26/MVA) lacks day 180 and Hamburg (rVSV) lacks day 365,
## define post-response at day 365 for Ad26/MVA and at day 180 for rVSV and placebo.

# df_clinical_all = df_clinical_all %>%
#   mutate(response_pre = ab_p_0,
#          response_post = ifelse(group == "Ad26MVA", ab_p_365, ab_p_180))

processed_data_path <- fs::path("data")

saveRDS(df_clinical_all,
        file = fs::path(processed_data_path, "df_clinical_all.rds"))

rm(list = ls())
