# R Script to preprocess files from the prevac trial transcriptomic substudy for analysis

# ---- Libraries ----
library(dplyr)
library(stringr)
library(readxl)
library(readr)
library(tidyr)
library(lubridate)

# ---- Specify paths to raw data ----
raw_data_path <- fs::path("data-raw", "prevac")
counts_path   <- fs::path(raw_data_path, "norm_counts_delete.RData")
clinical_path <- fs::path(raw_data_path, "clinical_PREVACImmuno.xlsx")
sample_path   <- fs::path(raw_data_path, "DM Prevac fichier Bordeaux.csv")
antibody_path <- fs::path(raw_data_path, "prevacup_db_reduite.rda")

# ---- Load raw data ----
# .RData / .rda files may contain multiple objects; load() returns their names
counts   <- get(load(counts_path))
clinical <- read_excel(clinical_path)

# Sample metadata (latin1 encoding due to special characters)
sample <- read_csv(sample_path,
                   locale = locale(encoding = "latin1"),
                   show_col_types = FALSE)

# Antibody measurements
antibody <- get(load(antibody_path))

# ---- Harmonize column names across datasets ----
# Ensure consistent identifiers for merging (PID and Timepoint)

sample <- sample %>%
  rename(PID = ParticipantID, Timepoint = TimePoint) %>%
  # Remove sampling date (not needed for analysis)
  dplyr::select(-`Sampling date`)

antibody <- antibody %>%
  rename(PID = SUBJID_N)


# ---- Reformat PID in antibody data (add missing hyphens) ----
# Original format: XXXYYYYZ  ->  XXX-YYYY-Z
antibody <- antibody %>%
  mutate(PID = str_c(str_sub(PID, 1, 3), "-", str_sub(PID, 4, 7), "-", str_sub(PID, 8, 8)))

# ---- Pivot antibody data to wide format ----
# One row per participant, one column per visit
antibody <- antibody %>%
  dplyr::select(PID, VISIT, TEST_MED_FINAL) %>%
  pivot_wider(names_from  = VISIT,
              values_from = TEST_MED_FINAL,
              id_cols     = PID)

# Prefix antibody columns with "Ab_"
colnames(antibody)[-1] <- paste0("Ab_", colnames(antibody)[-1])

# ---- Extract identifiers from counts rownames ----
# Rownames contain: "PID Date Timepoint"
counts_identifiers <- do.call(rbind, str_split(rownames(counts), pattern = " ")) %>%
  as.data.frame()

colnames(counts_identifiers) <- c("PID", "Date", "Timepoint")

# Harmonize date format (replace "/" by "-")
counts_identifiers$Date <- str_replace_all(counts_identifiers$Date,
                                           pattern = "/",
                                           replacement = "-")

# Combine identifiers with expression counts
counts_id <- data.frame(cbind(counts_identifiers, counts))

# ---- Merge all datasets ----
# Merge clinical, antibody, sample, and counts data into a single dataframe
prevac_merged <- clinical %>%
  # Merge with antibody data by PID
  inner_join(antibody, by = "PID") %>%
  # Merge with sample metadata
  inner_join(sample, by = "PID") %>%
  # Merge with counts identifiers + counts by PID and Timepoint
  inner_join(counts_id, by = c("PID", "Timepoint"))

# ---- Recode treatment arm into simplified categories ----
prevac_merged <- prevac_merged %>%
  mutate(
    arm_short = recode_factor(
      Arm,
      "Undiluted rVSVDG-ZEBOV-GP/Placebo boost" = "rVSV",
      "rHAd26/MVA-BN-Filo boost" = "Ad26MVA",
      "Undiluted rVSVDG-ZEBOV-GP/VSVDG-ZEBOV-GP boost" = "rVSV",
      "Pooled Placebo/Placebo boost" = "placebo"
    )
  )

# ---- Reorder columns for clarity ----
prevac_merged <- prevac_merged %>%
  dplyr::select(
    all_of(
      c(
        "PID",
        "Date",
        "Timepoint",
        "Arm",
        "arm_short",
        "Age at enrollement",
        "Sex",
        colnames(antibody)[-1]  # antibody columns
      )
    ),
    everything()  # keep remaining columns (e.g., counts) at the end
  )

# ---- Clean Timepoint strings ----
# Remove '+' characters from Timepoint labels
prevac_merged <- prevac_merged %>%
  mutate(Timepoint = str_replace_all(Timepoint, pattern = "\\+", ""))

# Remove samples (rows) with RNA Integrity number < 6
prevac_merged = prevac_merged %>% 
  filter(RIN > 6)

# Clean all names
prevac_merged = prevac_merged %>% 
  janitor::clean_names()

prevac_clinical = prevac_merged %>% 
  dplyr::select(-(a1cf:zzz3))

# ---- Specify path to save processed data ----
processed_data_path <- fs::path("data", "prevac")

# ---- Save merged and preprocessed dataset ----
saveRDS(prevac_merged,
        file = fs::path(processed_data_path, "prevac_merged.RData"))

saveRDS(prevac_clinical,
        file = fs::path(processed_data_path, "prevac_clinical.RData"))

rm(list = ls())
