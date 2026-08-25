# -------------------------------------------------------------------------
# R Script to preprocess files from the ebovac2 trial transcriptomic substudy
# Purpose: load and harmonize clinical, transcriptomic, treatment, and antibody data
# -------------------------------------------------------------------------

# ---- Libraries ----
library(janitor)   # clean column names
library(dplyr)     # data manipulation
library(stringr)   # string handling
library(readxl)    # read Excel files
library(readr)     # read RDS / CSV files
library(tidyr)     # pivoting / reshaping
library(lubridate) # date handling
library(haven)     # read SAS datasets
library(fs)        # file path handling

# ---- File paths ----
raw_data_path    <- fs::path("data-raw", "ebovac2")
clinical_path    <- fs::path(raw_data_path, "EBOVAC2_recap_final.xlsx")
counts_path      <- fs::path(raw_data_path, "raw_counts_trimmo.RDS")
treatment_path   <- fs::path(raw_data_path, "adsl.sas7bdat")
antibody_path    <- fs::path(raw_data_path, "adishum.sas7bdat")

# ---- Load raw data ----
# Clinical metadata
clinical <- read_excel(clinical_path) %>%
  as.data.frame()

# Transcriptomic counts
counts <- read_rds(counts_path) %>%
  as.data.frame()

# Treatment information
treatment <- read_sas(treatment_path) %>%
  as.data.frame()

# Antibody measurements
antibody <- read_sas(antibody_path)

# ---- Harmonize clinical identifiers ----
# Extract numeric subject ID from 'Subject' column to match treatment data
clinical <- clinical %>%
  mutate(Subject = as.numeric(substr(Subject, 6, 11)))

# ---- Keep only relevant columns ----
clinical <- clinical %>%
  dplyr::select(Sample, Subject, Timepoint) %>%
  rename(Sample = Sample,
         SUBJID = Subject,
         Time   = Timepoint)

clinical = clinical %>%
  mutate(SUBJID = as.numeric(SUBJID))

antibody = antibody %>%
  mutate(SUBJID = as.numeric(SUBJID))

treatment = treatment %>%
  mutate(SUBJID = as.numeric(SUBJID))

# ---- Recode Time variable to standardized format ----
clinical <- clinical %>%
  mutate(
    Time = recode_factor(
      Time,
      "J 0"     = "P+0",
      "P + 3h"  = "P+3H",
      "P + 1 J" = "P+1D",
      "P + 7 J" = "P+7D",
      "P + 29 J" = "P+29D",
      "P + 57 J" = "P+57D",
      "P + 85 j" = "P+85D",
      "B + 3h"  = "B+3H",
      "B + 1 J" = "B+1D",
      "B + 7 J" = "B+7D"
    )
  )

# ---- Merge clinical data with treatment information ----
clinical2 <- clinical %>%
  left_join(treatment, by = "SUBJID") %>%
  arrange(Sample)  # reorder by sample

# ---- Create simplified arm variable ----
clinical2 <- clinical2 %>%
  mutate(arm = recode_factor(
    as.factor(TRTSEQP),
    "Ad26, MVA"       = "active",
    "Placebo, Placebo" = "placebo"
  ))

# ---- Harmonize antibody identifiers ----
# Ensure numeric SUBJID to match clinical/treatment data
antibody <- antibody %>%
  mutate(SUBJID = as.numeric(SUBJID))

# ---- Engineer Time variable for antibody measurements ----
# Map planned visit days (VISITDY) to standardized time points
antibody <- antibody %>%
  mutate(
    Time = case_when(
      VISITDY %in% c(50, 78, 106) ~ "B+21D",
      VISITDY %in% c(29, 57, 85)  ~ "B+0D",
      VISITDY == 365              ~ "P+364D",
      VISITDY == 1                ~ "P+0D",
      TRUE                        ~ NA_character_
    )
  ) %>%
  mutate(Time = as.factor(Time))


# ---- Filter antibody measurements ----
# Keep only rows with valid Time values
antibody <- antibody %>%
  filter(!is.na(Time))

# Keep only patients with clinical data
clinical_pid <- unique(clinical2$SUBJID)
antibody <- antibody %>%
  filter(SUBJID %in% clinical_pid)

# Keep only relevant assays: ELISA at Dose 1 with positive analysis flag
antibody <- antibody %>%
  filter(PARAMCD == "ELISA", BASETYPE == "DOSE 1", ANL01FL == "Y")

# ---- Add log-transformed antibody values ----
antibody <- antibody %>%
  mutate(log10AVAL = log10(AVAL))

# ---- Pivot antibody data to wide format ----
antibody2 <- antibody %>%
  dplyr::select(SUBJID, Time, AVAL) %>%
  pivot_wider(names_from  = Time,
              values_from = AVAL,
              id_cols     = SUBJID)

# Add prefix to antibody columns
colnames(antibody2)[-1] <- paste0("Ab_", colnames(antibody2)[-1])

# ---- Merge clinical and antibody data ----
df_merged <- clinical2 %>%
  inner_join(antibody2, by = "SUBJID") %>%
  arrange(Sample) %>%
  mutate(
    RACE = recode(
      RACE,
      "WHITE" = "White",
      "OTHER" = "Unknown",
      "BLACK OR AFRICAN AMERICAN" = "Black or African American"
    ),
    ETHNIC = recode(ETHNIC,
                    "NOT HISPANIC OR LATINO" = "Other",
                    "HISPANIC OR LATINO" = "Hispanic or Latino")
  )

# ---- Prepare transcriptomic counts ----
# Remove first 3 columns (non-gene identifiers)
counts <- counts[, -c(1:3)]

# Add explicit sample identifier column
counts_raw <- counts %>%
  tibble::rownames_to_column(var = "Sample")

# Normalize counts: log2 CPM-like transformation
counts_norm <- apply(counts, 1, function(v) {
  log2((v + 0.5) / (sum(v) + 1) * 1e6)
}) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample")

# Ensure Sample is character
counts_norm$Sample <- as.numeric(counts_norm$Sample)
counts_raw$Sample <- as.numeric(counts_raw$Sample)

# ---- Merge counts with clinical + antibody data ----
ebovac_merged_raw <- df_merged %>%
  inner_join(counts_raw, by = "Sample")

ebovac_merged_norm <- df_merged %>%
  inner_join(counts_norm, by = "Sample")

# ---- Select and rename final columns ----
final_cols <- c(
  "Sample",
  "SUBJID",
  "Time",
  "AGE",
  "SEX",
  "RACE",
  "ETHNIC",
  "arm",
  "ACTARM",
  "BMICALC1",
  "Ab_P+0D",
  "Ab_P+364D",
  "Ab_B+0D",
  "Ab_B+21D"
)

ebovac_merged_raw <- ebovac_merged_raw %>%
  dplyr::select(all_of(final_cols), A1CF:ZZZ3) %>%
  rename(
    PID       = SUBJID,
    Age       = AGE,
    Sex       = SEX,
    Race      = RACE,
    Ethnicity = ETHNIC,
    Arm       = arm,
    BMI       = BMICALC1
  ) %>%
  janitor::clean_names()

ebovac_merged_norm <- ebovac_merged_norm %>%
  dplyr::select(all_of(final_cols), A1CF:ZZZ3) %>%
  rename(
    PID       = SUBJID,
    Age       = AGE,
    Sex       = SEX,
    Race      = RACE,
    Ethnicity = ETHNIC,
    Arm       = arm,
    BMI       = BMICALC1
  ) %>%
  janitor::clean_names()

ebovac2_clinical = ebovac_merged_norm %>%
  dplyr::select(-c(a1cf:zzz3))

# ---- Save processed datasets ----
processed_data_path <- fs::path("data", "ebovac2")

saveRDS(ebovac_merged_raw,
     file = fs::path(processed_data_path, "ebovac2_merged_raw.rds"))

saveRDS(
  ebovac_merged_norm,
  file = fs::path(processed_data_path, "ebovac2_merged_norm.rds")
)

saveRDS(
  ebovac2_clinical,
  file = fs::path(processed_data_path, "ebovac2_clinical.rds")
)

rm(list = ls())
