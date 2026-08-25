# File to pre-process the IS2 clinical characteristics data for use
## This includes some data engineering on columns (e.g. collapse "unknown" and "Not Specified" gender to same value)
## and also the creation of some data for downstream analysis purposes (e.g. the associating a colour with each vaccine)

# Load necessary packages
library(fs)
library(Biobase)
library(dplyr)
library(forcats)
library(stringr)

# Specify folder within folder root where the raw data lives
raw_data_folder = fs::path("data-raw", "is2")

# Use fs::path() to specify the data path robustly
p_load <- fs::path(raw_data_folder, "all_noNorm_eset.rds")

# Read in the rds file
all_noNorm_eset <- readRDS(p_load)

# Extract the clinical data as a dataframe
is2_clinical = all_noNorm_eset@phenoData@data %>%
  as.data.frame()

# Data engineering
# "collapsing" the values from certain important columns
# For gender and race, "Not Specified" and "Unknown" represent the same thing
# For ethnicity, the values "Not Hispanic or Latino" and "Other" represent the same thing (the only other values are
# "Not Specified" and "Hispanic or Latino").
# In addition, we round numerical time columns to avoid long recurring numbers
is2_clinical <- is2_clinical %>%
  mutate(
    gender    = fct_collapse(factor(gender), Unknown = c("Not Specified", "Unknown")),
    race      = fct_collapse(factor(race), Unknown = c("Not Specified", "Unknown")),
    ethnicity = fct_collapse(factor(ethnicity), Other   = c("Not Hispanic or Latino", "Other")),
    study_time_collected = round(study_time_collected, 2),
    time_post_last_vax = round(time_post_last_vax, 2)
  )

# Abbreviate the vaccine type column
is2_clinical$vaccine_type <- is2_clinical$vaccine_type %>%
  as.factor() %>%
  recode_factor(
    "Conjugate" = "CJ",
    "Inactivated" = "IN",
    "Inactivated/Recombinant protein" = "IN/RP",
    "Live virus" = "LV",
    "Polysaccharide" = "PS",
    "Recombinant Viral Vector" = "RVV",
    "Recombinant protein" = "RP"
  )

# Create vaccine name column by combining pathogen and vaccine type
is2_clinical <- is2_clinical %>%
  mutate(group_long = str_c(pathogen, " (", vaccine_type, ")"))

# Create a shortened vaccine name variable
is2_clinical$group_long_short <- recode(
  is2_clinical$group_long,
  "Ebola (RVV)"           = "Ebola (RVV)",
  "Yellow Fever (LV)"     = "Y.F. (LV)",
  "Smallpox (LV)"         = "Smallpox (LV)",
  "Tuberculosis (RVV)"    = "T.B. (RVV)",
  "Hepatitis A/B (IN/RP)" = "Hep. (IN/RP)",
  "Meningococcus (CJ)"    = "Men. (CJ)",
  "Meningococcus (PS)"    = "Men. (PS)",
  "Malaria (RP)"          = "Malaria (RP)",
  "Influenza (IN)"        = "Inf. (IN)",
  "Influenza (LV)"        = "Inf. (LV)",
  "Pneumococcus (PS)"     = "Pneumo. (PS)",
  "Varicella Zoster (LV)" = "Varicella (LV)",
  "HIV (RVV)"             = "HIV (RVV)"
)

is2_clinical = is2_clinical %>% 
  filter(group_long == "Influenza (IN)")

# Define an ordering for the vaccines (this is for later to make figures consistent)
conditions_order <- c(
  "Tuberculosis (RVV)",
  "Varicella Zoster (LV)",
  "Yellow Fever (LV)",
  "Ebola (RVV)",
  "Hepatitis A/B (IN/RP)",
  "HIV (RVV)",
  "Influenza (IN)",
  "Influenza (LV)",
  "Malaria (RP)",
  "Meningococcus (CJ)",
  "Meningococcus (PS)",
  "Pneumococcus (PS)",
  "Smallpox (LV)"
)

# Assign this order to the vaccine names
is2_clinical <- is2_clinical %>%
  mutate(group_long = factor(group_long, levels = conditions_order))

# Same for the shortened names
conditions_order_short <- c(
  "T.B. (RVV)",
  "Varicella (LV)",
  "Y.F. (LV)",
  "Ebola (RVV)",
  "Hep. (IN/RP)",
  "HIV (RVV)",
  "Inf. (IN)",
  "Inf. (LV)",
  "Malaria (RP)",
  "Men. (CJ)",
  "Men. (PS)",
  "Pneumo. (PS)",
  "Smallpox (LV)"
)


# Assign this order to the shortened vaccine names
is2_clinical <- is2_clinical %>%
  mutate(group_long_short = factor(group_long_short, levels = conditions_order_short))

# Define a colour for each vaccine
## This colour palette was chosen to maximise visual distinctiveness for 13 vaccines
## Using the "iwanthue" tool (https://medialab.github.io/iwanthue/)
color_palette_vaccine = c(
  "#b94a73",
  "#c6aa3c",
  "#6f71d9",
  "#64c46a",
  "#be62c2",
  "#7d973c",
  "#563382",
  "#4ea76e",
  "#bc69b0",
  "#33d4d1",
  "#bb4c41",
  "#6a87d3",
  "#b57736"
)

# Write a helper function to assign the colours
assign_color <- function(group_long) {
  return(color_palette_vaccine[match(is2_clinical$group_long,
                                     levels(is2_clinical$group_long))])
}

# Assign the colours to the vaccine names
is2_clinical$vaccine_colour <-
  assign_color(is2_clinical$group_long)

# Define an ordering for the studies (this is for later to make figures consistent)
study_descriptions = is2_clinical %>%
  dplyr::select(group_long, study_accession) %>%
  distinct() %>%
  arrange(group_long)

# There are three studies which include multiple vaccines
# SDY1260 and SDY1325 for Meningococcus CJ and PS,
# SDY269 for influenza IN and LV,
# SDY180 for influenza IN and Pneumococcus (PS)
# Define a new variable renaming these to "SDYXa" (CJ, IN) and "SDYXb" (PS, LV) where X is replaced with the appropriate study number
is2_clinical = is2_clinical %>%
  mutate(study_accession_unique = study_accession)

is2_clinical = is2_clinical %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1260" &
        group_long == "Meningococcus (CJ)",
      "SDY1260a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1260" &
        group_long == "Meningococcus (PS)",
      "SDY1260b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1325" &
        group_long == "Meningococcus (CJ)",
      "SDY1325a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY1325" &
        group_long == "Meningococcus (PS)",
      "SDY1325b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY269" &
        group_long == "Influenza (IN)",
      "SDY269a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY269" &
        group_long == "Influenza (LV)",
      "SDY269b",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY180" &
        group_long == "Influenza (IN)",
      "SDY180a",
      study_accession_unique
    )
  ) %>%
  mutate(
    study_accession_unique = ifelse(
      study_accession_unique == "SDY180" &
        group_long == "Pneumococcus (PS)",
      "SDY180b",
      study_accession_unique
    )
  )

# Redefine the study order
study_descriptions = is2_clinical %>%
  dplyr::select(group_long, study_accession_unique) %>%
  distinct() %>%
  arrange(group_long)

study_order <- study_descriptions$study_accession_unique

# Assign this order to the vaccine names
is2_clinical <- is2_clinical %>%
  mutate(study_accession_unique = factor(study_accession_unique, levels = study_order))

# Define a colour for each study
## This colour palette was chosen to represent 30 studies across 13 vaccines
## The study colours are chosen to be a gradient within each vaccine based on the previously chosen vaccine colours
## Using the "supercolorpalette" tool (https://supercolorpalette.com)
color_palette_study = c(
  "#b94a73", # Tuberculosis (RVV)
  "#c6aa3c", # Varicella Zoster (LV)
  "#F1D986",
  "#4042FB",  # Yellow Fever (LV)
  "#6077FC",
  "#889FFD",
  "#B4C4FD",
  "#64c46a", # Ebola (RVV)
  "#be62c2", # Hepatitis A/B (IN/RP)
  "#7d973c", # HIV (RVV)
  "#170931",# Influenza (IN)
  "#230C48",
  "#2E0D5F",
  "#3A0A77",
  "#460090",
  "#5100A5",
  "#5C00BA",
  "#6701D0",
  "#7301E5", #
  "#7E01FB",
  "#8632FF",
  "#8D4DFF",
  "#9562FF",
  "#9D74FF",
  "#A684FF",
  "#AF94FF",
  "#4ea76e", # Influenza (LV)
  "#bc69b0", # Malaria (RP)
  "#33d4d1", # Meningococcus (CJ)
  "#B2FBF8",
  "#bb4c41", # Meningococcus (PS)
  "#F39C90",
  "#6a87d3", # Pneumococcus (PS)
  "#b57736"  # Smallpox (LV)
)

# Write a helper function to assign the colours
assign_color <- function(study_accession_unique) {
  return(color_palette_study[match(
    is2_clinical$study_accession_unique,
    levels(is2_clinical$study_accession_unique)
  )])
}

# Assign the colours to the vaccine names
is2_clinical$study_colour <-
  assign_color(is2_clinical$study_accession_unique)

# Final dataframe to be saved has samples as rows and variables as columns
dim(is2_clinical)

is2_clinical = is2_clinical %>%
  dplyr::select(participant_id, study_time_collected, age_imputed, gender, race, ethnicity, study_accession, group_long) %>%
  distinct()

# Save processed dataframe

# Specify folder within folder root where the processed data lives
processed_data_folder = fs::path("data", "is2")

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "is2_clinical.rds")

# Save dataframe
saveRDS(is2_clinical, file = p_save)

rm(list = ls())
