# =============================================================================
# Script to describe antibody (NAb, HAI) and gene expression sample
# availability across timepoints for the influenza (IS2) studies.
# Styled after Descriptive-Ebolavirus.R (Chapter3-surrogacyTLS repo).
# =============================================================================

# ---- Libraries ----
library(tidyverse)   # loads dplyr, tidyr, ggplot2, stringr, etc.
library(fs)
library(patchwork)

# ---- Paths ----
processed_data_path        <- fs::path("data")
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# ---- Load data ----
# Per-assay (nAb / hai) antibody response availability at the EXACT collection
# day (not the +-7 day nominal-timepoint windowing used for the analysis
# proper), and before the chosen-assay-per-study collapsing used elsewhere in
# the pipeline
is2_immResp_by_assay <- readRDS(fs::path(processed_data_path, "is2", "is2_immResp_by_assay_exact.rds"))

# Harmonised clinical data (one row per GE sample) used for the GE panel
df_clinical_all <- readRDS(fs::path(processed_data_path, "df_clinical_all.rds"))

# ---- Shared ordering / labelling helpers ----

# Canonical influenza study order, matching analysis/descriptive/data_description.R
study_order <- c(
  "SDY61", "SDY269", "SDY270", "SDY180", "SDY56", "SDY67", "SDY1119",
  "SDY224", "SDY404", "SDY63", "SDY400", "SDY1276", "SDY520", "SDY80", "SDY640"
)
study_order <- intersect(study_order, unique(is2_immResp_by_assay$study_accession))

# =============================================================================
# Panels A/B: Antibody measurement availability (nAb, HAI) by study and EXACT
# collection day (not the +-7 day nominal timepoint used for the analysis
# proper, e.g. "day 28" there covers measurements from day 21-35). Restricted
# to day <= 35, matching the latest nominal timepoint (28 +- 7 days) used in
# the analysis.
# =============================================================================

max_day <- 35

count_by_assay <- function(assay_name) {
  df <- is2_immResp_by_assay %>%
    filter(
      assay == assay_name,
      study_accession %in% study_order,
      !is.na(response_mean),
      study_time_collected <= max_day
    ) %>%
    # Round to the nearest tenth of a day purely for a manageable, readable
    # x-axis (raw collection days can carry many decimal places) - this is
    # still the exact day, not the +-7 day nominal timepoint used elsewhere
    mutate(study_time_collected = round(study_time_collected, 1))

  day_order <- sort(unique(df$study_time_collected))

  df %>%
    mutate(
      timepoint = factor(study_time_collected, levels = day_order),
      study_accession = factor(study_accession, levels = rev(study_order))
    ) %>%
    group_by(study_accession, timepoint) %>%
    summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
    complete(study_accession, timepoint, fill = list(n_participants = 0))
}

df_counts_nab <- count_by_assay("nAb")
df_counts_hai <- count_by_assay("hai")

label_fun_day <- function(x) paste0("Day ", x)

# =============================================================================
# Panel C: Gene expression sample availability by study and timepoint
# =============================================================================

# Derived from the timepoints actually present for these studies (rather than
# a fixed list), so an unused timepoint (e.g. "P+3H", relevant to other
# vaccines' prime/boost designs but not sampled here) doesn't show up as an
# empty column. Restricted to <= 7 days post-vaccination; negative (baseline)
# timepoints are kept and correctly sorted before day 0.
max_ge_hours <- 7 * 24

# Extracts signed hours from a "P+XD"/"P+XH" label (X may be negative, e.g.
# "P+-7D" for day -7), so timepoints sort chronologically and filter correctly
time_to_hours <- function(x) {
  if (str_detect(x, "H$")) as.numeric(str_extract(x, "-?\\d+(?=H$)"))
  else as.numeric(str_extract(x, "-?\\d+(?=D$)")) * 24
}

ge_time_present <- df_clinical_all %>%
  filter(study_accession %in% study_order, !is.na(time)) %>%
  distinct(time) %>%
  pull(time) %>%
  as.character()

ge_time_present <- ge_time_present[vapply(ge_time_present, time_to_hours, numeric(1)) <= max_ge_hours]
ge_time_order <- ge_time_present[order(vapply(ge_time_present, time_to_hours, numeric(1)))]

df_counts_ge <- df_clinical_all %>%
  filter(study_accession %in% study_order, time %in% ge_time_order) %>%
  mutate(
    timepoint = factor(time, levels = ge_time_order),
    study_accession = factor(study_accession, levels = rev(study_order))
  ) %>%
  group_by(study_accession, timepoint) %>%
  summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
  complete(study_accession, timepoint, fill = list(n_participants = 0))

label_fun_ge <- function(x) {
  x <- as.character(x)
  ifelse(
    str_detect(x, "^P\\+\\d+H$"),
    paste0(str_remove(str_remove(x, "^P\\+"), "H"), " hours"),
    paste0("Day ", str_remove_all(str_remove(x, "^P\\+"), "D"))
  )
}

# =============================================================================
# Combined, harmonised figure
# =============================================================================

# All three panels share one colour scale so tile darkness is directly
# comparable, and one legend, collected by patchwork below.
shared_fill_max <- max(df_counts_nab$n_participants,
                       df_counts_hai$n_participants,
                       df_counts_ge$n_participants)

# Shared heatmap theme/style, reused for all three panels
heatmap_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
    axis.text.y     = element_text(size = 13, face = "bold"),
    axis.title.x    = element_text(size = 14, margin = margin(t = 10)),
    plot.title      = element_text(size = 17, face = "bold", hjust = 0.5),
    panel.grid      = element_blank(),
    plot.margin     = margin(15, 15, 15, 15)
  )

# Build a green heatmap tile + label layer, reused for all three panels, on a
# shared fill scale so the panels are visually comparable. Empty (0-count)
# cells are recoloured light grey with no numeric label, rather than sitting
# at the pale end of the green gradient.
heatmap_layers <- function(df, fill_max) {
  list(
    geom_tile(aes(fill = ifelse(n_participants == 0, NA_real_, n_participants)),
              color = "white", linewidth = 0.8),
    geom_text(
      aes(label = ifelse(n_participants > 0, n_participants, "")),
      size = 4.5, fontface = "bold",
      color = ifelse(df$n_participants > fill_max * 0.55, "white", "grey20")
    ),
    scale_fill_gradient(
      low = "#E8F5E9", high = "#237A21",
      limits = c(0, fill_max),
      na.value = "grey92",
      name = "Samples"
    )
  )
}

p1 <- ggplot(df_counts_nab, aes(x = timepoint, y = study_accession)) +
  heatmap_layers(df_counts_nab, shared_fill_max) +
  scale_x_discrete(labels = label_fun_day) +
  labs(x = "Timepoint", y = NULL, title = "Neutralising antibody (NAb) measurements") +
  heatmap_theme

p2 <- ggplot(df_counts_hai, aes(x = timepoint, y = study_accession)) +
  heatmap_layers(df_counts_hai, shared_fill_max) +
  scale_x_discrete(labels = label_fun_day) +
  labs(x = "Timepoint", y = NULL, title = "Haemagglutination inhibition (HAI) measurements") +
  heatmap_theme

p3 <- ggplot(df_counts_ge, aes(x = timepoint, y = study_accession)) +
  heatmap_layers(df_counts_ge, shared_fill_max) +
  scale_x_discrete(labels = label_fun_ge) +
  labs(x = "Timepoint", y = NULL, title = "Gene expression samples") +
  heatmap_theme

p_combined <- (p1 / p2 / p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Availability of Antibody and Gene Expression Measurements (Influenza Studies)",
    theme = theme(
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5)
    )
  )

p_combined

ggsave(
  filename = "influenza_sample_availability.pdf",
  path = descriptive_figures_folder,
  plot = p_combined,
  width = 35, height = 36, units = "cm"
)

rm(list = ls())
