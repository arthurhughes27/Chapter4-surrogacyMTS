# Descriptive summaries and figures for the HIPC influenza vaccine cohorts
library(tidyverse)

# Paths to processed data and output figures
processed_data_folder <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# Load merged gene expression dataset
hipc_merged_all_noNorm <- readRDS(fs::path(processed_data_folder, "df_merged_all.rds"))
raw_response_influenzain = readRDS(fs::path(
  processed_data_folder,
  "raw_response_influenzain_all_noNorm.rds"
))

# Gene columns present in the data with no missing values
gene_names <- hipc_merged_all_noNorm %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where( ~ !any(is.na(.)))) %>%
  colnames()

# Inclusion criteria: inactivated influenza vaccine, MFC available, observations
# at days 0-7 only, exactly one baseline row and at least one post-vaccination row
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(
    !is.na(ab_p_0),
    !is.na(ab_p_28),
    group_long == "Influenza (IN)",
    time %in% c("P+0D", "P+1D", "P+2D", "P+3D", "P+7D")
  ) %>%
  group_by(participant_id) %>%
  filter(sum(time == "P+0D") == 1,
         sum(time != "P+0D") > 0) %>%
  ungroup() %>%
  dplyr::select(
    participant_id,
    study_accession,
    arm_accession,
    age_imputed,
    time,
    ab_p_0,
    ab_p_28,
    gender,
    all_of(gene_names)
  ) %>%
  mutate(gender = droplevels(gender)) %>%
  arrange(participant_id)

# Canonical study order used consistently across all plots
study_order <- c(
  "SDY61",
  "SDY269",
  "SDY270",
  "SDY180",
  "SDY56",
  "SDY67",
  "SDY1119",
  "SDY224",
  "SDY404",
  "SDY63",
  "SDY400",
  "SDY1276",
  "SDY520",
  "SDY80",
  "SDY640"
)

# Apply canonical study ordering
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm_filtered %>%
  mutate(study_accession = factor(study_accession, levels = study_order))

# ===== Bubble plot: transcriptomic sample counts by study and timepoint =====

# Count observations per study and timepoint
counts <- hipc_merged_all_noNorm_filtered %>%
  filter(!is.na(time)) %>%
  group_by(study_accession, time) %>%
  summarise(n = n(), .groups = "drop")

# Convert timepoint to an ordered factor for correct x-axis ordering
time_levels <- counts %>%
  distinct(time) %>%
  arrange(as.numeric(time)) %>%
  pull(time) %>%
  as.character()

counts <- counts %>%
  mutate(
    time = factor(
      as.character(time),
      levels = time_levels,
      ordered = TRUE
    ),
    # Reverse study order so the first study appears at the top of the y-axis
    study_accession = factor(study_accession, levels = rev(levels(study_accession)))
  )

# Bubble size is proportional to raw count; labels show the count inside each bubble
p1 <- ggplot(counts, aes(x = time, y = study_accession)) +
  geom_point(
    aes(size = n),
    fill = "#5062FF",
    shape = 21,
    colour = "black",
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5,
    show.legend = TRUE
  ) +
  scale_fill_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x        = "Days post-vaccination",
    y        = "Study",
    title    = "Number of transcriptomic observations across time per study",
    subtitle = "Selected timepoints, healthy adults with baseline and post-vaccination data"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = 17),
    axis.text          = element_text(size = 15),
    plot.title         = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle      = element_text(size = 17, hjust = 0.5)
  )

p1

ggsave(
  filename = "tiv_bubble_plot.pdf",
  path     = descriptive_figures_folder,
  plot     = p1,
  width    = 30,
  height   = 26,
  units    = "cm"
)

rm(list = ls())
