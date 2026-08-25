# File to ensure the intended project structure is in place (as some of these folders are not uploaded to github)

# Ensure superfolders are created
dir.create("output", recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

dir.create("data-raw/prevac", recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/ebovac2", recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/hamburg", recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/is2", recursive = TRUE, showWarnings = FALSE)

dir.create("data/prevac", recursive = TRUE, showWarnings = FALSE)
dir.create("data/ebovac2", recursive = TRUE, showWarnings = FALSE)
dir.create("data/hamburg", recursive = TRUE, showWarnings = FALSE)
dir.create("data/is2", recursive = TRUE, showWarnings = FALSE)

# Output substructure
## Figures
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/application", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/application/main", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/application/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/simulation", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/simulation/main", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/simulation/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/descriptive", recursive = TRUE, showWarnings = FALSE)

# Results
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/application", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/application/main", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/application/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/simulation", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/simulation/main", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/simulation/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results/descriptive", recursive = TRUE, showWarnings = FALSE)

rm(list = ls())
