# Script to derive results for parametric simulation of trial-level surrogate effects
# In this script, we generate only invalid surrogates, in order to examine the false positive rate

global.seed = 08012025
set.seed(global.seed)

# Libraries
library(SurrogateRank)
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(tibble)
library(scales)
library(patchwork)

# Folder to store results and figures
simulation_figures_folder = fs::path("output", "figures", "simulation", "main")
simulation_results_folder = fs::path("output", "results", "simulation", "main")

# Number of independent markers to generate trial-level effects for
J = 100000

# Value of epsilon defining the validity region
epsilon <- 0.1

# Grid of number of studies to generate
M_grid <- c(3, 10, 25)

# Grid of number of samples within each study to generate
nm_grid = c(10, 50, 250)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon *
                      100)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon *
                     100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5)

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
seeds <- sample.int(n = .Machine$integer.max, size = nrow(results))

# Add structure to skip for loop if results already exist
output_file <- fs::path(simulation_results_folder,
                        "simulation_parametric_calibration_main.rds")

if (file.exists(output_file)) {
  message("Output already exists: skipping full simulation loop.")
  results <- readRDS(output_file)
  # optionally still load plots or downstream code
  skip_loop <- TRUE
} else {
  skip_loop <- FALSE
}

if (!skip_loop) {
  for (i in seq_len(nrow(results))) {
    M <- results$M[i]
    sample_sizes <- rep(results$nm[i], M)
    u_tau_max <- results$u_tau_max[i]
    u_nu_max <- results$u_nu_max[i]
    
    # --- DATA GENERATION ---
    data <- generate.example.data.highdim.multistudy(
      epsilon = epsilon,
      M = M,
      sample_sizes = sample_sizes,
      J = J,
      prop_valid = 0,
      # 0% valid surrogates
      u_tau_min = 0,
      u_tau_max = u_tau_max,
      u_nu_min = 0,
      u_nu_max = u_nu_max,
      prop_invalid_under = 0.5,
      invalid_at_boundary = TRUE,
      seed = seeds[i]
    )
    
    # Storage for p-values
    p_vals <- numeric(J)
    
    for (j in seq_len(J)) {
      # For each marker
      resj <- delta.reml.meta(
        delta = data$delta[, j],
        sd.delta = data$sd.delta[, j],
        epsilon = epsilon,
        alternative = "two.sided",
        test = "knha",
        meta.analysis.method = "RE"
      )
      p_vals[j] <- resj$results$p
    }
    
    # Derive FPR for each value of alpha
    fpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
    
    # Store results
    results_list[[i]] <- tibble(
      M = M,
      nm = results$nm[i],
      u_tau_max = u_tau_max,
      u_nu_max = u_nu_max,
      alpha = alpha_grid,
      fpr = fpr_vec
    )
  }
  # Extract results into a dataframe
  results <- bind_rows(results_list)
  
  # Save results
  saveRDS(
    results,
    file = fs::path(
      simulation_results_folder,
      "simulation_parametric_calibration_main.rds"
    )
  )
} else {
  results = readRDS(
    fs::path(
      simulation_results_folder,
      "simulation_parametric_calibration_main.rds"
    )
  )
  
}

# Plot 1: Grid of calibration plots for number of studies and sample sizes

fixed_u_tau_max = epsilon / 10
fixed_u_nu_max = epsilon * 100

# Keep only the fixed heterogeneity setting
results_plot <- results %>%
  filter(u_tau_max == fixed_u_tau_max, u_nu_max == fixed_u_nu_max) %>%
  mutate(
    # Pre-format labels as plotmath expressions
    M  = factor(
      M,
      levels = sort(unique(M)),
      labels = paste0("M == ", sort(unique(M)))
    ),
    nm = factor(
      nm,
      levels = sort(unique(nm)),
      labels = paste0("n[m] == ", sort(unique(nm)))
    )
  )

# Define plotting limits
min_val <- 1e-5
max_val <- 0.5

# Custom percent labels
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

# Gradient palette for p1: stronger blue gradient (avoids near-white)
p1_cols <- c("#9EC5E0", "#4D74B8", "#00158F")
names(p1_cols) <- levels(results_plot$nm)

# Upper plotting bound: just above 0.5, but not enough to change the label display
y_upper <- max(0.5, max(results_plot$fpr, na.rm = TRUE))

p1 <- ggplot(results_plot, aes(
  x = alpha,
  y = fpr,
  colour = nm,
  group = nm
)) +
  geom_ribbon(
    data = tibble(x = c(min_val, 0.5)),
    aes(x = x, ymin = x, ymax = 0.5),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(size = 3, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.6) +
  scale_x_log10(
    breaks = c(alpha_grid, 0.5),
    labels = fmt_percent
  ) +
  scale_y_log10(
    breaks = c(alpha_grid, 0.5),
    labels = fmt_percent
  ) +
  coord_cartesian(
    xlim = c(min_val, 0.5),
    ylim = c(min_val, y_upper)
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 2,
    alpha = 0.5
  ) +
  labs(
    title = "Calibration plots across number of studies and sample sizes",
    subtitle = expression(
      "Max between-study variability = " * epsilon / 10 *
        ", max within-study variability = " * 100 * epsilon
    ),
    x = "Nominal FPR",
    y = "Empirical FPR",
    colour = NULL
  ) +
  facet_grid(. ~ M, labeller = label_parsed) +
  scale_colour_manual(values = p1_cols, labels = scales::parse_format()) +
  guides(colour = guide_legend(override.aes = list(
    linewidth = 3,
    size = 6,
    alpha = 0.6
  ))) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(size = 55),
    plot.title = element_text(size = 100, hjust = 0.5),
    plot.subtitle = element_text(size = 60, hjust = 0.5),
    axis.title = element_text(size = 75),
    axis.text = element_text(size = 35),
    legend.text = element_text(size = 55),
    legend.key.width = unit(3, "cm")
  )

p1


# Plot 2: Grid of calibration plots for between- and within-study variability
fixed_M  = 10
fixed_nm = 10

tau_vals <- c(epsilon / 10, epsilon, epsilon * 10)
nu_vals  <- c(epsilon, epsilon * 10, epsilon * 100)

# Prepare data + create parsed labels
results_plot2 <- results %>%
  filter(M == fixed_M,
         nm == fixed_nm,
         u_tau_max %in% tau_vals,
         u_nu_max %in% nu_vals) %>%
  mutate(
    u_tau_lab = factor(
      case_when(
        dplyr::near(u_tau_max, epsilon / 10)  ~ "u[tau^2*','*max] == epsilon/10",
        dplyr::near(u_tau_max, epsilon)       ~ "u[tau^2*','*max] == epsilon",
        dplyr::near(u_tau_max, epsilon * 10)  ~ "u[tau^2*','*max] == 10*epsilon"
      ),
      levels = c(
        "u[tau^2*','*max] == epsilon/10",
        "u[tau^2*','*max] == epsilon",
        "u[tau^2*','*max] == 10*epsilon"
      )
    ),
    u_nu_lab = factor(
      case_when(
        dplyr::near(u_nu_max / fixed_nm, epsilon / 10)  ~ "u[nu*','*max]/n == epsilon/10",
        dplyr::near(u_nu_max / fixed_nm, epsilon)       ~ "u[nu*','*max]/n == epsilon",
        dplyr::near(u_nu_max / fixed_nm, epsilon * 10)  ~ "u[nu*','*max]/n == 10*epsilon"
      ),
      levels = c(
        "u[nu*','*max]/n == epsilon/10",
        "u[nu*','*max]/n == epsilon",
        "u[nu*','*max]/n == 10*epsilon"
      )
    )
  )

# Gradient palette for p2
p2_cols <- rev(c(
  # "#6954F2",
  "#7400E0",
  "#BA00CF",
  "#BD0087"#,
  # "#AB0039"
))
names(p2_cols) <- levels(results_plot2$u_nu_lab)

# Plot
p2 <- ggplot(results_plot2,
             aes(
               x = alpha,
               y = fpr,
               colour = u_nu_lab,
               group = u_nu_lab
             )) +
  geom_ribbon(
    data = tibble(x = c(min_val, max_val)),
    aes(x = x, ymin = x, ymax = max_val),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(size = 3, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.6) +
  scale_x_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  scale_y_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 2,
    alpha = 0.5
  ) +
  labs(
    title = "Calibration plots across between- and within-study variability",
    subtitle = paste0("Number of studies = number of samples = ", fixed_M),
    x = "Nominal FPR",
    y = "Empirical FPR",
    colour = NULL
  ) +
  facet_grid(. ~ u_tau_lab, labeller = label_parsed) +
  scale_colour_manual(values = p2_cols, labels = scales::parse_format()) +
  guides(colour = guide_legend(override.aes = list(
    linewidth = 3,
    size = 6,
    alpha = 0.6
  ))) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(size = 55),
    plot.title = element_text(size = 100, hjust = 0.5),
    plot.subtitle = element_text(size = 60, hjust = 0.5),
    axis.title = element_text(size = 75),
    axis.text = element_text(size = 35),
    legend.text = element_text(size = 55),
    legend.key.width = unit(3, "cm")
  )

p2

combined_plot <- (p1 / plot_spacer() / p2) +
  plot_layout(heights = c(1, 0.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 100),
    plot.tag.position = c(0.02, 0.98)
  )

combined_plot

ggsave(
  filename = "Figure1.pdf",
  path = simulation_figures_folder,
  plot = combined_plot,
  width = 150,
  height = 90,
  # adjust if needed
  units = "cm",
  limitsize = FALSE
)

rm(list = ls())
