# Script to investigate the impact of the meta-analysis model (random versus fixed effects)
# And the estimation procedure for the variance of the pooled effect on the test power

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
simulation_figures_folder = fs::path("output", "figures", "simulation", "supplementary")
simulation_results_folder = fs::path("output", "results", "simulation", "supplementary")

# Number of independent markers to generate trial-level effects for
J = 100000

# Value of epsilon defining the validity region
epsilon <- 0.2

# Grid of number of studies to generate
M_grid <- c(3, 10, 25)

# Grid of number of samples within each study to generate
nm_grid = c(50)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon / 10, epsilon, epsilon * 10)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon / 10, epsilon, epsilon * 10)

# Grid of alpha values for the power plot
alpha_grid <- c(0.05)

# Grid of meta-analysis model specification parameters
model_specification_grid = c("RE")

# Grid of meta-analysis estimation parameters
test_grid = c("knha")

# Grid of invalid surrogate means
valid_mean_discrete_grid = c(seq(-epsilon, epsilon, 0.01))

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  model_specification = model_specification_grid,
  test = test_grid,
  valid_mean_discrete = valid_mean_discrete_grid,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
seeds <- sample.int(n = .Machine$integer.max, size = nrow(results))

# Add structure to skip for loop if results already exist
output_file <- fs::path(simulation_results_folder,
                        "simulation_parametric_power_trueMean.rds")

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
    model_specification = results$model_specification[i]
    test = results$test[i]
    valid_mean_discrete = results$valid_mean_discrete[i]
    
    # --- DATA GENERATION ---
    data <- generate.example.data.highdim.multistudy(
      epsilon = epsilon,
      M = M,
      sample_sizes = sample_sizes,
      J = J,
      prop_valid = 1,
      u_tau_min = 0,
      u_tau_max = u_tau_max,
      u_nu_min = 0,
      u_nu_max = u_nu_max,
      prop_invalid_under = 0.5,
      invalid_at_boundary = FALSE,
      valid_mean_discrete = valid_mean_discrete,
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
        test = test,
        meta.analysis.method = model_specification
      )
      p_vals[j] <- resj$results$p
    }
    
    # Derive tpr for each value of alpha
    tpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
    
    # Store results
    results_list[[i]] <- tibble(
      M = M,
      nm = results$nm[i],
      u_tau_max = u_tau_max,
      u_nu_max = u_nu_max,
      test = test,
      meta.analysis.method = model_specification,
      valid_mean_discrete = valid_mean_discrete,
      alpha = alpha_grid,
      tpr = tpr_vec
    )
  }
  
  # Extract results into a dataframe
  results <- bind_rows(results_list)
  
  # Save results
  saveRDS(
    results,
    file = fs::path(
      simulation_results_folder,
      "simulation_parametric_power_trueMean.rds"
    )
  )
} else {
  results = readRDS(
    fs::path(
      simulation_results_folder,
      "simulation_parametric_power_trueMean.rds"
    )
  )
  
}

u_nu_max_fixed = epsilon / 10
u_tau_max_fixed = epsilon / 10

results = results %>%
  filter(u_nu_max == u_nu_max_fixed, u_tau_max == u_tau_max_fixed)

# Define breaks and labels in terms of epsilon
x_breaks <- c(-epsilon, -epsilon / 2, 0, epsilon / 2, epsilon)
x_labels <- c(
  expression(-epsilon),
  expression(-epsilon / 2),
  "0",
  expression(epsilon / 2),
  expression(epsilon)
)

p2_cols <- c("#FF7575", "#AB000E", "#570010")

p1 <- ggplot(results,
             aes(
               x = valid_mean_discrete,
               y = tpr,
               color = factor(M),
               group = M
             )) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.5) +
  labs(
    title = "Empirical power as a function of distribution mean",
    x = expression(mu),
    y = "Empirical Power",
    color = "Number of studies"
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_color_manual(values = p2_cols, name = "Number of studies") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p1

ggsave(
  filename = "power_validMean.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 30,
  height = 15,
  units = "cm"
)

rm(list = ls())
