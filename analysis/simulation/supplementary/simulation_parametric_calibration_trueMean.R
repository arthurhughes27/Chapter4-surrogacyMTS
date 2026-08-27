# Script to investigate the impact of the meta-analysis model (random versus fixed effects)
# And the estimation procedure for the variance of the pooled effect on the test calibration

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
M_grid <- c(25)

# Grid of number of samples within each study to generate
nm_grid = c(250)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon / 10, epsilon, epsilon * 10)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon / 100 , epsilon / 10, epsilon, epsilon * 10, epsilon *
                     100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.05)

# Grid of meta-analysis model specification parameters
model_specification_grid = c("RE")

# Grid of meta-analysis estimation parameters
test_grid = c("knha")

# Grid of invalid surrogate means
invalid_mean_discrete_grid = c(seq(epsilon, 1, 0.05))

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  model_specification = model_specification_grid,
  test = test_grid,
  invalid_mean_discrete = invalid_mean_discrete_grid,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
seeds <- sample.int(n = .Machine$integer.max, size = nrow(results))

# Add structure to skip for loop if results already exist
output_file <- fs::path(simulation_results_folder,
                        "simulation_parametric_calibration_trueMean.rds")

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
    invalid_mean_discrete = c(-results$invalid_mean_discrete[i],
                              results$invalid_mean_discrete[i])
    
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
      invalid_at_boundary = FALSE,
      invalid_mean_discrete = invalid_mean_discrete,
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
    
    # Derive FPR for each value of alpha
    fpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
    
    # Store results
    results_list[[i]] <- tibble(
      M = M,
      nm = results$nm[i],
      u_tau_max = u_tau_max,
      u_nu_max = u_nu_max,
      test = test,
      meta.analysis.method = model_specification,
      invalid_mean_discrete = invalid_mean_discrete,
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
      "simulation_parametric_calibration_trueMean.rds"
    )
  )
} else {
  results = readRDS(
    fs::path(
      simulation_results_folder,
      "simulation_parametric_calibration_trueMean.rds"
    )
  )
}

fixed_alpha <- 0.05
u_nu_max_fixed <- epsilon / 10

# Helper for legend labels relative to epsilon
rel_eps_text <- function(x, eps = epsilon) {
  x <- as.numeric(x)
  ratio <- x / eps
  
  out <- character(length(ratio))
  out[abs(ratio - 0.1) < 1e-8]  <- "epsilon/10"
  out[abs(ratio - 1)   < 1e-8]  <- "epsilon"
  out[abs(ratio - 10)  < 1e-8]  <- "10*epsilon"
  
  out[out == ""] <- as.character(x[out == ""])
  out
}

plot_df <- results %>%
  filter(alpha == fixed_alpha, u_nu_max == u_nu_max_fixed) %>%
  mutate(u_tau_max = factor(u_tau_max, levels = sort(unique(u_tau_max)))) %>%
  arrange(u_tau_max, invalid_mean_discrete)

ylim <- 0.075

# Split the data so lines do not connect across [-epsilon, epsilon]
plot_left <- plot_df %>%
  filter(invalid_mean_discrete <= -epsilon)

plot_right <- plot_df %>%
  filter(invalid_mean_discrete >= epsilon)

# Custom x-axis breaks and labels
x_breaks <- sort(unique(c(
  pretty(plot_df$invalid_mean_discrete, n = 6), -epsilon, epsilon
)))

x_labels <- function(x) {
  out <- character(length(x))
  out[abs(x + epsilon) < 1e-8] <- parse(text = "-epsilon")
  out[abs(x - epsilon) < 1e-8] <- parse(text = "epsilon")
  out[!(abs(x + epsilon) < 1e-8 |
          abs(x - epsilon) < 1e-8)] <- as.character(x[!(abs(x + epsilon) < 1e-8 |
                                                          abs(x - epsilon) < 1e-8)])
  out
}

p2_cols <- rev(c("#7BD64A", "#0B8500", "#003309"))

names(p2_cols) <- levels(plot_df$u_tau_max)

p1 <- ggplot(plot_df,
             aes(
               x = invalid_mean_discrete,
               y = fpr,
               color = u_tau_max,
               group = u_tau_max
             )) +
  geom_rect(
    inherit.aes = FALSE,
    aes(
      xmin = -epsilon,
      xmax = epsilon,
      ymin = -Inf,
      ymax = Inf,
      fill = "Valid region"
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_line(data = plot_left,
            linewidth = 1.2,
            alpha = 0.6) +
  geom_line(data = plot_right,
            linewidth = 1.2,
            alpha = 0.6) +
  geom_point(size = 2, alpha = 0.5) +
  geom_hline(
    yintercept = fixed_alpha,
    linetype = "dashed",
    color = "black",
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = c(-epsilon, epsilon),
    linetype = "dashed",
    color = "black",
    alpha = 0.7
  ) +
  scale_color_manual(
    values = p2_cols,
    name = expression(paste("Max between-study variance")),
    labels = function(x)
      parse(text = rel_eps_text(x))
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Valid region" = "#C1C1C1"),
    labels = c(expression(
      paste("Valid region: ", -epsilon, " to ", epsilon)
    ))
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = function(x) {
      lab <- rep("", length(x))
      lab[abs(x + epsilon) < 1e-8] <- "-epsilon"
      lab[abs(x - epsilon) < 1e-8] <- "epsilon"
      lab[!(abs(x + epsilon) < 1e-8 | abs(x - epsilon) < 1e-8)] <-
        as.character(x[!(abs(x + epsilon) < 1e-8 |
                           abs(x - epsilon) < 1e-8)])
      parse(text = lab)
    }
  ) +
  coord_cartesian(ylim = c(0, ylim)) +
  labs(x = expression(paste("True mean of invalid surrogates ", mu)),
       y = "Empirical FPR",
       title = "Empirical false positive rates across invalid surrogate mean") +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(size = 35, hjust = 0.5),
    panel.spacing = unit(3, "lines"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(size = 25),
    axis.text.x = element_text(hjust = 1),
    axis.title = element_text(size = 35),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2, override.aes = list(fill = NA, alpha = 1)))

p1

ggsave(
  filename = "fpr_invalidMean.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 45,
  height = 18,
  units = "cm"
)

rm(list = ls())
