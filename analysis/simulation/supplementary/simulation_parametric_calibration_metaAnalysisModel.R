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
epsilon <- 0.1

# Grid of number of studies to generate
M_grid <- c(3, 10, 25)

# Grid of number of samples within each study to generate
nm_grid = c(50)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon *
                      100)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon *
                     100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5)

# Grid of meta-analysis model specification parameters
model_specification_grid = c("RE", "FE")

# Grid of meta-analysis estimation parameters
test_grid = c("knha", "z")

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  model_specification = model_specification_grid,
  test = test_grid,
  stringsAsFactors = FALSE
)

results = results %>%
  filter(!(model_specification == "FE" & test == "knha"))

# Initialise list for storage
results_list <- vector("list", nrow(results))

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
seeds <- sample.int(n = .Machine$integer.max, size = nrow(results))

# Add structure to skip for loop if results already exist
output_file <- fs::path(
  simulation_results_folder,
  "simulation_parametric_calibration_metaAnalysisModel.rds"
)

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
      "simulation_parametric_calibration_metaAnalysisModel.rds"
    )
  )
  
} else {
  results = readRDS(
    fs::path(
      simulation_results_folder,
      "simulation_parametric_calibration_metaAnalysisModel.rds"
    )
  )
}

# Fixed scenario for this figure
fixed_alpha <- 0.05

# Helper: convert numeric values to plotmath text relative to epsilon
rel_eps_text <- function(x, eps = epsilon) {
  x <- as.numeric(x)
  tol <- sqrt(.Machine$double.eps)
  
  out <- rep(NA_character_, length(x))
  
  out[abs(x - eps / 100) < tol] <- "epsilon/100"
  out[abs(x - eps / 10)  < tol] <- "epsilon/10"
  out[abs(x - eps)       < tol] <- "epsilon"
  out[abs(x - 10 * eps)  < tol] <- "10*epsilon"
  out[abs(x - 100 * eps) < tol] <- "100*epsilon"
  
  out[is.na(out)] <- paste0(x[is.na(out)] / eps, "*epsilon")
  out
}

row_labs <- c("FE / z"     = "FE",
              "RE / z"     = "RE / conventional",
              "RE / knha"  = "RE / HKSJ")

plot_df <- results %>%
  filter(alpha == fixed_alpha) %>%
  mutate(
    model_row = paste0(meta.analysis.method, " / ", test),
    model_row = factor(model_row, levels = c("FE / z", "RE / z", "RE / knha")),
    M = factor(M, levels = sort(unique(M))),
    u_tau_max = factor(u_tau_max, levels = sort(unique(u_tau_max))),
    u_nu_max  = factor(u_nu_max, levels = sort(unique(u_nu_max)))
  )

ylim = max(plot_df$fpr)

p2_cols <- rev(c("#6954F2", "#7400E0", "#BA00CF", "#BD0087", "#AB0039"))

names(p2_cols) <- levels(plot_df$u_nu_max)

p1 <- ggplot(plot_df,
             aes(
               x = u_tau_max,
               y = fpr,
               color = u_nu_max,
               group = u_nu_max
             )) +
  geom_line(size = 1.2, alpha = 0.65) +
  geom_point(size = 3, alpha = 0.65) +
  facet_grid(
    rows = vars(model_row),
    cols = vars(M),
    labeller = labeller(
      model_row = as_labeller(row_labs),
      M = function(x)
        paste("M =", x)
    )
  ) +
  scale_x_discrete(
    labels = function(x)
      parse(text = rel_eps_text(x))
  ) +
  scale_color_manual(
    values = p2_cols,
    name = expression("Max within-study variance"),
    labels = function(x)
      parse(text = rel_eps_text(x))
  ) +
  ylim(0, ylim) +
  geom_hline(
    yintercept = fixed_alpha,
    linetype = "dashed",
    color = "black",
    alpha = 0.7
  ) +
  labs(
    x = expression("Maximum between-trial variance"),
    y = "Empirical FPR",
    title = "Empirical false positive rates across meta-analysis model specification and estimation procedures"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    panel.spacing = unit(3, "lines"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 35),
    legend.position = "right"
  )

p1

ggsave(
  filename = "fpr_between_trial_variance.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 50,
  height = 25,
  units = "cm"
)

rm(list = ls())
