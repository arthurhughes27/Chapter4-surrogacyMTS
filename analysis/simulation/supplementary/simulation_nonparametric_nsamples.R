# Simulation study with permuted dataset to examine the calibration of the approach
# R script to illustrate a meta-analysis with RISE-Meta

# Libraries
library(tidyverse)
library(SurrogateRank)
library(knitr)
library(kableExtra)
library(parallel)

# Set a global seed
global.seed = 08012025
set.seed(global.seed)

# Directory containing processed data files
processed_data_folder <- "data"

# Folder to store images
results_folder <- fs::path("output", "results", "simulation")
figures_folder <- fs::path("output", "figures", "simulation")
functions_folder <- fs::path("R")

# Source the function to generate a permuted dataset
source(fs::path(functions_folder, "generate_permuted_dataset.R"))

# Paths to real dataset
p_load_expr_all_noNorm <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

# Load data objects
hipc_merged_all_noNorm <- readRDS(p_load_expr_all_noNorm)

# Extract the names of the markers
gene_names <- hipc_merged_all_noNorm %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  colnames()

# Set the nominal significance level
alpha <- 0.05

# Set the post-vaccination timepoint of interest
# In this case, we choose day 7, because we have the most studies at this timepoint
tp <- c(7)

# Keep day 0 and day 7
timepoints_to_keep <- c(0, tp)

# Sample size limit
sample_size_limit <- 25

# Within-study sample size grid
sample_sizes <- seq(5, sample_size_limit, 5)

# Set a label for sample sizes for later plotting
size_labels <- ifelse(is.na(sample_sizes), "Full", as.character(sample_sizes))

# Filter the data to contain only relevant samples
df <- hipc_merged_all_noNorm %>%
  mutate(
    response_pre = 
      ab_p_0,
    response_post = ab_p_28
  ) %>%
  filter(
    !is.na(response_post),
    !is.na(response_pre),
    # Participants must have both pre and post response
    group_long == "Influenza (IN)",
    # Participants receive inactivated influenza vaccine
    time %in% timepoints_to_keep # Samples at chosen timepoints
  ) %>%
  group_by(participant_id) %>% # Check participants have both pre and post measurements
  filter(sum(time == "P+0D") == 1,
         sum(time == tp) == 1) %>%
  ungroup() %>%
  group_by(study_accession) %>%
  filter(length(unique(participant_id)) >= sample_size_limit) %>% # Filter out studies with less than 20 participants
  ungroup() %>%
  dplyr::select(
    # Select relevant columns
    participant_id,
    age_imputed,
    gender,
    race,
    study_accession,
    time,
    response_pre,
    response_post,
    all_of(gene_names)
  ) %>%
  arrange(participant_id) # Fix order

n_studies_max <- df$study_accession %>% # How many total studies
  unique() %>%
  length()

# Number-of-studies grid. In this script, take all studies
n_studies_grid <- c(n_studies_max)

# Set labels for later plotting
n_studies_labels <- ifelse(is.na(n_studies_grid), "All", as.character(n_studies_grid))

# Number of permutation replicates
J = 100000

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
perm_seeds <- sample.int(n = .Machine$integer.max, size = J)

p_save <- fs::path(results_folder, "simulation_nonparametric_nsamples.rds")

# If results already exist, skip the simulation entirely and go straight to plotting
if (file.exists(p_save)) {
  {
    cat("Results file found — loading and skipping simulation.\n")
    fpr_df <- readRDS(p_save)
  }
} else {
  # Flatten the nested grid into a single indexed list for checkpointing
  sim_grid <- expand.grid(
    j = seq_along(n_studies_grid),
    i = seq_along(sample_sizes),
    stringsAsFactors = FALSE
  )
  n_cells <- nrow(sim_grid)
  
  # Hash of all parameters that define the simulation — used to detect whether
  # a checkpoint was produced under different settings and should be discarded.
  sim_params <- list(
    J = J,
    alpha = alpha,
    epsilon_meta = 0.2,
    epsilon_study = 0.2,
    n_studies_grid = n_studies_grid,
    sample_sizes = sample_sizes,
    perm_seeds = perm_seeds
  )
  sim_hash <- digest::digest(sim_params)
  
  # Checkpoint path — allows the run to resume after a crash
  checkpoint_path <- fs::path(results_folder,
                              "simulation_nonparametric_nsamples_checkpoint.rds")
  
  # Save progress during the inner loop every N permutations
  checkpoint_every <- 100L
  
  if (file.exists(checkpoint_path)) {
    {
      checkpoint <- readRDS(checkpoint_path)
      if (!identical(checkpoint$sim_hash, sim_hash)) {
        {
          warning(
            "Checkpoint found but simulation parameters have changed — ",
            "starting fresh. Delete the checkpoint file manually if you want to suppress this message."
          )
          results <- vector("list", n_cells)
          start_row <- 1L
          start_b <- 1L
          current_fpr_vec <- NULL
        }
      } else {
        {
          results <- checkpoint$results
          start_row <- checkpoint$next_row
          start_b <- checkpoint$next_b
          current_fpr_vec <- checkpoint$current_fpr_vec
          cat(
            sprintf(
              "Resuming from checkpoint: grid cell %d / %d, permutation %d / %d\n",
              start_row,
              n_cells,
              start_b,
              J
            )
          )
        }
      }
    }
  } else {
    {
      results <- vector("list", n_cells)
      start_row <- 1L
      start_b <- 1L
      current_fpr_vec <- NULL
    }
  }
  
  for (row in seq_len(n_cells)) {
    # Skip already-completed cells when resuming
    if (row < start_row) {
      next
    }
    
    j <- sim_grid$j[row]
    i <- sim_grid$i[row]
    
    n_st <- n_studies_grid[j]
    n_st_lbl <- n_studies_labels[j]
    n_study <- sample_sizes[i]
    lbl <- size_labels[i]
    
    cat(
      sprintf(
        "\n==== Grid cell %d / %d | n_studies: %s | n_within_study: %s ====\n",
        row,
        n_cells,
        n_st_lbl,
        lbl
      )
    )
    
    # Resume the current row if the checkpoint stopped mid-row
    if (row == start_row &&
        !is.null(current_fpr_vec) && length(current_fpr_vec) == J) {
      {
        fpr_vec <- current_fpr_vec
      }
    } else {
      {
        fpr_vec <- rep(NA_real_, J)
      }
    }
    
    b_begin <- if (row == start_row) {
      start_b
    } else {
      1L
    }
    
    for (b in seq(from = b_begin, to = J)) {
      fpr_vec[b] <- tryCatch({
        df_perm_list <- generate_permuted_dataset(
          df = df,
          gene_names = gene_names,
          seed = perm_seeds[b],
          n_within_study = if (is.na(n_study)) {
            NULL
          } else {
            n_study
          },
          n_studies = if (is.na(n_st)) {
            NULL
          } else {
            n_st
          }
        )
        
        rise_result <- rise.screen.meta(
          df_perm_list$yone,
          df_perm_list$yzero,
          sone = df_perm_list$sone,
          szero = df_perm_list$szero,
          df_perm_list$studyone,
          df_perm_list$studyzero,
          alpha = alpha,
          epsilon.meta = 0.2,
          alternative = "two.sided",
          paired.all = TRUE,
          return.all.screen = TRUE,
          epsilon.study = 0.2,
          p.correction = "none",
          show.pooled.effect = TRUE,
          return.study.similarity.plot = FALSE,
          n.cores = parallel::detectCores(all.tests = FALSE, logical = TRUE) /
            2,
          return.fit.plot = FALSE,
          return.forest.plot = FALSE,
          normalise.weights = FALSE,
          weight.mode = "diff.epsilon",
          return.all.weights = FALSE,
          return.evaluate.results = FALSE
        )
        
        p_vals <- rise_result[["screening.metrics.meta"]]$p.unadjusted
        fpr_b <- mean(p_vals < alpha)
        
        #rm(df_perm_list, rise_result, p_vals)
        gc()
        fpr_b
        
      }, error = function(e) {
        {
          cat(sprintf(
            "  ERROR at permutation %d: %s\n",
            b,
            conditionMessage(e)
          ))
          NA_real_
        }
      })
      
      cat(sprintf("  Permutation %d / %d  |  FPR = %.4f\n", b, J, fpr_vec[b]))
      
      # Save checkpoint every checkpoint_every iterations, and always at the end
      if (b %% checkpoint_every == 0L || b == J) {
        {
          saveRDS(
            list(
              results = results,
              next_row = row,
              next_b = b + 1L,
              current_fpr_vec = fpr_vec,
              sim_hash = sim_hash
            ),
            checkpoint_path
          )
        }
      }
    }
    
    cat(
      sprintf(
        "Mean FPR (n_studies=%s, n_within_study=%s): %.4f  [%d NA(s)]\n",
        n_st_lbl,
        lbl,
        mean(fpr_vec, na.rm = TRUE),
        sum(is.na(fpr_vec))
      )
    )
    
    results[[row]] <- data.frame(fpr = fpr_vec,
                                 n_within_study = lbl,
                                 n_studies = n_st_lbl)
    
    # Save checkpoint after each completed cell so the run can resume after a crash
    saveRDS(
      list(
        results = results,
        next_row = row + 1L,
        next_b = 1L,
        current_fpr_vec = NULL,
        sim_hash = sim_hash
      ),
      checkpoint_path
    )
    
    current_fpr_vec <- NULL
    gc()
  }
  
  # Remove checkpoint now that the full run is complete
  if (file.exists(checkpoint_path)) {
    {
      file.remove(checkpoint_path)
    }
  }
  
  # Convert to long data frame with ordered factors
  fpr_df <- bind_rows(results)
  
  fpr_df <- fpr_df %>%
    mutate(n_within_study = 2 * as.numeric(n_within_study)) %>%
    mutate(n_within_study = factor(n_within_study, levels = 2 * sample_sizes))
  
  # Save dataframe
  saveRDS(fpr_df, file = p_save)
  
} # end if/else results file exists

# Boxplots for each number of studies
p1 <- ggplot(fpr_df, aes(x = n_within_study, y = fpr)) +
  geom_violin() +
  geom_hline(
    yintercept = alpha,
    colour = "black",
    linetype = "dashed",
    alpha = 0.8
  ) +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(
    x = "Number of observations within studies",
    y = "False positive rate",
    title = sprintf(
      "FPR distribution per within-study sample size over %d permutations",
      J
    )
  ) +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.title = element_text(size = 23))

p1

ggsave(
  filename = "WebFigure7.pdf",
  path = figures_folder,
  plot = p1,
  width = 40,
  height = 18,
  units = "cm"
)

rm(list = ls())
