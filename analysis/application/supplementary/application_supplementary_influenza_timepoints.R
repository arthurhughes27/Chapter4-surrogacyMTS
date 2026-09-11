# Supplementary analysis: same influenza (TIV) RISE-meta application as
# application_highDim.R, repeated for gene expression at days 2, 3 and 7
# post-vaccination (day 1 is the main analysis).
#
# Difference from the main application: if no significant markers are found
# at the screening stage on the (split) training data, no evaluation stage is
# possible. In that case, screening is re-run on the FULL (unsplit) data
# instead, and those results are interpreted as the final output for that
# timepoint - rather than reporting an evaluation on a held-out set that was
# never reached. See R/run_influenza_timepoint_supplementary.R.

# Libraries
library(tidyverse)
library(SurrogateRank)
library(parallel)

# Define global hyperparameters for analysis (identical to application_highDim.R,
# except tp is set per-timepoint by run_influenza_timepoint_supplementary())
hyperparameter_list = list(
  # Hyperparameters for data pre-processing
  screen.fraction = 0.66,
  # Fraction of data for screening
  seed = 10012025,
  # seed for random data splitting

  # Hyperparameters to define methodology
  meta.analysis.method = "RE",
  # meta analysis method (random or fixed effects)
  test = "knha",
  test.target = "ci",
  # method for variance estimation of pooled effect
  alternative = "two.sided",
  # form of alternative hypothesis
  epsilon.meta.mode = "user",
  # choice of how to define epsilon
  paired.all = TRUE,
  # paired mode
  paired.studies = NULL,
  # which studies are paired
  evaluate.weights = TRUE,
  # Whether to use weighting for evaluation stage

  # Numeric hyperparameters for testing procedure
  alpha = 0.05,
  # significance level
  power.want.s.study = NULL,
  # within-study power for epsilon
  epsilon.meta = 0.2,
  # fixed value for epsilon
  epsilon.study = 0.2,
  # epsilon for within-study testing
  p.correction = "BH",
  # multiplicity correction for p-values
  u.y.hyp = NULL,
  # hypothesised effect size on y
  weight.mode = "diff.epsilon",
  # How to weight surrogates in combination
  normalise.weights = TRUE,
  # normalise weights for the combination

  # Hyperparameters to define which objects to return
  return.all.screen = TRUE,
  show.pooled.effect = TRUE,
  return.study.similarity.plot = FALSE,
  return.forest.plot = TRUE,
  return.fit.plot = TRUE,
  return.evaluate.results = TRUE,
  return.screen.plot = TRUE,
  return.all.weights = FALSE,
  return.all.evaluate = FALSE,

  # Predictor transformation parameters
  aggregation_function = mean,
  # Function defining aggregation from gene to geneset level
  geneset_definition = "BTM",
  # Argument stating the definition of the genesets (options are BTM or BG3M)

  # Other hyperparameters
  n.cores = parallel::detectCores(all.tests = FALSE, logical = TRUE) / 2,
  screen.plot.topN = 20,
  screen.plot.point.estimate = F,

  # Graphical parameters
  screen.plot.width = 40,
  screen.plot.height = 23,
  forest.plot.width = 32,
  forest.plot.height = 15,
  fit.plot.width = 37,
  fit.plot.height = 20
)

# Load internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

# Paths to processed data and output figures
processed_data_folder <- fs::path("data")
application_figures_folder <- fs::path("output", "figures", "application", "supplementary")

# Load merged gene expression and GS_list gene set objects
df <- readRDS(fs::path(processed_data_folder, "df_merged_all.rds"))

GS_list <- readRDS(fs::path(
  processed_data_folder,
  paste0(hyperparameter_list$geneset_definition, "_processed.rds")
))

df_filtered = df %>%
  filter(group_long == "Influenza (IN)")

run_influenza_timepoint_supplementary(
  df_filtered = df_filtered,
  GS_list = GS_list,
  timepoints = c("P+2D", "P+3D", "P+7D"),
  hyperparameter_list = hyperparameter_list,
  application_figures_folder = application_figures_folder
)

rm(list = ls())
