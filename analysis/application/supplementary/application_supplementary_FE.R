# Script to run the main analysis for the paper : Meta-Analytic Evaluation of High-Dimensional Surrogate Markers : Application to vaccinology

# Libraries
library(tidyverse)
library(SurrogateRank)
library(parallel)

# Define global hyperparameters for analysis
hyperparameter_list = list(
  # Hyperparameters for data pre-processing
  tp = "P+1D",
  # Timepoint for gene expression
  screen.fraction = 0.66,
  # Fraction of data for screening
  seed = 10012025,
  # seed for random data splitting
  
  # Hyperparameters to define methodology
  meta.analysis.method = "FE",
  # meta analysis method (random or fixed effects)
  test = "knha",
  # method for variance estimation of pooled effect
  alternative = "two.sided",
  # form of alternative hypothesis
  epsilon.meta.mode = "mean.power",
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
  power.want.s.study = 0.8,
  # within-study power for epsilon
  epsilon.meta = NULL,
  # fixed value for epsilon
  epsilon.study = NULL,
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
  # returns all screening results
  show.pooled.effect = TRUE,
  # show pooled effect in forest plot?
  return.study.similarity.plot = FALSE,
  # return similarity of within-study analyses plot
  return.forest.plot = TRUE,
  # return forest plot for combined marker
  return.fit.plot = TRUE,
  # return fit plot for combined marker
  return.evaluate.results = TRUE,
  # return evaluation results for screening data
  return.screen.plot = TRUE,
  # return screening plot
  return.all.weights = FALSE,
  # return weights for all predictors
  return.all.evaluate = FALSE, 
  
  # Predictor transformation parameters
  aggregation_function = mean,
  # Function defining aggregation from gene to geneset level
  geneset_definition = "BTM",
  # Argument stating the definition of the genesets (options are BTM or BG3M)
  
  
  # Other hyperparameters
  n.cores = parallel::detectCores(all.tests = FALSE, logical = TRUE) / 2,
  # number of cores for parallel computing
  screen.plot.topN = 25,
  # how many predictors to plot
  screen.plot.point.estimate = F,
  
  # Graphical parameters
  screen.plot.width = 40,
  screen.plot.height = 23,
  forest.plot.width = 32,
  forest.plot.height = 15,
  fit.plot.width = 37,
  fit.plot.height = 20
)

file_name_tag = paste0(
  "_timepoint",
  hyperparameter_list$tp,
  "_method",
  hyperparameter_list$meta.analysis.method,
  "_test",
  hyperparameter_list$test,
  "_epsMode",
  hyperparameter_list$epsilon.meta.mode,
  ifelse(
    hyperparameter_list$epsilon.meta.mode == "user",
    paste0("_eps", hyperparameter_list$epsilon.meta),
    paste0("_power", hyperparameter_list$power.want.s.study)
  )
)

# Load internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

# Paths to processed data and output figures
processed_data_folder <- "data"
application_figures_folder <- fs::path("output", "figures", "application", "supplementary")

# Load merged gene expression and GS_list gene set objects
df <- readRDS(fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds"))
GS_list <- readRDS(fs::path(
  processed_data_folder,
  paste0(hyperparameter_list$geneset_definition, "_processed.rds")
))

preprocessed_data_list = preprocess_data(
  df = df,
  tp = hyperparameter_list$tp,
  screen.fraction = hyperparameter_list$screen.fraction,
  seed = hyperparameter_list$seed
)

preprocessed_data_list[["df.full"]]$participant_id %>% unique() %>% length()
preprocessed_data_list[["df.full"]]$study_accession %>% unique() %>% length()

preprocessed_data_list[["df.full"]] %>%
  dplyr::select(participant_id, study_accession) %>%
  distinct() %>%
  group_by(study_accession) %>%
  summarize(n = n())

df_train = preprocessed_data_list[["df.screen"]]
df_test = preprocessed_data_list[["df.evaluate"]]

predictor_names = df_train %>%
  dplyr::select(a1cf:zzz3) %>%
  colnames()

# ----- Screening on training data -----

train_inputs <- extract_rise_inputs(
  df = df_train,
  predictor_names = predictor_names,
  genesets = GS_list[["genesets"]],
  geneset_names = GS_list[["geneset.names.descriptions"]],
  aggregation_function = hyperparameter_list$aggregation_function
)

# Screen for surrogate markers across studies using BH-corrected meta-analysis
rise_screen_result <- rise.screen.meta(
  yone                         = train_inputs$yone,
  yzero                        = train_inputs$yzero,
  sone                         = train_inputs$sone,
  szero                        = train_inputs$szero,
  studyone                     = train_inputs$studyone,
  studyzero                    = train_inputs$studyzero,
  alpha                        = hyperparameter_list$alpha,
  epsilon.meta.mode            = hyperparameter_list$epsilon.meta.mode,
  power.want.s.study           = hyperparameter_list$power.want.s.study,
  epsilon.meta                 = hyperparameter_list$epsilon.meta,
  alternative                  = hyperparameter_list$alternative ,
  paired.all                   = hyperparameter_list$paired.all,
  return.all.screen            = hyperparameter_list$return.all.screen,
  epsilon.study                = hyperparameter_list$epsilon.study,
  p.correction                 = hyperparameter_list$p.correction,
  show.pooled.effect           = hyperparameter_list$show.pooled.effect,
  return.study.similarity.plot = hyperparameter_list$return.study.similarity.plot,
  test                         = hyperparameter_list$test,
  meta.analysis.method         = hyperparameter_list$meta.analysis.method,
  n.cores                      = hyperparameter_list$n.cores,
  screen.plot.topN             = hyperparameter_list$screen.plot.topN,
  screen.plot.point.estimate   = hyperparameter_list$screen.plot.point.estimate,
  return.evaluate.results      = hyperparameter_list$return.evaluate.results,
  return.fit.plot              = hyperparameter_list$return.fit.plot,
  return.forest.plot           = hyperparameter_list$return.forest.plot,
  normalise.weights            = hyperparameter_list$normalise.weights,
  return.screen.plot           = hyperparameter_list$return.screen.plot,
  weight.mode                  = hyperparameter_list$weight.mode,
  return.all.weights           = hyperparameter_list$return.all.weights,
  paired.studies               = hyperparameter_list$paired.studies,
  u.y.hyp                      = hyperparameter_list$u.y.hyp
)

screen_output = extract_rise_outputs(screen_result = rise_screen_result)

# LaTeX table formatting the significant results of the analysis
screen_output$screen_table

# Extract and show the graphics for the screening stage
screen_plot_1 = screen_output$screen_plot
screen_forest_1 = screen_output$screen_forest
screen_fit_1 = screen_output$screen_fit

screen_plot_1
screen_forest_1
screen_fit_1

# Save these graphics with an informative name
ggsave(
  filename = paste0("WebFigure9.pdf"),
  path     = application_figures_folder,
  plot     = screen_plot_1,
  width    = hyperparameter_list$screen.plot.width,
  height   = hyperparameter_list$screen.plot.height,
  units    = "cm"
)

# ----- Evaluation on test data -----
test_inputs <- extract_rise_inputs(
  df_test,
  predictor_names = predictor_names,
  genesets = GS_list[["genesets"]],
  geneset_names = GS_list[["geneset.names.descriptions"]],
  aggregation_function = hyperparameter_list$aggregation_function
)

# Evaluate significant markers from screening on held-out test data
rise_evaluation_result <- rise.evaluate.meta(
  yone                 = test_inputs$yone,
  yzero                = test_inputs$yzero,
  sone                 = test_inputs$sone,
  szero                = test_inputs$szero,
  studyone             = test_inputs$studyone,
  studyzero            = test_inputs$studyzero,
  screening.weights    = rise_screen_result[["screening.weights"]],
  markers              = rise_screen_result[["significant.markers"]],
  alpha                = hyperparameter_list$alpha,
  epsilon.meta         = hyperparameter_list$epsilon.meta,
  alternative          = hyperparameter_list$alternative,
  paired.all           = hyperparameter_list$paired.all,
  epsilon.study        = hyperparameter_list$epsilon.study,
  p.correction         = hyperparameter_list$p.correction,
  show.pooled.effect   = hyperparameter_list$show.pooled.effect,
  test                 = hyperparameter_list$test,
  epsilon.meta.mode    = hyperparameter_list$epsilon.meta.mode,
  power.want.s.study   = hyperparameter_list$power.want.s.study,
  meta.analysis.method = hyperparameter_list$meta.analysis.method,
  return.fit.plot      = hyperparameter_list$return.fit.plot,
  return.all.evaluate  = hyperparameter_list$return.all.evaluate,
  return.forest.plot   = hyperparameter_list$return.forest.plot,
  weight.mode          = hyperparameter_list$weight.mode,
  evaluate.weights     = hyperparameter_list$evaluate.weights,
  paired.studies       = hyperparameter_list$paired.studies,
  n.cores              = hyperparameter_list$n.cores,
  u.y.hyp              = hyperparameter_list$u.y.hyp
)

evaluation_output = extract_rise_outputs(evaluation_result = rise_evaluation_result)

evaluation_output$evaluation_table

# Extract and show the graphics for the screening stage
evaluation_forest_1 = evaluation_output$evaluation_forest
evaluation_fit_1 = evaluation_output$evaluation_fit

evaluation_forest_1
evaluation_fit_1

ggsave(
  filename = paste0("WebFigure10.pdf"),
  path     = application_figures_folder,
  plot     = evaluation_forest_1,
  width    = hyperparameter_list$forest.plot.width,
  height   = hyperparameter_list$forest.plot.height,
  units    = "cm"
)

rm(list = ls())
