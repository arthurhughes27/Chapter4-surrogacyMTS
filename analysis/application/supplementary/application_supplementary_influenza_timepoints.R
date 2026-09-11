# Supplementary analysis: same influenza (TIV) RISE-meta application as
# application_highDim.R, repeated for gene expression at days 2, 3 and 7
# post-vaccination (day 1 is the main analysis).
#
# Difference from the main application: if no significant markers are found
# at the screening stage on the (split) training data, no evaluation stage is
# possible. In that case, screening is re-run on the FULL (unsplit) data
# instead, and those results are interpreted as the final output for that
# timepoint - rather than reporting an evaluation on a held-out set that was
# never reached.

# Libraries
library(tidyverse)
library(SurrogateRank)
library(parallel)

# Define global hyperparameters for analysis (identical to application_highDim.R,
# except tp is set per-timepoint below)
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

# Runs the RISE-meta screening stage with the given inputs, using the shared
# hyperparameters above (tp-independent, reused for both the split-data and
# full-data screening calls below)
run_screen <- function(inputs, hp) {
  rise.screen.meta(
    yone                          = inputs$yone,
    yzero                         = inputs$yzero,
    sone                          = inputs$sone,
    szero                         = inputs$szero,
    studyone                      = inputs$studyone,
    studyzero                     = inputs$studyzero,
    alpha                         = hp$alpha,
    epsilon.meta.mode             = hp$epsilon.meta.mode,
    power.want.s.study            = hp$power.want.s.study,
    epsilon.meta                  = hp$epsilon.meta,
    alternative                   = hp$alternative,
    paired.all                    = hp$paired.all,
    return.all.screen             = hp$return.all.screen,
    epsilon.study                 = hp$epsilon.study,
    p.correction                  = hp$p.correction,
    show.pooled.effect            = hp$show.pooled.effect,
    return.study.similarity.plot  = hp$return.study.similarity.plot,
    test                          = hp$test,
    test.target                   = hp$test.target,
    meta.analysis.method          = hp$meta.analysis.method,
    n.cores                       = hp$n.cores,
    screen.plot.topN              = hp$screen.plot.topN,
    screen.plot.point.estimate    = hp$screen.plot.point.estimate,
    return.evaluate.results       = hp$return.evaluate.results,
    return.fit.plot               = hp$return.fit.plot,
    return.forest.plot            = hp$return.forest.plot,
    normalise.weights             = hp$normalise.weights,
    return.screen.plot            = hp$return.screen.plot,
    weight.mode                   = hp$weight.mode,
    return.all.weights            = hp$return.all.weights,
    paired.studies                = hp$paired.studies,
    u.y.hyp                       = hp$u.y.hyp
  )
}

timepoints <- c("P+2D", "P+3D", "P+7D")

for (tp in timepoints) {
  tp_tag <- str_remove(tp, "^P\\+")  # e.g. "2D", "3D", "7D"

  preprocessed_data_list = preprocess_data(
    df = df_filtered,
    tp = tp,
    screen.fraction = hyperparameter_list$screen.fraction,
    seed = hyperparameter_list$seed
  )

  preprocessed_data_list[["df.full"]] %>%
    dplyr::select(participant_id, study_accession) %>%
    distinct() %>%
    group_by(study_accession) %>%
    summarize(n = n())

  df_train = preprocessed_data_list[["df.screen"]]
  df_test  = preprocessed_data_list[["df.evaluate"]]

  predictor_names = df_train %>%
    dplyr::select(a1cf:zzz3) %>%
    colnames()

  # ----- Screening on (split) training data -----

  train_inputs <- extract_rise_inputs(
    df = df_train,
    predictor_names = predictor_names,
    genesets = GS_list[["genesets"]],
    geneset_names = GS_list[["geneset.names.descriptions"]],
    aggregation_function = hyperparameter_list$aggregation_function
  )

  rise_screen_result <- run_screen(train_inputs, hyperparameter_list)

  n_significant <- length(rise_screen_result[["significant.markers"]])

  if (n_significant > 0) {
    # ----- Significant markers found: proceed exactly as in the main
    # application, evaluating on the held-out test data -----
    screen_output = extract_rise_outputs(screen_result = rise_screen_result)

    ggsave(
      filename = paste0("risemeta-tiv-screening-day", tp_tag, ".pdf"),
      path     = application_figures_folder,
      plot     = screen_output$screen_plot,
      width    = hyperparameter_list$screen.plot.width,
      height   = hyperparameter_list$screen.plot.height,
      units    = "cm"
    )

    test_inputs <- extract_rise_inputs(
      df_test,
      predictor_names = predictor_names,
      genesets = GS_list[["genesets"]],
      geneset_names = GS_list[["geneset.names.descriptions"]],
      aggregation_function = hyperparameter_list$aggregation_function
    )

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
      test.target          = hyperparameter_list$test.target,
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

    ggsave(
      filename = paste0("risemeta-tiv-evaluation-day", tp_tag, ".pdf"),
      path     = application_figures_folder,
      plot     = evaluation_output$evaluation_forest,
      width    = hyperparameter_list$forest.plot.width,
      height   = hyperparameter_list$forest.plot.height,
      units    = "cm"
    )

  } else {
    # ----- No significant markers on the split training data: an evaluation
    # stage isn't possible. Re-run screening on the FULL (unsplit) data
    # instead, and interpret those results as the final output for this
    # timepoint. -----
    preprocessed_data_list_full = preprocess_data(
      df = df_filtered,
      tp = tp,
      screen.fraction = 1,
      seed = hyperparameter_list$seed
    )

    df_full = preprocessed_data_list_full[["df.screen"]]

    full_inputs <- extract_rise_inputs(
      df = df_full,
      predictor_names = predictor_names,
      genesets = GS_list[["genesets"]],
      geneset_names = GS_list[["geneset.names.descriptions"]],
      aggregation_function = hyperparameter_list$aggregation_function
    )

    rise_screen_result_full <- run_screen(full_inputs, hyperparameter_list)

    screen_output_full = extract_rise_outputs(screen_result = rise_screen_result_full)

    ggsave(
      filename = paste0("risemeta-tiv-screening-fulldata-day", tp_tag, ".pdf"),
      path     = application_figures_folder,
      plot     = screen_output_full$screen_plot,
      width    = hyperparameter_list$screen.plot.width,
      height   = hyperparameter_list$screen.plot.height,
      units    = "cm"
    )
  }
}

rm(list = ls())
