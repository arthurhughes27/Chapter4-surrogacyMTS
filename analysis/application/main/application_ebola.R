# Script to run the main analysis for the Ebola vaccine studies, using the same
# RISE-meta screening/evaluation methodology as the influenza (TIV) application
# (see application_highDim.R), adapted to the Ebola studies' trial designs:
#
# - ebovac2 (Ad26.ZEBOV/MVA-BN-Filo) has no placebo transcriptomic arm available for
#   this analysis, so it is analysed in a PAIRED (pre- vs. post-vaccination,
#   within-participant) manner, exactly as in the influenza application.
# - prevac (rVSV-ZEBOV and Ad26.ZEBOV/MVA-BN-Filo arms) is analysed UNPAIRED: each
#   active vaccine arm is compared against the prevac placebo arm at the same
#   post-vaccination timepoint.
# - hamburg (rVSV) is explicitly excluded.
#
# Each study + vaccine arm combination (e.g. "ebovac2-Ad26MVA", "prevac-rVSV",
# "prevac-Ad26MVA") is treated as its own "study" in the meta-analysis.

# Libraries
library(tidyverse)
library(SurrogateRank)
library(parallel)

# Define global hyperparameters for analysis
hyperparameter_list = list(
  # Hyperparameters for data pre-processing
  tp = "P+7D",
  # Timepoint for gene expression (earliest post-dose timepoint common to both
  # ebovac2 and prevac)
  response_pre_col = "ab_p_0",
  response_post_col = "ab_p_365",
  # Pre/post immune response columns (baseline vs. 1-year post-prime antibody level,
  # the only "far" post-vaccination readout with non-missing data in both studies)
  screen.fraction = 1,
  # Fraction of data for screening
  seed = 10012025,
  # seed for random data splitting

  # Hyperparameters to define methodology
  meta.analysis.method = "RE",
  # meta analysis method (random or fixed effects)
  test = "knha",
  # method for variance estimation of pooled effect
  alternative = "two.sided",
  # form of alternative hypothesis
  epsilon.meta.mode = "mean.power",
  # choice of how to define epsilon
  paired.all = FALSE,
  # paired mode (studies not listed in paired.studies are treated as unpaired)
  paired.studies = c("ebovac2-Ad26MVA"),
  # ebovac2 is analysed paired (pre/post); prevac arms are analysed unpaired
  # (active vs. placebo)
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
  geneset_definition = "BTM",

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
application_figures_folder <- fs::path("output", "figures", "application", "main")

# Load merged gene expression and GS_list gene set objects
df <- readRDS(fs::path(processed_data_folder, "df_merged_all.rds"))

GS_list <- readRDS(fs::path(
  processed_data_folder,
  paste0(hyperparameter_list$geneset_definition, "_processed.rds")
))

# ---- Restrict to the Ebola studies, explicitly excluding hamburg ----
ebola_studies <- c("ebovac2", "prevac")

df_ebola <- df %>%
  filter(study_accession %in% ebola_studies,
         study_accession != "hamburg")

stopifnot(!("hamburg" %in% df_ebola$study_accession))

# Gene columns present with no missing values across the whole Ebola population,
# used consistently across both the paired (ebovac2) and unpaired (prevac) study
# units so their surrogate marker matrices can be combined
gene_names_ebola <- df_ebola %>%
  dplyr::select(a1cf:zzz3) %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  colnames()

prevac_arms    <- c("prevac-rVSV", "prevac-Ad26MVA")
prevac_control <- "prevac-placebo"

# ---- Paired study unit: ebovac2 (pre- vs. post-vaccination) ----
df_ebovac2 <- df_ebola %>%
  filter(study_vaccine == "ebovac2-Ad26MVA")

preprocessed_ebovac2 <- preprocess_data(
  df = df_ebovac2,
  tp = hyperparameter_list$tp,
  study_unit_col = "study_vaccine",
  response_pre_col = hyperparameter_list$response_pre_col,
  response_post_col = hyperparameter_list$response_post_col,
  screen.fraction = hyperparameter_list$screen.fraction,
  seed = hyperparameter_list$seed
)

# ---- Unpaired study units: prevac vaccine arms vs. prevac placebo ----
# Restricted to the post-vaccination timepoint only; no baseline pairing is needed
# for an active-vs-placebo comparison.
set.seed(hyperparameter_list$seed)

df_prevac_tp <- df_ebola %>%
  filter(
    study_vaccine %in% c(prevac_arms, prevac_control),
    time == hyperparameter_list$tp,
    !is.na(.data[[hyperparameter_list$response_post_col]])
  ) %>%
  mutate(response_post = .data[[hyperparameter_list$response_post_col]]) %>%
  dplyr::select(participant_id, study_vaccine, response_post, all_of(gene_names_ebola)) %>%
  arrange(participant_id) %>%
  group_by(study_vaccine) %>%
  filter(n_distinct(participant_id) > 5) %>%
  ungroup()

# Split each arm into screening/evaluation sets, matching the participant-level
# train/test splitting logic used for the paired studies in preprocess_data()
split_by_arm <- function(df, screen.fraction) {
  if (screen.fraction == 1) {
    return(list(train = df, test = df[0, ]))
  } else if (screen.fraction == 0) {
    return(list(train = df[0, ], test = df))
  }
  train_ids <- df %>%
    distinct(study_vaccine, participant_id) %>%
    group_by(study_vaccine) %>%
    slice_sample(prop = screen.fraction) %>%
    ungroup()
  list(
    train = df %>% semi_join(train_ids, by = c("study_vaccine", "participant_id")),
    test  = df %>% anti_join(train_ids, by = c("study_vaccine", "participant_id"))
  )
}

prevac_split <- split_by_arm(df_prevac_tp, hyperparameter_list$screen.fraction)

# Combine the paired ebovac2 inputs with the unpaired prevac-arm-vs-placebo inputs
# for a given split (screening/training or evaluation/test)
build_ebola_inputs <- function(df_ebovac2_split, df_prevac_split_df) {
  input_list <- list()

  input_list[["ebovac2-Ad26MVA"]] <- extract_rise_inputs(
    df = df_ebovac2_split,
    predictor_names = gene_names_ebola,
    genesets = GS_list[["genesets"]],
    geneset_names = GS_list[["geneset.names.descriptions"]],
    aggregation_function = hyperparameter_list$aggregation_function
  )

  df_control <- df_prevac_split_df %>% filter(study_vaccine == prevac_control)

  for (arm in prevac_arms) {
    df_active <- df_prevac_split_df %>% filter(study_vaccine == arm)
    if (nrow(df_active) == 0 || nrow(df_control) == 0) next

    input_list[[arm]] <- extract_rise_inputs_unpaired(
      df_active = df_active,
      df_control = df_control,
      study_label = arm,
      predictor_names = gene_names_ebola,
      genesets = GS_list[["genesets"]],
      geneset_names = GS_list[["geneset.names.descriptions"]],
      aggregation_function = hyperparameter_list$aggregation_function
    )
  }

  combine_rise_inputs(input_list)
}

# Sample sizes per study unit, for reporting
bind_rows(
  preprocessed_ebovac2[["df.full"]] %>%
    dplyr::select(participant_id, study_accession) %>%
    distinct(),
  df_prevac_tp %>%
    dplyr::select(participant_id, study_accession = study_vaccine) %>%
    distinct()
) %>%
  group_by(study_accession) %>%
  summarize(n = n())

# ----- Screening on training data -----

train_inputs <- build_ebola_inputs(
  df_ebovac2_split = preprocessed_ebovac2[["df.screen"]],
  df_prevac_split_df = prevac_split$train
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
  alternative                  = hyperparameter_list$alternative,
  paired.all                   = hyperparameter_list$paired.all,
  paired.studies               = hyperparameter_list$paired.studies,
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

ggsave(
  filename = "risemeta_ebola_screening.pdf",
  path     = application_figures_folder,
  plot     = screen_plot_1,
  width    = hyperparameter_list$screen.plot.width,
  height   = hyperparameter_list$screen.plot.height,
  units    = "cm"
)

# ----- Evaluation on test data -----

test_inputs <- build_ebola_inputs(
  df_ebovac2_split = preprocessed_ebovac2[["df.evaluate"]],
  df_prevac_split_df = prevac_split$test
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
  paired.studies       = hyperparameter_list$paired.studies,
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
  n.cores              = hyperparameter_list$n.cores,
  u.y.hyp              = hyperparameter_list$u.y.hyp
)

evaluation_output = extract_rise_outputs(evaluation_result = rise_evaluation_result)

evaluation_output$evaluation_table

evaluation_forest_1 = evaluation_output$evaluation_forest
evaluation_fit_1 = evaluation_output$evaluation_fit

evaluation_forest_1
evaluation_fit_1

ggsave(
  filename = "risemeta_ebola_evaluation.pdf",
  path     = application_figures_folder,
  plot     = evaluation_forest_1,
  width    = hyperparameter_list$forest.plot.width,
  height   = hyperparameter_list$forest.plot.height,
  units    = "cm"
)

# rm(list = ls())
