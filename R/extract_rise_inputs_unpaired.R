# Extract RISE inputs for an UNPAIRED (two-arm, active vs. control) comparison at a
# single timepoint, e.g. an active vaccine arm compared against a placebo arm within
# the same trial. Unlike extract_rise_inputs() (which splits a single arm's data into
# pre/post rows from the same participants), df_active and df_control here are
# independent sets of participants, both already filtered to the timepoint of interest.
extract_rise_inputs_unpaired <- function(df_active,
                                         df_control,
                                         study_label,
                                         predictor_names,
                                         genesets,
                                         geneset_names,
                                         aggregation_function) {
  sone_raw  <- df_active  %>% dplyr::select(any_of(predictor_names))
  szero_raw <- df_control %>% dplyr::select(any_of(predictor_names))
  yone      <- df_active  %>% pull(response_post)
  yzero     <- df_control %>% pull(response_post)

  if (!is.null(genesets)) {
    sone = aggregate_to_geneset(
      df = sone_raw,
      genesets = genesets,
      geneset_names = geneset_names,
      FUN = aggregation_function
    )
    szero = aggregate_to_geneset(
      df = szero_raw,
      genesets = genesets,
      geneset_names = geneset_names,
      FUN = aggregation_function
    )
  } else {
    sone = sone_raw
    szero = szero_raw
  }

  # Both arms are tagged with the same study_label so rise.screen.meta/rise.evaluate.meta
  # treat this active-vs-control contrast as a single "study" in the meta-analysis
  studyone  <- rep(study_label, nrow(df_active))
  studyzero <- rep(study_label, nrow(df_control))

  list(
    "yone" = yone,
    "yzero" = yzero,
    "sone" = sone,
    "szero" = szero,
    "studyone" = studyone,
    "studyzero" = studyzero
  )
}

# Combine several extract_rise_inputs()/extract_rise_inputs_unpaired() results (each a
# list with yone/yzero/sone/szero/studyone/studyzero) into a single set of inputs
# spanning all study units, suitable for rise.screen.meta()/rise.evaluate.meta().
combine_rise_inputs <- function(input_list) {
  list(
    "yone"      = do.call(c, lapply(input_list, `[[`, "yone")),
    "yzero"     = do.call(c, lapply(input_list, `[[`, "yzero")),
    "sone"      = do.call(rbind, lapply(input_list, `[[`, "sone")),
    "szero"     = do.call(rbind, lapply(input_list, `[[`, "szero")),
    "studyone"  = do.call(c, lapply(input_list, `[[`, "studyone")),
    "studyzero" = do.call(c, lapply(input_list, `[[`, "studyzero"))
  )
}
