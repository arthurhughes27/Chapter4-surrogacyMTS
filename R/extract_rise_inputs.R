# Extract RISE inputs (response vectors, gene-set matrices, study labels)
# from a data frame containing both baseline and post-vaccination rows
extract_rise_inputs <- function(df,
                                predictor_names,
                                genesets,
                                geneset_names,
                                aggregation_function) {
  sone_raw  <- df %>% filter(time != "P+0D")  %>% dplyr::select(any_of(predictor_names))
  szero_raw <- df %>% filter(time == "P+0D") %>% dplyr::select(any_of(predictor_names))
  yone      = df %>% filter(time != "P+0D")  %>% pull(response_post)
  yzero     = df %>% filter(time == "P+0D") %>% pull(response_pre)
  
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
  studyone  = df %>% filter(time != "P+0D")  %>% pull(study_accession)
  studyzero = df %>% filter(time == "P+0D") %>% pull(study_accession)
  
  res = list(
    "yone" = yone,
    "yzero" = yzero,
    "sone" = sone,
    "szero" = szero,
    "studyone" = studyone,
    "studyzero" = studyzero
  )
  
  return(res)
}
