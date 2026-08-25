# Function to pre-process the IS2 dataframe to extract valid participants for a given timepoint,
# and performing sample splitting if desired
preprocess_data = function(df,
                           tp,
                           screen.fraction = 1,
                           seed = 12345) {
  set.seed(seed)
  
  # Gene columns present in the data with no missing values
  gene_names <- df %>%
    dplyr::select(a1cf:zzz3) %>%
    dplyr::select(where( ~ !any(is.na(.)))) %>%
    colnames()
  
  # Define the timepoints to keep
  timepoints_to_keep <- c("P+0D", tp)
  
  # Define the studies to use neutralising antibody response
  nab_studies = df %>%
    dplyr::select(study_accession,
                  ab_p_0) %>%
    filter(!is.na(ab_p_0)) %>%
    pull(study_accession) %>%
    unique()
  
  # Filter only inactivated influenza vaccine participants
  df_filtered <- df  %>%
    filter(group_long == "Influenza (IN)")
  
  # Filter for only the timepoints we care about
  df_filtered <- df_filtered %>%
    filter(time %in% timepoints_to_keep)
  
  # Define pre and post vaccination responses (nAb or HAI)
  df_filtered <- df_filtered %>%
    mutate(
      response_pre = ab_p_0,
      response_post = ab_p_28
    )
  
  # Filter any participants lacking immune responses
  df_filtered = df_filtered %>%
    filter(!is.na(response_pre), !is.na(response_post))
  
  
  # Remove participants lacking expression data at at least one timepoint
  df_filtered = df_filtered %>%
    group_by(participant_id) %>%
    filter(sum(time == "P+0D") == 1,
           sum(time == tp) == 1)  %>%
    ungroup()
  
  # Remove studies without at least 5 participants
  df_filtered = df_filtered %>%
    group_by(study_accession) %>%
    filter(length(unique(participant_id)) > 5) %>%
    ungroup()
  
  # Select relevant columns to keep
  df_filtered = df_filtered %>%
    dplyr::select(
      participant_id,
      study_accession,
      time,
      response_pre,
      response_post,
      all_of(gene_names)
    ) %>%
    arrange(participant_id)
  
  # Subsampling if desired
  if (!(screen.fraction %in% c(0, 1))) {
    # Sample 66% of participants per study for training; remainder becomes test set
    train_indices <- df_filtered %>%
      distinct(study_accession, participant_id) %>%
      group_by(study_accession) %>%
      slice_sample(prop = screen.fraction) %>%
      ungroup()
    
    df_train <- df_filtered %>%
      semi_join(train_indices, by = c("study_accession", "participant_id"))
    
    df_test <- df_filtered %>%
      anti_join(train_indices, by = c("study_accession", "participant_id"))
  } else if (screen.fraction == 1) {
    df_train <- df_filtered
    
    df_test = NULL
  } else if (screen.fraction == 0) {
    df_train <- NULL
    
    df_test = df_filtered
  }
  
  res = list(
    "df.full" = df_filtered,
    "df.screen" = df_train,
    "df.evaluate" = df_test
  )
  
  return(res)
}
