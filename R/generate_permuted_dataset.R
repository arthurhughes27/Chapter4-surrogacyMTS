#' Generate a permuted dataset by randomly reassigning surrogate subject rows within each study.
#'
#' For each study, subject rows of the surrogate matrices (pre and post together)
#' are permuted (same permutation applied to pre and post), thereby reassigning
#' surrogate trajectories to different response trajectories while preserving
#' between-marker correlation and the pre/post pairing for each surrogate.
#' The ordering of individuals (rows) is retained.
#'
#' @param df A preprocessed data frame with exactly one pre-timepoint
#'   (time == "P+0D") and one post-timepoint (time != "P+0D")
#'   row per participant. Must contain columns: participant_id, time,
#'   study_accession, response_pre, response_post, and all columns in gene_names.
#' @param gene_names Character vector of gene/biomarker column names.
#' @param seed numeric seed for reproducibility.
#' @param n_within_study numeric value giving a fixed sample size to subsample from each study.
#'   Defaults to \code{NULL}, in which case the original study sizes are used.
#' @param n_studies integer giving the number of studies to randomly sample before permutation.
#'   Defaults to \code{NULL}, in which case all studies are retained.
#'
#' @return A named list with elements yone, yzero, sone, szero, studyone, studyzero.
generate_permuted_dataset <- function(df,
                                      gene_names,
                                      seed,
                                      n_within_study = NULL,
                                      n_studies = NULL) {
  library(tidyverse)
  
  set.seed(seed)
  
  # optional subsampling of studies
  if (!is.null(n_studies)) {
    all_studies    <- unique(df$study_accession)
    sampled_studies <- sample(all_studies, size = min(n_studies, length(all_studies)))
    df <- df %>% filter(study_accession %in% sampled_studies)
  }
  
  # optional subsampling of participants within each study
  if (!is.null(n_within_study)) {
    # subsample within each study
    keep_ids <- df %>%
      distinct(study_accession, participant_id) %>%
      group_by(study_accession) %>%
      group_modify( ~ .x %>% slice_sample(n = min(n_within_study, nrow(.x)))) %>%
      ungroup() %>%
      pull(participant_id)
    
    df <- df %>% filter(participant_id %in% keep_ids)
  }
  
  # arrange so pre/post rows align by study and participant
  pre_rows  <- df %>% filter(time == "P+0D")  %>% arrange(study_accession, participant_id)
  post_rows <- df %>% filter(time != "P+0D") %>% arrange(study_accession, participant_id)
  
  pre_genes  <- as.matrix(pre_rows[, gene_names])
  post_genes <- as.matrix(post_rows[, gene_names])
  
  # permute rows within each study
  studies <- unique(pre_rows$study_accession)
  for (st in studies) {
    idx <- which(pre_rows$study_accession == st)
    perm <- sample(seq_along(idx))
    pre_genes[idx, ]  <- pre_genes[idx[perm], , drop = FALSE]
    post_genes[idx, ] <- post_genes[idx[perm], , drop = FALSE]
  }
  
  # responses unchanged
  yzero <- pre_rows$response_pre
  yone  <- post_rows$response_post
  
  res = list(
    yone      = yone,
    yzero     = yzero,
    sone      = as.data.frame(post_genes),
    szero     = as.data.frame(pre_genes),
    studyone  = post_rows$study_accession,
    studyzero = pre_rows$study_accession
  )
  
  return(res)
}
