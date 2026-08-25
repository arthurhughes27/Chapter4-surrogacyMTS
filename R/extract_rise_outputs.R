# Function to extract outputs from RISE-Meta screen and evaluation functions
extract_rise_outputs = function(screen_result = NULL,
                                evaluation_result = NULL) {
  library(knitr)
  library(kableExtra)
  
  if (!is.null(screen_result)) {
    # Extract per-marker screening metrics from result
    screening_metrics <- screen_result[["screening.metrics.meta"]]
    
    n_significant = length(screen_result[["significant.markers"]])
    
    epsilon = mean(screen_result[["screening.metrics.study"]][["epsilon"]])
    
    # Format significant markers into a publication-ready table
    screening_metrics_select <- screening_metrics %>%
      arrange(p.unadjusted) %>%
      mutate(mu.delta.ci = paste0(
        round(mu.delta, 3),
        " (",
        round(ci.delta.lower, 3),
        ", ",
        round(ci.delta.upper, 3),
        ")"
      )) %>%
      dplyr::select(marker, mu.delta.ci, p.unadjusted, p.adjusted) %>%
      filter(p.adjusted < 0.05)
    
    # Render LaTeX table of significant screening markers
    screen_table =
      kable(
        screening_metrics_select,
        format   = "latex",
        booktabs = TRUE,
        caption  = "Screening Metrics Table"
      ) %>%
      kable_styling(latex_options = "hold_position") %>%
      row_spec(0, bold = TRUE) %>%
      column_spec(1:ncol(screening_metrics_select))
    
    forest_plot = screen_result[["gamma.s.plot"]]$forest.plot
    
    fit_plot = screen_result[["gamma.s.plot"]]$fit.plot
    
    screen_plot = screen_result[["gamma.s.plot"]]$screen.plot
    
    res = list(
      "n_significant" = n_significant,
      "epsilon" = epsilon,
      "screen_table" = screen_table,
      "screen_plot" = screen_plot,
      "screen_forest" = forest_plot,
      "screen_fit" = fit_plot
    )
  } else if (!is.null(evaluation_result)) {
    # Extract per-marker screening metrics from result
    evaluation_metrics <- evaluation_result[["evaluation.metrics.meta"]] %>%
      distinct()
    
    # Format significant markers into a publication-ready table
    evaluation_metrics_select <- evaluation_metrics %>%
      arrange(p) %>%
      mutate(mu.delta.ci = paste0(
        round(mu.delta, 3),
        " (",
        round(ci.delta.lower, 3),
        ", ",
        round(ci.delta.upper, 3),
        ")"
      )) %>%
      dplyr::select(marker, mu.delta.ci, p)  %>%
      distinct()
    
    # Render LaTeX table of significant screening markers
    evaluation_table =
      kable(
        evaluation_metrics_select,
        format   = "latex",
        booktabs = TRUE,
        caption  = "Screening Metrics Table"
      ) %>%
      kable_styling(latex_options = "hold_position") %>%
      row_spec(0, bold = TRUE) %>%
      column_spec(1:ncol(evaluation_metrics_select))
    
    
    fit_plot = evaluation_result[["gamma.s.plot"]][["fit.plot"]]
    
    forest_plot = evaluation_result[["gamma.s.plot"]][["forest.plot"]]
    
    res = list(
      "evaluation_table" = evaluation_table,
      "evaluation_forest" = forest_plot,
      "evaluation_fit" = fit_plot
    )
    
    return(res)
  }
}
