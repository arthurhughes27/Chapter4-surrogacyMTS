# Script to produce a LaTeX table summarising the cell type (PBMC or whole
# blood) used for transcriptomic profiling in each IS2 influenza study.

library(tidyverse)
library(knitr)
library(kableExtra)

processed_data_path <- fs::path("data")
tables_path <- fs::path("output", "tables", "descriptive", "main")

is2_cell_type <- readRDS(fs::path(processed_data_path, "is2", "is2_cell_type.rds"))

cell_type_table <- is2_cell_type %>%
  rename(Study = study_accession, `Cell type` = cell_type)

cell_type_latex <- kable(
  cell_type_table,
  format = "latex",
  booktabs = TRUE,
  caption = "Cell type used for transcriptomic profiling (PBMC or whole blood) in each IS2 influenza study."
) %>%
  kable_styling(latex_options = "hold_position") %>%
  row_spec(0, bold = TRUE)

writeLines(cell_type_latex, fs::path(tables_path, "is2_cell_type_table.tex"))

rm(list = ls())
