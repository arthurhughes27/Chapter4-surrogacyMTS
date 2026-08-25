# File to pre-process the Blood Gen 3 Modules for downstream use
## This includes removing empty genesets and genesets with no genes in our data

# Load necessary packages
library(fs)
library(dplyr)
library(stringr)
library(GSA)
library(readr)
library(readxl)

# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_BG3M_raw <- fs::path(raw_data_folder, "Suppl_File_1_BIOINF.xls")
p_load_expression <- fs::path(processed_data_folder, "df_merged_all.rds")

# Load data
BG3M_raw = read_excel(p_load_BG3M_raw)
all_noNorm_eset <- readRDS(p_load_expression)

# First we must preprocess the BG3M_raw raw data into a .gmt format

# Member genes 
genesets <- strsplit(BG3M_raw$`Member genes`, ",\\s*")

# Trim spaces and remove duplicates/empties
genesets <- lapply(genesets, function(x) unique(trimws(x[x != ""])))

BG3M <- list(
  genesets = genesets,
  geneset.names = BG3M_raw$ID,
  geneset.descriptions = BG3M_raw$`Module title`
)

# First take the expression data
expr_data = all_noNorm_eset %>% 
  dplyr::select(a1cf:zzz3)

# Now extract the gene names from the expression matrix - we are only interested in keeping these
gene_names = colnames(expr_data) %>%
  tolower()

# Now in the geneset data - make the geneset names lowercase
BG3M[["genesets"]] = lapply(
  BG3M[["genesets"]],
  FUN = function(x) {
    tolower(x)
  }
)

# Remove any empty elements
BG3M[["genesets"]]  = lapply(
  BG3M[["genesets"]],
  FUN = function(x) {
    x[x != ""]
  }
)

# Use only genesets which contain at least one of the commonly aligned genes in IS2
## Find the indices of the genesets containing at least one relevant gene
gs_membership <-
  sapply(BG3M[["genesets"]], function(gs)
    any(gs %in% gene_names))

# Filter the BG3M object to remove genesets without IS2 member
BG3M[["genesets"]] <- BG3M[["genesets"]][gs_membership]

BG3M[["geneset.names"]] <- BG3M[["geneset.names"]][gs_membership]

BG3M[["geneset.descriptions"]] <-
  BG3M[["geneset.descriptions"]][gs_membership]

# Rename the column names.descriptions to names
BG3M[["geneset.names.descriptions"]] = paste0(BG3M[["geneset.descriptions"]], " (", BG3M[["geneset.names"]], ")")

names(BG3M)[names(BG3M) == "Cluster"] = "geneset.aggregates"

names(BG3M[["genesets"]]) = NULL

# Remove TBA genesets
BG3M_filter <- !grepl("TBD", BG3M[["geneset.names.descriptions"]])
BG3M[["geneset.names.descriptions"]] <- BG3M[["geneset.names.descriptions"]][BG3M_filter]
BG3M[["geneset.names"]] <- BG3M[["geneset.names"]][BG3M_filter]
BG3M[["geneset.descriptions"]] <- BG3M[["geneset.descriptions"]][BG3M_filter]
BG3M[["geneset.aggregates"]] <- BG3M[["geneset.aggregates"]][BG3M_filter]
BG3M[["genesets"]] <- BG3M[["genesets"]][BG3M_filter]

# Save the processed geneset object

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "BG3M_processed.rds")

# Save dataframe
saveRDS(BG3M, file = p_save)

rm(list = ls())

