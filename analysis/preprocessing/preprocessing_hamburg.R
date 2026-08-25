# -------------------------------------------------------------------------
# R Script to preprocess files from the hamburg rVSV trial transcriptomic substudy
# Purpose: load and harmonize transcriptomic and clinical data for analysis
# -------------------------------------------------------------------------

# ---- Libraries ----
library(dplyr)
library(biomaRt)
library(edgeR)
library(tidyr)
library(ggplot2)
library(fs)

# ---- File paths ----
raw_data_path   <- fs::path("data-raw", "hamburg")
counts_path     <- fs::path(raw_data_path, "GSE97590_counts.tsv")
clinical_path   <- fs::path(raw_data_path, "rVSV_Hamburg_rVSV_clinical_anonym.txt")

# ---- Load raw data ----
# Transcriptomic counts
hamburg_expr <- read.table(counts_path, sep = "\t", header = TRUE)

# Clinical metadata
hamburg_clinical <- read.delim(clinical_path)

# ---- Engineer clinical identifiers ----
# Create a pid column by extracting the part before "_" in RNA column IDs
hamburg_clinical <- hamburg_clinical %>%
  mutate(pid = sub("_.*", "", Anonym_RNA_columnID))

# ---- Link clinical identifiers to transcriptomic identifiers ----
# Keep only unique, non-missing pid associations
link <- hamburg_clinical %>%
  dplyr::select(id_num, pid) %>%
  filter(!is.na(pid)) %>%
  distinct()

# ---- Merge clinical metadata with pid link ----
# Remove original pid.x and RNA column ID, keep consistent pid
# Drop patients without valid pid
hamburg_clinical <- hamburg_clinical %>%
  left_join(link, by = "id_num") %>%
  dplyr::select(-pid.x, -Anonym_RNA_columnID) %>%
  rename(pid = pid.y) %>%
  filter(!is.na(pid)) %>%
  relocate(pid, .after = id_num)

# ---- Pivot clinical data to wide format ----
# Goal: create one column per antibody type and day
hamburg_clinical <- hamburg_clinical %>%
  pivot_wider(
    names_from  = c(AB_type, day),
    values_from = GP_specific_AB,
    names_sep   = "_",
    names_prefix = "GP_specific_AB_"
  ) %>%
  # Remove columns with NA measurements
  dplyr::select(-c(GP_specific_AB_NA_1, GP_specific_AB_NA_3))

# ---- Prepare transcriptomic expression data ----
# Transpose so that samples are rows and genes are columns
hamburg_expr <- t(hamburg_expr) %>%
  as.data.frame()

# Add full identifier column to link with clinical data
hamburg_expr <- hamburg_expr %>%
  mutate(
    pid_full = rownames(hamburg_expr),
    # Extract pid (sample ID) and timepoint from the identifier
    pid = sub("_.*", "", pid_full),
    time = sub(".*_", "", pid_full)
  ) %>%
  # Move identifier columns to the front
  relocate(pid_full, pid, time)

# ---- Merge clinical and transcriptomic data ----
hamburg_merged <- merge(
  x = hamburg_clinical,
  y = hamburg_expr,
  by = "pid"
)

# ---- Separate counts and clinical metadata ----
# First 19 columns correspond to clinical data / antibody measurements
hamburg_clinical <- hamburg_merged[, 1:19]
hamburg_counts   <- hamburg_merged[, -c(1:19)]

# ---- Create edgeR DGEList and normalize counts ----
dge_list <- DGEList(t(hamburg_counts))

# Normalize library sizes using TMM
dge <- calcNormFactors(dge_list, method = "TMM")

# Compute log2 CPM with TMM normalization
counts_log2_cpm_tmm <- cpm(dge, log = TRUE) %>%
  t() %>%
  as.data.frame()

# ---- Convert Ensembl gene IDs to official HGNC symbols ----

# Get all gene column names from normalized counts
ensembl_gene_ID <- colnames(counts_log2_cpm_tmm)

# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve mapping from Ensembl IDs to HGNC symbols
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_gene_ID,
  mart       = ensembl
)

# Keep only genes with existing HGNC symbols
ensembl_existing <- gene_conversion %>%
  filter(hgnc_symbol != "") %>%
  pull(ensembl_gene_id)

# Subset counts to only genes with valid HGNC symbols
counts_log2_cpm_tmm_hgnc <- counts_log2_cpm_tmm %>%
  dplyr::select(any_of(ensembl_existing))

# Create named vector to map Ensembl IDs to HGNC symbols
id_to_symbol <- setNames(gene_conversion$hgnc_symbol, gene_conversion$ensembl_gene_id)

# Update column names to HGNC symbols
colnames(counts_log2_cpm_tmm_hgnc) <- id_to_symbol[colnames(counts_log2_cpm_tmm_hgnc)]

# Order columns alphabetically
counts_log2_cpm_tmm_hgnc <- counts_log2_cpm_tmm_hgnc[, order(colnames(counts_log2_cpm_tmm_hgnc))]

# ---- Combine clinical metadata with normalized expression data ----
hamburg_norm_final <- bind_cols(hamburg_clinical, counts_log2_cpm_tmm_hgnc) %>%
  janitor::clean_names()  # make column names consistent

hamburg_clinical = hamburg_norm_final %>% 
  dplyr::select(-(a1cf:zzz3))

# ---- Save processed dataset ----
processed_data_path <- fs::path("data", "hamburg")

# Save as RDS for single object storage
saveRDS(hamburg_norm_final, file = fs::path(processed_data_path, "hamburg_merged_norm.rds"))

saveRDS(hamburg_clinical, file = fs::path(processed_data_path, "hamburg_clinical.rds"))


# ---- Clean environment ----
rm(list = ls())
