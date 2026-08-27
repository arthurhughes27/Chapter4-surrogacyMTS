library(tidyverse)
library(MVN)
library(furrr)

set.seed(18042025)

processed_data_path <- fs::path("data")
figure_path <- fs::path("output", "figures", "application")
results_path = fs::path("output", "results", "descriptive")
p_load_merged_all = fs::path(file = fs::path(processed_data_path, "df_merged_all.rds"))

df_merged_all = readRDS(p_load_merged_all)

# Identify participants with entries at BOTH P+0D and P+1D
participants_both_times <- df_merged_all %>%
  filter(time %in% c("P+0D", "P+1D")) %>%
  distinct(participant_id, time) %>%
  count(participant_id) %>%
  filter(n == 2) %>%
  pull(participant_id)

df <- df_merged_all %>%
  filter(
    study_accession == "SDY1276",
    !is.na(ab_p_28),!is.na(ab_p_0),
    participant_id %in% participants_both_times
  ) %>%
  arrange(participant_id) %>%
  dplyr::select(where( ~ !any(is.na(.))))

# How many input genes
n_inputs = df %>%
  dplyr::select(a1cf:zzz3) %>%
  ncol()

data = df %>%
  filter(time == "P+0D") %>%
  dplyr::select(ab_p_0, a1cf:zzz3)

gene_cols <- df %>%
  filter(time == "P+0D") %>%
  dplyr::select(a1cf:zzz3) %>%
  names()

gene_cols = gene_cols[1:100]

data <- df %>%
  filter(time == "P+0D") %>%
  dplyr::select(ab_p_0, all_of(gene_cols))

# Function to run hz's test for one (response, gene) pair
run_hz_bvn <- function(gene_name, data, response = "ab_p_0",
                            B = 1000, bootstrap = F) {
  pair_data <- data %>% dplyr::select(all_of(c(response, gene_name)))

  res <- tryCatch(
    hz(data = pair_data, B = B, cores = 1, bootstrap = F),
    error = function(e) NULL
  )

  tibble(
    gene      = gene_name,
    statistic = if (is.null(res)) NA_real_ else res$Statistic,
    p_value   = if (is.null(res)) NA_real_ else res$p.value,
    error     = is.null(res)
  )
}

# Parallelize across genes (adjust workers to your machine)
plan(multisession, workers = 8)

set.seed(18042025)
hz_results <- future_map_dfr(
  gene_cols,
  ~ run_hz_bvn(.x, data = data, B = 5000, bootstrap = F),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)

# Flag which pairs reject bivariate normality at alpha = 0.05
hz_results <- hz_results %>%
  mutate(bvn_rejected = p_value < 0.05,
         p_adjusted = p.adjust(p_value, method = "BH")) %>%
  mutate(bvn_adjusted_rejected = p_adjusted < 0.05)

sum(hz_results$bvn_adjusted_rejected)/nrow(hz_results)

saveRDS(hz_results, fs::path(results_path, "hz_bvn_results.rds"))
