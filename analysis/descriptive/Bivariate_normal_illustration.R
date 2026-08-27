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

# Function to run a battery of MVN bivariate normality tests for one
# (response, gene) pair. All tests are run for a given gene within the same
# worker call, reusing the same parallel (one-task-per-gene) structure instead
# of parallelizing separately per test.
# The energy (E-statistic) test is deliberately excluded here: unlike the
# others it has no closed-form/asymptotic p-value and always requires
# bootstrap resampling, which would be far more expensive to run across
# thousands of genes.
run_bvn_tests <- function(gene_name, data, response = "ab_p_0") {
  pair_data <- data %>% dplyr::select(all_of(c(response, gene_name)))

  safe_test <- function(fn, ...) {
    tryCatch(fn(data = pair_data, ...), error = function(e) NULL)
  }

  single_row <- function(res, test_label) {
    tibble(
      gene      = gene_name,
      test      = test_label,
      statistic = if (is.null(res)) NA_real_ else res$Statistic,
      p_value   = if (is.null(res)) NA_real_ else res$p.value,
      error     = is.null(res)
    )
  }

  res_hz  <- safe_test(hz,             bootstrap = FALSE)
  res_hw  <- safe_test(hw,             bootstrap = FALSE)
  res_roy <- safe_test(royston)
  res_dh  <- safe_test(doornik_hansen)
  res_mar <- safe_test(mardia,         bootstrap = FALSE)

  bind_rows(
    single_row(res_hz,  "Henze-Zirkler"),
    single_row(res_hw,  "Henze-Wagner"),
    single_row(res_roy, "Royston"),
    single_row(res_dh,  "Doornik-Hansen"),
    if (is.null(res_mar)) {
      bind_rows(
        single_row(NULL, "Mardia Skewness"),
        single_row(NULL, "Mardia Kurtosis")
      )
    } else {
      tibble(
        gene      = gene_name,
        test      = res_mar$Test,
        statistic = res_mar$Statistic,
        p_value   = res_mar$p.value,
        error     = FALSE
      )
    }
  )
}

# Parallelize across genes (adjust workers to your machine)
plan(multisession, workers = 8)

set.seed(18042025)
bvn_results <- future_map_dfr(
  gene_cols,
  ~ run_bvn_tests(.x, data = data),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)

# Flag which pairs reject bivariate normality at alpha = 0.05, with BH
# adjustment applied within each test (different tests are not pooled
# together for multiplicity correction)
bvn_results <- bvn_results %>%
  group_by(test) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(bvn_rejected = p_value < 0.05,
         bvn_adjusted_rejected = p_adjusted < 0.05)

# Proportion of gene pairs rejecting bivariate normality, per test
bvn_results %>%
  group_by(test) %>%
  summarise(
    prop_rejected_raw      = mean(bvn_rejected, na.rm = TRUE),
    prop_rejected_adjusted = mean(bvn_adjusted_rejected, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(bvn_results, fs::path(results_path, "bvn_results.rds"))
