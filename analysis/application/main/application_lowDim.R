library(Surrogate)
library(SurrogateRank)
library(ggplot2)
library(dplyr)
library(patchwork)
library(fs)
library(boot)

application_figures_folder <- fs::path("output", "figures", "application", "main")

ccc_fun <- function(x, y) {
  if (length(x) < 2 || length(y) < 2)
    return(NA_real_)
  if (var(x, na.rm = TRUE) == 0 ||
      var(y, na.rm = TRUE) == 0)
    return(NA_real_)
  2 * cov(x, y, use = "complete.obs") /
    (var(x, na.rm = TRUE) +
       var(y, na.rm = TRUE) +
       (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))^2)
}

make_legend_breaks <- function(n_vals) {
  n_vals <- n_vals[!is.na(n_vals)]   # <-- remove NA values
  
  if (length(n_vals) == 0) {
    return(list(breaks = NULL, labels = NULL))
  }
  
  if (length(unique(n_vals)) == 1) {
    legend_breaks <- unique(n_vals)
    legend_labels <- as.character(unique(n_vals))
  } else {
    min_val <- min(n_vals)
    max_val <- max(n_vals)
    median_val <- median(n_vals)
    mid_val <- round(median_val)
    
    if (mid_val <= min_val)
      mid_val <- min_val + 1
    if (mid_val >= max_val)
      mid_val <- max_val - 1
    
    legend_breaks <- c(min_val, mid_val, max_val)
    legend_labels <- as.character(legend_breaks)
  }
  
  list(breaks = legend_breaks, labels = legend_labels)
}

make_size_scale <- function(n_vals) {
  leg <- make_legend_breaks(n_vals)
  scale_size_continuous(
    range = c(5, 20),
    breaks = leg$breaks,
    labels = leg$labels
  )
}

# ============================================================
# ARMD
# ============================================================

data(ARMD)

ok_centers_ARMD <- ARMD %>%
  group_by(Center, Treat) %>%
  summarise(n_patients = n(), .groups = "drop") %>%
  group_by(Center) %>%
  filter(all(n_patients >= 2), sum(n_patients) >= 5) %>%
  pull(Center)

dat_ARMD <- ARMD %>%
  filter(Center %in% ok_centers_ARMD)

dat_ARMD$Id %>% unique() %>% length()
dat_ARMD$Center %>% unique() %>% length()

fit_full_ARMD <- BifixedContCont(
  Dataset  = dat_ARMD,
  Surr     = Diff24,
  True     = Diff52,
  Treat    = Treat,
  Trial.ID = Center,
  Pat.ID   = Id,
  Model    = "Full",
  Weighted = TRUE
)

trial_R2_ARMD <- fit_full_ARMD$Trial.R2$`R2 Trial`

cat("ARMD full-model trial-level R^2:", trial_R2_ARMD, "\n")

trial_effects_ARMD <- fit_full_ARMD$Results.Stage.1 %>%
  transmute(
    Center = as.factor(Trial),
    n      = Obs.per.trial,
    u.s    = Treatment.S,
    u.y    = Treatment.T
  )

trial_ccc_ARMD <- ccc_fun(trial_effects_ARMD$u.s, trial_effects_ARMD$u.y)

# Bootstrap CI for CCC from the bivariate joint modelling trial-level effects
ccc_stat_ARMD <- function(data, indices) {
  d <- data[indices, ]
  ccc_fun(d$u.s, d$u.y)
}

set.seed(123)
trial_ccc_ARMD_boot <- boot(data = trial_effects_ARMD,
                            statistic = ccc_stat_ARMD,
                            R = 2000)

trial_ccc_ARMD_CI <- boot.ci(trial_ccc_ARMD_boot, type = c("bca"))

trial_ccc_ARMD_CI_bca <- c(lower = trial_ccc_ARMD_CI$bca[1, 4], upper = trial_ccc_ARMD_CI$bca[1, 5])


plot.min.global_ARMD <- min(trial_effects_ARMD$u.s, trial_effects_ARMD$u.y) - 2

jointModel_plot_ARMD <- trial_effects_ARMD %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(
    aes(size = n),
    shape = 21,
    alpha = 0.5,
    stroke = 1,
    fill = "#BCBCBC"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "#FF2128",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5
  ) +
  annotate(
    "text",
    x = plot.min.global_ARMD,
    y = max(trial_effects_ARMD$u.y, na.rm = TRUE),
    label = paste0(
      "Trial R2 = ",
      round(trial_R2_ARMD, 2),
      " [",
      round(fit_full_ARMD$Trial.R2$`CI lower limit`, 2),
      ", ",
      round(fit_full_ARMD$Trial.R2$`CI upper limit`, 2),
      "]",
      "\nCCC = ",
      round(trial_ccc_ARMD, 2),
      " [",
      round(trial_ccc_ARMD_CI_bca[1], 2),
      ", ",
      round(trial_ccc_ARMD_CI_bca[2], 2),
      "]"
    ),
    hjust = 0,
    vjust = 1,
    color = "red",
    size = 8
  ) +
  make_size_scale(c(trial_effects_ARMD$n, NA)) +
  scale_x_continuous(limits = c(
    plot.min.global_ARMD,
    max(trial_effects_ARMD$u.s, na.rm = TRUE) + 2
  ),
  expand = c(0, 0)) +
  scale_y_continuous(limits = c(
    plot.min.global_ARMD,
    max(trial_effects_ARMD$u.y, na.rm = TRUE) + 2
  ),
  expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Treatment effect on CVA at 6 months", y = "Treatment effect on CVA at 12 months", size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "right"
  )

yone_ARMD       <- dat_ARMD$Diff52[dat_ARMD$Treat ==  1]
yzero_ARMD      <- dat_ARMD$Diff52[dat_ARMD$Treat == -1]
sone_ARMD       <- dat_ARMD$Diff24[dat_ARMD$Treat ==  1]
szero_ARMD      <- dat_ARMD$Diff24[dat_ARMD$Treat == -1]
Centerone_ARMD  <- as.character(dat_ARMD$Center[dat_ARMD$Treat ==  1])
Centerzero_ARMD <- as.character(dat_ARMD$Center[dat_ARMD$Treat == -1])

rise_fit_ARMD <- rise.screen.meta(
  yone = yone_ARMD,
  yzero = yzero_ARMD,
  sone = sone_ARMD,
  szero = szero_ARMD,
  studyone = Centerone_ARMD,
  studyzero = Centerzero_ARMD,
  alpha = 0.05,
  epsilon.study = 0.1,
  epsilon.meta = 0.1,
  p.correction = "none",
  return.study.similarity.plot = FALSE,
  paired.all = FALSE,
  test = "knha",
  meta.analysis.method = "RE"
)

gamma_df_ARMD <- rise_fit_ARMD[["screening.metrics.study"]]

R2_rise_w_ARMD <- summary(lm(u.y ~ u.s, data = gamma_df_ARMD, weights = n))$r.squared
rise_ccc_ARMD  <- ccc_fun(gamma_df_ARMD$u.s, gamma_df_ARMD$u.y)

# Bootstrap CI for weighted R^2 from RISE-Meta
r2_rise_w_stat_ARMD <- function(data, indices) {
  d <- data[indices, ]
  summary(lm(u.y ~ u.s, data = d, weights = n))$r.squared
}

set.seed(123)
R2_rise_w_ARMD_boot <- boot(data = gamma_df_ARMD,
                            statistic = r2_rise_w_stat_ARMD,
                            R = 2000)

R2_rise_w_ARMD_CI <- boot.ci(R2_rise_w_ARMD_boot, type = c("bca"))

R2_rise_w_ARMD_CI_bca <- c(lower = R2_rise_w_ARMD_CI$bca[1, 4], upper = R2_rise_w_ARMD_CI$bca[1, 5])

# Bootstrap CI for CCC from RISE-Meta trial-level effects
rise_ccc_stat_ARMD <- function(data, indices) {
  d <- data[indices, ]
  ccc_fun(d$u.s, d$u.y)
}

set.seed(123)
rise_ccc_ARMD_boot <- boot(data = gamma_df_ARMD,
                           statistic = rise_ccc_stat_ARMD,
                           R = 2000)

rise_ccc_ARMD_CI <- boot.ci(rise_ccc_ARMD_boot, type = c("bca"))

rise_ccc_ARMD_CI_bca <- c(lower = rise_ccc_ARMD_CI$bca[1, 4], upper = rise_ccc_ARMD_CI$bca[1, 5])

riseMeta_plot_ARMD <- gamma_df_ARMD %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(
    aes(size = n),
    shape = 21,
    alpha = 0.5,
    stroke = 1,
    fill = "#BCBCBC"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "#FF2128",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5
  ) +
  annotate(
    "text",
    x = -0.1,
    y = 1.05,
    label = paste0(
      "Trial R2 = ",
      round(R2_rise_w_ARMD, 2),
      " [",
      round(R2_rise_w_ARMD_CI_bca["lower"], 2),
      ", ",
      round(R2_rise_w_ARMD_CI_bca["upper"], 2),
      "]",
      "\nCCC = ",
      round(rise_ccc_ARMD, 2),
      " [",
      round(rise_ccc_ARMD_CI_bca["lower"], 2),
      ", ",
      round(rise_ccc_ARMD_CI_bca["upper"], 2),
      "]"
    ),
    hjust = 0,
    vjust = 1,
    color = "red",
    size = 8
  ) +
  make_size_scale(c(gamma_df_ARMD$n, NA)) +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Treatment effect on CVA at 6 months", y = "Treatment effect on CVA at 12 months", size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "right"
  )

# Shared legend breaks within ARMD only
armd_n_all <- c(trial_effects_ARMD$n, gamma_df_ARMD$n)
armd_size_scale <- make_size_scale(armd_n_all)

jointModel_plot_mod_ARMD <- jointModel_plot_ARMD +
  scale_size_continuous(
    range = c(5, 20),
    breaks = make_legend_breaks(armd_n_all)$breaks,
    labels = make_legend_breaks(armd_n_all)$labels
  ) +
  labs(title = "Bivariate Joint Modelling")

riseMeta_plot_mod_ARMD <- riseMeta_plot_ARMD +
  scale_size_continuous(
    range = c(5, 20),
    breaks = make_legend_breaks(armd_n_all)$breaks,
    labels = make_legend_breaks(armd_n_all)$labels
  ) +
  labs(title = "RISE-Meta") +
  theme(axis.title.y = element_blank())

combined_plot_ARMD <-
  (
    jointModel_plot_mod_ARMD + plot_spacer() + riseMeta_plot_mod_ARMD +
      plot_layout(guides = "collect", widths = c(4, 0.5, 4))
  ) +
  plot_annotation(title = "A) Age-related macular degeneration dataset",
                  theme = theme(plot.title = element_text(
                    size = 40, face = "bold", hjust = 0.5
                  )))

# ============================================================
# Ovarian cancer
# ============================================================

data(Ovarian)

ok_centers_Ovarian <- Ovarian %>%
  group_by(Center, Treat) %>%
  summarise(n_patients = n(), .groups = "drop") %>%
  group_by(Center) %>%
  filter(all(n_patients >= 2), sum(n_patients) >= 5) %>%
  pull(Center)

dat_Ovarian <- Ovarian %>%
  filter(Center %in% ok_centers_Ovarian) %>%
  mutate(Surv = log(Surv), Pfs = log(Pfs))

dat_Ovarian$Patient %>% unique() %>% length()
dat_Ovarian$Center %>% unique() %>% length()

fit_full_Ovarian <- BifixedContCont(
  Dataset  = dat_Ovarian,
  Surr     = Pfs,
  True     = Surv,
  Treat    = Treat,
  Trial.ID = Center,
  Pat.ID   = Patient,
  Model    = "Full",
  Weighted = TRUE
)

trial_R2_Ovarian <- fit_full_Ovarian$Trial.R2$`R2 Trial`

cat("Ovarian full-model trial-level R^2:", trial_R2_Ovarian, "\n")

trial_effects_Ovarian <- fit_full_Ovarian$Results.Stage.1 %>%
  transmute(
    Center = as.factor(Trial),
    n      = Obs.per.trial,
    u.s    = Treatment.S,
    u.y    = Treatment.T
  )

trial_ccc_Ovarian <- ccc_fun(trial_effects_Ovarian$u.s, trial_effects_Ovarian$u.y)

# Bootstrap CI for CCC from the bivariate joint modelling trial-level effects
ccc_stat_Ovarian <- function(data, indices) {
  d <- data[indices, ]
  ccc_fun(d$u.s, d$u.y)
}

set.seed(123)
trial_ccc_Ovarian_boot <- boot(data = trial_effects_Ovarian,
                               statistic = ccc_stat_Ovarian,
                               R = 2000)

trial_ccc_Ovarian_CI <- boot.ci(trial_ccc_Ovarian_boot, type = c("bca"))

trial_ccc_Ovarian_CI_bca <- c(lower = trial_ccc_Ovarian_CI$bca[1, 4], upper = trial_ccc_Ovarian_CI$bca[1, 5])

plot.min.global_Ovarian <- min(trial_effects_Ovarian$u.s, trial_effects_Ovarian$u.y) - 0.2

jointModel_plot_Ovarian <- trial_effects_Ovarian %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(
    aes(size = n),
    shape = 21,
    alpha = 0.5,
    stroke = 1,
    fill = "#BCBCBC"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "#FF2128",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5
  ) +
  annotate(
    "text",
    x = plot.min.global_Ovarian,
    y = max(trial_effects_Ovarian$u.y, na.rm = TRUE),
    label = paste0(
      "Trial R2 = ",
      round(trial_R2_Ovarian, 2),
      " [",
      round(fit_full_Ovarian$Trial.R2$`CI lower limit`, 2),
      ", ",
      round(fit_full_Ovarian$Trial.R2$`CI upper limit`, 2),
      "]",
      "\nCCC = ",
      round(trial_ccc_Ovarian, 2),
      " [",
      round(trial_ccc_Ovarian_CI_bca["lower"], 2),
      ", ",
      round(trial_ccc_Ovarian_CI_bca["upper"], 2),
      "]"
    ),
    hjust = 0,
    vjust = 1,
    color = "red",
    size = 8
  ) +
  make_size_scale(c(trial_effects_Ovarian$n, NA)) +
  scale_x_continuous(limits = c(
    plot.min.global_Ovarian,
    max(trial_effects_Ovarian$u.s, na.rm = TRUE) + 0.2
  ),
  expand = c(0, 0)) +
  scale_y_continuous(limits = c(
    plot.min.global_Ovarian,
    max(trial_effects_Ovarian$u.y, na.rm = TRUE) + 0.2
  ),
  expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Treatment effect on log(progression-free survival)", y = "Treatment effect on log(overall survival)", size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "right"
  )

yone_Ovarian       <- dat_Ovarian$Surv[dat_Ovarian$Treat == 1]
yzero_Ovarian      <- dat_Ovarian$Surv[dat_Ovarian$Treat == 0]
sone_Ovarian       <- dat_Ovarian$Pfs[dat_Ovarian$Treat == 1]
szero_Ovarian      <- dat_Ovarian$Pfs[dat_Ovarian$Treat == 0]
Centerone_Ovarian  <- as.character(dat_Ovarian$Center[dat_Ovarian$Treat == 1])
Centerzero_Ovarian <- as.character(dat_Ovarian$Center[dat_Ovarian$Treat == 0])

rise_fit_Ovarian <- rise.screen.meta(
  yone = yone_Ovarian,
  yzero = yzero_Ovarian,
  sone = sone_Ovarian,
  szero = szero_Ovarian,
  studyone = Centerone_Ovarian,
  studyzero = Centerzero_Ovarian,
  alpha = 0.05,
  epsilon.study = 0.1,
  epsilon.meta = 0.1,
  p.correction = "none",
  return.study.similarity.plot = FALSE,
  paired.all = FALSE,
  test = "knha",
  meta.analysis.method = "RE"
)

gamma_df_Ovarian <- rise_fit_Ovarian[["screening.metrics.study"]]

R2_rise_w_Ovarian <- summary(lm(u.y ~ u.s, data = gamma_df_Ovarian, weights = n))$r.squared
rise_ccc_Ovarian  <- ccc_fun(gamma_df_Ovarian$u.s, gamma_df_Ovarian$u.y)

# Bootstrap CI for weighted R^2 from RISE-Meta
r2_rise_w_stat_Ovarian <- function(data, indices) {
  d <- data[indices, ]
  summary(lm(u.y ~ u.s, data = d, weights = n))$r.squared
}

set.seed(123)
R2_rise_w_Ovarian_boot <- boot(data = gamma_df_Ovarian,
                               statistic = r2_rise_w_stat_Ovarian,
                               R = 2000)

R2_rise_w_Ovarian_CI <- boot.ci(R2_rise_w_Ovarian_boot, type = c("bca"))

R2_rise_w_Ovarian_CI_bca <- c(lower = R2_rise_w_Ovarian_CI$bca[1, 4], upper = R2_rise_w_Ovarian_CI$bca[1, 5])

# Bootstrap CI for CCC from RISE-Meta trial-level effects
rise_ccc_stat_Ovarian <- function(data, indices) {
  d <- data[indices, ]
  ccc_fun(d$u.s, d$u.y)
}

set.seed(123)
rise_ccc_Ovarian_boot <- boot(data = gamma_df_Ovarian,
                              statistic = rise_ccc_stat_Ovarian,
                              R = 2000)

rise_ccc_Ovarian_CI <- boot.ci(rise_ccc_Ovarian_boot, type = c("bca"))

rise_ccc_Ovarian_CI_bca <- c(lower = rise_ccc_Ovarian_CI$bca[1, 4], upper = rise_ccc_Ovarian_CI$bca[1, 5])

riseMeta_plot_Ovarian <- gamma_df_Ovarian %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(
    aes(size = n),
    shape = 21,
    alpha = 0.5,
    stroke = 1,
    fill = "#BCBCBC"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "#FF2128",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5
  ) +
  annotate(
    "text",
    x = -0.1,
    y = 1.05,
    label = paste0(
      "Trial R2 = ",
      round(R2_rise_w_Ovarian, 2),
      " [",
      round(R2_rise_w_Ovarian_CI_bca["lower"], 2),
      ", ",
      round(R2_rise_w_Ovarian_CI_bca["upper"], 2),
      "]",
      "\nCCC = ",
      round(rise_ccc_Ovarian, 2),
      " [",
      round(rise_ccc_Ovarian_CI_bca["lower"], 2),
      ", ",
      round(rise_ccc_Ovarian_CI_bca["upper"], 2),
      "]"
    ),
    hjust = 0,
    vjust = 1,
    color = "red",
    size = 8
  ) +
  make_size_scale(c(gamma_df_Ovarian$n, NA)) +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "Treatment effect on log(progression-free survival)", y = "Treatment effect on log(overall survival)", size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "right"
  )

# Shared legend breaks within Ovarian only
ov_n_all <- c(trial_effects_Ovarian$n, gamma_df_Ovarian$n)
ov_size_breaks <- make_legend_breaks(ov_n_all)

jointModel_plot_mod_Ovarian <- jointModel_plot_Ovarian +
  scale_size_continuous(
    range = c(5, 20),
    breaks = ov_size_breaks$breaks,
    labels = ov_size_breaks$labels
  ) +
  labs(title = "Bivariate Joint Modelling")

riseMeta_plot_mod_Ovarian <- riseMeta_plot_Ovarian +
  scale_size_continuous(
    range = c(5, 20),
    breaks = ov_size_breaks$breaks,
    labels = ov_size_breaks$labels
  ) +
  labs(title = "RISE-Meta") +
  theme(axis.title.y = element_blank())

combined_plot_Ovarian <-
  (
    jointModel_plot_mod_Ovarian + plot_spacer() + riseMeta_plot_mod_Ovarian +
      plot_layout(guides = "collect", widths = c(4, 0.5, 4))
  ) +
  plot_annotation(title = "B) Ovarian cancer dataset",
                  theme = theme(plot.title = element_text(
                    size = 40, face = "bold", hjust = 0.5
                  )))

# ============================================================
# Overall combined plot
# ============================================================

overall_combined_plot <-
  wrap_elements(full = combined_plot_ARMD) /
  plot_spacer() /
  wrap_elements(full = combined_plot_Ovarian) +
  plot_layout(heights = c(1, 0.06, 1))

overall_combined_plot

# ============================================================
# Save plots
# ============================================================

ggsave(
  filename = "risemeta-lowdim.pdf",
  plot = overall_combined_plot,
  path  = application_figures_folder,
  width = 40,
  height = 42,
  units = "cm"
)

rm(list = ls())
