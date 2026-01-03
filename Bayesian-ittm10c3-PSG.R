# ======== Bayesian BLMM for PSG data ========
# Author: Teng Gao
# Primary outcomes analysis
#Version:2025.12.13

# ========== Step 1. Load Libraries ==========
setwd("D:/R/Rdata/PHD-main work/Part2_FMT/PSG/PSGmice10C3/WASO")

library(mice)
library(brms)
library(bayesplot)
library(ggplot2)
library(posterior)
library(dplyr)
library(loo)
library(tidyr)
library(readr)
library(broom.mixed)
library(cowplot)
library(ggforce)
library(scales)


# Step 1. Main modle
#Read Data
    data <- read.csv("WASO.csv")
names(data) <- make.names(names(data))  # Clean names

# Convert variable types
data$Group    <- factor(data$Group, levels = c(0, 1))
data$subgroup <- factor(data$subgroup, levels = c(0, 1), labels = c("A", "B"))
data$center   <- factor(data$center)

#  Define covariates
covariates <- c("Sex", "BMI", "Age")


formula_covars <- paste(covariates, collapse = " + ")

# Multiple Imputation
mids_data <- mice(data[, c("post", "baseline", "Group", "subgroup", "center", covariates)], m = 10, seed = 2025)
completed_data_list <- complete(mids_data, action = "all")

#  Estimate Prior SD from baseline
baseline_sd <- sd(data$baseline[data$Group == 0], na.rm = TRUE)
adaptive_prior <- c(
  prior(normal(0, baseline_sd), class = "b"),
  prior(cauchy(0, 1), class = "sd")
)

#  Build formulas
formula_all <- as.formula(paste("post ~ Group + subgroup + baseline +(1|center)"))
formula_sub <- as.formula(paste("post ~ subgroup + baseline +(1|center)"))

# Step 8. Fit Main Model (fit_all)
fit_all <- brm_multiple(
  formula = formula_all,
  data = completed_data_list,
  family = gaussian(),
  prior = c(
    prior(normal(0, 57), class = "b"),
    prior(student_t(3, 0, 0.5), class = "sd")
  ),
  iter = 12000,
  warmup = 6000,
  chains = 4,
  seed = 2025,
  init = 0,
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

# Fit Stratified Subgroup Models
fit_group0 <- brm_multiple(
  formula = formula_sub,
  data = lapply(completed_data_list, function(df) df[df$Group == 0, ]),
  family = gaussian(),
  seed = 2025,
  init = 0,
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

fit_group1 <- brm_multiple(
  formula = formula_sub,
  data = lapply(completed_data_list, function(df) df[df$Group == 1, ]),
  family = gaussian(),
  seed = 2025,
  init = 0,
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

#  Posterior Extraction
p_all <- as_draws_df(fit_all)$b_Group1
p_subgroup <- as_draws_df(fit_all)$b_subgroupB
p_g0  <- as_draws_df(fit_group0)$b_subgroupB
p_g1  <- as_draws_df(fit_group1)$b_subgroupB

#  Posterior Summary Table
s_all    <- summary(fit_all)$fixed["Group1", ]
s_group0 <- summary(fit_group0)$fixed["subgroupB", ]
s_group1 <- summary(fit_group1)$fixed["subgroupB", ]

p_all_prob <- mean(p_all > 0)
p_g0_prob  <- mean(p_g0  > 0)
p_g1_prob  <- mean(p_g1  > 0)

model_summary <- data.frame(
  Model           = c("Overall sample (Group1 vs 0)", "Group 0 (B vs A)", "Group 1 (B vs A)"),
  Posterior_Mean = c(s_all["Estimate"], s_group0["Estimate"], s_group1["Estimate"]),
  CI_lower       = c(s_all["l-95% CI"], s_group0["l-95% CI"], s_group1["l-95% CI"]),
  CI_upper       = c(s_all["u-95% CI"], s_group0["u-95% CI"], s_group1["u-95% CI"]),
  P_gt_0         = c(p_all_prob, p_g0_prob, p_g1_prob),
  Rhat           = c(s_all["Rhat"], s_group0["Rhat"], s_group1["Rhat"]),
  ESS_bulk       = c(s_all["Bulk_ESS"], s_group0["Bulk_ESS"], s_group1["Bulk_ESS"]),
  ESS_tail       = c(s_all["Tail_ESS"], s_group0["Tail_ESS"], s_group1["Tail_ESS"])
)

write.csv(model_summary, "Posterior_Model_Summary_WASO.csv", row.names = FALSE)



# Step 2.VIF diagnosis----------------------------------------
# 1. Prepare data (exclude non-numeric and grouping variables)
vif_data <- data %>%
  select(all_of(covariates)) %>%   # Select only covariates
  select(where(is.numeric))       # Keep only numeric variables

# 2. Check if numeric variables exist
if (ncol(vif_data) == 0) {
  stop("No numeric covariates available for VIF analysis! Please check covariates list.")
}

# 3. Build temporary linear model (for VIF diagnostics only)
vif_model <- lm(as.formula(paste("post ~", paste(colnames(vif_data), collapse = " + "))), 
                data = data)

# 4. Calculate VIF and output results
vif_results <- car::vif(vif_model)
vif_summary <- data.frame(
  Variable = names(vif_results),
  VIF = round(vif_results, 2),
  Problem = ifelse(vif_results > 5, "High VIF (>5)", 
                   ifelse(vif_results > 3, "Moderate", "OK"))
)

# 5. Print high VIF warnings
high_vif_vars <- vif_summary %>% filter(VIF > 5)
if (nrow(high_vif_vars) > 0) {
  warning("High collinearity detected:\n", 
          paste(high_vif_vars$Variable, " (VIF=", high_vif_vars$VIF, ")", sep = "", collapse = "\n"))
} else {
  message("All covariates have VIF ≤5 - no serious collinearity issues")
}

# 6. Save results to CSV
write.csv(vif_summary, "VIF_diagnosis_WASO.csv", row.names = FALSE)
print(vif_summary)



# Step 3: Posterior Density Plot for Group/Subgroup Effects on SL
# Combine posterior samples into a long dataframe for plotting
df_post <- rbind(
  data.frame(value = p_all, group = "Group1 vs 0 (All)"),
  data.frame(value = p_g0,  group = "Subgroup B vs A (Group 0)"),
  data.frame(value = p_g1,  group = "Subgroup B vs A (Group 1)")
)

# Plot posterior densities
p <- ggplot(df_post, aes(x = value, fill = group, color = group)) +
  geom_density(alpha = 0.3, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Distributions of R Effects by Group/Subgroup",
    x = "Effect Estimate (ΔSL in minutes)", 
    y = "Posterior Density"
  ) +
  scale_fill_manual(values = c("#1f77b4", "#E69F00", "#7570b3")) +
  scale_color_manual(values = c("#1f77b4", "#E69F00", "#7570b3")) +
  theme(legend.title = element_blank())

# Save the plot
ggsave("Posterior_Comparison_Group_Subgroup.pdf", plot = p, width = 8, height = 5)

# (selected)Step . LOO-CV 模型比较与图示
formula_m1 <- as.formula(paste("post ~ Group + subgroup +", formula_covars, "+ (1|center)"))
formula_m2 <- as.formula(paste("post ~ Group + subgroup + baseline +", formula_covars, "+ (1|center)"))  
formula_m3 <- as.formula(paste("post ~ Group + subgroup + baseline +", formula_covars))
formula_m4 <- formula_all  # previously defined as main model

fit_m1 <- brm_multiple(formula_m1, data = completed_data_list, family = gaussian(),
                       prior = c(
                         prior(normal(0, 36), class = "b"),
                         prior(cauchy(0, 1), class = "sd")
                       ), iter = 8000,
                       warmup = 4000,
                       chains = 4, seed = 2025, control = list(adapt_delta = 0.95))
fit_m2 <- brm_multiple(formula_m2, data = completed_data_list, family = gaussian(),
                       prior = c(
                         prior(normal(0, 36), class = "b"),
                         prior(cauchy(0, 1), class = "sd")
                       ), iter = 8000,
                       warmup = 4000,
                       chains = 4, seed = 2025, control = list(adapt_delta = 0.95))
fit_m3 <- brm_multiple(formula = formula_m3, data = completed_data_list, family = gaussian(),
                       prior = prior(normal(0, 36), class = "b"),  # 只保留 b 的先验
                       iter = 8000,
                       warmup = 4000,
                       chains = 4, seed = 2025, control = list(adapt_delta = 0.95))
fit_m4 <- fit_all  # main model used



loo1 <- loo(fit_m1)
loo2 <- loo(fit_m2)
loo3 <- loo(fit_m3)
loo4 <- loo(fit_m4)

print(loo_compare(loo1, loo2, loo3, loo4))

# 保存 LOO 
loo_df <- tibble(
  model = c("fit_m1", "fit_m2 ", "fit_m3", "fit_m4(Main)"),
  elpd_loo = c(loo1$estimates["elpd_loo", "Estimate"],
               loo2$estimates["elpd_loo", "Estimate"],
               loo3$estimates["elpd_loo", "Estimate"],
               loo4$estimates["elpd_loo", "Estimate"]),
  se_elpd = c(loo1$estimates["elpd_loo", "R"],
              loo2$estimates["elpd_loo", "R"],
              loo3$estimates["elpd_loo", "R"],
              loo4$estimates["elpd_loo", "R"])
)

loo_df_sorted <- loo_df |> arrange(desc(elpd_loo))

ggplot(loo_df_sorted, aes(x = reorder(model, elpd_loo), y = elpd_loo)) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  geom_errorbar(aes(ymin = elpd_loo - se_elpd, ymax = elpd_loo + se_elpd), width = 0.05) +
  coord_flip() +
  labs(title = "Model Comparison via LOO-CV", x = NULL, y = "ELPD (higher = better)") +
  theme_minimal(base_size = 14)
ggsave("LOO_Comparison.pdf", width = 5, height = 3)

# Step 4: Export Fixed Effects Table for Main Model
fixed_effects <- summary(fit_all)$fixed
fixed_df <- as.data.frame(fixed_effects)
fixed_df$Predictor <- rownames(fixed_df)
fixed_df <- fixed_df[, c("Predictor", "Estimate", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS")]

param_names <- rownames(fixed_df)

draws <- as_draws_df(fit_all)
probs_all <- sapply(param_names, function(p) {
  colname <- paste0("b_", p)
  if (colname %in% colnames(draws)) {
    mean(draws[[colname]] > 0)
  } else {
    NA
  }
})
probs_all[which(param_names == "Intercept")] <- NA

fixed_df$P_gt_0 <- round(probs_all, 3)

write.csv(fixed_df, "Fixed_Effects_Estimates.csv", row.names = FALSE)

# Step 14 Posterior Predictive Checks and Trace Plots
# For each model
model_list <- list(fit_all = fit_all, fit_group0 = fit_group0, fit_group1 = fit_group1)

for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  
  # Posterior Predictive Check
  p_ppc <- pp_check(model, ndraws = 100) +
    ggtitle(paste("Posterior Predictive Check:", model_name)) +
    theme_minimal(base_size = 13)
  ggsave(paste0("PPC_SE", model_name, ".pdf"), p_ppc, width = 6, height = 4)
  
  # Trace Plot
  draw_df <- as_draws_df(model)
  pars_to_plot <- grep("b_", colnames(draw_df), value = TRUE)[1:min(4, length(grep("b_", colnames(draw_df), value = TRUE)))]
  
  p_trace <- mcmc_trace(as_draws_df(fit_all), pars = c("b_Group1", "b_subgroupB", "b_Intercept"))+
    ggtitle(paste("Trace Plot:", model_name))
  ggsave(paste0("Trace_SE", model_name, ".pdf"), p_trace, width = 8, height = 6)
}

# Step 5. Sensitivity Analysis with Stricter Priors
# 1. Fit stricter prior model
fit_strict <- brm_multiple(
  formula = formula_all,
  data = completed_data_list,
  family = gaussian(),
  prior = c(
    prior(normal(0,28), class = "b"),  # User-defined stricter prior
    prior(cauchy(0, 1), class = "sd")
  ),
  seed = 2025,
  control = list(adapt_delta = 0.95)
)

# 2. Extract posteriors
post_default <- as_draws_df(fit_all)$b_Group1
post_strict  <- as_draws_df(fit_strict)$b_Group1

# 3. Combine
df_post <- bind_rows(
  tibble(value = post_default, model = "Default Prior"),
  tibble(value = post_strict,  model = "Stricter Prior")
)

# 4. Compute posterior summary stats
summary_stats <- df_post %>%
  group_by(model) %>%
  summarise(
    Mean = mean(value),
    Lower = quantile(value, 0.025),
    Upper = quantile(value, 0.975),
    Prob = mean(value > 0)
  )

# 5. Construct label text 
summary_stats <- summary_stats %>%
  mutate(
    label = sprintf(
      "%s\nMean = %.2f\n95%% CrI = [%.2f, %.2f]\nP(Δ > 0) = %.2f",
      model, Mean, Lower, Upper, Prob
    ),
    x = ifelse(model == "Default Prior", Mean - 0.3, Mean +0.25),
    y = ifelse(model == "Default Prior",10, 10)
  )


# 6. Plot
summary_stats <- df_post %>%
  group_by(model) %>%
  summarise(
    mean_val = mean(value),
    median_val = median(value),
    ci_lower = quantile(value, 0.025),
    ci_upper = quantile(value, 0.975),
    p_gt_0 = mean(value > 0)
  ) %>%
  mutate(
    label = sprintf(
      "Mean: %.3f\nMedian: %.3f\n95%% CrI: [%.3f, %.3f]\nP(>0) = %.1f%%",
      mean_val, median_val, ci_lower, ci_upper, p_gt_0 * 100
    ),
    x_pos = ifelse(model == "Default Prior", 
                   min(df_post$value) + diff(range(df_post$value)) * 0.7,
                   min(df_post$value) + diff(range(df_post$value)) * 0.05),
    y_pos = max(density(df_post$value)$y) * 0.8
  )


library(ggplot2)
library(ggridges)
library(scales)

p <- ggplot(df_post, aes(x = value, fill = model, color = model)) +
 
  geom_density(alpha = 0.2, size = 1, adjust = 1.5) +
  stat_density(geom = "area", position = "identity", alpha = 0.3, adjust = 1.5) +
  geom_vline(data = summary_stats, 
             aes(xintercept = median_val, color = model), 
             linetype = "dashed", size = 0.8, alpha = 0.7) +
 
  geom_segment(data = summary_stats,
               aes(x = ci_lower, xend = ci_upper, y = 0, yend = 0, color = model),
               size = 3, alpha = 0.5, lineend = "round") +

  geom_vline(xintercept = 0, linetype = "solid", color = "grey30", 
             size = 0.8, alpha = 0.7) +

  geom_label(data = summary_stats,
             aes(x = x_pos, y = y_pos, label = label, fill = model),
             color = "white", size = 3.5, hjust = 0.5, vjust = 0.5,
             alpha = 0.9, label.size = 0, 
             label.padding = unit(0.4, "lines")) +

  theme_minimal(base_size = 14) +
  labs(
    title = "Prior Sensitivity Analysis: Group 1 Effect",
    subtitle = "Comparison of posterior distributions under different prior specifications",
    x = "Effect Estimate (Group 1 vs Group 0)",
    y = "Posterior Density",
    caption = "Vertical dashed lines: posterior medians\nThick horizontal lines: 95% credible intervals"
  ) +
  
  scale_fill_manual(
    values = c("Default Prior" = "#3B9AB2", "Stricter Prior" = "#E1AF00"),
    name = "Prior Specification"
  ) +
  scale_color_manual(
    values = c("Default Prior" = "#3B9AB2", "Stricter Prior" = "#E1AF00"),
    name = "Prior Specification"
  ) +

  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, 
                              margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12,
                                 margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0, color = "grey50", size = 10,
                                margin = margin(t = 10)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  
  annotate("text", x = 0, y = Inf, label = "No effect", 
           vjust = 2, hjust = 0.5, color = "grey40", size = 4, 
           fontface = "italic") +
  
  scale_x_continuous(
    expand = expansion(mult = 0.05)
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  )

print(p)

ggsave("Posterior_Group1_Prior_Comparison_Enhanced.pdf", 
       plot = p, width = 12, height = 8, dpi = 300)


p_simple <- ggplot(df_post, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 1) +
  scale_fill_manual(values = c("Default Prior" = "#1f77b4", "Stricter Prior" = "#E69F00")) +
  labs(
    title = "Prior Sensitivity Analysis",
    subtitle = "Group 1 Effect on R",
    x = "Effect Estimate",
    y = "Density",
    fill = "Prior"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

ggsave("Posterior_Group1_Prior_Comparison_Simple.pdf", 
       plot = p_simple, width = 10, height = 6)

stat_table <- summary_stats %>%
  select(model, mean_val, median_val, ci_lower, ci_upper, p_gt_0) %>%
  rename(
    Model = model,
    Mean = mean_val,
    Median = median_val,
    `CI Lower` = ci_lower,
    `CI Upper` = ci_upper,
    `P(>0)` = p_gt_0
  ) %>%
  mutate(
    Mean = round(Mean, 3),
    Median = round(Median, 3),
    `CI Lower` = round(`CI Lower`, 3),
    `CI Upper` = round(`CI Upper`, 3),
    `P(>0)` = round(`P(>0)` * 100, 1)
  )

write.csv(stat_table, "Prior_Sensitivity_Statistics.csv", row.names = FALSE)

cat("\n")
cat(paste0(rep("=", 60), collapse = ""))
cat("\n")
cat("PRIOR SENSITIVITY ANALYSIS RESULTS\n")
cat(paste0(rep("=", 60), collapse = ""))
cat("\n\n")

for (i in 1:nrow(stat_table)) {
  cat(sprintf("%s:\n", stat_table$Model[i]))
  cat(sprintf("  Mean effect: %.3f\n", stat_table$Mean[i]))
  cat(sprintf("  95%% CrI: [%.3f, %.3f]\n", stat_table$`CI Lower`[i], stat_table$`CI Upper`[i]))
  cat(sprintf("  Probability of positive effect: %.1f%%\n", stat_table$`P(>0)`[i]))
  cat("\n")
}

diff_means <- stat_table$Mean[1] - stat_table$Mean[2]
cat(sprintf("Difference in means: %.3f\n", diff_means))

if (abs(diff_means) < 0.05) {
  cat("Interpretation: Results are robust to prior specification (difference < 0.05)\n")
} else if (abs(diff_means) < 0.1) {
  cat("Interpretation: Moderate sensitivity to prior specification\n")
} else {
  cat("Interpretation: High sensitivity to prior specification\n")
}

cat("\n")
cat(paste0(rep("-", 60), collapse = ""))
cat("\n")
cat("Files created:\n")
cat("1. Posterior_Group1_Prior_Comparison_Enhanced_R.pdf/png - Enhanced plot\n")
cat("2. Posterior_Group1_Prior_Comparison_Simple_R.pdf - Simplified plot\n")
cat("3. Prior_Sensitivity_Statistics_R.csv - Statistical summary\n")
cat(paste0(rep("-", 60), collapse = ""))
cat("\n")


ggsave("Posterior_Group1_Prior_Comparison.pdf", plot = p, width = 8, height = 5)




# Step 6. Posterior Threshold Probability Plot (SE % Increase)
draws <- as_draws_df(fit_all)
baseline_mean <- mean(data$baseline[data$Group == 1], na.rm = TRUE)
delta_pct <- (draws$b_Group1 / baseline_mean) * 100
thresholds <- seq(0, 100, by = 1)
threshold_probs <- sapply(thresholds, function(th) mean(delta_pct > th))
threshold_df <- data.frame(Threshold = thresholds, Posterior_Prob = threshold_probs)

# Save plot
p <- ggplot(threshold_df, aes(x = Threshold, y = Posterior_Prob)) +
  geom_line(size = 1.2, color = "#1b9e77") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "grey40") +
  annotate("text", x = 10, y = 0.9, label = "≥10% increase\nClinically meaningful", 
           hjust = -0.1, size = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Posterior Probability of ArI Increase after FMT",
    x = "Percentage Increase in ArI (Δ%)",
    y = "Posterior Probability P(Δ > x%)"
  ) +
  theme_minimal(base_size = 14)
ggsave("R_Threshold_Posterior_Probability.pdf", plot = p, width = 8, height = 5)
write.csv(threshold_df, "R_Percent_Increase_Threshold.csv", row.names = FALSE)

# Step 16.1: Compute posterior probabilities for threshold-based reductions in SL (Sleep Latency)

# Extract posterior draws from the Bayesian model (fit_all)
draws <- as_draws_df(fit_all)

# Calculate the baseline mean SL in the treatment group (Group == 1)
baseline_mean <- mean(data$baseline[data$Group == 1], na.rm = TRUE)

# Compute percentage change in SL relative to baseline
# A negative delta_pct indicates SL reduction, which is considered improvement
delta_pct <- (draws$b_Group1 / baseline_mean) * 100

# Define a sequence of percentage reduction thresholds (e.g., 1% to 100%)
thresholds <- seq(1, 100, by = 1)

# For each threshold, calculate the posterior probability that SL decreased by at least x%
threshold_probs <- sapply(thresholds, function(th) mean(delta_pct < -th))

# Combine threshold and probability into a data frame
threshold_df <- data.frame(Threshold = thresholds, Posterior_Prob = threshold_probs)

# Generate a line plot showing the probability of achieving at least x% reduction in SL
p <- ggplot(threshold_df, aes(x = Threshold, y = Posterior_Prob)) +
  geom_line(size = 1.2, color = "#d95f02") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "grey40") +
  annotate("text", x = 10, y = 0.9, label = "≥10% reduction\nClinically meaningful", 
           hjust = -0.1, size = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Posterior Probability of ArI Reduction after FMT",
    x = "Percentage Reduction in SL (Δ%)",
    y = "Posterior Probability P(Δ < -x%)"
  ) +
  theme_minimal(base_size = 14)

# Save the figure as PDF
ggsave("Threshold_Posterior_Probability.pdf", plot = p, width = 8, height = 5)

# Export the threshold-probability table for reference or supplementary material
write.csv(threshold_df, "Percent_Reduction_Threshold.csv", row.names = FALSE)





# Step 7.---- Posterior draws of variance components ----
draws_all <- as_draws_df(fit_all)

# sd parameters from brms: sd_center__Intercept and sigma
stopifnot(all(c("sd_center__Intercept", "sigma") %in% names(draws_all)))

tau  <- draws_all$sd_center__Intercept           # posterior draws of SD for center intercepts
sig  <- draws_all$sigma                          # posterior draws of residual SD
tau2 <- tau^2
sig2 <- sig^2

# ICC as heterogeneity proportion (I^2 analog for LMM)
ICC  <- tau2 / (tau2 + sig2)

# Summaries (posterior mean and 95% CrI)
summ <- function(x) c(mean = mean(x), l95 = quantile(x, 0.025), u95 = quantile(x, 0.975))
tau_summ  <- summ(tau2)
sig_summ  <- summ(sig2)
icc_summ  <- summ(ICC)

heterogeneity_summary <- data.frame(
  Parameter = c("Between-center variance (tau^2)", "Residual variance (sigma^2)", "ICC (I^2 analog)"),
  Mean      = c(tau_summ["mean"], sig_summ["mean"], icc_summ["mean"]),
  L95CrI    = c(tau_summ["l95"],  sig_summ["l95"],  icc_summ["l95"]),
  U95CrI    = c(tau_summ["u95"],  sig_summ["u95"],  icc_summ["u95"])
)

write.csv(heterogeneity_summary, "center_heterogeneity_I2_summary.cREMsv", row.names = FALSE)


# ---- Center-specific random intercepts (conditional modes) ----
# ranef() returns posterior summaries for random effects by group level

# --- Inputs assumed: brms model `fit_all` already fitted ---

library(dplyr)
library(ggplot2)
library(brms)
library(posterior)  # for as_draws_df

# 1) Pick up data
build_center_forest_df <- function(fit, group = "center", coef_name = "Intercept") {
  re <- ranef(fit, summary = TRUE)[[group]]
  dn <- dimnames(re)
  
  stat_dim <- which(sapply(dn, function(x) all(c("Estimate","Est.Error","Q2.5","Q97.5") %in% x)))
  coef_dim <- which(sapply(dn, function(x) coef_name %in% x))
  grp_dim  <- setdiff(1:3, c(stat_dim, coef_dim))
  
  get_vals <- function(stat) {
    idx <- list(NULL,NULL,NULL)
    idx[[grp_dim]]  <- seq_len(dim(re)[grp_dim])
    idx[[stat_dim]] <- which(dn[[stat_dim]] == stat)
    idx[[coef_dim]] <- which(dn[[coef_dim]] == coef_name)
    as.numeric(do.call(`[`, c(list(re), idx)))
  }
  
  centers <- dn[[grp_dim]]
  out <- data.frame(
    center   = centers,
    estimate = get_vals("Estimate"),
    lower    = get_vals("Q2.5"),
    upper    = get_vals("Q97.5"),
    stringsAsFactors = FALSE
  )
  out
}

forest_df <- build_center_forest_df(fit_all, group = "center", coef_name = "Intercept") %>%
  arrange(estimate) %>%
  mutate(center = factor(center, levels = center),
         label_ci = sprintf("%.3f  (%.3f to %.3f)", estimate, lower, upper))

# 2) ICC ≈ tau^2 / (tau^2 + sigma^2)
draws <- as_draws_df(fit_all)
tau   <- draws$sd_center__Intercept       # SD of center intercepts
sig   <- draws$sigma                      # residual SD
tau2  <- tau^2;  sig2 <- sig^2
ICC   <- tau2 / (tau2 + sig2)

icc_mean <- mean(ICC)
icc_l95  <- quantile(ICC, 0.025)
icc_u95  <- quantile(ICC, 0.975)

# 3) 
x_rng <- range(c(forest_df$lower, forest_df$upper), na.rm = TRUE)
pad   <- diff(x_rng) * 0.25
x_min <- x_rng[1] - 0.05 * diff(x_rng)
x_max <- x_rng[2] + pad

# ==== Add sample size per center (robust & clean) ====

center_n <- fit_all$data %>%
  dplyr::count(center, name = "n_center") %>%
  dplyr::mutate(
    center = as.character(center),
    n_center = as.integer(n_center)
  )

forest_df <- forest_df %>%
  dplyr::mutate(center = as.character(center)) %>%
  dplyr::left_join(center_n, by = "center") %>%
  dplyr::mutate(
    n_center = ifelse(is.na(n_center), 0L, n_center),
    center_label = sprintf("%s (n=%d)", center, n_center)
  ) %>%
  dplyr::arrange(estimate) %>%
  dplyr::mutate(center_label = factor(center_label, levels = center_label))

p_icc <- ggplot(forest_df, aes(x = estimate, y = center_label)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, size = 0.5) +
  geom_point(size = 1.8, shape = 16) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text(aes(x = x_rng[2] + pad * 0.05, label = label_ci),
            hjust = 0, size = 3.2) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  labs(
    title = "Center-level effects from random-intercept model",
    subtitle = "Points: posterior mean; bars: 95% credible intervals",
    x = "Center random intercept (identity scale)",
    y = "Center",
    caption = sprintf("ICC (I^2 analog) = %.3f [%.3f, %.3f]", icc_mean, icc_l95, icc_u95)
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.title.y    = element_text(margin = margin(r = 8)),
    axis.title.x    = element_text(margin = margin(t = 6)),
    axis.text.y     = element_text(size = 9),
    axis.line       = element_line(linewidth = 0.4),
    axis.ticks      = element_line(linewidth = 0.3),
    plot.margin     = margin(t = 10, r = 80, b = 10, l = 10),
    legend.position = "none"
  )

# Display forest plot in R session
print(p_icc)

ggsave(
  filename = "forest_center_random_intercept_with.pdf",
  plot     = p_icc,
  width    = 7,
  height   = 6,
  device   = cairo_pdf
)


