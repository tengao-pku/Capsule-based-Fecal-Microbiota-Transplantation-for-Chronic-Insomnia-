setwd("D:/R/Rdata/PHD-main work/Part2_FMT/PSG/PSGmice10C3/SE")# ---- Packages ----
library(tidyverse)
library(readr)
library(mice)
library(emmeans)

# ============================================================
# Helper: Rubin pooling for scalar estimates (robust to NA)
# ============================================================
pool_scalar_rubin <- function(est_vec, se_vec, conf_level = 0.95) {
  ok <- is.finite(est_vec) & is.finite(se_vec)
  Q <- est_vec[ok]
  U <- se_vec[ok]^2
  m <- length(Q)
  stopifnot(m >= 2)
  
  qbar <- mean(Q)
  ubar <- mean(U)
  b    <- var(Q)
  tvar <- ubar + (1 + 1/m) * b
  se   <- sqrt(tvar)
  
  if (b < .Machine$double.eps || ubar < .Machine$double.eps) {
    df <- Inf
  } else {
    r  <- ((1 + 1/m) * b) / ubar
    df <- (m - 1) * (1 + 1/r)^2
  }
  
  alpha <- 1 - conf_level
  tcrit <- if (is.finite(df)) qt(1 - alpha/2, df) else qnorm(1 - alpha/2)
  
  ci_low  <- qbar - tcrit * se
  ci_high <- qbar + tcrit * se
  
  tval <- qbar / se
  pval <- if (is.finite(df)) 2 * pt(abs(tval), df = df, lower.tail = FALSE) else 2 * pnorm(abs(tval), lower.tail = FALSE)
  
  data.frame(
    estimate  = qbar,
    std.error = se,
    df        = df,
    statistic = tval,
    p.value   = pval,
    CI_low    = ci_low,
    CI_high   = ci_high,
    m_used    = m
  )
}

# ============================================================
# 1) Read data
# ============================================================
file <- "SE.csv"
dat0 <- read_csv(file, show_col_types = FALSE)

# ============================================================
# 2) Factor coding (centre as factor for dummy variables)
# ============================================================
dat0 <- dat0 %>%
  mutate(
    center   = factor(center),
    Group    = factor(Group, levels = c(0, 1), labels = c("nonFMT", "FMT")),
    subgroup = factor(subgroup, levels = c(0, 1), labels = c("Placebo", "Synbiotic"))
  )

# Missingness overview
miss_tbl <- sapply(dat0[, c("post","baseline","Group","subgroup","center")], function(x) sum(is.na(x)))
print(miss_tbl)

# ============================================================
# 3) Multiple Imputation (MICE)
#   - Following Lancet logic: MI under MAR using baseline + design variables as predictors
#   - NOTE: Lancet example used n=3000 imputations; that is usually heavy.
#         If you truly want 3000, set m = 3000 below.
# ============================================================
vars_for_mi <- c("post", "baseline", "Group", "subgroup", "center")
dat_mi <- dat0 %>% select(all_of(vars_for_mi))

ini  <- mice(dat_mi, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

meth[c("baseline", "Group", "subgroup", "center")] <- ""
meth["post"] <- "pmm"

pred[,] <- 0
pred["post", c("baseline","Group","subgroup","center")] <- 1

set.seed(20251220)
imp <- mice(
  dat_mi,
  m = 20,            
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  printFlag = TRUE
)

# ============================================================
# 4) Primary analysis: baseline ANCOVA with centre fixed (dummy variables)
#   Equivalent to "stratified for study centre using dummy variables"
# ============================================================
fit_mi <- with(
  imp,
  lm(post ~ Group * subgroup + baseline + center)  # <-- centre fixed effect
)

fit_list <- fit_mi$analyses
m_total  <- length(fit_list)

# Pool ANCOVA coefficients (Rubin via mice::pool)
pool_fit <- pool(fit_mi)

pool_sum <- summary(pool_fit)
print(pool_sum)
write.csv(pool_sum, "ANCOVA_MI_pooled_fixed_effects.csv", row.names = FALSE)

pool_ci <- summary(pool_fit, conf.int = TRUE, conf.level = 0.95)
print(pool_ci)
write.csv(pool_ci, "ANCOVA_MI_pooled_fixed_effects_CI95.csv", row.names = FALSE)

# ============================================================
# 5) Primary treatment contrast: FMT - nonFMT marginal over subgroup (Lancet-style "adjusted difference")
#   Rationale: if Group*subgroup exists, do marginal contrast to represent overall treatment effect.
# ============================================================
extract_marginal_group_contrast <- function(model) {
  emm_g <- emmeans(model, ~ Group)  # averaged over subgroup and centre distribution in design matrix
  con   <- contrast(emm_g, method = list("FMT - nonFMT" = c(-1, 1)))
  dfc   <- as.data.frame(con)
  dfc %>% select(contrast, estimate, SE, df)
}

mg_list <- lapply(fit_list, extract_marginal_group_contrast)

mg_long <- bind_rows(lapply(seq_along(mg_list), function(i) {
  d <- mg_list[[i]]
  d$imp <- i
  d
}))

mg_pooled <- pool_scalar_rubin(mg_long$estimate, mg_long$SE, conf_level = 0.95)

# Effect size (Lancet): adjusted difference / SD at baseline (total sample)
sd_baseline_total <- sd(dat0$baseline, na.rm = TRUE)
effect_size <- mg_pooled$estimate / sd_baseline_total

treat_effect_report <- data.frame(
  Outcome     = "TST",
  Contrast    = "FMT - nonFMT (marginal over subgroup; pooled via Rubin)",
  Estimate    = mg_pooled$estimate,
  SE          = mg_pooled$std.error,
  df          = mg_pooled$df,
  t           = mg_pooled$statistic,
  p           = mg_pooled$p.value,
  CI_low      = mg_pooled$CI_low,
  CI_high     = mg_pooled$CI_high,
  EffectSize  = effect_size,
  Baseline_SD = sd_baseline_total,
  m_used      = mg_pooled$m_used
)

print(treat_effect_report)
write.csv(treat_effect_report, "ANCOVA_MI_treatment_effect_marginal_group_contrast.csv", row.names = FALSE)

# ============================================================
# 6) Within-arm subgroup contrasts (Synbiotic - Placebo within each Group), pooled by Rubin
# ============================================================
extract_within_arm_contrasts <- function(model) {
  emm_within <- emmeans(model, ~ subgroup | Group)
  con <- contrast(emm_within, method = list("Synbiotic - Placebo" = c(-1, 1)))
  as.data.frame(con) %>% select(Group, contrast, estimate, SE, df)
}

all_contrasts <- lapply(fit_list, extract_within_arm_contrasts)

contr_long <- bind_rows(lapply(seq_along(all_contrasts), function(i) {
  d <- all_contrasts[[i]]
  d$imp <- i
  d
}))

pool_one_arm <- function(df_arm, arm_label) {
  pooled <- pool_scalar_rubin(df_arm$estimate, df_arm$SE, conf_level = 0.95)
  data.frame(
    Outcome  = "TST",
    Arm      = arm_label,
    Contrast = unique(df_arm$contrast),
    Estimate = pooled$estimate,
    SE       = pooled$std.error,
    df       = pooled$df,
    t        = pooled$statistic,
    p        = pooled$p.value,
    CI_low   = pooled$CI_low,
    CI_high  = pooled$CI_high,
    m_used   = pooled$m_used
  )
}

within_arm_pooled <- bind_rows(
  pool_one_arm(contr_long %>% filter(Group == "FMT"),    "FMT"),
  pool_one_arm(contr_long %>% filter(Group == "nonFMT"), "nonFMT")
)

print(within_arm_pooled)
write.csv(within_arm_pooled, "ANCOVA_MI_within_arm_subgroup_contrasts.csv", row.names = FALSE)

# ============================================================
# 7) ANCOVA assumption checks (Lancet-style residual diagnostics)
#   - Inspect skewness/kurtosis of residuals
#   - Regress absolute residuals on fitted values (heteroscedasticity check)
#   We do this on one representative completed dataset (imp #1) for reporting.
# ============================================================
skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s < .Machine$double.eps) return(NA_real_)
  mean(((x - m)/s)^3)
}
kurtosis_excess <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s < .Machine$double.eps) return(NA_real_)
  mean(((x - m)/s)^4) - 3
}

dat_imp1 <- complete(imp, action = 1)
fit_imp1 <- lm(post ~ Group * subgroup + baseline + center, data = dat_imp1)

resid1  <- residuals(fit_imp1)
fitted1 <- fitted(fit_imp1)

diag_tbl <- data.frame(
  resid_skewness = skewness(resid1),
  resid_kurtosis_excess = kurtosis_excess(resid1)
)

# regression of absolute residuals on predicted values
abs_resid_lm <- lm(abs(resid1) ~ fitted1)
abs_resid_sum <- summary(abs_resid_lm)

diag_out <- list(
  residual_moments = diag_tbl,
  abs_resid_on_fitted = data.frame(
    Estimate = coef(abs_resid_lm)[2],
    StdError = abs_resid_sum$coefficients[2,2],
    t        = abs_resid_sum$coefficients[2,3],
    p        = abs_resid_sum$coefficients[2,4],
    R2       = abs_resid_sum$r.squared
  )
)

print(diag_out$residual_moments)
print(diag_out$abs_resid_on_fitted)
write.csv(diag_out$residual_moments, "ANCOVA_assumption_residual_moments_imp1.csv", row.names = FALSE)
write.csv(diag_out$abs_resid_on_fitted, "ANCOVA_assumption_absresid_on_fitted_imp1.csv", row.names = FALSE)

# ============================================================
# 8) Sensitivity analysis: Complete-case ANCOVA (same model specification)
# ============================================================
dat_cc <- dat0 %>%
  filter(!is.na(post), !is.na(baseline), !is.na(Group), !is.na(subgroup), !is.na(center))

fit_cc <- lm(post ~ Group * subgroup + baseline + center, data = dat_cc)
cc_sum <- summary(fit_cc)
print(cc_sum)

cc_fixed <- broom::tidy(fit_cc, conf.int = TRUE)
write.csv(cc_fixed, "ANCOVA_complete_case_fixed_effects.csv", row.names = FALSE)

# Complete-case within-arm contrasts
emm_within_cc <- emmeans(fit_cc, ~ subgroup | Group)
within_cc <- as.data.frame(contrast(emm_within_cc, method = list("Synbiotic - Placebo" = c(-1, 1))))
write.csv(within_cc, "ANCOVA_complete_case_within_arm_subgroup_contrasts.csv", row.names = FALSE)

cat("\nDONE.\n",
    "- Primary: MI + baseline ANCOVA with centre fixed (dummy variables)\n",
    "- Primary contrast: FMT - nonFMT marginal over subgroup pooled via Rubin\n",
    "- Effect size: adjusted difference / SD at baseline (total sample)\n",
    "- ANCOVA assumption checks: residual skewness/kurtosis + abs(resid)~fitted\n",
    "- Sensitivity: complete-case ANCOVA (+ within-arm contrasts)\n", sep = "")
