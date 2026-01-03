setwd("D:/R/Rdata/PHD-main work/Part2_FMT/Scale data/1218scale")

# ---- Packages ----
library(tidyverse)
library(lme4)
library(lmerTest)   # p-values (Satterthwaite)
library(emmeans)
library(patchwork)

# ---- Load data (keep original column names) ----
# 你上传的文件也可用绝对路径："/mnt/data/ff65bf67-0f3f-4280-b11e-af2905690ffa.csv"
df <- read.csv("Analysis data for FMT-P20708.csv", check.names = FALSE)

# =========================================================
# 1) Reshape wide -> long (PSQI/ISI/SAS/SDS/GSRS across 0/1/2/3/6 months)
# =========================================================
scales <- c("PSQI","ISI","SAS","SDS","GSRS")

long_df <- df %>%
  select(
    ID, Group, Center, Sex, BMI, Age, Subgroup,
    matches(paste0("^(", paste(scales, collapse="|"), ")"))
  ) %>%
  pivot_longer(
    cols = matches(paste0("^(", paste(scales, collapse="|"), ")")),
    names_to = "ScaleMonth",
    values_to = "Score"
  ) %>%
  mutate(
    Scale = case_when(
      str_detect(ScaleMonth, "^PSQI") ~ "PSQI",
      str_detect(ScaleMonth, "^ISI")  ~ "ISI",
      str_detect(ScaleMonth, "^SAS")  ~ "SAS",
      str_detect(ScaleMonth, "^SDS")  ~ "SDS",
      str_detect(ScaleMonth, "^GSRS") ~ "GSRS",
      TRUE ~ NA_character_
    ),
    # 提取月份数字（0/1/2/3/6）
    MonthNum = str_extract(ScaleMonth, "\\d+") %>% as.integer(),
    Month = factor(
      MonthNum,
      levels = c(0, 1, 2, 3, 6),
      labels = c("Baseline", "1M", "2M", "3M", "6M")
    ),
    ID = factor(ID),
    Center = factor(Center),
    Group = factor(Group, levels = c(0, 1), labels = c("non-FMT", "FMT")),
    Subgroup = factor(Subgroup),
    Sex = factor(Sex)
  ) %>%
  filter(!is.na(Scale), !is.na(Month), !is.na(Score))

# Quick sanity check
table(long_df$Scale, long_df$Month, useNA = "ifany")

# =========================================================
# 2) Fit LMM per scale + extract EMMs (Group within Month) + model p-values
# =========================================================
fit_one_scale <- function(dat_scale) {
  # LMM: repeated measures by subject, adjust clustering by center
  # Month as categorical (most conservative, reviewer-friendly)
  fit <- lmer(Score ~ Group * Month + (1 | ID) + (1 | Center),
              data = dat_scale, REML = FALSE)
  
  # EMMs: model-estimated marginal means per Group at each Month
  emm <- emmeans(fit, ~ Group | Month)
  emm_df <- as.data.frame(emm) %>%
    mutate(
      Scale = unique(dat_scale$Scale),
      Month = factor(Month, levels = levels(long_df$Month))
    )
  
  # Key tests: Group main effect, Month, interaction
  an <- anova(fit) %>% as.data.frame()
  an$Effect <- rownames(an); rownames(an) <- NULL
  an <- an %>% mutate(Scale = unique(dat_scale$Scale))
  
  list(fit = fit, emm = emm_df, anova = an)
}

# Run all scales
res_list <- long_df %>%
  split(.$Scale) %>%
  lapply(fit_one_scale)

emm_all <- bind_rows(lapply(res_list, `[[`, "emm"))
anova_all <- bind_rows(lapply(res_list, `[[`, "anova"))

# Save statistics tables
write.csv(emm_all, "LMM_EMMeans_byScale_GroupMonth.csv", row.names = FALSE)
write.csv(anova_all, "LMM_ANOVA_byScale.csv", row.names = FALSE)

# Also save fixed effects tables (optional but useful)
fixed_all <- bind_rows(lapply(names(res_list), function(s){
  fit <- res_list[[s]]$fit
  fx <- summary(fit)$coefficients %>% as.data.frame()
  fx$Term <- rownames(fx); rownames(fx) <- NULL
  fx$Scale <- s
  fx
}))
write.csv(fixed_all, "LMM_FixedEffects_byScale.csv", row.names = FALSE)

# =========================================================
# 3) Plot settings (Lancet Psychiatry-like red/blue)
# =========================================================
col_blue <- "#1F5AA6"   # vivid blue
col_red  <- "#D11F2D"   # vivid red

group_cols <- c("non-FMT" = col_blue, "FMT" = col_red)

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )

plot_scale <- function(emm_df, title_text = NULL, ylab_text = "Score") {
  ggplot(emm_df, aes(x = Month, y = emmean, group = Group, color = Group)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.6) +
    scale_color_manual(values = group_cols) +
    labs(
      title = title_text %||% unique(emm_df$Scale),
      subtitle = "Model-estimated marginal means ± 95% CI",
      y = ylab_text
    ) +
    base_theme
}

# =========================================================
# 4) Build requested panels:
#    - Panel A: PSQI + ISI stacked
#    - Panel B: SAS + SDS + GSRS stacked
# =========================================================
emm_psqi <- emm_all %>% filter(Scale == "PSQI")
emm_isi  <- emm_all %>% filter(Scale == "ISI")
emm_sas  <- emm_all %>% filter(Scale == "SAS")
emm_sds  <- emm_all %>% filter(Scale == "SDS")
emm_gsrs <- emm_all %>% filter(Scale == "GSRS")

p_psqi <- plot_scale(emm_psqi, title_text = "PSQI")
p_isi  <- plot_scale(emm_isi,  title_text = "ISI")

p_sas  <- plot_scale(emm_sas,  title_text = "SAS")
p_sds  <- plot_scale(emm_sds,  title_text = "SDS")
p_gsrs <- plot_scale(emm_gsrs, title_text = "GSRS")

# Stack: PSQI + ISI (one figure)
fig_sleep <- p_psqi / p_isi + plot_layout(guides = "collect")
ggsave("Fig_Scales_PSQI_ISI_LMM_EMM.pdf", fig_sleep, width = 7.2, height = 7.6)

# Stack: SAS + SDS + GSRS (one figure)
fig_others <- p_sas / p_sds / p_gsrs + plot_layout(guides = "collect")
ggsave("Fig_Scales_SAS_SDS_GSRS_LMM_EMM.pdf", fig_others, width = 7.2, height = 10.5)

# Optional: also save PNGs for drafts
ggsave("Fig_Scales_PSQI_ISI_LMM_EMM.png", fig_sleep, width = 7.2, height = 7.6, dpi = 300)
ggsave("Fig_Scales_SAS_SDS_GSRS_LMM_EMM.png", fig_others, width = 7.2, height = 10.5, dpi = 300)

cat("\nDONE.\nOutputs:\n",
    "- LMM_ANOVA_byScale.csv (Group, Month, Group×Month tests)\n",
    "- LMM_FixedEffects_byScale.csv (coef table)\n",
    "- LMM_EMMeans_byScale_GroupMonth.csv (EMMs for plotting/report)\n",
    "- Fig_Scales_PSQI_ISI_LMM_EMM.pdf/png\n",
    "- Fig_Scales_SAS_SDS_GSRS_LMM_EMM.pdf/png\n")


# =========================================================
# Plot code: Lancet Psy red/blue, tight y-limits, square-like combined figure
# Requires: emm_all from emmeans (columns: Scale, Month, Group, emmean, lower.CL, upper.CL)
# =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ---- Lancet Psychiatry–style colors ----
col_blue <- "#003A8F"   # deep blue
col_red  <- "#C41E3A"   # crimson red
group_cols <- c("non-FMT" = col_blue, "FMT" = col_red)

# ---- Helper: compute tight y-limits from CI with padding ----
tight_ylim <- function(df, pad_frac = 0.08) {
  y_min <- min(df$lower.CL, na.rm = TRUE)
  y_max <- max(df$upper.CL, na.rm = TRUE)
  rng <- y_max - y_min
  if (!is.finite(rng) || rng <= 0) rng <- max(abs(c(y_min, y_max)), 1)
  pad <- rng * pad_frac
  c(y_min - pad, y_max + pad)
}

# ---- Theme (clean, publication) ----
base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4)
  )

# ---- Plot function ----
# baseline_deemph: make baseline points slightly transparent (optional)
plot_scale_tight <- function(emm_df, title_text = NULL, ylab_text = "Score",
                             pad_frac = 0.08, baseline_deemph = TRUE) {
  # ensure Month order is correct
  if (!is.factor(emm_df$Month)) {
    emm_df$Month <- factor(emm_df$Month, levels = c("Baseline","1M","2M","3M","6M"))
  } else {
    emm_df$Month <- factor(emm_df$Month, levels = c("Baseline","1M","2M","3M","6M"))
  }
  
  # baseline deemphasis
  emm_df <- emm_df %>% mutate(is_baseline = (Month == "Baseline"))
  
  yl <- tight_ylim(emm_df, pad_frac = pad_frac)
  
  ggplot(emm_df, aes(x = Month, y = emmean, group = Group, color = Group)) +
    geom_line(linewidth = 1.2) +
    geom_point(aes(alpha = ifelse(baseline_deemph & is_baseline, 0.65, 1)),
               size = 2.4, shape = 16) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL,
                      alpha = ifelse(baseline_deemph & is_baseline, 0.75, 1)),
                  width = 0.08, linewidth = 0.6) +
    scale_color_manual(values = group_cols) +
    scale_alpha_identity() +
    coord_cartesian(ylim = yl, clip = "off") +
    labs(
      title = title_text %||% unique(emm_df$Scale),
      subtitle = "Model-estimated marginal means (95% CI)",
      y = ylab_text
    ) +
    base_theme
}

# =========================================================
# Build PSQI + ISI stacked figure
# =========================================================

# --- assuming you already have emm_all from earlier code ---
# emm_psqi <- emm_all %>% filter(Scale == "PSQI")
# emm_isi  <- emm_all %>% filter(Scale == "ISI")

p_psqi <- plot_scale_tight(emm_psqi, title_text = "PSQI", pad_frac = 0.07)
p_isi  <- plot_scale_tight(emm_isi,  title_text = "ISI",  pad_frac = 0.07)

# Collect legend once; keep panels compact
fig_sleep <- (p_psqi / p_isi) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")

# ---- Save: make overall figure more "square" ----
# Adjust width/height: 7.2 x 7.2 gives square; legend on right increases width demand slightly
ggsave("Fig_Scales_PSQI_ISI_LMM_EMM_tightSquare.pdf",
       fig_sleep, width = 7.8, height = 7.2)

ggsave("Fig_Scales_PSQI_ISI_LMM_EMM_tightSquare.png",
       fig_sleep, width = 7.8, height = 7.2, dpi = 300)

print(fig_sleep)
contrast_all <- bind_rows(lapply(res_list, function(x){
  emmeans(x$fit, ~ Group | Month) %>%
    contrast(method = "revpairwise") %>%   # FMT - nonFMT
    as.data.frame()
})) %>%
  mutate(
    Scale = rep(names(res_list),
                each = nrow(.) / length(names(res_list)))
  )

write.csv(contrast_all, "LMM_GroupDiff_byScale_Month.csv", row.names = FALSE)
