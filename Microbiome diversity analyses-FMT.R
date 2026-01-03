# ======== Alpha & Beta Diversity ========
# Author: Teng Gao
#Version:2025.08.01
# 
# -- Packages --
library(vegan)        # diversity, vegdist, adonis2
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)       # ggboxplot
library(zCompositions)  # zero replacement for CLR
library(compositions)   # clr transform

# -- Paths --
setwd("D:/R/Rdata/FMT/16s")  

# -- Load --
asv_data <- read.csv("nonFMT.csv", header = TRUE, check.names = FALSE)
stopifnot(colnames(asv_data)[1] == "Group")
asv_data$Group <- factor(trimws(as.character(asv_data$Group)),
                         levels = c("B","1"))  

# non-rarefied count tables
otu <- as.data.frame(lapply(asv_data[ , -1, drop = FALSE], as.numeric))
rownames(otu) <- paste0("S", seq_len(nrow(otu)))
total_reads <- rowSums(otu)
asv_data$total_reads <- total_reads

# ========= Alpha diversity =========
rel_abund <- sweep(otu, 1, rowSums(otu), FUN = "/")
rel_abund[is.na(rel_abund)] <- 0

Observed <- rowSums(otu > 0)
Shannon  <- diversity(rel_abund, index = "shannon")
Simpson  <- diversity(rel_abund, index = "simpson")

alpha_df <- data.frame(
  Sample      = rownames(otu),
  Group       = asv_data$Group,
  total_reads = total_reads,
  Observed    = Observed,
  Shannon     = Shannon,
  Simpson     = Simpson
)
#  Chao1：
Chao1 <- apply(otu, 1, function(x){
  S_obs <- sum(x > 0); F1 <- sum(x == 1); F2 <- sum(x == 2)
  if (F2 == 0) return(S_obs) else return(S_obs + (F1^2)/(2*F2))
})
alpha_df$Chao1 <- Chao1

# save alpha summary
write.csv(alpha_df, "alpha_diversity_summaryFMT.csv", row.names = FALSE)

# draw figures
alpha_long <- melt(alpha_df,
                   id.vars = c("Sample","Group"),
                   measure.vars = c("Observed","Shannon","Simpson","Chao1"),
                   variable.name = "Index", value.name = "Diversity")


jama_palette <- c("B" = "black", "1" = "#E69F00")#E69F00#1f77b4


box_plot <- function(index) {
  ggboxplot(
    subset(alpha_long, Index == index),
    x = "Group", y = "Diversity",
    color = "Group", palette = jama_palette,
    add = "jitter", title = index
  ) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    theme_minimal(base_size = 12)
}

# output
pdf("nonAlpha_nonrarefied_B_vs_1.pdf", width = 8, height = 6)
print(box_plot("Observed"))
print(box_plot("Shannon"))
print(box_plot("Simpson"))
print(box_plot("Chao1"))
dev.off()



# ========= Beta diversity（Aitchison/CLR） =========

otu_czm <- cmultRepl(as.matrix(otu), label = 0, method = "CZM")
clr_mat <- t(apply(otu_czm, 1, compositions::clr))
dist_ait <- dist(clr_mat, method = "euclidean")

# PCoA
pcoa <- cmdscale(dist_ait, k = 3, eig = TRUE)
coords <- as.data.frame(pcoa$points)
colnames(coords) <- c("PC1","PC2","PC3")
coords$Group <- asv_data$Group
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 2)

# PERMANOVA
adon <- adonis2(dist_ait ~ Group + log10(total_reads),
                data = asv_data, permutations = 999, by = "margin")
cat("\n== PERMANOVA (Aitchison) with total_reads covariate ==\n")
print(adon)

# 2D PCoA（PC1 vs PC2）
p_beta <- ggplot(coords, aes(PC1, PC2, color = Group)) +
  stat_ellipse(geom = "path", size = 0.9, alpha = 0.9) +
  geom_point(alpha = 0.9, size = 2.6) +
  scale_color_manual(values = jama_palette) +
  labs(
    title = sprintf("PCoA (Aitchison)\nPERMANOVA R²(Group)=%.3f, p=%s",
                    adon$R2[1], format.pval(adon$`Pr(>F)`[1], digits = 3)),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
ggsave("nonPCoA_Aitchison_B_vs_1.pdf", p_beta, width = 7.5, height = 6)

# =========Bray–Curtis  =========
bc <- vegdist(sweep(otu, 1, rowSums(otu), "/"), method = "bray")
adon_bc <- adonis2(bc ~ Group + log10(total_reads), data = asv_data, permutations = 999, by = "margin")
cat("\n== PERMANOVA (Bray–Curtis; relative abundance) ==\n")
print(adon_bc)


# ---- Packages ----
library(vegan)          # adonis2, diversity, vegdist
library(ggplot2)
library(dplyr)
library(zCompositions)  # cmultRepl
library(compositions)   # clr

# ---- Paths ----
setwd("D:/R/Rdata/FMT/16s")   

# ---- Load data ----
asv_data <- read.csv("nonFMT.csv", header = TRUE, check.names = FALSE)
stopifnot(colnames(asv_data)[1] == "Group")
asv_data$Group <- factor(trimws(as.character(asv_data$Group)), levels = c("B","1"))

# Count table (non-rarefied)
otu <- as.data.frame(lapply(asv_data[, -1, drop = FALSE], as.numeric))
rownames(otu) <- paste0("S", seq_len(nrow(otu)))
asv_data$total_reads <- rowSums(otu)

# ============================================================
# 1) Feature prevalence filter  (recommended to stabilize CLR)
#    keep ASVs present in >= 5% samples (you can set 0.01–0.10)
# ============================================================
prev_thr <- 0.05
prev <- colSums(otu > 0) / nrow(otu)
keep_feat <- prev >= prev_thr
otu_f <- otu[, keep_feat, drop = FALSE]

# Remove all-zero samples after filtering (rare, but safer)
keep_samp <- rowSums(otu_f) > 0
otu_f <- otu_f[keep_samp, , drop = FALSE]
asv_data_f <- asv_data[keep_samp, , drop = FALSE]

cat(sprintf("\n[INFO] Prevalence >= %.0f%% kept %d/%d features; kept %d/%d samples.\n",
            prev_thr*100, ncol(otu_f), ncol(otu), nrow(otu_f), nrow(otu)))

# ============================================================
# 2) Zero replacement (no deletion) -> 3) CLR transform
# ============================================================
otu_czm <- zCompositions::cmultRepl(
  as.matrix(otu_f),
  label    = 0,
  method   = "CZM",
  z.warning = FALSE,
  z.delete  = FALSE   # <-- critical: do NOT delete rows/cols
)

clr_mat  <- t(apply(otu_czm, 1, compositions::clr))
dist_ait <- dist(clr_mat, method = "euclidean")  # Aitchison

# ============================================================
# 4) PCoA + coords aligned with filtered metadata
# ============================================================
pcoa <- cmdscale(dist_ait, k = 3, eig = TRUE)
coords <- as.data.frame(pcoa$points)
stopifnot(nrow(coords) == nrow(asv_data_f))
colnames(coords) <- c("PC1","PC2","PC3")
coords$Group <- asv_data_f$Group
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 2)

# ============================================================
# 5) PERMANOVA with sequencing depth as covariate
# ============================================================
adon <- adonis2(
  dist_ait ~ Group + log10(total_reads),
  data = asv_data_f,
  permutations = 999,
  by = "margin"
)
cat("\n== PERMANOVA (Aitchison) with total_reads covariate ==\n")
print(adon)

# Save PERMANOVA table
adon_tbl <- as.data.frame(adon)
adon_tbl$Term <- rownames(adon_tbl); rownames(adon_tbl) <- NULL
adon_tbl <- adon_tbl[, c("Term", setdiff(names(adon_tbl), "Term"))]
write.csv(adon_tbl, "PERMANOVA_Aitchison_prevalenceFiltered.csv", row.names = FALSE)

# ============================================================
# 6) Plot PCoA (PC1 vs PC2) — with R2 & P in title
# ============================================================
jama_palette <- c("B" = "black", "1" = "#1f77b4")
p_beta <- ggplot(coords, aes(PC1, PC2, color = Group)) +
  stat_ellipse(geom = "path", size = 0.9, alpha = 0.9) +
  geom_point(alpha = 0.9, size = 2.6) +
  scale_color_manual(values = jama_palette) +
  labs(
    title = sprintf("PCoA (Aitchison; prevalence >= %.0f%%)\nPERMANOVA R²(Group)=%.3f, p=%s",
                    prev_thr*100, adon$R2[1], format.pval(adon$`Pr(>F)`[1], digits = 3)),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
ggsave("nonPCoA_Aitchison_prevalenceFiltered.pdf", p_beta, width = 7.5, height = 6)

# ============================================================
# 7) Sensitivity: Bray–Curtis on relative abundance (same filter)
# ============================================================
rel_abund_f <- sweep(otu_f, 1, rowSums(otu_f), "/"); rel_abund_f[is.na(rel_abund_f)] <- 0
bc <- vegdist(rel_abund_f, method = "bray")
adon_bc <- adonis2(bc ~ Group + log10(total_reads), data = asv_data_f, permutations = 999, by = "margin")
cat("\n== PERMANOVA (Bray–Curtis; filtered; relative abundance) ==\n")
print(adon_bc)

# Save BC PERMANOVA
bc_tbl <- as.data.frame(adon_bc)
bc_tbl$Term <- rownames(bc_tbl); rownames(bc_tbl) <- NULL
bc_tbl <- bc_tbl[, c("Term", setdiff(names(bc_tbl), "Term"))]
write.csv(bc_tbl, "PERMANOVA_Bray_prevalenceFiltered.csv", row.names = FALSE)
