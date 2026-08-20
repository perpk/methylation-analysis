rm(list = ls())
gc(full = TRUE)

library(tidyverse)

m_values_combined <- readRDS("/workspace/results/consolidated_m_values.rds")
targets_combined <- readRDS("/workspace/results/consolidated_targets.rds")

targets_combined$Cohort <- as.factor(targets_combined$Cohort)
targets_combined$Sample_Group <- as.factor(targets_combined$Sample_Group)
targets_combined$Sex <- as.factor(targets_combined$Sex)

design <- model.matrix(~ Sample_Group + Sex + CD4T + CD8T + NK + Bcell + Mono + Gran, data = targets_combined)

library(sva)

m_values_combat_corr <- ComBat(
    dat = m_values_combined,
    batch = targets_combined$Cohort,
    mod = design,
    par.prior = TRUE,
    ref.batch = "PEG1",
    prior.plots = FALSE
)

pca_prev <- prcomp(t(m_values_combined), center = TRUE, scale. = TRUE)
pca_df_before <- data.frame(pca_prev$x[, 1:2], Cohort = targets_combined$Cohort, Group = targets_combined$Sample_Group)

plot_before_combat <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = Cohort)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(
        title = "PCA of M-Values Before ComBat Harmonization"
    )
print(plot_before_combat)

ggsave(plot_before_combat, filename = "/workspace/results/pca_before_combat.png", width = 8, height = 6, dpi = 300)

pca_after <- prcomp(t(m_values_combat_corr), center = TRUE, scale. = TRUE)
pca_df_after <- data.frame(pca_after$x[, 1:2], Cohort = targets_combined$Cohort, Group = targets_combined$Sample_Group)

plot_after_combat <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Cohort)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(
        title = "PCA of M-Values After ComBat Harmonization"
    )
print(plot_after_combat)
ggsave(plot_after_combat, filename = "/workspace/results/pca_after_combat.png", width = 8, height = 6, dpi = 300)


m_values_combat_corr_sgpd <- ComBat(
    dat = m_values_combined,
    batch = targets_combined$Cohort,
    mod = design,
    par.prior = TRUE,
    ref.batch = "SGPD",
    prior.plots = FALSE
)
saveRDS(m_values_combat_corr_sgpd, "/workspace/results/combat_corrected_m_values_sgpd.rds")

m_values_combat_corr_ppmi <- ComBat(
    dat = m_values_combined,
    batch = targets_combined$Cohort,
    mod = design,
    par.prior = TRUE,
    ref.batch = "PPMI",
    prior.plots = FALSE
)
saveRDS(m_values_combat_corr_ppmi, "/workspace/results/combat_corrected_m_values_ppmi.rds")
