rm(list = ls())
gc(full = TRUE)

library(tidyverse)

setClass(
    "DataSources",
    slots = list(
        targets = "ANY",
        m_values = "ANY",
        pca_df = "ANY"
    )
)

data_sources <- list(
    sgpd = new("DataSources",
        targets = "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/processed/GSE145361_harmonized_targets.rds",
        m_values = "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/results/m_values_bmiq_no_outliers.rds",
        pca_df = "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/processed/GSE145361_pca_df_enriched_meta.rds"
    ),
    peg1 = new("DataSources",
        targets = "/root/workspace/methyl-pipe-out/peg1_harmonized_targets.rds",
        m_values = "/root/workspace/methyl-pipe-out/peg1_m_values_bmiq_no_outliers.rds",
        pca_df = "/root/workspace/methyl-pipe-out/peg1_pca_df_with_outliers.rds"
    ),
    ppmi = new("DataSources",
        targets = "/root/workspace/methyl-pipe-out/ppmi_harmonized_targets.rds",
        m_values = "/root/workspace/methyl-pipe-out/ppmi_m_values_bmiq_no_outliers.rds",
        pca_df = "/root/workspace/methyl-pipe-out/ppmi_pca_df_with_outliers.rds"
    )
)

m_values_sgpd <- readRDS(data_sources$sgpd@m_values)
targets_sgpd <- readRDS(data_sources$sgpd@targets)
pca_df_sgpd <- readRDS(data_sources$sgpd@pca_df)

m_values_peg1 <- readRDS(data_sources$peg1@m_values)
targets_peg1 <- readRDS(data_sources$peg1@targets)
pca_df_peg1 <- readRDS(data_sources$peg1@pca_df)

m_values_ppmi <- readRDS(data_sources$ppmi@m_values)
targets_ppmi <- readRDS(data_sources$ppmi@targets)
pca_df_ppmi <- readRDS(data_sources$ppmi@pca_df)

m_values_ppmi %>% head()
m_values_peg1 %>% head()
m_values_sgpd %>% head()

common_cpgs <- Reduce(intersect, list(rownames(m_values_sgpd), rownames(m_values_peg1), rownames(m_values_ppmi)))
common_cpgs %>% length()

m_values_sgpd_common <- m_values_sgpd[common_cpgs, ]
m_values_peg1_common <- m_values_peg1[common_cpgs, ]
m_values_ppmi_common <- m_values_ppmi[common_cpgs, ]

# check whether there are duplicated sample IDs among the cohorts (unlikely but better safe than sorry).
common_samples <- Reduce(intersect, list(colnames(m_values_sgpd_common), colnames(m_values_peg1_common), colnames(m_values_ppmi_common)))
common_samples %>% length()

m_values_combined <- cbind(
    m_values_sgpd_common[common_cpgs, ],
    m_values_peg1_common[common_cpgs, ],
    m_values_ppmi_common[common_cpgs, ]
)

targets_sgpd %>% head()

targets_sgpd_reduced <- targets_sgpd[, !colnames(targets_sgpd) %in% "ScanDate"]
targets_sgpd_reduced$Age <- NA
targets_sgpd_reduced$Sex <- targets_sgpd_reduced$`gender:ch1`

targets_peg1 %>% head()
targets_peg1$Age <- targets_peg1$`age:ch1`
targets_peg1 %>% pull(Age) %>% head()
targets_peg1_reduced <- targets_peg1[, !colnames(targets_peg1) %in% c("age:ch1", "ScanDate")]

targets_ppmi %>% head()
targets_ppmi$Age <- targets_ppmi$ENROLL_AGE %>% round()
targets_ppmi$Sex <- targets_ppmi$SEX
targets_ppmi_reduced <- targets_ppmi[, !colnames(targets_ppmi) %in% c("ENROLL_AGE", "ScanDate", "SEX")]

targets_common_cols <- Reduce(intersect, list(colnames(targets_sgpd_reduced), colnames(targets_peg1_reduced), colnames(targets_ppmi_reduced)))
targets_common_cols %>% length()
targets_common_cols

targets_ppmi_reduced <- targets_ppmi_reduced[, targets_common_cols]
targets_peg1_reduced <- targets_peg1_reduced[, targets_common_cols]
targets_sgpd_reduced <- targets_sgpd_reduced[, targets_common_cols]

targets_ppmi_reduced$Cohort <- "PPMI"
targets_peg1_reduced$Cohort <- "PEG1"
targets_sgpd_reduced$Cohort <- "SGPD"

targets_combined <- rbind(targets_sgpd_reduced, targets_peg1_reduced, targets_ppmi_reduced)

targets_combined %>% head()

saveRDS(m_values_combined, file = "/workspace/results/consolidated_m_values.rds")
saveRDS(targets_combined, file = "/workspace/results/consolidated_targets.rds")