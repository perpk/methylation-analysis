rm(list = ls())
gc(full = TRUE)


library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(FlowSorted.Blood.EPIC)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE145361_20260714_080452"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets_gender_mismatch <- readRDS(file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"))

beta_matrix <- readRDS(file.path(project_context$paths$results, "beta_matrix.rds"))

pca_df <- readRDS(file.path(project_context$paths$results, "pca_df.rds"))

pca_scores <- pca_df[, 1:2]
center <- colMeans(pca_scores)
cov_matrix <- cov(pca_scores)
distances <- mahalanobis(pca_scores, center, cov_matrix)
outlier_threshold <- qchisq(0.975, df = 2)
outliers <- which(distances > outlier_threshold)

outlier_samples <- rownames(pca_scores)[outliers]
pca_df$Is_Outlier <- rownames(pca_df) %in% outlier_samples

beta_matrix_no_outliers <- beta_matrix[, !colnames(beta_matrix) %in% outlier_samples]

sum(beta_matrix_no_outliers == 0)
saveRDS(beta_matrix_no_outliers, file.path(project_context$paths$results, "beta_matrix_no_outliers.rds"))

source("R/results_container.R")
source("R/apply_BMIQ.R")
res <- apply_BMIQ(project_context,
  beta_matrix = beta_matrix_no_outliers,
  plot = FALSE
)

m_bmiq_no_outliers <- readRDS(file.path(project_context$paths$results, "m_values_bmiq.rds"))
beta_bmiq_no_outliers <- readRDS(file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))

targets_gender_mismatch$Basename %>% head()

rownames(targets_gender_mismatch) <- targets_gender_mismatch$Basename %>% str_remove(paste0("GSE145361_RAW", "/"))

pd_samples <- rownames(targets_gender_mismatch[targets_gender_mismatch$Sample_Group == "PD", ])
hc_samples <- rownames(targets_gender_mismatch[targets_gender_mismatch$Sample_Group == "Control", ])

colnames_pd <- which(colnames(beta_bmiq_no_outliers) %in% pd_samples)
colnames_hc <- which(colnames(beta_bmiq_no_outliers) %in% hc_samples)

mean_beta_pd <- rowMeans(beta_bmiq_no_outliers[, colnames_pd], na.rm = TRUE)
mean_beta_hc <- rowMeans(beta_bmiq_no_outliers[, colnames_hc], na.rm = TRUE)
delta_beta <- mean_beta_pd - mean_beta_hc

beta_means <- data.frame(
  mean_beta_pd,
  mean_beta_hc,
  delta_beta
)

write.csv(beta_means, file.path(project_context$paths$results, "beta_means.csv"))

targets <- readRDS(file.path(project_context$paths$processed, "targets.rds"))

source("R/extract_methyl_set.R")
source("R/intermediate_data_proxy.R")
res_extract_methyl_set <- intermediate_data_proxy(
  extract_methyl_set, project_context,
  targets = targets
)

source("R/cell_cnt_estimate.R")
res_cell_cnt_estimate <- intermediate_data_proxy(
  cell_cnt_estimate,
  project_context,
  rg_set = res_extract_methyl_set$rg_set_container@object,
  targets = targets
)
saveRDS(res_cell_cnt_estimate$targets_container@object, file.path(project_context$paths$qc, "targets_s_mismatch_cells.rds"))

source("R/plot_cell_proportions.R")
plot_cell_proportions(
  context = project_context,
  targets = res_cell_cnt_estimate$targets_container@object,
  targets_file = res_cell_cnt_estimate$targets_container@filename,
  cell_types = getCellTypesForPlatform(project_context$platform)
)
