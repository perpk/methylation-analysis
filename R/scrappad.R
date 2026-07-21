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
library(tidyverse)

source("R/progress_mgr.R")
source("R/project_context.R")
source("R/intermediate_data_proxy.R")
source("R/results_container.R")

project_to_load <- "ppmi_20260721_075730"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "EPIC"
data_folder <- "ppmi"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)
mode <- results_mode()$memory_only

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
project_context$mode <- mode
project_context$mode

beta_before <- readRDS("/root/workspace/methyl-pipe-out/ppmi_20260721_075730/results/beta_matrix.rds")
beta_bmiq <- readRDS("/root/workspace/methyl-pipe-out/ppmi_20260721_075730/results/beta_matrix_bmiq.rds")

project_context$paths$processed
source("R/apply_BMIQ.R")
plot_BMIQ(
  context = project_context,
  beta_before = beta_before,
  beta_after = beta_bmiq
) 

# Principal Component Analysis
source("R/principal_component_analysis.R")
source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = "ppmi")
col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- var_mappings$gender_var
col_map[["Age"]] <- var_mappings$age_var
keys <- c("Sample_Group", "Gender", "Age")

m_values_bmiq <- readRDS(file.path(project_context$paths$results, "m_values_bmiq.rds"))
targets <- readRDS(file.path(project_context$paths$processed, "targets_remove_mismatch.rds"))

m_bmiq_container <- new("ResultsContainer", filename = file.path(project_context$paths$results, "m_values_bmiq.rds"), object = m_values_bmiq, future = NULL)
targets_container <- new("ResultsContainer", filename = file.path(project_context$paths$processed, "targets_remove_mismatch.rds"), object = targets, future = NULL)  

bmiq_res <- list(
  m_bmiq_container = m_bmiq_container
)

res_bio_gender_mismatch <- list(
  targets_container = targets_container
)

res_pca <- intermediate_data_proxy(
  principal_component_analysis,
  auto_clean = FALSE,
  project_context,
  m_values = bmiq_res$m_bmiq_container@object,
  m_values_filename = bmiq_res$m_bmiq_container@filename,
  targets = res_bio_gender_mismatch$targets_container@object,
  targets_filename = res_bio_gender_mismatch$targets_container@filename,
  col_maps = col_map,
  keys = keys,
  npc = 5
)

saveRDS(res_pca$pca_df_container@object, file.path(project_context$paths$results, "pca_df.rds"))

source("R/plot_PCA.R")
pca_vars <- c("Sample_Group" = "By Diagnosis", "Gender" = "By Biological Gender", "Age" = "By Age")
convert_f <- function(df) {
  (transform(df, Age = as.numeric(Age)))
}
plot_PCA(
  context = project_context,
  pca_df = res_pca$pca_df_container@object,
  pca_vars = pca_vars,
  convert_fun = convert_f,
  continuously_scaled = c("Age"),
  "pca_plot_",
  pairplot_color_by = "Age"
)

res_pca$pca_df_container@object %>% head()

beta_bmiq <- readRDS(file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))
beta_bmiq_container <- new("ResultsContainer", filename = file.path(project_context$paths$results, "beta_matrix_bmiq.rds"), object = beta_bmiq, future = NULL)
bmiq_res <- list(
  beta_bmiq_container = beta_bmiq_container,
  m_bmiq_container = m_bmiq_container
)

source("R/outlier_analysis.R")
res_outlier <- intermediate_data_proxy(
  outlier_analysis,
  auto_clean = FALSE,
  project_context,
  pca = res_pca$pca_df_container@object,
  pca_filename = res_pca$pca_df_container@filename,
  targets = res_bio_gender_mismatch$targets_container@object,
  targets_filename = res_bio_gender_mismatch$targets_container@filename,
  beta_bmiq = bmiq_res$beta_bmiq_container@object,
  beta_bmiq_filename = bmiq_res$beta_bmiq_container@filename
)

saveRDS(res_outlier$pca_outliers_container@object, res_outlier$pca_outliers_container@filename)

source("R/outlier_remove_redo_BMIQ.R")
outlier_removed_bmiq_res <- outlier_remove_redo_BMIQ(
  context = project_context,
  pca = res_outlier$pca_outliers_container@object,
  beta_matrix = bmiq_res$beta_bmiq_container@object,
  pca_filename = res_outlier$pca_outliers_container@filename,
  beta_matrix_filename = bmiq_res$beta_bmiq_container@filename
)

beta_matrix <- outlier_removed_bmiq_res$beta_bmiq_container@object
targets <- res_bio_gender_mismatch$targets_container@object

  library(dplyr)
  library(stringr)

rownames(targets) <- targets$Basename %>% str_remove(paste0("ppmi/", "Project120_IDATS_n524final_toLONI_030718", "/"))
rownames(targets) <- targets$Basename %>% str_remove("ppmi")

rownames(targets) %>% head()

pd_samples <- rownames(targets[targets$Sample_Group == "PD", ])
hc_samples <- rownames(targets[targets$Sample_Group == "Control", ])

colnames_pd <- which(colnames(beta_matrix) %in% pd_samples)
colnames_hc <- which(colnames(beta_matrix) %in% hc_samples)

mean_beta_pd <- rowMeans(beta_matrix[, colnames_pd], na.rm = TRUE)
mean_beta_hc <- rowMeans(beta_matrix[, colnames_hc], na.rm = TRUE)
delta_beta <- mean_beta_pd - mean_beta_hc

beta_means <- data.frame(
  mean_beta_pd,
  mean_beta_hc,
  delta_beta
)

write.csv(beta_means, file.path(project_context$paths$results, "beta_means.csv"))

# --

rg_set <- readRDS(file.path(project_context$paths$raw_data, "rg_set.rds"))
targets <- readRDS(file.path(project_context$paths$processed, "targets.rds"))


source("R/cell_cnt_estimate.R")
res_cell_cnt_estimate <- intermediate_data_proxy(
  cell_cnt_estimate,
  auto_clean = FALSE,
  project_context,
  rg_set = rg_set,
  targets = targets
)

if (project_context$mode == results_mode()$memory_only) {
  print(paste("Saving targets after cell count estimation to", file.path(project_context$paths$processed, "targets_after_cell_count_estimation.rds")))
  saveRDS(res_cell_cnt_estimate$targets_container@object, file.path(project_context$paths$processed, "targets_after_cell_count_estimation.rds"))
}

# Plot Cell proportion per cohort and write results to files for each cell type
source("R/plot_cell_proportions.R")
plot_cell_proportions(
  context = project_context,
  targets = res_cell_cnt_estimate$targets_container@object,
  cell_types = getCellTypesForPlatform(project_context$platform)
)
getCellTypesForPlatform(project_context$platform)
project_context$platform

res_cell_cnt_estimate$targets_container@object %>% head()  
