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

targets <- readRDS(file.path(project_context$paths$processed, "targets.rds"))
targets_gender_mismatch <- readRDS(file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"))

targets_after_qc <- readRDS(file.path(project_context$paths$processed, "targets_after_qc.rds"))

dim(targets)
dim(targets_gender_mismatch)
dim(targets_after_qc)



m_bmiq <- readRDS(file.path(project_context$paths$results, "m_values_bmiq.rds"))
beta_bmiq <- readRDS(file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))

source("R/principal_component_analysis.R")
col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- var_mapping$gender_var
keys <- c("Sample_Group", "Gender")

pca <- prcomp(t(m_bmiq), center = TRUE, scale. = FALSE)

pca_df <- data.frame(matrix(NA, nrow = ncol(m_bmiq), ncol = 5))

rownames(pca_df) <- colnames(m_bmiq)
colnames(pca_df) <- sapply(1:5, function(pca_df) {
  (paste0("PC", pca_df))
})

for (index in 1:5) {
  pca_df[[paste0("PC", index)]] <- pca$x[, index]
}

source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = "GSE145361")

col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- var_mappings$gender_var
col_map[["Age"]] <- var_mappings$age_var

keys <- c("Sample_Group", "Gender", "Age")
for (k in keys) {
  if (!is.null(col_map[[k]])) {
    pca_df[[k]] <- targets_gender_mismatch[[col_map[[k]]]]
  }
}


pca_df_filepath <- file.path(project_context$paths$results, "pca_df.rds")
saveRDS(pca_df, pca_df_filepath)

source("R/plot_PCA.R")
pca_vars <- c("Sample_Group" = "By Diagnosis", "Gender" = "By Biological Gender", "Age" = "By Age")

# GSE145361 doesn't have Age included in the source meta information
# convert_f <- function(df) {
#   (transform(df, Age = as.numeric(Age)))
# }
plot_PCA(
  context = project_context,
  pca_df = pca_df,
  pca_results_rds_filename = pca_df_filepath,
  pca_vars = pca_vars,
  convert_fun = NULL,
  continuously_scaled = c("Age"),
  "pca_plot_",
  pairplot_color_by = "Age"
)

source("R/outlier_analysis.R")
res_outlier <- intermediate_data_proxy(
  outlier_analysis,
  project_context,
  pca = pca_df,
  pca_filename = pca_df_filepath,
  targets = targets_gender_mismatch,
  targets_filename = file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"),
  beta_bmiq = m_bmiq,
  beta_bmiq_filename = file.path(project_context$paths$results, "m_values_bmiq.rds")
)

source("R/outlier_remove_redo_BMIQ.R")
outlier_removed_bmiq_res <- outlier_remove_redo_BMIQ(
  context = project_context,
  pca = res_outlier$pca_outliers_container@object,
  beta_matrix = beta_bmiq,
  pca_filename = res_outlier$pca_outliers_container@filename,
  beta_matrix_filename = NULL
)
