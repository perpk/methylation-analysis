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

beta_matrix <- readRDS(file.path(project_context$paths$results, "beta_matrix.rds"))
m_matrix <- readRDS(file.path(project_context$paths$results, "m_matrix.rds"))

source("R/intermediate_data_proxy.R")
source("R/progress_mgr.R")
source("R/project_context.R")

source("R/apply_BMIQ.R")
bmiq_res <- intermediate_data_proxy(
  apply_BMIQ,
  project_context,
  beta_matrix = beta_matrix,
  beta_matrix_file = file.path(project_context$paths$results, "beta_matrix.rds"),
  plot = TRUE
)

m_bmiq <- readRDS(file.path(project_context$paths$results, "m_values_bmiq.rds"))
beta_bmiq <- readRDS(file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))

m_bmiq %>%
  is.na() %>%
  sum()

beta_bmiq %>%
  is.na() %>%
  sum()

library(ggplot2)
library(patchwork)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

plot_data <- data.frame(
  Beta = c(as.vector(beta_matrix), as.vector(beta_bmiq)),
  Method = rep(c("Before BMIQ", "After BMIQ"),
    each = length(beta_matrix)
  ),
  ProbeType = rep(.get_probe_types(rownames(beta_matrix), IlluminaHumanMethylation450kanno.ilmn12.hg19),
    times = ncol(beta_matrix) * 2
  )
)

p1 <- ggplot(plot_data, aes(x = Beta, color = Method)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Overall Beta Distribution",
    x = "Beta Value", y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

p2 <- ggplot(plot_data, aes(x = Beta, color = Method)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ProbeType, ncol = 1) +
  labs(
    title = "Beta Distribution by Probe Type",
    x = "Beta Value", y = "Density"
  ) +
  theme_minimal()

combined_plot <- p1 / p2 +
  plot_annotation(title = "BMIQ Normalization Effect")

bmiq_plot_file <- file.path(project_context$paths$processed, "BMIQ_normalization_comparison.png")
# ggsave fails for this dataset...
ggsave(bmiq_plot_file,
  plot = combined_plot,
  width = 10, height = 8, dpi = 300
)

source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = "GSE145361")

source("R/principal_component_analysis.R")
col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- var_mappings$gender_var
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

pca_df <- readRDS(pca_df_filepath)
head(m_bmiq)
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

is.null(outlier_removed_bmiq_res$beta_bmiq_container@object)
"package:wateRmelon" %in% search()

beta_no_outliers <- readRDS(file.path(project_context$paths$results, "beta_matrix_no_outliers.rds"))
head(beta_no_outliers)

source("R/results_container.R")
res_outlier <- intermediate_data_proxy(
  outlier_analysis,
  project_context,
  pca = pca_df,
  pca_filename = NULL,
  targets = targets_gender_mismatch,
  beta_bmiq = m_bmiq
)

source("R/outlier_remove_redo_BMIQ.R")
outlier_removed_bmiq_res <- outlier_remove_redo_BMIQ(
  context = project_context,
  pca = res_outlier$pca_outliers_container@object,
  beta_matrix = m_bmiq,
  pca_filename = res_outlier$pca_outliers_container@filename,
  beta_matrix_filename = file.path(project_context$paths$results, "m_values_bmiq.rds")
)
