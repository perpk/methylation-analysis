rm(list = ls())
gc(full = TRUE)

library(arrow)
library(dplyr)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "ppmi_20260721_075730"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "EPIC"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
targets <- readRDS(file.path(project_context$paths$processed, "ppmi_harmonized_targets.rds"))
m_values <- readRDS(file.path(project_context$paths$results, "ppmi_combat_m_values.rds"))

if (all(rownames(targets) %in% colnames(m_values))) {
    combined <- cbind(targets, t(m_values))
    dim(combined)
    write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path(project_context$paths$results, "ppmi_data.parquet"))
}
dim(targets)
dim(m_values)
ncol(m_values)
colnames(m_values)
rownames(targets)

sum(targets$Sample_Name == colnames(m_values))
targets$Basename

sum(colnames(m_values) %in% rownames(targets))
all(rownames(targets) %in% colnames(m_values))
