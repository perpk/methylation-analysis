rm(list = ls())
gc(full = TRUE)

library(arrow)
library(dplyr)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE111629_20260722_083339"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
targets <- readRDS(file.path(project_context$paths$processed, "GSE111629_harmonized_targets.rds"))
m_values <- readRDS(file.path(project_context$paths$results, "GSE111629_combat_m_values.rds"))

if (all(rownames(targets) == colnames(m_values))) {
    combined <- cbind(targets, t(m_values))
    dim(combined)
    write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path(project_context$paths$results, "GSE111629_data.parquet"))
}

dim(m_values)
ncol(m_values)
