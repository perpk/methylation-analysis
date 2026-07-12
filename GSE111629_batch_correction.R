rm(list = ls())
gc(full = TRUE)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE111629_20260709_222156"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "GSE111629_harmonized_targets.rds"))
m_values <- readRDS(file.path(project_context$paths$processed, "GSE111629_harmonized_m_values.rds"))

library(sva)

run_combat <- function(m_values, meta_df, batch_colname, mod_colnames) {
    library(sva)
    batch <- meta_df[[batch_colname]]
    mod_matrix <- model.matrix(~ as.factor(Sample_Group), data = meta_df)
    combat_m_values <- ComBat(dat = m_values, batch = batch, mod = mod_matrix)
    return(combat_m_values)
}

combat_m_values <- run_combat(
    m_values = m_values,
    meta_df = targets,
    batch_colname = "ScanDate",
    mod_colnames = "Sample_Group"
)
saveRDS(combat_m_values, file.path(project_context$paths$results, "GSE111629_combat_m_values.rds"))
