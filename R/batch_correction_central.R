args <- commandArgs(trailingOnly = TRUE)

projects <- args[1]
rootDir <- args[2]

print(paste("Project Name:", projects))
print(paste("Root Directory:", rootDir))

library(dplyr)
library(stringr)

project_list <- str_split(projects, ",")[[1]]

m_values_template <- "_harmonized_m_values.rds"
targets_template <- "_harmonized_targets.rds"

source("R/run_combat.R")

for (project_name in project_list) {
    print(paste("Processing project:", project_name))
    
    m_values_loc <- paste0(project_name, m_values_template)
    targets_loc <- paste0(project_name, targets_template)
    
    print(paste("Loading target DataFrame from:", paste0(rootDir, targets_loc)))
    target_df <- readRDS(paste0(rootDir, targets_loc))
    print(paste("Target DataFrame loaded with dimensions:", dim(target_df)[1], "rows and", dim(target_df)[2], "columns"))
    m_values <- readRDS(paste0(rootDir, m_values_loc))
    print(paste("M-values loaded with dimensions:", dim(m_values)[1], "rows and", dim(m_values)[2], "columns"))
    print(paste("Running ComBat batch correction for project:", project_name))
    combat_m_values <- run_combat(
        m_values = m_values,
        meta_df = target_df,
        batch_colname = "ScanDate",
        mod_colnames = "Sample_Group")
    print(paste("ComBat batch correction completed for project:", project_name))
    print(paste("Saving batch-corrected m-values for project:", project_name))
    saveRDS(combat_m_values, paste0(rootDir, project_name, "_combat_m_values.rds"))
}

print("All done!")


