rm(list = ls())
gc(full = TRUE)

library(arrow)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

project_name <- "PPMI"
root_dir <- "/Volumes/Elements/vastai/combat/"
target_df_loc <- "PPMI_harmonized_targets.rds"
m_values_loc <- "PPMI_combat_m_values.rds"
harmonize_targets <- as.logical("TRUE")
pattern_harm <- "PPMI_RAW/"

harmonize_meta <- function(m_values, meta, pattern) {
    rownames(meta) <- gsub(meta$Basename, pattern = pattern, replacement = "")
    return(meta[colnames(m_values),])
}
target_df <- readRDS(paste0(root_dir, target_df_loc))
m_values <- readRDS(paste0(root_dir, m_values_loc))

target_df %>% rownames() %>% head()
m_values %>% colnames() %>% head()

if (harmonize_targets) {
    target_df <- harmonize_meta(m_values, target_df, pattern_harm)
}

# PPMI had also SWEDD cohort samples, which we want to exclude from the final dataset.
target_df <- target_df %>% filter(Sample_Group != "SWEDD")
dim(target_df)
target_df %>% head()

dim(m_values)

common <- intersect(rownames(target_df), colnames(m_values))
length(common)

target_df <- target_df[common,]
dim(target_df)

target_df %>% rownames() %>% head()
m_values %>% colnames() %>% head()

dim(target_df)
dim(m_values)

diff <- setdiff(rownames(target_df), colnames(m_values))

all(rownames(target_df) == colnames(m_values))
if (all(rownames(target_df) == colnames(m_values))) {
    combined <- cbind(target_df, t(m_values))
    dim(combined)
    write_parquet(combined, write_statistics=FALSE, use_dictionary=FALSE, file.path(root_dir, paste0(project_name, "_data.parquet")))
} else {
    stop("Row names of meta do not match column names of m_values.")
}
