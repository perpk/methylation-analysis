library(arrow)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
root_dir <- args[2]
target_df_loc <- args[3]
m_values_loc <- args[4]
harmonize_targets <- as.logical(args[5])
pattern_harm <- args[6]

harmonize_meta <- function(m_values, meta, pattern) {
    rownames(meta) <- gsub(meta$Basename, pattern = pattern, replacement = "")
    meta <- meta[rownames(meta) %in% colnames(m_values), ]
    return(meta)
}
target_df <- readRDS(paste0(root_dir, target_df_loc))
m_values <- readRDS(paste0(root_dir, m_values_loc))

if (harmonize_targets) {
    target_df <- harmonize_meta(m_values, target_df, pattern_harm)
}
    
if (all(rownames(target_df) == colnames(m_values))) {
    combined <- cbind(target_df, t(m_values))
    write_parquet(combined, write_statistics=FALSE, use_dictionary=FALSE, file.path(root_dir, paste0(project_name, "_data.parquet")))
} else {
    stop("Row names of meta do not match column names of m_values.")
}
