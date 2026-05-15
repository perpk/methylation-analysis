args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
root_dir <- args[2]
target_df_loc <- args[3]
m_values_loc <- args[4]
idat_folder_loc <- args[5]
extract_sentrix_id_from_basename <- as.logical(args[6])

print(paste("Project Name:", project_name))
print(paste("Root Directory:", root_dir))
print(paste("Target DataFrame Location:", target_df_loc))
print(paste("M-values Location:", m_values_loc))
print(paste("IDAT Folder Location:", idat_folder_loc))
print(paste("Extract Sentrix ID from Basename:", extract_sentrix_id_from_basename))
print(paste("Loading target DataFrame from:", paste0(root_dir, target_df_loc)))
target_df <- readRDS(paste0(root_dir, target_df_loc))
print(paste("Target DataFrame loaded with dimensions:", dim(target_df)[1], "rows and", dim(target_df)[2], "columns"))

if (extract_sentrix_id_from_basename) {
    target_df$Sentrix_ID <- sapply(target_df$Basename, function(x) {
        basename(x) %>% str_extract("\\d{10}_R\\d{2}C\\d{2}")
    })
}

print(paste("Loading M-values from:", paste0(root_dir, m_values_loc)))
m_values <- readRDS(paste0(root_dir, m_values_loc))
print(paste("M-values loaded with dimensions:", dim(m_values)[1], "rows and", dim(m_values)[2], "columns"))

source("R/extract_scandate_from_idat.R")
print(paste("Extracting scan dates from IDAT files in:", paste0(root_dir, idat_folder_loc)))
scan_dates <- extract_scandate_from_idat(
    file_path=paste0(root_dir, idat_folder_loc)
)
print(paste("Scan dates extracted with dimensions:", dim(scan_dates)[1], "rows and", dim(scan_dates)[2], "columns"))

source("R/enrich_meta_with_batch.R")
print("Enriching target DataFrame with scan date information...")
enriched_target_df <- enricht_meta_with_batch(
    meta_df = target_df, 
    batch_df = scan_dates, 
    x_colname = "Sentrix_ID", 
    y_colname = "SentrixID"
)
print(paste("Enriched target DataFrame dimensions:", dim(enriched_target_df)[1], "rows and", dim(enriched_target_df)[2], "columns"))

source("R/run_combat.R")
print("Running ComBat for batch correction...")
combat_m_values <- run_combat(
    m_values = m_values,
    meta_df = enriched_target_df,
    batch_colname = "ScanDate",
    mod_colnames = c("Sample_Group")
)
print("ComBat batch correction completed.")

print(paste("Saving ComBat-corrected M-values to:", paste0(root_dir, paste0(project_name, "_combat_m_values.rds"))))
saveRDS(
    combat_m_values, 
    paste0(root_dir, 
        paste0(project_name, "_combat_m_values.rds")
    )
)
print("All done!")


