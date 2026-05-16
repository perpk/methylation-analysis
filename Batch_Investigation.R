extract_scandate_from_idat <- function(file_path=NULL, idat_ptn="^.+\\.idat(\\.gz)?$") {
    library(illuminaio)
    library(stringr)
    library(dplyr)
    results_df = data.frame(
        SentrixID = character(),
        ScanDate = character()
    )  
    all_idat_files <- list.files(file_path, pattern = idat_ptn, full.names = TRUE)
    print(paste("Found", length(all_idat_files), "IDAT files in the specified directory."))
    for (idat_file in all_idat_files) {
        idat_data <- readIDAT(idat_file)
        run_metadata <- idat_data$RunInfo
        scan_row <- which(run_metadata[, "BlockType"] == "Scan")[1]
        scan_date_string <- run_metadata[scan_row, "RunTime"]
        scan_date <- as.POSIXct(scan_date_string, format="%m/%d/%Y") %>% as.character() %>% str_extract("^(\\d{4}-\\d{2})")
        sentrix_id <- paste0(idat_data$Barcode, "_", idat_data$Unknowns$MostlyA)
        results_df <- rbind(results_df, data.frame(SentrixID = sentrix_id, ScanDate = scan_date))
    }
    return (results_df)
}
root_data_folder <- "/Volumes/Elements/methylation-analysis/"
root_dir <- "/Volumes/Elements/vastai/gse111629/GSE111629_20260515_083253"
m_values_loc <- "/results/m_values_bmiq.rds"
target_df_loc <- "/processed/targets_remove_mismatch.rds"

scan_dates_gse <- extract_scandate_from_idat(file_path=paste0(root_data_folder, "GSE111629_RAW/"))

m_values <- readRDS(paste0(root_dir, m_values_loc))

target_df <- readRDS(paste0(root_dir, target_df_loc))

target_df$Sample_Name <- target_df$Basename %>% str_extract("GSM\\d+_\\d{10}_R\\d{2}C\\d{2}")
head(target_df$Sample_Name)

target_df <- target_df[target_df$Sample_Name %in% colnames(m_values), ]
dim(target_df)

enriched_meta <- merge(
  x = target_df,
  y = scan_dates_gse,
  by.x = "Sentrix_ID",
  by.y = "SentrixID",
  all.x = FALSE
)

head(enriched_meta)
dim(enriched_meta)

enriched_meta$Sample_Name <- enriched_meta$Basename %>% str_extract("GSM\\d+_\\d{10}_R\\d{2}C\\d{2}")
head(enriched_meta$Sample_Name)
head(colnames(m_values))

head(enriched_meta$Sample_Name)

common_samples <- intersect(colnames(m_values), enriched_meta$Sample_Name)
length(common_samples)
head(common_samples)

enriched_meta_filtered <- enriched_meta[enriched_meta$Sample_Name %in% common_samples, ]
dim(enriched_meta_filtered)
