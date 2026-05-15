idat_dir = "./ppmi/Project120_IDATS_n524final_toLONI_030718/"
all_idat_files <- list.files(idat_dir, pattern = "\\.idat$", full.names = TRUE)

idat_dir_gse = "./GSE111629_RAW/"
all_gse_idat_files <- list.files(idat_dir_gse, pattern = "^GSM.*\\.gz$", full.names = TRUE)
head(all_gse_idat_files)

results_df = data.frame(
    SentrixID = character(),
    ScanDate = character()
)

library(illuminaio)

idat <- readIDAT(all_idat_files[1])

idat_gse <- readIDAT(gzfile(all_gse_idat_files[1]))
all_gse_idat_files[1]


str(idat)



library(illuminaio)
for (idat_file in all_idat_files) {
    idat_data <- readIDAT(idat_file)
    run_metadata <- idat_data$RunInfo
    scan_row <- which(run_metadata[, "BlockType"] == "Scan")[1]
    scan_date_string <- run_metadata[scan_row, "RunTime"]
    scan_date <- as.POSIXct(scan_date_string, format="%m/%d/%Y")
    
    sentrix_id <- paste0(idat_data$Barcode, "_", idat_data$Unknowns$MostlyA)

    results_df <- rbind(results_df, data.frame(SentrixID = sentrix_id, ScanDate = scan_date))
}

print(results_df)

write.csv(results_df, "./ppmi/ppmi_scan_dates.csv", row.names = FALSE)

head(ppmi_meth_120_meta)

ppmi_meth_120_dates <- merge(
  x = ppmi_meth_120_meta,
  y = results_df,
  by.x = "Basename",
  by.y = "SentrixID",
  all.x = FALSE
)

dim(ppmi_meth_120_meta)
dim(ppmi_meth_120_dates)
View(ppmi_meth_120_dates)
View(ppmi_meth_120_meta)
