extract_scandate_from_idat <- function(file_path=NULL, idat_ptn="^.+\\.idat(\\.gz)?$") {
    library(illuminaio)
    library(stringr)
    library(dplyr)
    library(lubridate)
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
        scan_date <- NA
        # if (is.vector(scan_date_string)) {
        #     print(paste("Scan date string for file", idat_file, "is a vector. Attempting to parse the first element."))
        #     tmp_d_str <- scan_date_string[1]
        #     scan_date <- dmy_hms(tmp_d_str) %>% as.character() %>% str_extract("^(\\d{4}-\\d{2})")
        #     print("successfully parsed scan date string using lubridate's dmy_hms function.")
        # } else {
            print(paste("Scan date string for file", idat_file, "is not a vector. Attempting to parse directly."))
            scan_date <- as.POSIXct(scan_date_string, format="%m/%d/%Y") %>% as.character() %>% str_extract("^(\\d{4}-\\d{2})")
            print("successfully parsed scan date string using as.POSIXct function.")
        # }
        if (is.na(scan_date)) {
            warning(paste("Could not extract scan date from file:", idat_file))
        }
        sentrix_id <- paste0(idat_data$Barcode, "_", idat_data$Unknowns$MostlyA)
        results_df <- rbind(results_df, data.frame(SentrixID = sentrix_id, ScanDate = scan_date))
    }
    return (results_df)
}

# root_data_folder <- "/Volumes/Elements/methylation-analysis/"
# scan_dates_ppmi <- extract_scandate_from_idat(file_path=paste0(root_data_folder, "ppmi/Project120_IDATS_n524final_toLONI_030718/"))
# scan_dates_gse <- extract_scandate_from_idat(file_path=paste0(root_data_folder, "GSE111629_RAW/"))
# scan_dates_gse2 <- extract_scandate_from_idat(file_path=paste0(root_data_folder, "GSE145361_RAW/"))
# head(scan_dates_gse)

