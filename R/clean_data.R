clean_data <- function(samples_to_remove = NULL, data_to_clean = NULL) {
  names_list <- names(data_to_clean)
  results <- lapply(names_list, function(name) {
    print(paste("Cleaning Dataset", name))
    data <- readRDS(data_to_clean[[name]]$filepath)
    if (!is.null(data_to_clean[[name]]$clean_rows) &&
      data_to_clean[[name]]$clean_rows == TRUE) {
      by_col <- data_to_clean[[name]]$by_col
      if (!is.null(by_col)) {
        s <- data[[by_col]]
      } else {
        s <- rownames(data)
      }
      indices <- which(s %in% samples_to_remove)
      clean_data <- data[-indices, ]
    } else {
      s <- colnames(data)
      indices <- which(s %in% samples_to_remove)
      clean_data <- data[,-indices]
    }
    rm(data)
    saveRDS(clean_data, data_to_clean[[name]]$clean_filename)
    rm(clean_data)
    gc(full=T)
  })
}
