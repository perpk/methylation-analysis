principal_component_analysis <- function(context = NULL, m_values_filename = "m_values_bmiq.rds", targets_filename = "targets_remove_mismatch.rds", col_maps, keys, npc) {
  prog <- .create_progress_manager(4)

  targets_filepath <- file.path(context$paths$processed, targets_filename)
  prog$update(1, paste("Reading targets from ", targets_filepath))
  targets <- readRDS(targets_filepath)

  m_values_filepath <- file.path(context$paths$results, m_values_filename)
  prog$update(2, paste("Reading m-values matrix from", m_values_filepath))
  cat(paste("Reading m-values from", m_values_filepath))
  m_values <- readRDS(m_values_filepath)

  prog$update(3, "Running PCA")
  pca <- prcomp(t(m_values), center = TRUE, scale. = FALSE)

  pca_df <- data.frame(matrix(NA, nrow = ncol(m_values), ncol = npc))
  rownames(pca_df) <- colnames(m_values)
  colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
  })
  for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
  }

  for (k in keys) {
    if (!is.null(col_maps[[k]])) {
      pca_df[[k]] <- targets[[col_maps[[k]]]]
    }
  }

  pca_df_filepath <- file.path(context$paths$results, "pca_df.rds")
  prog$update(4, paste("Writing PCA results under", pca_df_filepath))
  saveRDS(pca_df, pca_df_filepath)

  prog$complete()
  rm(list = ls())
  gc(full = TRUE)
}
