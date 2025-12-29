principal_component_analysis <- function(context = NULL, beta_matrix_filename = "beta_matrix_bmiq.rds", targets_filename = "targets_remove_mismatch.rds", col_maps, keys, npc) {

  prog <- .create_progress_manager(4)

  targets_filepath <- file.path(context$paths$processed, targets_filename)
  prog$update(1, paste("Reading targets from ", targets_filepath))
  targets <- readRDS(targets_filepath)

  beta_val_filepath <- file.path(context$paths$results, beta_matrix_filename)
  prog$update(2, paste("Reading beta matrix from", beta_val_filepath))
  cat(paste("Reading beta-values from", beta_val_filepath))
  beta_matrix <- readRDS(beta_val_filepath)

  prog$update(3, "Running PCA")
  pca <- prcomp(t(beta_matrix))

  pca_df <- data.frame(matrix(NA, nrow = ncol(beta_matrix), ncol = npc))
  rownames(pca_df) <- colnames(beta_matrix)
  colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
  })
  for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[,index]
  }

  for (k in keys) {
    pca_df[[k]] <- targets[[col_maps[[k]]]]
  }

  pca_df_filepath <- file.path(context$paths$results, "pca_df.rds")
  prog$update(4, paste("Writing PCA results under", pca_df_filepath))
  saveRDS(pca_df, pca_df_filepath)

  prog$complete()
  rm(list = ls())
  gc(full = T)
}
