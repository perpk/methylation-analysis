principal_component_analysis <- function(context = NULL, col_maps, keys, npc) {
  targets_filepath <- file.path(context$paths$qc, "targets_clean.rds")
  cat(paste("\nReading targets from", targets_filepath, "\n"))
  targets <- readRDS(targets_filepath)

  beta_val_filepath <- file.path(context$paths$processed, "beta_matrix.rds")
  cat(paste("\nReading beta-values from", beta_val_filepath, "\n"))
  beta_matrix <- readRDS(beta_val_filepath)

  cat("\nRunning PCA\n")
  pca <- prcomp(t(beta_matrix))

  pca_df <- data.frame(matrix(NA, nrow = ncol(beta_matrix[1:100,]), ncol = npc))
  colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
  })
  for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[,index]
  }

  for (k in keys) {
    pca_df[[k]] <- targets[[col_maps[[k]]]]
  }

  pca_df_filepath <- file.path(context$paths$results, "pca_df_before_cleanup.rds")
  saveRDS(pca_df, pca_df_filepath)
  (pca_df)
}
