source("R/apply_BMIQ.R")

outlier_remove_redo_BMIQ <- function(context = NULL) {
    pca_filename <- "pca_df_with_outliers.rds"
    targets_filename <- "targets_remove_mismatch.rds"
    beta_matrix_filename <- "beta_matrix.rds"

    prog <- .create_progress_manager(4)
    prog$update(1, "Reading PCA results with outlier information")
    pca <- readRDS(file.path(context$paths$results, pca_filename))

    prog$update(2, "Reading beta-matrix")
    beta_matrix <- readRDS(file.path(context$paths$results, beta_matrix_filename))

    prog$update(3, "Removing outliers from beta-matrix")
    outlier_samples <- rownames(pca)[pca$Is_Outlier]
    beta_matrix_no_outliers <- beta_matrix[, !colnames(beta_matrix) %in% outlier_samples]

    saveRDS(beta_matrix_no_outliers, file.path(context$paths$results, "beta_matrix_no_outliers.rds"))

    prog$update(4, "Re-running BMIQ normalization without outliers")
    beta_matrix_no_outliers <- apply_BMIQ(context, "beta_matrix_no_outliers.rds")
}
