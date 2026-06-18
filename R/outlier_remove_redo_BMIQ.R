source("R/apply_BMIQ.R")

outlier_remove_redo_BMIQ <- function(
    context = NULL,
    pca = NULL,
    beta_matrix = NULL,
    pca_filename = NULL,
    beta_matrix_filename = NULL
) {

    prog <- .create_progress_manager(4)
    prog$update(1, "Reading PCA results with outlier information")
    if (is.null(pca)) {
        pca <- readRDS(pca_filename)
    }

    prog$update(2, "Reading beta-matrix")
    if (is.null(beta_matrix)) {
        beta_matrix <- readRDS(beta_matrix_filename)
    }

    prog$update(3, "Removing outliers from beta-matrix")
    outlier_samples <- rownames(pca)[pca$Is_Outlier]
    beta_matrix_no_outliers <- beta_matrix[, !colnames(beta_matrix) %in% outlier_samples]

    beta_matrix_no_outliers_filepath <- file.path(context$paths$results, "beta_matrix_no_outliers.rds")
    if (context$mode == results_mode()$disk_only || context$mode == results_mode()$disk_and_memory) {
        prog$update(3, paste("Saving beta-matrix without outliers to", beta_matrix_no_outliers_filepath))
    }
    saveRDS(beta_matrix_no_outliers, beta_matrix_no_outliers_filepath)

    prog$update(4, "Re-running BMIQ normalization without outliers")
    return(apply_BMIQ(context, beta_matrix = beta_matrix_no_outliers, beta_matrix_filename = beta_matrix_no_outliers_filepath))
}
