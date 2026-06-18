outlier_analysis <- function(
  context = NULL, 
  pca = NULL,
  targets = NULL,
  beta_bmiq = NULL,
  pca_filename = NULL, 
  targets_filename = NULL, 
  beta_bmiq_filename = NULL,
  sample_metadata = NULL
) {
  if (is.null(pca)) {
    pca <- readRDS(file.path(context$paths$results, pca_filename))
  }
  if (is.null(targets)) {
    targets <- readRDS(file.path(context$paths$processed, targets_filename))
  }

  prog <- .create_progress_manager(4)

  prog$update(1, "Detecting Outliers")
  pca_scores <- pca[, 1:2]
  center <- colMeans(pca_scores)
  cov_matrix <- cov(pca_scores)
  distances <- mahalanobis(pca_scores, center, cov_matrix)
  outlier_threshold <- qchisq(0.975, df = 2)
  outliers <- which(distances > outlier_threshold)
  print("PCA outliers detected:")
  print(rownames(pca_scores)[outliers])

  outlier_metadata <- targets[outliers, ]
  print("Outlier sample metadata:")
  print(outlier_metadata[, sample_metadata])
  print(table(outlier_metadata$Sample_Group))

  if (is.null(beta_bmiq)) {
    beta_bmiq <- readRDS(file.path(context$paths$results, beta_bmiq_filename))
  }

  prog$update(2, "Fetching Outlier beta ranges")
  outlier_betas <- beta_bmiq[, outliers]
  print("Outlier beta value ranges:")
  print(apply(outlier_betas, 2, function(x) c(mean = mean(x), sd = sd(x))))

  non_outlier_betas <- beta_bmiq[, -outliers]
  print("Non-outlier beta value ranges:")
  print(c(mean = mean(non_outlier_betas), sd = sd(non_outlier_betas)))

  outlier_samples <- rownames(pca_scores)[outliers]
  pca$Is_Outlier <- rownames(pca) %in% outlier_samples

  prog$update(3, "Creating Outlier Scatterplot")
  plot_outliers <- ggplot(pca, aes(x = PC1, y = PC2, color = Is_Outlier, shape = Sample_Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    ggtitle("PCA Showing Outlier Positions")

  filename <- file.path(context$paths$plots, "pca_outliers.png")

  ggsave(
    filename = filename,
    plot = plot_outliers,
    width = 8,
    height = 6,
    dpi = 300
  )
  print(paste("Saved:", filename))

  prog$update(4, "Creating Outlier Log")
  outlier_cases_logfile <- file.path(context$paths$results, "outlier_log.csv")
  outlier_cases <- pca[pca$Is_Outlier, ]
  write.csv(outlier_cases, file = outlier_cases_logfile)
  prog$complete()

  pca_outliers_container <- new("ResultsContainer", filename = file.path(context$paths$results, "pca_df_with_outliers.rds"), object = pca, future = NULL) 
}
