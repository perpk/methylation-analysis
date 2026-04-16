rm(list = ls())
gc(full = TRUE)

rootDir <- "/Volumes/Elements/vastai/ppmi/ppmi_20260406"

pca_df <- readRDS(file.path(rootDir, "local/pca_df.rds"))

pca_scores <- pca_df[, 1:2]
center <- colMeans(pca_scores)
cov_matrix <- cov(pca_scores)
distances <- mahalanobis(pca_scores, center, cov_matrix)
outlier_threshold <- qchisq(0.975, df = 2)
outliers <- which(distances > outlier_threshold)

print("PCA outliers detected:")
print(rownames(pca_scores)[outliers])

targets <- readRDS(file.path(rootDir, "targets.rds"))

outlier_metadata <- targets[outliers, ]
print(table(outlier_metadata$Sample_Group))

beta_matrix <- readRDS(file.path(rootDir, "results/beta_matrix_bmiq.rds"))

outlier_betas <- beta_matrix[, outliers]
non_outlier_betas <- beta_matrix[, -outliers]
outlier_samples <- rownames(pca_scores)[outliers]
pca_df$Is_Outlier <- rownames(pca_df) %in% outlier_samples



library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Is_Outlier, shape = Sample_Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    ggtitle("PCA Showing Outlier Positions")
