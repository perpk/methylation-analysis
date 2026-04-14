#####################################################
# Outlier removal and BMIQ normalization redo
# This script performs outlier detection and removal based on PCA, followed by BMIQ normalization on
# the remaining samples. It saves the cleaned beta matrix and the corresponding M-values.
# Files needed:
#   pca_df.rds
#   targets.rds
#   beta_matrix_bmiq.rds
# Outputs:
#   pca_outliers.png
#   outlier_log.csv
#   beta_matrix_no_outliers.rds
#   beta_matrix_bmiq_no_outliers.rds
#   m_matrix_bmiq_no_outliers.rds
#####################################################

rootDir = ""
platform = "450K"

pca_df <- readRDS(file.path(rootDir, "pca_df.rds"))

pca_scores <- pca_df[, 1:2]
center <- colMeans(pca_scores)
cov_matrix <- cov(pca_scores)
distances <- mahalanobis(pca_scores, center, cov_matrix)
outlier_threshold <- qchisq(0.975, df = 2)
outliers <- which(distances > outlier_threshold)

targets <- readRDS(file.path(rootDir, "targets.rds"))
outlier_metadata <- targets[outliers, ]

beta_matrix <- readRDS(file.path(rootDir, "beta_matrix_bmiq.rds"))

outlier_betas <- beta_matrix[, outliers]
non_outlier_betas <- beta_matrix[, -outliers]
outlier_samples <- rownames(pca_scores)[outliers]
pca_df$Is_Outlier <- rownames(pca_df) %in% outlier_samples

library(ggplot2)
plot_outliers <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Is_Outlier, shape = Sample_Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    ggtitle("PCA Showing Outlier Positions")

filename <- file.path(rootDir, "pca_outliers.png")

ggsave(
    filename = filename,
    plot = plot_outliers,
    width = 8,
    height = 6,
    dpi = 300
)
outlier_cases_logfile <- file.path(rootDir, "outlier_log.csv")
outlier_cases <- pca_df[pca_df$Is_Outlier, ]
write.csv(outlier_cases, file = outlier_cases_logfile)

beta_matrix <- beta_matrix[, -outliers]
saveRDS(beta_matrix, file.path(rootDir, "beta_matrix_no_outliers.rds"))

source("R/apply_BMIQ.R")

print("Beta-Mixture quantile normalization (BMIQ)")

beta_matrix <- readRDS(file.path(rootDir, "beta_matrix_no_outliers.rds"))
design.v <- .get_probe_design_vector(rownames(beta_matrix), platform)

print(paste("Type I probes:", sum(design.v == 1)))
print(paste("Type II probes:", sum(design.v == 2)))

beta_bmiq <- matrix(NA, nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))

rownames(beta_bmiq) <- rownames(beta_matrix)
colnames(beta_bmiq) <- colnames(beta_matrix)

for (i in 1:ncol(beta_matrix)) {
    # Extract beta values for this sample
    beta.v <- beta_matrix[, i]

    # Apply BMIQ to this sample
    tryCatch(
        {
            bmiq_result <- BMIQ(
                beta.v = beta.v,
                design.v = design.v,
                nfit = 10000,
                plots = FALSE,
                pri = TRUE
            )

            # Store the normalized values
            beta_bmiq[, i] <- bmiq_result$nbeta
        },
        error = function(e) {
            warning(paste("BMIQ failed for sample", colnames(beta_matrix)[i], ":", e$message))
            beta_bmiq[, i] <- beta.v # Keep original values if BMIQ fails
        }
    )
}

m_bmiq <- log2(beta_bmiq / (1 - beta_bmiq))
m_bmiq[is.infinite(m_bmiq)] <- NA

saveRDS(beta_bmiq, file.path(rootDir, "beta_matrix_bmiq_no_outliers.rds"))
saveRDS(m_bmiq, file.path(rootDir, "m_matrix_bmiq_no_outliers.rds"))


write.csv(beta_bmiq, file.path(rootDir, "beta_matrix_bmiq_no_outliers.csv"))
write.csv(m_bmiq, file.path(rootDir, "m_matrix_bmiq_no_outliers.csv"))
