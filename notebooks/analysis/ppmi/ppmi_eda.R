root_data_dir <- "/Volumes/Elements/vastai/ppmi/ppmi_20260406"

targets <- readRDS(file.path(root_data_dir, "/results/targets_remove_mismatch.rds"))

beta_matrix <- readRDS(file.path(root_data_dir, "results/beta_matrix_bmiq.rds"))

dim(targets)
dim(beta_matrix)

head(targets)

## ____
npc <- 6

pca <- prcomp(t(beta_matrix))
pca_df <- data.frame(matrix(NA, nrow = ncol(beta_matrix), ncol = npc))

rownames(pca_df) <- colnames(beta_matrix)
colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
})
for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}

colnames(targets)

pca_df[["Sample_Group"]] <- targets[["Sample_Group"]]
pca_df[["SEX"]] <- targets[["SEX"]]
pca_df[["Age_Group"]] <- targets[["Age_Group"]]
pca_df[["Array"]] <- targets[["Array"]]
pca_df[["Slide"]] <- targets[["Slide"]]

saveRDS(pca_df, file.path(root_data_dir, "local/pca_df.rds"))

library(ggplot2)
library(gridExtra)

pca_df <- pca_df[pca_df[["Sample_Group"]] %in% c("PD", "Control"), ]

ggplot(pca_df, aes(x = PC1, y = PC2, color = SEX)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
        title = "PCA of PPMI Methylation Data",
        x = "Principal Component 1",
        y = "Principal Component 2"
    ) +
    theme_minimal()
