rootDir <- "/Volumes/Elements/vastai/gse111629"

targets <- readRDS(file.path(rootDir, "qc/targets_s_mismatch_cells.rds"))

library(sva)
library(dplyr)
library(stringr)

targets$Sex <- targets$"gender:ch1"

targets <- targets %>%
    mutate(
        Age_Group =
            case_when(
                targets$"age:ch1" >= 30 & targets$"age:ch1" < 50 ~ "30-49",
                targets$"age:ch1" >= 50 & targets$"age:ch1" < 60 ~ "50-59",
                targets$"age:ch1" >= 60 & targets$"age:ch1" < 70 ~ "60-69",
                targets$"age:ch1" >= 70 & targets$"age:ch1" < 80 ~ "70-79",
                targets$"age:ch1" >= 80 ~ "80+",
            )
    )

targets <- targets %>%
    mutate(
        Sentrix_ID = str_extract(source_name_ch1, "[0-9]{10}")
    )

targets <- targets %>%
    mutate(
        Sentrix_Pos = str_extract(source_name_ch1, "R\\d{2}C\\d{2}")
    )

View(targets)

beta_vals <- readRDS(file.path(rootDir, "results/beta_matrix_bmiq.rds"))
targets <- targets[match(colnames(beta_vals), rownames(targets)), ]
m_vals <- readRDS(file.path(rootDir, "results/m_values_bmiq.rds"))
dim(targets)

library(ggplot2)

pca <- prcomp(t(m_vals))
pca_df <- data.frame(matrix(NA, nrow = ncol(beta_vals), ncol = 6))
rownames(pca_df) <- colnames(beta_vals)
colnames(pca_df) <- sapply(1:6, function(pca_df) {
    (paste0("PC", pca_df))
})
for (index in 1:6) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}
pca_df[["Sample_Group"]] <- targets[["Sample_Group"]]
pca_df[["Sex"]] <- targets[["Sex"]]
pca_df[["Age_Group"]] <- targets[["Age_Group"]]
pca_df[["Sentrix_ID"]] <- targets[["Sentrix_ID"]]
pca_df[["Sentrix_Pos"]] <- targets[["Sentrix_Pos"]]

saveRDS(pca_df, file.path(rootDir, "local/pca_df.rds"))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Age_Group)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
        title = "PCA of PPMI Methylation Data",
        x = "Principal Component 1",
        y = "Principal Component 2"
    ) +
    theme_minimal()

library(patchwork)
colnames(pca_df)

unique(pca_df$Sentrix_ID)

pca_df[pca_df$Sentrix_ID == "9721367028", ]

n_pcs <- 6
pc_columns <- grep("^PC\\d+$", colnames(pca_df), value = TRUE)
pc_columns <- pc_columns[1:min(n_pcs, length(pc_columns))]
combos <- expand.grid(x = pc_columns, y = pc_columns)
View(targets)
drivers <- c("Sentrix_ID", "Sample_Group", "Sex", "Age_Group", "Sentrix_Pos", "CD4T", "CD8T", "Bcell")

for (driver in drivers) {
    color_by <- driver

    plots <- list()
    plot_idx <- 1
    for (i in 1:nrow(combos)) {
        x_var <- as.character(combos$x[i])
        y_var <- as.character(combos$y[i])

        if (x_var == y_var) {
            # Diagonal: Density plot
            p <- ggplot(pca_df, aes(x = .data[[x_var]])) +
                geom_density(fill = "steelblue", alpha = 0.5) +
                labs(x = x_var, y = "Density") +
                theme_minimal() +
                theme(axis.text = element_text(size = 6))
        } else {
            # Off-diagonal: Scatter plot
            if (!is.null(color_by) && color_by %in% colnames(pca_df)) {
                p <- ggplot(pca_df, aes(
                    x = .data[[x_var]], y = .data[[y_var]],
                    color = as.factor(.data[[color_by]]), shape = Sample_Group
                )) +
                    geom_point(size = 0.5, alpha = 0.5) +
                    scale_color_viridis_d(name = color_by)
            } else {
                p <- ggplot(pca_df, aes(x = .data[[x_var]], y = .data[[y_var]], color = Sample_Group)) +
                    geom_point(size = 0.5, alpha = 0.5)
            }
            legendPos <- "right"
            if (driver == "Sentrix_ID" || driver == "Sentrix_Pos") {
                legendPos <- "none"
            }
            p <- p +
                labs(x = x_var, y = y_var) +
                theme_minimal() +
                theme(
                    axis.text = element_text(size = 6),
                    legend.position = legendPos
                )
        }
        plots[[plot_idx]] <- p
        plot_idx <- plot_idx + 1
    }

    pairplot <- wrap_plots(plots, ncol = length(pc_columns)) +
        plot_annotation(title = paste("PCA Pairplot - First", length(pc_columns), "PCs"))

    ggsave(file.path(rootDir, paste0("local/pca_pairplot_", driver, ".png")), pairplot, width = 3 * length(pc_columns), height = 3 * length(pc_columns))
}

# m_vals <- readRDS(file.path(rootDir, "results/m_values_bmiq.rds"))
# mod <- model.matrix(~ Sample_Group + Age_Group + Sex + CD4T + Bcell, data = targets)
# m_combat <- ComBat(dat = m_vals, batch = )
