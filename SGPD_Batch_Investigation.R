rm(list = ls())
gc(full = TRUE)

library(tidyverse)
library(minfi)

data_dir <- "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/"

m_values_bmiq_no_outliers <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_m_values.rds"))
targets <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_targets.rds"))

qc_dmr <- function(m_values, design) {
    library(limma)

    # until Bcell λ=1.007
    # with Bcell and scandate at the end λ=1.055

    fit <- lmFit(m_values, design)
    fit <- eBayes(fit)
    stats <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    p_values <- stats$P.Value
    chisq <- qchisq(1 - p_values, 1)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

    print(paste("Genomic Inflation Factor (Lambda):", round(lambda, 3)))

    expected_p <- ppoints(length(p_values))
    observed_p <- sort(p_values)

    png(file.path(data_dir, "plots/GSE145361_QQ_plot.png"), width = 800, height = 800)

    plot(-log10(expected_p), -log10(observed_p),
        main = paste("QQ Plot of PD Methylation\nLambda =", round(lambda, 3)),
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        pch = 16, cex = 0.5, col = "darkblue"
    )

    abline(0, 1, col = "red", lwd = 2)
    dev.off()
}

targets$Slide <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets$Array <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])

table(targets$ScanDate, targets$Sample_Group)
table(targets$ScanDate, targets$Array)
table(targets$ScanDate, targets$Slide)

targets$Sentrix_ID %>% head()

targets$Sex <- targets$`gender:ch1`

targets %>% head()

# design <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array, data = targets)
design <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Array, data = targets)
qc_dmr(m_values_bmiq_no_outliers, design)

# 0.732, 0.738, 0.792, 0.972, 0.987, 1.029,
# Up to Mono = 1.029
# With NK = 1.028
# With Gran = 0.987
# Deteriorates with ScanDate: λ = 0.643
# With Array = 1.007

# R03C02, etc are array positions. The other part is the slide
pca <- prcomp(t(m_values_bmiq_no_outliers))

pca_df <- data.frame(matrix(NA, nrow = ncol(m_values_bmiq_no_outliers), ncol = 10))
colnames(pca_df) <- sapply(1:10, function(i) paste0("PC", i))
rownames(pca_df) <- colnames(m_values_bmiq_no_outliers)
for (index in 1:10) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}

targets_for_pca <- targets[, c("Sample_Group", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "Sex", "ScanDate", "Sample_Name", "Slide", "Array")]

pca_df_enriched_meta <- merge(
    x = pca_df,
    y = targets_for_pca,
    by.x = "row.names",
    by.y = "row.names",
    all = FALSE
)

saveRDS(
    pca_df_enriched_meta,
    file.path(data_dir, "processed/GSE145361_pca_df_enriched_meta.rds")
)

library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)

pca_pairplot <- function(pca_df_with_dates, color_by = NULL, size_factor = 3, n_pcs = 6, prefix = "GSE145361_pca_pairplot_") {
    pc_columns <- grep("^PC\\d+$", colnames(pca_df_with_dates), value = TRUE)
    pc_columns <- pc_columns[1:min(n_pcs, length(pc_columns))]
    combos <- expand.grid(x = pc_columns, y = pc_columns)
    plots <- list()
    plot_idx <- 1

    for (i in 1:nrow(combos)) {
        x_var <- as.character(combos$x[i])
        y_var <- as.character(combos$y[i])

        if (x_var == y_var) {
            # Diagonal: Density plot
            p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]])) +
                geom_density(fill = "steelblue", alpha = 0.5) +
                labs(x = x_var, y = "Density") +
                theme_minimal() +
                theme(axis.text = element_text(size = 6))
        } else {
            # Check if color_by is valid and exists
            if (!is.null(color_by) && color_by %in% colnames(pca_df_with_dates)) {
                # BRANCH A: Continuous / Numeric
                if (is.numeric(pca_df_with_dates[[color_by]])) {
                    p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[color_by]])) +
                        geom_point(size = size_factor, alpha = 0.8) +
                        scale_color_viridis_c(name = color_by) +
                        theme_minimal()

                    # BRANCH B: Categorical / Discrete
                } else {
                    p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]], y = .data[[y_var]], color = as.factor(.data[[color_by]]))) +
                        geom_point(size = size_factor, alpha = 0.5) +
                        scale_color_viridis_d(name = color_by) +
                        theme_minimal()
                }

                # BRANCH C: No color_by provided
            } else {
                p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]], y = .data[[y_var]])) +
                    geom_point(size = size_factor, alpha = 0.5, color = "steelblue") +
                    theme_minimal()
            }

            # Apply common labels and themes
            p <- p +
                labs(x = x_var, y = y_var) +
                theme(
                    axis.text = element_text(size = 6),
                    legend.position = "right"
                )
        }
        plots[[plot_idx]] <- p
        plot_idx <- plot_idx + 1
    }

    pplot <- wrap_plots(plots, ncol = length(pc_columns)) +
        plot_annotation(title = paste("PCA Pairplot - First", length(pc_columns), "PCs"))

    ggsave(file.path(data_dir, "plots", paste0(prefix, color_by, ".png")),
        pplot,
        width = size_factor * length(pc_columns),
        height = size_factor * length(pc_columns),
        dpi = 300,
        limitsize = FALSE
    )
}

pca_pairplot(pca_df_enriched_meta, color_by = "ScanDate", n_pcs = 5)
# pca_pairplot(pca_df_enriched_meta, color_by = "Slide_num", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Array", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Sex", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Sample_Group", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "CD8T", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "CD4T", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "NK", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Bcell", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Mono", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Gran", n_pcs = 5)


# Remove batch effects via limma
dd <- model.matrix(~Sample_Group, data = targets)
# design <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Array, data = targets)
sgpd_covariates <- as.matrix(targets[, c("CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran")])
sex_numeric <- as.numeric(as.factor(targets$Sex)) - 1
sgpd_covariates <- cbind(sgpd_covariates, Sex = sex_numeric)
m_matrix_clean <- removeBatchEffect(
    x = m_values_bmiq_no_outliers,
    batch = targets$Array,
    covariates = sgpd_covariates,
    design = dd
)

saveRDS(
    m_matrix_clean,
    file.path(data_dir, "processed/GSE145361_harmonized_m_values_cleaned.rds")
)

library(arrow)

all(rownames(targets) %in% colnames(m_matrix_clean))

all(colnames(m_matrix_clean) %in% rownames(targets))

combined <- cbind(targets, t(m_matrix_clean))
dim(combined)
write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path(data_dir, "processed/GSE145361_data_test.parquet"))

library(lumi)
beta_matrix_clean <- m2beta(m_matrix_clean)

saveRDS(
    beta_matrix_clean,
    file.path(data_dir, "processed/GSE145361_harmonized_beta_values_cleaned.rds")
)


m_matrix_clean <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_m_values_cleaned.rds"))

pca <- prcomp(t(m_matrix_clean))

pca_df_clean <- data.frame(matrix(NA, nrow = ncol(m_matrix_clean), ncol = 10))
colnames(pca_df_clean) <- sapply(1:10, function(i) paste0("PC", i))
rownames(pca_df_clean) <- colnames(m_matrix_clean)
for (index in 1:10) {
    pca_df_clean[[paste0("PC", index)]] <- pca$x[, index]
}

targets_for_pca <- targets[, c("Sample_Group", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "Sex", "ScanDate", "Sample_Name", "Slide", "Array")]

pca_df_enriched_meta_clean <- merge(
    x = pca_df_clean,
    y = targets_for_pca,
    by.x = "row.names",
    by.y = "row.names",
    all = FALSE
)

library(ggplot2)
library(stringr)
library(patchwork)

pca_pairplot(pca_df_enriched_meta_clean, color_by = "ScanDate", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
# pca_pairplot(pca_df_enriched_meta, color_by = "Slide_num", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Array", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Sex", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Sample_Group", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "CD8T", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "CD4T", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "NK", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Bcell", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Mono", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")
pca_pairplot(pca_df_enriched_meta_clean, color_by = "Gran", n_pcs = 5, prefix = "GSE145361_pca_pairplot_cleaned_")

saveRDS(
    m_matrix_clean,
    file.path(data_dir, "processed/GSE145361_harmonized_m_values_cleaned.rds")
)


library(arrow)

all(rownames(targets) %in% colnames(m_matrix_clean))

all(colnames(m_matrix_clean) %in% rownames(targets))

combined <- cbind(targets, t(m_matrix_clean))
dim(combined)
write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path(data_dir, "processed/GSE145361_data.parquet"))

library(lumi)
beta_matrix_clean <- m2beta(m_matrix_clean)

saveRDS(
    beta_matrix_clean,
    file.path(data_dir, "processed/GSE145361_harmonized_beta_values_cleaned.rds")
)
