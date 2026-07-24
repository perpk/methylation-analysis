rm(list = ls())
gc(full = TRUE)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE111629_20260722_083339"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "GSE111629_harmonized_targets.rds"))
m_values <- readRDS(file.path(project_context$paths$processed, "GSE111629_harmonized_m_values.rds"))

library(sva)

run_combat <- function(m_values, meta_df, batch_colname, mod_colnames) {
    library(sva)
    batch <- meta_df[[batch_colname]]
    mod_matrix <- model.matrix(~ as.factor(Sample_Group), data = meta_df)
    combat_m_values <- ComBat(dat = m_values, batch = batch, mod = mod_matrix)
    return(combat_m_values)
}

combat_m_values <- run_combat(
    m_values = m_values,
    meta_df = targets,
    batch_colname = "ScanDate",
    mod_colnames = "Sample_Group"
)
saveRDS(combat_m_values, file.path(project_context$paths$results, "GSE111629_combat_m_values.rds"))

# ==== re-do PCA

pca <- prcomp(t(combat_m_values))

pca_df <- data.frame(matrix(NA, nrow = ncol(combat_m_values), ncol = 10))
colnames(pca_df) <- sapply(1:10, function(i) paste0("PC", i))
rownames(pca_df) <- colnames(combat_m_values)
for (index in 1:10) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}
split_names <- strsplit(rownames(targets), "_")
targets$Slide <- sapply(split_names, function(x) x[2])
targets$Array <- sapply(split_names, function(x) x[3])
targets %>% head()

targets_for_pca <- targets[, c("Sample_Group", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "Sex", "Age_Group", "ScanDate", "Sample_Name", "Slide", "Array")]

pca_df_enriched_meta <- merge(
    x = pca_df,
    y = targets_for_pca,
    by.x = "row.names",
    by.y = "row.names",
    all = FALSE
)

saveRDS(
    pca_df_enriched_meta,
    file.path(project_context$paths$processed, "GSE111629_pca_df_enriched_meta_after_combat.rds")
)

library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)

pca_pairplot <- function(pca_df_with_dates, color_by = NULL, size_factor = 3, n_pcs = 6) {
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

    ggsave(file.path(project_context$paths$plots, paste0("GSE111629_pca_pairplot_after_combat", color_by, ".png")),
        pplot,
        width = size_factor * length(pc_columns),
        height = size_factor * length(pc_columns),
        dpi = 300,
        limitsize = FALSE
    )
}

pca_df_enriched_meta$Slide_num <- as.numeric(pca_df_enriched_meta$Slide)

pca_df_enriched_meta$Slide_num

pca_pairplot(pca_df_enriched_meta, color_by = "ScanDate", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Slide_num", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Array", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Sex", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Age_Group", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Sample_Group", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "CD8T", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "CD4T", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "NK", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Bcell", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Mono", n_pcs = 5)
pca_pairplot(pca_df_enriched_meta, color_by = "Gran", n_pcs = 5)
