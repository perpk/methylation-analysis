### Differential Methylation Analysis

library(dplyr)
library(stringr)

project_to_load <- "GSE111629_20251226_102044"
data_folder <- "GSE111629_RAW"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

source("R/project_context.R")
context <- .load_methylation_project(
    base_dir = "/Volumes/Elements/methyl-pipe-out",
    project_id = project_to_load,
    cohorts = cohorts,
    platform = "450K"
)

targets <- readRDS(file.path(context$paths$qc, "targets_s_mismatch_cells.rds"))
pcs <- readRDS(file.path(context$paths$results, "pca_df.rds"))
View(pcs)
View(targets)

pcs$Sample_Name <- pcs %>%
    tibble::rownames_to_column("Sample_Name") %>%
    mutate(Sample_Name = str_extract(Sample_Name, "^[^_]+")) %>%
    pull(Sample_Name)

targets <- merge(targets, pcs[c("Sample_Name", "PC1", "PC2", "PC3", "PC4", "PC5")],
    by.x = "Sample_Name", by.y = "Sample_Name", all.x = TRUE
)

saveRDS(targets, file.path(context$paths$qc, "targets_s_mismatch_cells_pcs.rds"))

design <- model.matrix(
    ~ 0 + Sample_Group + CD4T + Bcell + PC1 + PC2,
    data = targets
)

# targets$`age:ch1` <- as.numeric(targets$`age:ch1`)
# targets$`gender:ch1` <- as.character(targets$`gender:ch1`)
# names(targets)[names(targets) == "gender:ch1"] <- "Sex"
# names(targets)[names(targets) == "age:ch1"] <- "Age"

# # Have a look at correlation between covariates, just for kicks
# design_1 <- model.matrix(
#   ~ 0 + Sample_Group + PC1 + PC2 + PC3 + PC4 + PC5 + Age + Sex,
#   data = targets
# )

cor_matrix <- cor(design[, -1])
heatmap(cor_matrix,
    main = "Covariate Correlation Matrix",
    margins = c(10, 10)
)

library(limma)

m_values <- readRDS(file.path(context$paths$results, "m_values_bmiq.rds"))

fit <- lmFit(m_values, design)
cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number = Inf, coef = "Parkinsons_vs_Control")
print("Sum of statistically significant differentially methylated probes")
sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.3)

hist(results[results$adj.P.Val < 0.05, ]$logFC, breaks = 50, main = "LogFC Distribution of Significant DMPs", xlab = "LogFC")

results$ProbeID <- rownames(results)
ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_sub <- ann_450k[rownames(results), c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
annotated_results <- merge(results, ann_sub, by.x = "ProbeID", by.y = "row.names")
annotated_results <- annotated_results[order(annotated_results$adj.P.Val, -annotated_results$logFC), ]

library(ggplot2)
volcano_plot <- ggplot(annotated_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(logFC) > 0.3))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed") +
    ggtitle("Volcano Plot of Differential Methylation") +
    theme_minimal()
print(volcano_plot)

### PCA on HC only

beta <- readRDS(file.path(context$paths$results, "beta_matrix_bmiq.rds"))

colnames(beta)

targets$Sample_Id <- targets$Basename %>%
    str_remove(paste0(data_folder, "/"))

hc_samples <- targets %>%
    filter(Sample_Group == "Control") %>%
    pull(Sample_Id)

beta_hc <- beta[, colnames(beta) %in% hc_samples]
dim(beta_hc)

pca_hc <- prcomp(t(beta_hc))

beta_pd <- beta[, !colnames(beta) %in% hc_samples]

beta_all <- cbind(beta_hc, beta_pd)
all_projected <- predict(pca_hc, newdata = t(beta_all))
original_controls <- pca_hc$x
head(original_controls)

projected_controls <- all_projected[rownames(pca_hc$x), 1:10]
cor(original_controls[, 1], projected_controls[, 1])
summary(pca_hc)$importance[, 1:5]

View(projected_controls)

keep <- which(!colnames(targets) %in% c("PC1", "PC2", "PC3", "PC4", "PC5"))
targets <- targets[, keep]

targets <- merge(
    targets,
    data.frame(
        Sample_Id = rownames(all_projected),
        all_projected[, 1:5]
    ),
    by.x = "Sample_Id",
    by.y = "Sample_Id",
    all.x = TRUE
)
View(targets)

design <- model.matrix(
    ~ 0 + Sample_Group + CD4T + Bcell + PC1 + PC2 + PC3,
    data = targets
)

cor_matrix <- cor(design[, -1])
heatmap(cor_matrix,
    main = "Covariate Correlation Matrix",
    margins = c(10, 10)
)

library(limma)

m_values <- readRDS(file.path(context$paths$results, "m_values_bmiq.rds"))

fit <- lmFit(m_values, design)
cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number = Inf, coef = "Parkinsons_vs_Control")
print("Sum of statistically significant differentially methylated probes")
sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.2)

hist(results[results$adj.P.Val < 0.05, ]$logFC, breaks = 50, main = "LogFC Distribution of Significant DMPs", xlab = "LogFC")

results$ProbeID <- rownames(results)
ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_sub <- ann_450k[rownames(results), c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
annotated_results <- merge(results, ann_sub, by.x = "ProbeID", by.y = "row.names")
annotated_results <- annotated_results[order(annotated_results$adj.P.Val, -annotated_results$logFC), ]

library(ggplot2)
volcano_plot <- ggplot(annotated_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(logFC) > 0.2))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
    ggtitle("Volcano Plot of Differential Methylation") +
    theme_minimal()
print(volcano_plot)


main <- function() {}
