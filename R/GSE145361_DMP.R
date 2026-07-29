rm(list = ls())
gc(full = TRUE)

parent_path <- "/root/"
sub_path <- "/data/"

m_values <- readRDS(
    file.path(
        paste0(
            parent_path, sub_path,
            "GSE145361_combat_m_values.rds"
        )
    )
)

targets <- readRDS(
    file.path(
        paste0(
            parent_path, sub_path,
            "GSE145361_harmonized_targets.rds"
        )
    )
)

library(sva)

targets$Sex <- targets$`gender:ch1`

mod <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK, data = targets)
mod0 <- model.matrix(~ Sex + CD8T + CD4T + Bcell + Mono + NK, data = targets)

n.sv <- num.sv(m_values, mod, method = "leek")
svobj <- sva(m_values, mod, mod0, n.sv = n.sv)

design_sva <- cbind(mod, svobj$sv)
colnames(design_sva) <- make.names(colnames(design_sva))

fit <- lmFit(m_values, design_sva)
fit2 <- eBayes(fit)

results <- topTable(fit2, coef = "Sample_GroupPD", number = Inf, adjust.method = "BH")

results %>% summary()
table(results$adj.P.Val < 0.05)


pca <- prcomp(t(m_values), scale. = TRUE)
pcs <- pca$x[, 1:5]
colnames(pcs) <- paste0("PC", 1:5)
design_pc <- model.matrix(
    ~ Sample_Group + pcs,
    data = targets
)

fit <- lmFit(m_values, design_pc)
fit2 <- eBayes(fit)

results <- topTable(fit2, coef = "Sample_GroupPD", number = Inf, adjust.method = "BH")
table(results$adj.P.Val < 0.05)

results %>% head()


library(ChAMP)

library(lumi)
beta_matrix <- m2beta(m_values)
myDMP <- champ.DMP(beta = beta_matrix, pheno = targets$Sample_Group, adjPVal = 0.05)
table(myDMP$PD)


probe_sd <- apply(m_values, 1, sd)
variable_probes <- names(probe_sd[probe_sd > quantile(probe_sd, 0.80)])

m_filtered <- m_values[variable_probes, ]

design <- model.matrix(
    ~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK,
    data = targets
)

fit <- lmFit(m_filtered, design)
fit2 <- eBayes(fit)

results <- topTable(fit2, coef = "Sample_GroupPD", number = Inf, adjust.method = "BH")
results %>% dim()
table(results$adj.P.Val < 0.05)
# ///


library(ggplot2)
library(tidyverse)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


is.factor(targets$Sample_Group)
targets$Sample_Group <- as.factor(targets$Sample_Group)

is.factor(targets$Sex)
targets$Sex <- as.factor(targets$`gender:ch1`)

is.factor(targets$ScanDate)
targets$ScanDate <- as.factor(targets$ScanDate)

design <- model.matrix(
    ~ 0 + Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran,
    data = targets
)

targets %>% head()

fit <- lmFit(m_values, design)
cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# rm(m_values)
# gc(full = TRUE)

## fit2 <- eBayes(fit)
results <- topTable(fit2, number = Inf, coef = "Parkinsons_vs_Control")

results$ProbeID <- rownames(results)
ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_sub <- ann_450k[rownames(results), c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
annotated_results <- merge(results, ann_sub, by.x = "ProbeID", by.y = "row.names")
annotated_results <- annotated_results[order(annotated_results$adj.P.Val, -annotated_results$logFC), ]

results[order(results$adj.P.Val, -results$logFC), ] %>% head()

head(annotated_results)

results %>% summary()

write.csv(annotated_results, file.path(paste0(parent_path, sub_path), "GSE145361_dmp_annotated_results.csv"))

volcano_plot <- ggplot(annotated_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(logFC) > 0.15))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed") +
    ggtitle("Volcano Plot of Differential Methylation") +
    theme_minimal()

ggsave(
    filename = file.path(paste0(parent_path, sub_path), "GSE145361_dmp_volcano_plot.png"),
    plot = volcano_plot,
    dpi = 300
)

beta_means <- read.csv(file.path(paste0(parent_path, sub_path), "GSE145361_beta_means.csv"), row.names = 1)

beta_means %>%
    rownames() %>%
    head()

annotated_results %>%
    rownames() %>%
    head()

rownames(annotated_results) <- annotated_results$ProbeID

annotated_with_beta <- merge(annotated_results, beta_means, by.x = "ProbeID", by.y = "row.names")
annotated_with_beta <- annotated_with_beta[order(annotated_with_beta$adj.P.Val, -abs(annotated_with_beta$logFC), -abs(annotated_with_beta$delta_beta)), ]

annotated_with_beta %>%
    head()

max(annotated_with_beta$delta_beta)
min(annotated_with_beta$delta_beta)

write.csv(annotated_with_beta, file.path(paste0(parent_path, sub_path), "GSE145361_dmp_annotated_with_beta_results.csv"))

annotated_with_beta_filtered_delta_beta_sign <- annotated_with_beta[(annotated_with_beta$adj.P.Val < 0.05 & abs(annotated_with_beta$delta_beta) > 0.01), ]

write.csv(annotated_with_beta_filtered_delta_beta_sign, file.path(paste0(parent_path, sub_path), "GSE145361_dmp_annotated_with_beta_filtered_delta_beta_sign_results.csv"))


annotated_with_beta_filtered_delta_beta_sign %>% head()

volcano_plot_delta_beta <- ggplot(annotated_with_beta, aes(x = delta_beta, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(delta_beta) > 0.01))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed") +
    ggtitle("Volcano Plot of Differential Methylation (Delta Beta)") +
    theme_minimal()

ggsave(
    filename = file.path(paste0(parent_path, sub_path), "GSE145361_dmp_volcano_plot_delta_beta.png"),
    plot = volcano_plot_delta_beta,
    dpi = 300
)
