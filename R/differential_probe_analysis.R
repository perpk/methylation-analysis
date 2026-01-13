library(ggplot2)
library(tidyverse)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

differential_probe_analysis <- function(
  context = NULL,
  beta_value_filename = "beta_matrix_bmiq.rds",
  m_value_filename = "m_values_bmiq.rds",
  targets_filename = "targets_s_mismatch_cells.rds",
  pca_cnt = 3
) {
  tryCatch({
    prog <- .create_progress_manager(5)

    beta_filepath <- file.path(context$paths$results, beta_value_filename)
    prog$update(1, paste("Reading beta-values from", beta_filepath))
    beta_matrix <- readRDS(beta_filepath)

    targets_filepath <- file.path(context$paths$qc, targets_filename)
    prog$update(2, paste("Reading targets from", targets_filepath))
    targets <- readRDS(targets_filepath)

    prog$update(3, "Calculating HC PC's and projecting Cases to establish unbiased covariates")
    targets_pcas <- .principal_component_covariates(context, beta_matrix, targets, pca_cnt)

    saveRDS(targets_pcas, file.path(context$paths$qc, "targets_pcas.rds"))

    design <- model.matrix(
      formula(paste("~ ", context$design_formula)),
      data = targets_pcas
    )

    heatmap_filepath <- file.path(context$paths$plots, "heatmap_covariate_correlation.png")
    cor_matrix <- cor(design[, -1])
    prog$update(4, paste("Creating covariate correlation heatmap and saving under", heatmap_filepath))
    png(heatmap_filepath)
    heatmap(cor_matrix,
      main = "Covariate Correlation Matrix",
      margins = c(10, 10)
    )
    dev.off()

    m_value_filepath <- file.path(context$paths$results, m_value_filename)
    prog$update(5, "Performing differential probe methylation analysis")
    m_values <- readRDS(m_value_filepath)

    fit <- lmFit(m_values, design)
    cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    results <- topTable(fit2, number = Inf, coef = "Parkinsons_vs_Control")

    write.csv(results, file.path(context$paths$results, "dmp_results.csv"))

    results$ProbeID <- rownames(results)
    ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann_sub <- ann_450k[rownames(results), c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
    annotated_results <- merge(results, ann_sub, by.x = "ProbeID", by.y = "row.names")
    annotated_results <- annotated_results[order(annotated_results$adj.P.Val, -annotated_results$logFC), ]

    write.csv(annotated_results, file.path(context$paths$results, "dmp_annotated_results.csv"))

    volcano_plot <- ggplot(annotated_results, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(alpha = 0.6, aes(color = (adj.P.Val < 0.05 & abs(logFC) > 0.2))) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
      ggtitle("Volcano Plot of Differential Methylation") +
      theme_minimal()

    ggsave(
      filename = file.path(context$paths$plots, "dmp_volcano_plot.png"),
      plot = volcano_plot,
      dpi = 300
    )

    prog$complete()
  }, finally = {
    rm(list = ls())
    gc(full = TRUE)
  })
}

.principal_component_covariates <- function(context, beta, targets, pca_cnt) {
  targets$Sample_Id <- gsub(".*/", "", targets$Basename)

  hc_samples <- targets %>%
    filter(Sample_Group == "Control") %>%
    pull(Sample_Id)

  beta_hc <- beta[, colnames(beta) %in% hc_samples]

  pca_hc <- prcomp(t(beta_hc))
  beta_pd <- beta[, !colnames(beta) %in% hc_samples]
  beta_all <- cbind(beta_hc, beta_pd)

  all_projected <- predict(pca_hc, newdata = t(beta_all))
  targets_pcas <- merge(
    targets,
    data.frame(
      Sample_Id = rownames(all_projected),
      all_projected[, 1:pca_cnt]
    ),
    by.x = "Sample_Id",
    by.y = "Sample_Id",
    all.x = TRUE
  )
  (targets_pcas)
}
