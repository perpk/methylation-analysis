plot_cell_proportions <- function(context = NULL,
                                  targets_file = "targets_s_mismatch_cells.rd",
                                  cell_types = c("CD8T", "CD4T", "NK", "Bcell", "Mono")) {
  targets <- readRDS(file.path(context$paths$qc, "targets_s_mismatch_cells.rds"))

  library(ggpubr)
  cell_results <- list()
  for (cell_type in cell_types) {
    test_result <- wilcox.test(targets[[cell_type]] ~ targets$Sample_Group)

    cell_results[[cell_type]] <- data.frame(
      CellType = cell_type,
      p_value = test_result$p.value,
      pd_median = median(targets[[cell_type]][targets$Sample_Group == "PD"], na.rm = TRUE),
      hc_median = median(targets[[cell_type]][targets$Sample_Group == "Control"], na.rm = TRUE),
      pd_mean = mean(targets[[cell_type]][targets$Sample_Group == "PD"], na.rm = TRUE),
      hc_mean = mean(targets[[cell_type]][targets$Sample_Group == "Control"], na.rm = TRUE),
      effect_size = median(targets[[cell_type]][targets$Sample_Group == "PD"], na.rm = TRUE) -
        median(targets[[cell_type]][targets$Sample_Group == "Control"], na.rm = TRUE)
    )
  }

  results <- do.call(rbind, cell_results)
  results$fdr <- p.adjust(results$p_value, method = "fdr")
  results$bonf <- p.adjust(results$p_value, method = "bonferroni")

  results <- results[order(results$fdr), ]
  write.csv(results, file.path(context$paths$results, "cell_statistic_sign.csv"))

  plots <- list()

  for(ct in cell_types) {
    current <- results[results$CellType == ct, ]

    p <- .create_celltype_plot(
      targets = targets,
      cell_type = ct,
      fdr_p = current$fdr,
      bonf_p = current$bonf,
      raw_p = current$p_value
    )

    ggsave(
      filename = file.path(context$paths$plots, paste0(ct, "_proportions.png")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )

    plots[[ct]] <- p
  }

  library(patchwork)
  multi_plot <- wrap_plots(plots, ncol = 3) +
    plot_annotation(title = "Blood Cell Type Proportions in Parkinson's Disease",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

  ggsave(
    filename = file.path(context$paths$plots, "all_celltypes_multipanel.png"),
    plot = multi_plot,
    width = 16,
    height = 12,
    dpi = 300
  )

  rm(list = ls())
  gc(full = T)
}

.create_celltype_plot <- function(targets, cell_type, fdr_p, bonf_p, raw_p) {
  display_p <- fdr_p
  display_label <- paste0("p[FDR] = ", fdr_p)

  if (bonf_p < 0.05) {
    display_label <- paste0("p[FDR] = ", fdr_p,
                            "\np[Bonf] = ", bonf_p)
  }

  p <- ggplot(targets, aes(x = Sample_Group, y = .data[[cell_type]],
                           fill = Sample_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +

    annotate("text",
             x = 1.5,
             y = max(targets[[cell_type]], na.rm = TRUE) * 1.05,
             label = display_label,
             parse = TRUE,
             size = 4,
             fontface = "bold") +

    geom_signif(comparisons = list(c("Control", "PD")),
                annotations = ifelse(fdr_p < 0.001, "***",
                                     ifelse(fdr_p < 0.01, "**",
                                            ifelse(fdr_p < 0.05, "*", "ns"))),
                y_position = max(targets[[cell_type]], na.rm = TRUE) * 1.1,
                tip_length = 0.01) +

    labs(
      title = paste(cell_type, "Proportions in Parkinson's Disease"),
      subtitle = paste("FDR =", fdr_p,
                       ifelse(bonf_p < 0.05, paste0("(Bonferroni = ", bonf_p, ")"), "")),
      y = "Estimated Proportion",
      x = "",
      caption = ifelse(fdr_p < 0.05,
                       "Significant after FDR correction",
                       "Not significant after multiple testing correction")
    ) +
    scale_fill_manual(values = c("Control" = "#2E86AB", "PD" = "#A23B72")) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray40"))

  return(p)
}
