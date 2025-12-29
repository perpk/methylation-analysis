cell_types <- c("CD8T", "CD4T", "NK", "Bcell", "Mono")
cell_results <- list()

targets <- readRDS(file.path(project_context$paths$qc, "targets_s_mismatch_cells.rds"))

library(ggpubr)

for(cell_type in cell_types) {
  test_result <- wilcox.test(targets[[cell_type]] ~ targets$Sample_Group)

  cell_results[[cell_type]] <- data.frame(
    CellType = cell_type,
    p_value = test_result$p.value,
    W_statistic = test_result$statistic
  )

  p <- ggplot(targets, aes(x = Sample_Group, y = .data[[cell_type]], fill = Sample_Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    labs(title = paste(cell_type, "Proportions"),
         y = "Proportion", x = "") +
    theme_minimal()
  print(p)
}
