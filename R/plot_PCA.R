plot_PCA <- function(context = NULL, pca_results_rds_filename = NULL, pca_vars = NULL,
                     convert_fun = NULL, continuously_scaled = c(), pca_output_name = "pca_plot_", npc=5, create_pairplot=TRUE, pairplot_color_by=NULL) {

  # Load required libraries
  library(ggplot2)
  library(gridExtra)

  prog <- .create_progress_manager(1)

  # 1. Read PCA data
  m_filepath <- file.path(context$paths$results, pca_results_rds_filename)
  print(paste("Reading PCA from file", m_filepath))
  pca_df <- readRDS(m_filepath)


  # 2. Apply conversion function if provided (FIXED)
  if (!is.null(convert_fun)) {
    pca_df <- convert_fun(pca_df)  # Function must RETURN the modified dataframe
  }

  # Debug: Check pca_df structure
  print("pca_df columns:")
  print(colnames(pca_df))
  print("pca_df dimensions:")
  print(dim(pca_df))

  # 3. Check which pca_vars actually exist in pca_df
  available_vars <- names(pca_vars)[names(pca_vars) %in% colnames(pca_df)]
  if (length(available_vars) < length(pca_vars)) {
    missing_vars <- setdiff(names(pca_vars), colnames(pca_df))
    warning(paste("Some variables not found in pca_df:", paste(missing_vars, collapse = ", ")))
  }

  plots <- list()

  prog$update(1, "Creating PCA plots")

  # 4. Create plots only for available variables
  plot_index <- 1
  for (i in seq_along(available_vars)) {
    var_name <- available_vars[i]           # Column name in pca_df ("Sample_Group", "Gender", "Age")
    plot_title <- pca_vars[var_name]        # Plot title ("By Diagnosis", etc.)

    print(paste("Creating plot:", var_name, "->", plot_title))

    if (var_name %in% continuously_scaled) {
      # Continuous variable - use viridis color scale
      plots[[plot_index]] <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[var_name]])) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_viridis_c(name = var_name) +
        labs(title = plot_title,
             x = "Principal Component 1",
             y = "Principal Component 2") +
        theme_minimal()
    } else {
      # Categorical variable
      plots[[plot_index]] <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                                color = as.factor(.data[[var_name]]))) +
        geom_point(size = 3, alpha = 0.8) +
        labs(title = plot_title,
             x = "Principal Component 1",
             y = "Principal Component 2",
             color = var_name) +
        theme_minimal() +
        theme(legend.position = "right")
    }

    plot_index <- plot_index + 1
  }

  # 5. Save individual plots
  for (i in seq_along(plots)) {
    var_name <- available_vars[i]
    filename <- file.path(context$paths$results,
                          paste0(pca_output_name, var_name, ".png"))

    ggsave(filename = filename,
           plot = plots[[i]],
           width = 8,
           height = 6,
           dpi = 300)
    print(paste("Saved:", filename))
  }

  # 6. Save combined plot if we have any plots
  if (length(plots) > 0) {
    combined <- grid.arrange(grobs = plots, ncol = min(2, length(plots)))
    combined_filename <- file.path(context$paths$results,
                                   "pca_all_variables_grid.png")

    ggsave(filename = combined_filename,
           plot = combined,
           width = min(2, length(plots)) * 8,
           height = ceiling(length(plots)/2) * 6,
           dpi = 300)
    print(paste("Saved combined plot:", combined_filename))
  } else {
    print("No plots were created - check if variables exist in pca_df")
  }

  if (create_pairplot && npc > 1) {
    print("Creating PCA pairplot...")

    # Choose which variable to color by (first available by default)
    color_var <- if (!is.null(pairplot_color_by)) {
      pairplot_color_by
    } else if (length(available_vars) > 0) {
      available_vars[1]
    } else {
      NULL
    }

    # Create pairplot (using GGally version)
    pairplot <- .add_pca_pairplot(
      pca_df = pca_df,
      n_pcs = npc,
      color_by = color_var,
      context = context,
      output_name = paste0(pca_output_name, "pairplot")
    )
  }

  prog$complete()
  return(invisible(list(individual_plots = plots, pairplot = if(exists("pairplot")) pairplot else NULL)))
}

.add_pca_pairplot <- function(pca_df, n_pcs = 4, color_by = NULL,
                                    context = NULL, output_name = "pca_pairplot") {

  library(ggplot2)
  library(patchwork)

  # Select first n PCs
  pc_columns <- grep("^PC\\d+$", colnames(pca_df), value = TRUE)
  pc_columns <- pc_columns[1:min(n_pcs, length(pc_columns))]

  if (length(pc_columns) < 2) {
    warning("Need at least 2 PC columns for pairplot")
    return(NULL)
  }

  # Generate all combinations
  combos <- expand.grid(x = pc_columns, y = pc_columns)

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
        p <- ggplot(pca_df, aes(x = .data[[x_var]], y = .data[[y_var]],
                                color = as.factor(.data[[color_by]]))) +
          geom_point(size = 0.5, alpha = 0.5) +
          scale_color_viridis_d(name = color_by)
      } else {
        p <- ggplot(pca_df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
          geom_point(size = 0.5, alpha = 0.5, color = "steelblue")
      }
      p <- p +
        labs(x = x_var, y = y_var) +
        theme_minimal() +
        theme(axis.text = element_text(size = 6),
              legend.position = "none")
    }

    plots[[plot_idx]] <- p
    plot_idx <- plot_idx + 1
  }

  # Arrange in grid
  pairplot <- wrap_plots(plots, ncol = length(pc_columns)) +
    plot_annotation(title = paste("PCA Pairplot - First", length(pc_columns), "PCs"))

  # Save if context provided
  if (!is.null(context)) {
    filename <- file.path(context$paths$results, paste0(output_name, ".png"))
    ggsave(filename, pairplot,
           width = 3 * length(pc_columns),
           height = 3 * length(pc_columns),
           dpi = 300)
    print(paste("Saved pairplot to:", filename))
  }

  return(pairplot)
}
