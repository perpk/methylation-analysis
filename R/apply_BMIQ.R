library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

apply_BMIQ <- function(context = NULL,
                       beta_matrix_file = "beta_matrix_noob.rds",
                       plot = TRUE) {
  platform <- NULL
  if (context$platform == "450K") {
    platform <- IlluminaHumanMethylation450kanno.ilmn12.hg19
  } else if (context$platform == "EPIC") {
    platform <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
  } else {
    stop("Unsupported platform. Please specify '450K' or 'EPIC'.")
  }

  prog <- .create_progress_manager(2)

  beta_matrix_filepath <- file.path(context$paths$results, beta_matrix_file)

  prog$update(1, paste("Reading beta-matrix file from", beta_matrix_filepath))
  beta_matrix <- readRDS(beta_matrix_filepath)

  prog$update(2, paste("Creating probe design vector for", context$platform))
  design.v <- .get_probe_design_vector(rownames(beta_matrix), platform)

  print(paste("Type I probes:", sum(design.v == 1)))
  print(paste("Type II probes:", sum(design.v == 2)))

  beta_bmiq <- matrix(NA, nrow = nrow(beta_matrix), ncol = ncol(beta_matrix))
  rownames(beta_bmiq) <- rownames(beta_matrix)
  colnames(beta_bmiq) <- colnames(beta_matrix)

  prog$complete()

  prog <- .create_progress_manager(ncol(beta_matrix))
  for (i in 1:ncol(beta_matrix)) {
    prog$update(i, paste("Processing sample", i, "of", ncol(beta_matrix)))

    # Extract beta values for this sample
    beta.v <- beta_matrix[, i]

    # Apply BMIQ to this sample
    tryCatch(
      {
        bmiq_result <- BMIQ(
          beta.v = beta.v,
          design.v = design.v,
          nfit = 10000,
          plots = FALSE,
          pri = TRUE
        )

        # Store the normalized values
        beta_bmiq[, i] <- bmiq_result$nbeta
      },
      error = function(e) {
        warning(paste("BMIQ failed for sample", colnames(beta_matrix)[i], ":", e$message))
        beta_bmiq[, i] <- beta.v # Keep original values if BMIQ fails
      }
    )
  }
  prog$complete()

  prog <- .create_progress_manager(2)

  prog$update(1, "Calculating M-Values based on BMIQ normalized betas")
  m_bmiq <- log2(beta_bmiq / (1 - beta_bmiq))
  m_bmiq[is.infinite(m_bmiq)] <- NA

  prog$update(2, "Saving normalized beta and m-values")
  beta_bmiq_filepath <- file.path(context$paths$results, "beta_matrix_bmiq.rds")
  saveRDS(beta_bmiq, beta_bmiq_filepath)

  m_bmiq_filepath <- file.path(context$paths$results, "m_values_bmiq.rds")
  saveRDS(m_bmiq, m_bmiq_filepath)
  if (plot) {
    print("Plot beta distribution before and after normalization")
    .plot_BMIQ_comparison(beta_matrix, beta_bmiq, platform)
  }
  prog$complete()
  rm(list = ls())
  gc(full = TRUE)
}

.plot_BMIQ_comparison <- function(beta_before, beta_after, platform) {
  library(ggplot2)
  library(patchwork)

  plot_data <- data.frame(
    Beta = c(as.vector(beta_before), as.vector(beta_after)),
    Method = rep(c("Before BMIQ", "After BMIQ"),
      each = length(beta_before)
    ),
    ProbeType = rep(.get_probe_types(rownames(beta_before), platform),
      times = ncol(beta_before) * 2
    )
  )

  p1 <- ggplot(plot_data, aes(x = Beta, color = Method)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Overall Beta Distribution",
      x = "Beta Value", y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "top")

  p2 <- ggplot(plot_data, aes(x = Beta, color = Method)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ProbeType, ncol = 1) +
    labs(
      title = "Beta Distribution by Probe Type",
      x = "Beta Value", y = "Density"
    ) +
    theme_minimal()

  combined_plot <- p1 / p2 +
    plot_annotation(title = "BMIQ Normalization Effect")

  print(combined_plot)

  bmiq_plot_file <- file.path(context$paths$processed, "BMIQ_normalization_comparison.png")
  ggsave(bmiq_plot_file,
    plot = combined_plot,
    width = 10, height = 8, dpi = 300
  )
}

.get_probe_types <- function(probe_names, platform) {
  anno <- getAnnotation(platform)

  probe_types <- anno[probe_names, "Type"]

  probe_types <- factor(probe_types,
    levels = c("I", "II"),
    labels = c("Type I", "Type II")
  )

  return(probe_types)
}

.get_probe_design_vector <- function(probe_names, platform) {
  # Returns: 1 for Type I probes, 2 for Type II probes
  anno <- getAnnotation(platform)

  # Match probes to annotation (handle missing probes)
  matched_probes <- intersect(probe_names, rownames(anno))

  if (length(matched_probes) < length(probe_names)) {
    warning(paste(
      length(probe_names) - length(matched_probes),
      "probes not found in annotation"
    ))
  }

  # Get probe types
  probe_types <- rep(NA, length(probe_names))
  names(probe_types) <- probe_names

  # Type I = 1, Type II = 2 (BMIQ convention)
  probe_types[matched_probes] <- ifelse(anno[matched_probes, "Type"] == "I", 1, 2)

  # Check for NAs
  if (any(is.na(probe_types))) {
    warning("Some probes could not be classified. Using alternative method...")

    # Alternative: Use probe name patterns
    # Type I: typically start with cg and have specific Infinium I design
    # Type II: typically start with ch or some cg with Infinium II design
    # This is less reliable but works as fallback
    na_probes <- probe_names[is.na(probe_types)]
    probe_types[na_probes] <- ifelse(grepl("^ch", na_probes), 2, 1)
  }

  return(probe_types)
}
