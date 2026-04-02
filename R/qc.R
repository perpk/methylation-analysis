library(minfi)
library(minfiData)

qc <- function(context = NULL,
               targets = NULL,
               rg_set_filename = "rg_set.rds",
               methyl_set_filename = "methyl_set.rds",
               qc_threshold = 10.5,
               bisulfite_sd_multiplier = 2,
               bisulfite_absolute_min = NULL) {
  prog <- .create_progress_manager(8) # Increased to 8 steps

  methyl_set_filepath <- file.path(context$paths$raw_data, methyl_set_filename)
  prog$update(1, paste("Reading methyl set from ", methyl_set_filepath))
  methyl_set <- readRDS(methyl_set_filepath)
  rg_set_filepath <- file.path(context$paths$raw_data, rg_set_filename)
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, "Performing QC")

  qc <- getQC(methyl_set)
  saveRDS(qc, file.path(context$paths$qc, "qc_metrics.rds"))

  prog$update(3, "Creating QC plot")
  qc_plot_path <- file.path(context$paths$plots, "qc_intensity_plot.png")

  png(qc_plot_path, width = 800, height = 600)
  plotQC(qc)
  dev.off()

  prog$update(4, "Checking bisulfite conversion controls")
  bisulfite_results <- .check_bisulfite_conversion(rg_set, bisulfite_sd_multiplier, bisulfite_absolute_min, context, prog)
  bisulfite_failed_samples <- bisulfite_results$failed_samples
  bisulfite_thresholds <- bisulfite_results$thresholds

  prog$update(5, "Detecting failed samples")
  intensity_threshold <- .determine_qc_threshold(qc, qc_threshold, context, prog)

  # Get intensity-based failures
  intensity_bad_samples <- .identify_failed_samples(qc, intensity_threshold, prog)

  # Combine with bisulfite failures
  all_bad_indices <- unique(c(intensity_bad_samples$indices, bisulfite_failed_samples$indices))
  all_bad_names <- unique(c(intensity_bad_samples$names, bisulfite_failed_samples$names))

  # Create combined bad_samples structure
  bad_samples <- list(
    indices = all_bad_indices,
    names = all_bad_names,
    intensity_failures = intensity_bad_samples$indices,
    bisulfite_failures = bisulfite_failed_samples$indices
  )

  prog$update(6, "Filtering samples")

  if (length(bad_samples$indices) > 0) {
    rg_set_clean <- rg_set[, -bad_samples$indices, drop = FALSE]
    targets_clean <- targets[-bad_samples$indices, , drop = FALSE]
    methyl_set_clean <- methyl_set[, -bad_samples$indices, drop = FALSE]

    cat(paste(
      "Removing", length(bad_samples$indices),
      "failed sample(s):\n"
    ))
    cat(paste(bad_samples$names, collapse = ", "))
    cat("\n")
  } else {
    rg_set_clean <- rg_set
    targets_clean <- targets
    methyl_set_clean <- methyl_set
    cat("No samples failed QC\n")
  }

  rm(list = c("targets", "methyl_set"))
  gc(full = T)

  prog$update(7, "Saving results")

  saveRDS(
    rg_set_clean,
    file.path(context$paths$qc, "rg_set_clean.rds")
  )
  saveRDS(
    targets_clean,
    file.path(context$paths$qc, "targets_clean.rds")
  )
  saveRDS(
    methyl_set_clean,
    file.path(context$paths$qc, "methyl_set_clean.rds")
  )

  # Save bisulfite thresholds for documentation
  saveRDS(
    bisulfite_thresholds,
    file.path(context$paths$logs, "bisulfite_thresholds.rds")
  )

  qc_log <- .create_qc_log_with_bisulfite(rg_set, bad_samples, qc, intensity_threshold, bisulfite_thresholds, bisulfite_failed_samples)
  write.csv(qc_log,
    file.path(context$paths$logs, "sample_removal_log.csv"),
    row.names = FALSE
  )

  summary <- .create_qc_summary_with_bisulfite(rg_set, rg_set_clean, qc_log, intensity_threshold, bisulfite_thresholds)
  saveRDS(
    summary,
    file.path(context$paths$qc, "qc_summary.rds")
  )


  rm(list = setdiff(ls(), c("prog", "qc", "qc_log", "qc_summary", "rg_set_filepath", "qc_plot_path", "targets_clean", "rg_set_clean", "context")))
  gc(full = T)

  prog$update(8, "Creating Density Plot...")
  density_plot_path <- file.path(context$paths$plots, "density_plot.png")

  png(density_plot_path, width = 800, height = 600)
  densityPlot(rg_set_clean, sampGroups = targets_clean$Sample_Group, main = "Beta Value Distribution")
  dev.off()

  prog$complete()
}

.check_bisulfite_conversion <- function(rg_set, sd_multiplier, absolute_min, context, prog) {
  # Get addresses for both conversion control types
  ctrls_i <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION I")
  ctrls_ii <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION II")

  # Calculate average intensities per sample
  # BC I uses Green channel for the 'Converted' signal
  bc1_scores <- colMeans(getGreen(rg_set)[ctrls_i, , drop = FALSE], na.rm = TRUE)
  # BC II uses Red channel for the 'Converted' signal
  bc2_scores <- colMeans(getRed(rg_set)[ctrls_ii, , drop = FALSE], na.rm = TRUE)

  # Calculate thresholds
  bc1_mean <- mean(bc1_scores, na.rm = TRUE)
  bc1_sd <- sd(bc1_scores, na.rm = TRUE)
  bc2_mean <- mean(bc2_scores, na.rm = TRUE)
  bc2_sd <- sd(bc2_scores, na.rm = TRUE)

  bc1_threshold_stat <- bc1_mean - (sd_multiplier * bc1_sd)
  bc2_threshold_stat <- bc2_mean - (sd_multiplier * bc2_sd)

  # Apply absolute minimum floor if provided
  if (!is.null(absolute_min)) {
    bc1_threshold <- max(bc1_threshold_stat, absolute_min)
    bc2_threshold <- max(bc2_threshold_stat, absolute_min)
  } else {
    bc1_threshold <- bc1_threshold_stat
    bc2_threshold <- bc2_threshold_stat
  }

  # Identify failures
  fail_i_indices <- which(bc1_scores < bc1_threshold)
  fail_ii_indices <- which(bc2_scores < bc2_threshold)
  all_fail_indices <- unique(c(fail_i_indices, fail_ii_indices))

  # Get sample names
  sample_names <- colnames(rg_set)
  fail_names <- sample_names[all_fail_indices]

  # Store results
  failed_samples <- list(
    indices = all_fail_indices,
    names = fail_names,
    bc1_failures = sample_names[fail_i_indices],
    bc2_failures = sample_names[fail_ii_indices]
  )

  thresholds <- list(
    bc1 = list(
      threshold = bc1_threshold,
      threshold_statistical = bc1_threshold_stat,
      mean = bc1_mean,
      sd = bc1_sd,
      multiplier = sd_multiplier
    ),
    bc2 = list(
      threshold = bc2_threshold,
      threshold_statistical = bc2_threshold_stat,
      mean = bc2_mean,
      sd = bc2_sd,
      multiplier = sd_multiplier
    ),
    absolute_min = absolute_min,
    date = Sys.time()
  )

  bisulfite_plot_path <- file.path(context$paths$plots, "bisulfite_conversion_plot.png")
  png(bisulfite_plot_path, width = 1000, height = 600)
  par(mfrow = c(1, 2))

  # BC I plot
  plot(bc1_scores,
    ylab = "BC I Intensity (Green Channel)",
    xlab = "Sample Index",
    main = "Bisulfite Conversion I",
    pch = 19,
    col = ifelse(bc1_scores < bc1_threshold, "red", "black")
  )
  abline(h = bc1_threshold, col = "red", lty = 2)
  legend("topright",
    legend = c("Pass", "Fail", "Threshold"),
    col = c("black", "red", "red"), pch = c(19, 19, NA), lty = c(NA, NA, 2)
  )

  # BC II plot
  plot(bc2_scores,
    ylab = "BC II Intensity (Red Channel)",
    xlab = "Sample Index",
    main = "Bisulfite Conversion II",
    pch = 19,
    col = ifelse(bc2_scores < bc2_threshold, "red", "black")
  )
  abline(h = bc2_threshold, col = "red", lty = 2)
  legend("topright",
    legend = c("Pass", "Fail", "Threshold"),
    col = c("black", "red", "red"), pch = c(19, 19, NA), lty = c(NA, NA, 2)
  )

  dev.off()
  par(mfrow = c(1, 1))

  cat(sprintf("Bisulfite conversion QC: %d sample(s) failed\n", length(fail_names)))
  if (length(fail_names) > 0) {
    cat(paste("  Failed:", paste(fail_names, collapse = ", "), "\n"))
  }

  return(list(failed_samples = failed_samples, thresholds = thresholds))
}

.identify_failed_samples <- function(qc, threshold, prog) {
  bad_m <- qc$mMed < threshold
  bad_u <- qc$uMed < threshold
  bad_indices <- which(bad_m | bad_u)

  bad_names <- if (!is.null(rownames(qc))) {
    rownames(qc)[bad_indices]
  } else if (!is.null(names(qc$mMed))) {
    names(qc$mMed)[bad_indices]
  } else {
    paste0("Sample_", bad_indices)
  }

  result <- list(
    indices = bad_indices,
    names = bad_names,
    m_medians = qc$mMed[bad_indices],
    u_medians = qc$uMed[bad_indices],
    m_failed = bad_m[bad_indices],
    u_failed = bad_u[bad_indices]
  )

  return(result)
}

.determine_qc_threshold <- function(qc, qc_threshold, context, prog) {
  if (is.character(qc_threshold) && qc_threshold == "auto") {
    cat("Auto-detecting QC threshold...\n")

    # Simple auto-detection: Use median - 2*MAD TODO consinder Hampel to change factor to 3.
    m_threshold <- median(qc$mMed) - 2 * mad(qc$mMed)
    u_threshold <- median(qc$uMed) - 2 * mad(qc$uMed)

    threshold <- min(m_threshold, u_threshold)
    threshold <- max(threshold, 10.5) # Don't go below reasonable minimum

    cat(sprintf("Auto-detected threshold: %.2f", threshold))

    auto_info <- list(
      method = "median_minus_2MAD",
      m_median = median(qc$mMed),
      m_mad = mad(qc$mMed),
      u_median = median(qc$uMed),
      u_mad = mad(qc$uMed),
      threshold = threshold,
      date = Sys.time()
    )

    saveRDS(
      auto_info,
      file.path(context$paths$logs, "qc_auto_threshold.rds")
    )
  } else {
    threshold <- as.numeric(qc_threshold)
    cat(sprintf("Using specified threshold: %.2f", threshold))
  }

  return(threshold)
}

.create_qc_log_with_bisulfite <- function(rg_set, bad_samples, qc, intensity_threshold, bisulfite_thresholds, bisulfite_failed_samples) {
  if (length(bad_samples$indices) == 0) {
    return(data.frame(
      Sample_Name = character(0),
      Reason = character(0),
      mMed_Intensity = numeric(0),
      uMed_Intensity = numeric(0),
      BC1_Intensity = numeric(0),
      BC2_Intensity = numeric(0),
      Intensity_Threshold = numeric(0),
      BC1_Threshold = numeric(0),
      BC2_Threshold = numeric(0),
      Date_Removed = character(0)
    ))
  }

  # Get bisulfite scores for all samples
  ctrls_i <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION I")
  ctrls_ii <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION II")
  bc1_all <- colMeans(getGreen(rg_set)[ctrls_i, , drop = FALSE], na.rm = TRUE)
  bc2_all <- colMeans(getRed(rg_set)[ctrls_ii, , drop = FALSE], na.rm = TRUE)

  sample_names <- colnames(rg_set)

  # Build log for each failed sample
  log_entries <- list()

  for (idx in bad_samples$indices) {
    sample_name <- sample_names[idx]
    reasons <- c()

    # Check intensity failure
    if (idx %in% bad_samples$intensity_failures) {
      if (!is.null(qc$mMed[idx]) && qc$mMed[idx] < intensity_threshold) {
        reasons <- c(reasons, "Low methylated intensity")
      }
      if (!is.null(qc$uMed[idx]) && qc$uMed[idx] < intensity_threshold) {
        reasons <- c(reasons, "Low unmethylated intensity")
      }
    }

    # Check bisulfite failure
    if (idx %in% bisulfite_failed_samples$indices) {
      if (idx %in% which(bc1_all < bisulfite_thresholds$bc1$threshold)) {
        reasons <- c(reasons, "Poor BC I conversion")
      }
      if (idx %in% which(bc2_all < bisulfite_thresholds$bc2$threshold)) {
        reasons <- c(reasons, "Poor BC II conversion")
      }
    }

    log_entries[[length(log_entries) + 1]] <- data.frame(
      Sample_Name = sample_name,
      Reason = paste(reasons, collapse = " & "),
      mMed_Intensity = if (!is.null(qc$mMed[idx])) qc$mMed[idx] else NA,
      uMed_Intensity = if (!is.null(qc$uMed[idx])) qc$uMed[idx] else NA,
      BC1_Intensity = bc1_all[idx],
      BC2_Intensity = bc2_all[idx],
      Intensity_Threshold = intensity_threshold,
      BC1_Threshold = bisulfite_thresholds$bc1$threshold,
      BC2_Threshold = bisulfite_thresholds$bc2$threshold,
      Date_Removed = Sys.Date(),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, log_entries)
}

.create_qc_summary_with_bisulfite <- function(rg_set, rg_set_clean, qc_log, intensity_threshold, bisulfite_thresholds) {
  list(
    date = Sys.time(),
    total_samples_initial = ncol(rg_set),
    total_samples_after_qc = ncol(rg_set_clean),
    samples_removed = nrow(qc_log),
    intensity_threshold = intensity_threshold,
    bisulfite_bc1_threshold = bisulfite_thresholds$bc1$threshold,
    bisulfite_bc2_threshold = bisulfite_thresholds$bc2$threshold,
    bisulfite_bc1_mean = bisulfite_thresholds$bc1$mean,
    bisulfite_bc2_mean = bisulfite_thresholds$bc2$mean,
    bisulfite_sd_multiplier = bisulfite_thresholds$bc1$multiplier,
    removal_rate = if (ncol(rg_set) > 0) round(nrow(qc_log) / ncol(rg_set) * 100, 2) else 0,
    removed_samples = if (nrow(qc_log) > 0) qc_log$Sample_Name else "None",
    qc_method = "Median intensity and bisulfite conversion filtering (BC I and BC II)"
  )
}

.determine_qc_threshold <- function(qc, qc_threshold, context, prog) {
  if (is.character(qc_threshold) && qc_threshold == "auto") {
    cat("Auto-detecting QC threshold...\n")

    # Simple auto-detection: Use median - 2*MAD
    m_threshold <- median(qc$mMed) - 2 * mad(qc$mMed)
    u_threshold <- median(qc$uMed) - 2 * mad(qc$uMed)

    threshold <- min(m_threshold, u_threshold)
    threshold <- max(threshold, 10.5) # Don't go below reasonable minimum

    cat(sprintf("Auto-detected threshold: %.2f", threshold))

    auto_info <- list(
      method = "median_minus_2MAD",
      m_median = median(qc$mMed),
      m_mad = mad(qc$mMed),
      u_median = median(qc$uMed),
      u_mad = mad(qc$uMed),
      threshold = threshold,
      date = Sys.time()
    )

    saveRDS(
      auto_info,
      file.path(context$paths$logs, "qc_auto_threshold.rds")
    )
  } else {
    threshold <- as.numeric(qc_threshold)
    cat(sprintf("Using specified threshold: %.2f", threshold))
  }

  return(threshold)
}

.identify_failed_samples <- function(qc, threshold, prog) {
  bad_m <- qc$mMed < threshold
  bad_u <- qc$uMed < threshold
  bad_indices <- which(bad_m | bad_u)

  bad_names <- if (!is.null(rownames(qc))) {
    rownames(qc)[bad_indices]
  } else if (!is.null(names(qc$mMed))) {
    names(qc$mMed)[bad_indices]
  } else {
    paste0("Sample_", bad_indices)
  }

  result <- list(
    indices = bad_indices,
    names = bad_names,
    m_medians = qc$mMed[bad_indices],
    u_medians = qc$uMed[bad_indices],
    m_failed = bad_m[bad_indices],
    u_failed = bad_u[bad_indices]
  )

  return(result)
}

.create_qc_log <- function(rg_set, bad_samples, qc, threshold) {
  if (length(bad_samples$indices) == 0) {
    return(data.frame(
      Sample_Name = character(0),
      Reason = character(0),
      mMed_Intensity = numeric(0),
      uMed_Intensity = numeric(0),
      Threshold_Used = numeric(0),
      Date_Removed = character(0)
    ))
  }

  library(dplyr)
  data.frame(
    Sample_Name = bad_samples$names,
    Reason = paste(
      ifelse(bad_samples$m_failed, "Low methylated intensity", ""),
      ifelse(bad_samples$u_failed, "Low unmethylated intensity", ""),
      sep = ifelse(bad_samples$m_failed & bad_samples$u_failed, " & ", "")
    ) %>% trimws(),
    mMed_Intensity = bad_samples$m_medians,
    uMed_Intensity = bad_samples$u_medians,
    Threshold_Used = threshold,
    Date_Removed = Sys.Date(),
    stringsAsFactors = FALSE
  )
}

.create_qc_summary <- function(rg_set, rg_set_clean, qc_log, threshold) {
  list(
    date = Sys.time(),
    total_samples_initial = ncol(rg_set),
    total_samples_after_qc = ncol(rg_set_clean),
    samples_removed = nrow(qc_log),
    qc_threshold = threshold,
    removal_rate = round(nrow(qc_log) / ncol(rg_set) * 100, 2),
    removed_samples = if (nrow(qc_log) > 0) qc_log$Sample_Name else "None",
    qc_method = "Median intensity filtering"
  )
}
