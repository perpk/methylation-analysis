qc <- function(context = NULL,
               targets = NULL,
               rg_set_filename = "rg_set.rds",
               methyl_set_filename = "methyl_set.rds",
               qc_threshold = 10.5) {

  prog <- .create_progress_manager(7)

  methyl_set_filepath <- file.path(context$paths$raw_data, methyl_set_filename)
  rg_set_filepath <- file.path(context$paths$raw_data, rg_set_filename)
  prog$update(1, paste("Reading methyl set from ", methyl_set_filepath))
  methyl_set <- readRDS(methyl_set_filepath)
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, "Performing QC")

  qc <- getQC(methyl_set)
  saveRDS(qc, file.path(context$paths$qc, "qc_metrics.rds"))

  prog$update(3, "Creating QC plot")
  qc_plot_path <- file.path(context$paths$plots, "qc_intensity_plot.png")

  png(qc_plot_path, width = 800, height = 600)
  plotQC(qc)
  dev.off()

  prog$update(4, "Detecting failed samples")
  threshold <- .determine_qc_threshold(qc, qc_threshold, context, prog)


  prog$update(5, "Filtering samples")

  bad_samples <- .identify_failed_samples(qc, threshold, prog)

  if (length(bad_samples$indices) > 0) {
    rg_set_clean <- rg_set[, -bad_samples$indices, drop = FALSE]
    targets_clean <- targets[-bad_samples$indices, , drop = FALSE]
    methyl_set_clean <- methyl_set[, -bad_samples$indices, drop = FALSE]

    cat(paste("Removing", length(bad_samples$indices),
              "failed sample(s):"))
    cat(paste(bad_samples$names, collapse = ", "))
  } else {
    rg_set_clean <- rg_set
    targets_clean <- targets
    methyl_set_clean <- methyl_set
    cat("No samples failed QC\n")
  }

  rm(list = c("targets", "methyl_set"))
  gc(full = T)

  prog$update(6, "Saving results")

  saveRDS(rg_set_clean,
          file.path(context$paths$qc, "rg_set_clean.rds"))
  saveRDS(targets_clean,
          file.path(context$paths$qc, "targets_clean.rds"))
  saveRDS(methyl_set_clean,
          file.path(context$paths$qc, "methyl_set_clean.rds"))

  qc_log <- .create_qc_log(rg_set, bad_samples, qc, threshold)
  write.csv(qc_log,
            file.path(context$paths$logs, "sample_removal_log.csv"),
            row.names = FALSE)

  summary <- .create_qc_summary(rg_set, rg_set_clean, qc_log, threshold)
  saveRDS(summary,
          file.path(context$paths$qc, "qc_summary.rds"))


  rm(list = setdiff(ls(), c("prog", "qc", "qc_log", "qc_summary", "rg_set_filepath", "qc_plot_path", "targets_clean", "rg_set_clean", "context")))
  gc(full=T)

  prog$update(7, "Creating Density Plot...")
  density_plot_path <- file.path(context$paths$plots, "density_plot.png")

  png(density_plot_path, width = 800, height = 600)
  densityPlot(rg_set_clean, sampGroups = targets_clean$Sample_Group, main = "Beta Value Distribution")
  dev.off()

  prog$complete()
}

.determine_qc_threshold <- function(qc, qc_threshold, context, prog) {

  if (is.character(qc_threshold) && qc_threshold == "auto") {
    cat("Auto-detecting QC threshold...\n")

    # Simple auto-detection: Use median - 2*MAD
    m_threshold <- median(qc$mMed) - 2 * mad(qc$mMed)
    u_threshold <- median(qc$uMed) - 2 * mad(qc$uMed)

    threshold <- min(m_threshold, u_threshold)
    threshold <- max(threshold, 8)  # Don't go below reasonable minimum

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

    saveRDS(auto_info,
            file.path(context$paths$logs, "qc_auto_threshold.rds"))

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
