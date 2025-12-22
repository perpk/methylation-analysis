qc <- function(targets,
               context = NULL,
               progress = TRUE,
               qc_threshold = 10.5,
               force = FALSE) {

  prog <- .create_progress_manager(5)

  if (is.null(context)) {
    context <- .project_context$get()
  }

  prog$update(1, "Reading array data")

  rg_set_path <- file.path(context$paths$raw_data, "rg_set.rds")

  if (!force && file.exists(rg_set_path)) {
    cat("Loading existing RGChannelSet...\n")
    rg_set <- readRDS(rg_set_path)
  } else {
    tryCatch({
      cat("Reading methylation array data...\n")
      rg_set <- read.metharray.exp(targets = targets, extended = TRUE)

      cat("Saving RGChannelSet...\n")
      saveRDS(rg_set, rg_set_path)
    }, error = function(e) {
      stop("Failed to read methylation array data: ", e$message)
    })
  }

  prog$update(2, "Performing QC")

  tryCatch({
    cat("Computing raw intensities...")
    raw_data <- preprocessRaw(rg_set)
    qc <- getQC(raw_data)

    saveRDS(qc, file.path(context$paths$qc, "qc_metrics.rds"))

    cat("\nCreating QC plot...\n")
    qc_plot_path <- file.path(context$paths$plots, "qc_intensity_plot.png")

    png(qc_plot_path, width = 800, height = 600)
    plotQC(qc)
    dev.off()

  }, error = function(e) {
    stop("Failed during QC computation: ", e$message)
  })

  prog$update(3, "Detecting failed samples")

  threshold <- .determine_qc_threshold(qc, qc_threshold, context, prog)

  prog$update(4, "Filtering samples")

  tryCatch({
    bad_samples <- .identify_failed_samples(qc, threshold, prog)

    if (length(bad_samples$indices) > 0) {
      rg_set_clean <- rg_set[, -bad_samples$indices, drop = FALSE]
      targets_clean <- targets[-bad_samples$indices, , drop = FALSE]

      cat(paste("Removing", length(bad_samples$indices),
                         "failed sample(s):"))
      cat(paste(bad_samples$names, collapse = ", "))
    } else {
      rg_set_clean <- rg_set
      targets_clean <- targets
      cat("No samples failed QC\n")
    }

  }, error = function(e) {
    stop("Failed during sample filtering: ", e$message)
  })

  prog$update(5, "Saving results")

  tryCatch({
    saveRDS(rg_set_clean,
            file.path(context$paths$qc, "rg_set_clean.rds"))
    saveRDS(targets_clean,
            file.path(context$paths$qc, "targets_clean.rds"))

    qc_log <- .create_qc_log(rg_set, bad_samples, qc, threshold)
    write.csv(qc_log,
              file.path(context$paths$logs, "sample_removal_log.csv"),
              row.names = FALSE)

    summary <- .create_qc_summary(rg_set, rg_set_clean, qc_log, threshold)
    saveRDS(summary,
            file.path(context$paths$qc, "qc_summary.rds"))

  }, error = function(e) {
    stop("Failed to save QC results: ", e$message)
  })

  rm(list = setdiff(ls(), c("prog", "qc_metrics", "qc_log", "qc_summary", "rg_set_path", "qc_plot_path", "targets_clean", "rg_set_clean", "context")))
  gc(full=T)

  cat("Creating Density Plot...\n")
  density_plot_path <- file.path(context$paths$plots, "density_plot.png")

  png(density_plot_path, width = 800, height = 600)
  densityPlot(rg_set_clean, sampGroups = targets$Sample_Group, main = "Beta Value Distribution")
  dev.off()

  prog$complete()

  list(
    targets_clean = targets_clean,
    qc_metrics = qc,
    qc_log = qc_log,
    qc_summary = summary,
    paths = list(
      rg_set = rg_set_path,
      qc_plot = qc_plot_path,
      cleaned_data = file.path(context$paths$qc, "rg_set_clean.rds")
    )
  )
}

#' Determine QC threshold (auto or manual)
#' @keywords internal
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

#' Identify failed samples based on QC metrics
#' @keywords internal
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

#' Create QC log data frame
#' @keywords internal
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

#' Create QC summary
#' @keywords internal
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
