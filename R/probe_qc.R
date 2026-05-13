library(minfi)

probe_qc <- function(
  context = NULL,
  rg_set_filename = "rg_set.rds",
  det_p_threshold = 0.01,
  max_failed_samples = 0.05
) {
  tryCatch({
    prog <- .create_progress_manager(5)

    rg_set_filepath <- file.path(context$paths$raw, rg_set_filename)
    prog$update(1, paste("Reading cleaned RG set from", rg_set_filepath))
    rg_set <- readRDS(rg_set_filepath)

    prog$update(2, "Calculating detection p-values")
    det_p <- detectionP(rg_set)

    prog$update(3, "Identifying failed probes")
    failed <- det_p > det_p_threshold
    prop_failed <- rowMeans(failed, na.rm = TRUE)

    keep_probes <- prop_failed < max_failed_samples

    cat("\n=== Probe QC Results ===\n")
    cat("Total probes:", length(keep_probes), "\n")
    cat("Probes kept:", sum(keep_probes), "\n")
    cat("Probes removed:", sum(!keep_probes), "\n")
    cat("Percentage removed:", round(mean(!keep_probes) * 100, 2), "%\n")

    removed_indices <- which(!keep_probes)
    removed_probes <- rownames(rg_set)[removed_indices]
    failure_rates <- prop_failed[removed_indices] * 100

    if (length(removed_probes) != length(failure_rates)) {
      cat("WARNING: Length mismatch! Debug info:\n")
      cat("  removed_probes length:", length(removed_probes), "\n")
      cat("  failure_rates length:", length(failure_rates), "\n")
      cat("  removed_indices length:", length(removed_indices), "\n")
      cat("  prop_failed length:", length(prop_failed), "\n")
    }

    prog$update(5, "Creating removal records")
    removed_df <- data.frame(
      probe_id = removed_probes,
      failure_rate_percent = failure_rates,
      reason = ifelse(failure_rates >= max_failed_samples * 100,
        paste0("Failed in >= ", max_failed_samples * 100, "% samples"),
        "Other QC failure"
      ),
      stringsAsFactors = FALSE
    )

    saveRDS(
      removed_df,
      file.path(context$paths$qc, "removed_probes_detection_p.rds")
    )

    qc_summary <- list(
      date = Sys.time(),
      det_p_threshold = det_p_threshold,
      max_failed_samples = max_failed_samples,
      total_probes_initial = nrow(rg_set),
      total_probes_after = nrow(rg_set_filtered),
      probes_removed = sum(!keep_probes),
      removal_rate = round(mean(!keep_probes) * 100, 2),
      median_failure_rate_all = median(prop_failed * 100, na.rm = TRUE),
      median_failure_rate_removed = median(failure_rates, na.rm = TRUE),
      worst_probe_failure_rate = max(prop_failed * 100, na.rm = TRUE)
    )

    saveRDS(qc_summary, file.path(context$paths$qc, "probe_qc_summary.rds"))

    cat("\n=== Probe QC Complete ===\n")
    cat("Summary:\n")
    print(data.frame(
      Metric = names(qc_summary),
      Value = unlist(qc_summary)
    ))
  }, error = function(e) {
    cat("Error in probe_qc:", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
    stop(e)
  }, finally = {
    rm(list = ls())
    gc(full = T)
  })
}
