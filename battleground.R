project_to_load <- "GSE111629_20251226_102044"
project_location <- "/Volumes/Elements/methyl-pipe-out"
project_name <- "GSE111629"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)

source("R/project_context.R")
context <- .load_methylation_project(project_location, project_to_load, platform = "450k", cohorts = cohorts)

rg_set_filename <- "rg_set_clean.rds"
methyl_set_filename <- "methyl_set_clean.rds"
det_p_threshold <- 0.01
max_failed_samples <- 0.05

rg_set_filepath <- file.path(context$paths$qc, rg_set_filename)
rg_set <- readRDS(rg_set_filepath)

# methyl_set_filepath <- file.path(context$paths$qc, methyl_set_filename)
# methyl_set <- readRDS(methyl_set_filepath)

library(minfi)
det_p <- detectionP(rg_set)
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

removed_df <- data.frame(
  probe_id = removed_probes,
  failure_rate_percent = failure_rates,
  reason = ifelse(failure_rates >= max_failed_samples * 100,
    paste0("Failed in >= ", max_failed_samples * 100, "% samples"),
    "Other QC failure"
  ),
  stringsAsFactors = FALSE
)

rg_set_filtered <- rg_set[keep_probes, ]

methyl_set_filepath <- file.path(context$paths$qc, methyl_set_filename)
methyl_set <- readRDS(methyl_set_filepath)

methyl_set_filtered <- methyl_set[keep_probes, ]
