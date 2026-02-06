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
targets_filename <- "targets_remove_mismatch.rds"
det_p_threshold <- 0.01
max_failed_samples <- 0.05

rg_set_path <- file.path(context$paths$processed, rg_set_filename)
targets_path <- file.path(context$paths$processed, targets_filename)

print("Estimating cell counts")

prog <- .create_progress_manager(3)

prog$update(1, "Reading files")
rg_set <- readRDS(rg_set_path)
targets <- readRDS(targets_path)

prog$update(2, "Estimating cell counts")
cell_counts <- estimateCellCounts(rg_set, compositeType = "Blood", probeSelect = "IDOL")
targets <- cbind(targets, cell_counts)

prog$update(3, "Saving results in targets as metadata for downstream analysis")
targets_cell_types_path <- file.path(context$paths$qc, "targets_s_mismatch_cells.rds")
saveRDS(targets, targets_cell_types_path)
rm(list = ls())
gc(full = T)
