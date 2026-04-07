cell_cnt_estimate <- function(context = NULL, rg_set_filename = "rg_set.rds", targets_filename = "targets.rds") {
  rg_set_path <- file.path(context$paths$raw_data, rg_set_filename)
  targets_path <- file.path(context$paths$raw_data, targets_filename)

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
}
