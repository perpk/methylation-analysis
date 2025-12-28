filter_rg_set <- function(context = NULL, methyl_set_filename = "methyl_set_removed_cross_reactive.rds", rg_set_filename = "rg_set_remove_mismatch.rds") {
  prog <- .create_progress_manager(5)

  rg_set_filepath <- file.path(context$paths$processed, rg_set_filename)
  methyl_set_filepath <- file.path(context$paths$processed, methyl_set_filename)

  prog$update(1, paste("Reading rg_set from", rg_set_filepath))
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, paste("Reading methyl set from ", methyl_set_filepath))
  methyl_set <- readRDS(methyl_set_filepath)

  prog$update(3, "Retrieving annotations for rg_set")
  annotation <- getAnnotation(rg_set)

  prog$update(4, "Filtering rg_set with methyl_set's probes")
  filtered_probes <- rownames(methyl_set)
  rm(methyl_set)
  filtered_rg_set <- rg_set[annotation$Name %in% filtered_probes, ]

  filtered_rg_set_filepath <- file.path(context$paths$processed, "rg_set_filtered.rds")
  prog$update(5, paste("Saving filtered rg_set under", filtered_rg_set_filepath))
  saveRDS(filtered_rg_set, filtered_rg_set_filepath)

  prog$complete()
  rm(filtered_rg_set)
  gc(full = T)
}
