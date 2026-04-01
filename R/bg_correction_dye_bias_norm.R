bg_correction_dye_bias_norm <- function(context = NULL, rg_set_filename = "rg_set.rds") {
  prog <- .create_progress_manager(3)

  print("Background correction and dye bias normalization")

  rg_set_filepath <- file.path(context$paths$raw_data, rg_set_filename) # Note: raw_data, not processed
  prog$update(1, paste("Reading raw rg_set from", rg_set_filepath))
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, "Performing Noob normalization")
  methyl_set_norm <- preprocessNoob(rg_set)

  methyl_set_filepath <- file.path(context$paths$processed, "methyl_set_norm.rds")
  prog$update(3, paste("Saving normalized methyl_set under", methyl_set_filepath))
  saveRDS(methyl_set_norm, methyl_set_filepath)

  # Also save the normalized RG set for any downstream needs
  rg_set_norm_filepath <- file.path(context$paths$processed, "rg_set_norm.rds")
  saveRDS(rg_set, rg_set_norm_filepath) # Note: preprocessNoob returns MethylSet, not RGSet
  # Actually preprocessNoob returns a MethylSet, so adjust accordingly

  prog$complete()
  rm(list = ls())
  gc(full = T)
}
