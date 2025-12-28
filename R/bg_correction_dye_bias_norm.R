bg_correction_dye_bias_norm <- function(context = NULL, rg_set_filename = "rg_set_filtered.rds") {
  prog <- .create_progress_manager(4)

  rg_set_filepath <- file.path(context$paths$processed, rg_set_filename)
  prog$update(1, paste("Reading rg_set from", rg_set_filepath))
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, "Performing normal-exponential out-of-band correction with dye-bias normalization")
  methyl_set_norm <- preprocessNoob(rg_set)

  methyl_set_filepath <- file.path(context$paths$processed, "methyl_set_norm.rds")
  prog$update(3, paste("Saving corrected methyl_set under", methyl_set_filepath))
  saveRDS(methyl_set_norm, methyl_set_filepath)

  beta_matrix_noob_filepath <- file.path(context$paths$results, "beta_matrix_noob.rds")
  prog$update(4, paste("Extracting beta-matrix and saving under", beta_matrix_noob_filepath))
  beta_matrix_noob <- getBeta(methyl_set_norm)
  saveRDS(beta_matrix_noob, beta_matrix_noob_filepath)

  prog$complete()
  rm(list = ls())
  gc(full = T)
}
