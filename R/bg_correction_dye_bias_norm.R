bg_correction_dye_bias_norm <- function(context = NULL, rg_set = NULL, rg_set_filename = NULL) {
  prog <- .create_progress_manager(3)

  print("Background correction and dye bias normalization")

  prog$update(1, "Reading raw rg_set")
  if (is.null(rg_set)) {
    rg_set <- readRDS(rg_set_filename)
  }

  prog$update(2, "Performing Noob normalization")
  methyl_set_norm <- preprocessFunnorm(rg_set, nPCs = 2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE)

  methyl_set_filepath <- file.path(context$paths$processed, "methyl_set_norm.rds")  
  rg_set_norm_filepath <- file.path(context$paths$processed, "rg_set_norm.rds")
  
  methyl_set_container <- new("ResultsContainer", filename = methyl_set_filepath, object = methyl_set_norm, future = NULL)
  rg_set_container <- new("ResultsContainer", filename = rg_set_norm_filepath, object = rg_set, future = NULL)

  prog$complete()

  return(list(methyl_set_container, rg_set_container))
}
