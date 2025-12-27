extract_methyl_set <- function(context = NULL, targets = NULL) {
  prog <- .create_progress_manager(2)

  prog$update(1, "Reading raw intensity data from .idat files")
  rg_set <- read.metharray.exp(targets = targets, extended = TRUE)

  rg_set_filepath <- file.path(context$paths$raw_data, "rg_set.rds")
  saveRDS(rg_set, rg_set_filepath)
  print(paste("Saved rg_set under", rg_set_filepath))

  prog$update(2, "Reading raw intensity data from .idat files")
  methyl_set <- preprocessRaw(rg_set)

  methyl_set_filepath <- file.path(context$paths$raw_data, "methyl_set.rds")
  saveRDS(methyl_set, methyl_set_filepath)
  print(paste("Saved methyl_set under, methyl_set_filepath"))

  prog$complete()
}
