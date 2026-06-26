source("R/results_container.R")

extract_methyl_set <- function(context = NULL, read_chunkwise = FALSE, targets = NULL) {
  prog <- .create_progress_manager(2)

  prog$update(1, "Reading raw intensity data from .idat files")
  
  if (read_chunkwise) {
    rg_set <- .read_in_chunks(targets)
  } else {
    rg_set <- read.metharray.exp(targets = targets, extended = TRUE)
  }
  
  rg_set <- read.metharray.exp(targets = targets, extended = TRUE)

  prog$update(2, "Reading raw intensity data from .idat files")
  methyl_set <- preprocessRaw(rg_set)

  rg_set_container <- new("ResultsContainer", filename = file.path(context$paths$raw_data, "rg_set.rds"), object = rg_set, future = NULL) 
  m_set_container <- new("ResultsContainer", filename = file.path(context$paths$raw_data, "methyl_set.rds"), object = methyl_set, future = NULL)
  results <- list(rg_set_container = rg_set_container, m_set_container = m_set_container)
  prog$complete()
  return (results)
}

.read_in_chunks <- function(targets, chunk_size = 10) {
  target_chunks <- split(targets, ceiling(seq_len(nrow(targets)) / chunk_size))

  rg_set_list <- list()
  for (i in seq_along(target_chunks)) {
    message("Processing chunk ", i, " of ", length(target_chunks))
    rg_set_chunk <- read.metharray.exp(targets = target_chunks[[i]], extended = TRUE, force = TRUE)
    rg_set_list[[i]] <- rg_set_chunk
    rm(rg_set_chunk)
    gc(full = TRUE)
  }
  (Reduce(combineArrays, rg_set_list))
}
