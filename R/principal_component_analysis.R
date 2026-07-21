principal_component_analysis <- function(
  context = NULL,
  m_values = NULL,
  targets = NULL,
  m_values_filename = NULL,
  targets_filename = NULL, 
  col_maps, keys, npc
) {
  prog <- .create_progress_manager(4)

  if (is.null(targets)) {
    targets <- readRDS(targets_filename)
  }
  if (is.null(m_values)) {
    m_values <- readRDS(m_values_filename)
  }

  prog$update(3, "Running PCA")
  pca <- prcomp(t(m_values), center = TRUE, scale. = FALSE)

  pca_df <- data.frame(matrix(NA, nrow = ncol(m_values), ncol = npc))
  rownames(pca_df) <- colnames(m_values)
  colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
  })
  for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
  }

  for (k in keys) {
    if (!is.null(col_maps[[k]])) {
      pca_df[[k]] <- targets[[col_maps[[k]]]]
    }
  }

  pca_df_filepath <- file.path(context$paths$results, "pca_df.rds")
  
  prog$complete()

  pca_df_container <- new("ResultsContainer", filename = pca_df_filepath, object = pca_df, future = NULL)
  return(
    list(
      pca_df_container = pca_df_container
    )
  )
}
