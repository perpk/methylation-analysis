normalize_extract_methylset <- function(context = NULL, force=F) {
  if (is.null(context)) {
    context <- .project_context$get()
  }

  cat("Pre-processing and background correction\n")

  m_filepath <- file.path(context$paths$normalized, "methyl_set.rds")
  if (!force && file.exists(m_filepath)) {
    cat(paste("Read methylset from", m_filepath))
    methyl_set <- readRDS(m_filepath)
  } else {
    cat("Extract methylset\n")
    rg_set_path <- file.path(context$paths$qc, "rg_set_clean.rds")
    rg_set_clean <- readRDS(rg_set_path)
    methyl_set <- preprocessNoob(rg_set_clean)
    saveRDS(methyl_set, m_filepath)
    cat(paste("Saved methylset under\n", m_filepath))
    rm(rg_set_clean)
  }

  cat("\nExtract beta-values\n")
  b_val_filepath <- file.path(context$paths$processed, "beta_matrix.rds")
  beta_matrix <- getBeta(methyl_set)
  saveRDS(beta_matrix, b_val_filepath)
  cat(paste("Saved beta-values under", b_val_filepath, "\n"))

  cat("\nExtract m-values\n")
  m_val_filepath <- file.path(context$paths$processed, "m_values.rds")
  m_values <- getM(methyl_set)
  saveRDS(m_values, m_val_filepath)
  cat(paste("Saved m-values under", m_val_filepath, "\n"))

  list(
    m_values = m_values,
    beta_matrix = beta_matrix,
    paths = list(
        m_values = m_val_filepath,
        beta_values = b_val_filepath
      )
  )
}
