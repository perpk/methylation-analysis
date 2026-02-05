bg_correction_dye_bias_norm <- function(context = NULL, rg_set_filename = "rg_set_filtered.rds") {
  prog <- .create_progress_manager(4)

  print("bg correction dye bias norm")

  rg_set_filepath <- file.path(context$paths$processed, rg_set_filename)
  prog$update(1, paste("Reading rg_set from", rg_set_filepath))
  rg_set <- readRDS(rg_set_filepath)

  prog$update(2, "Performing normal-exponential out-of-band correction with dye-bias normalization")
  methyl_set_norm <- preprocessNoob(rg_set)


  gr_set <- ratioConvert(methyl_set_norm, type = "Illumina")
  gr_set <- mapToGenome(gr_set)

  is_sex_chr <- seqnames(gr_set) %in% c("chrX", "chrY")

  cat("Sex chromosome probes found after normalization:", sum(is_sex_chr), "\n")

  if (any(is_sex_chr)) {
    gr_set <- gr_set[!is_sex_chr, ]
    cat("Probes after sex chromosome removal:", nrow(gr_set), "\n")
  }

  methyl_set_filepath <- file.path(context$paths$processed, "methyl_set_norm.rds")
  prog$update(3, paste("Saving corrected methyl_set under", methyl_set_filepath))
  saveRDS(gr_set, methyl_set_filepath)

  beta_matrix_noob_filepath <- file.path(context$paths$results, "beta_matrix_noob.rds")
  prog$update(4, paste("Extracting beta-matrix and saving under", beta_matrix_noob_filepath))
  beta_matrix_noob <- getBeta(gr_set)
  saveRDS(beta_matrix_noob, beta_matrix_noob_filepath)

  prog$complete()
  rm(list = ls())
  gc(full = T)
}
