remove_sex_chromosome_probes <- function(
  context = NULL, 
  methyl_set = NULL,
  methyl_set_filename = NULL
) {
  print("remove sex chromosome probes")
  if (is.null(methyl_set)) {
    methyl_set <- readRDS(methyl_set_filename)
  }
  autosomal_probes <- !(seqnames(methyl_set) %in% c("chrX", "chrY"))
  methyl_set_final <- methyl_set[autosomal_probes, ]
  print(paste("Probes after removing sex chromosomes:", nrow(methyl_set_final)))

  mismatch_container <- new("ResultsContainer", filename = file.path(context$paths$processed, "methyl_set_filtered_chrom.rds"), object = methyl_set_final, future = NULL)

  return(
    list(
      mismatch_container = mismatch_container
    )
  )
}
