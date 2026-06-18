remove_snp <- function(
  context = NULL, 
  methyl_set = NULL, 
  methyl_set_file = NULL) {

  prog <- .create_progress_manager(4)
  print("remove single nucleotide polymorphisms")
  ann <- NULL
  if (context$platform == "450k") {
    prog$update(1, "Fetching annotation for Illummina 450k array")
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
  if (context$platform == "EPIC") {
    prog$update(1, "Fetching annotation for Illummina EPIC array")
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }

  prog$update(2, "Reading methyl set data")
  if (is.null(methyl_set)) {
    methyl_set <- readRDS(methyl_set_file)
  }

  print(paste("Original probes in dataset:", nrow(getBeta(methyl_set))))

  prog$update(3, "Map probes to genomic regions")
  methyl_set_genomic <- mapToGenome(methyl_set)

  prog$update(4, "Drop loci from methyl set that contain known SNPs")
  methyl_set_removed_snp <- dropLociWithSnps(methyl_set_genomic, snps = c("CpG", "SBE"), maf = 0.01)

  print(paste("Probes after removing single nucleotide polymorphisms:", nrow(methyl_set_removed_snp)))

  prog$complete()

  methyl_set_removed_snps_path <- file.path(context$paths$processed, "methyl_set_removed_snps.rds")
  methyl_set_removed_snps_container <- new("ResultsContainer", filename = methyl_set_removed_snps_path, object = methyl_set_removed_snp, future = NULL) 

  return(
    list(
      methyl_set_removed_snps_container = methyl_set_removed_snps_container
    )
  )
}
