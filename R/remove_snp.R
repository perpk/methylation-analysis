remove_snp <- function(context = NULL, methyl_set_file = "methyl_set_remove_mismatch.rds") {

  prog <- .create_progress_manager(5)

  ann <- NULL
  if (context$platform == "450k") {
    prog$update(1, "Fetching annotation for Illummina 450k array")
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
  if (context$platform == "EPIC") {
    prog$update(1, "Fetching annotation for Illummina EPIC array")
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }

  methyl_set_path <- file.path(context$paths$processed, methyl_set_file)
  prog$update(2, paste("Reading methyl set from", methyl_set_path))
  methyl_set <- readRDS(methyl_set_path)

  print(paste("Original probes in dataset:", nrow(getBeta(methyl_set))))

  prog$update(3, "Map probes to genomic regions")
  methyl_set_genomic <- mapToGenome(methyl_set)

  prog$update(4, "Drop loci from methyl set that contain known SNPs")
  methyl_set_removed_snp <- dropLociWithSnps(methyl_set_genomic, snps = c("CpG", "SBE"), maf=0.01)

  print(paste("Probes after removing single nucleotide polymorphisms:", nrow(methyl_set_removed_snp)))

  methyl_set_removed_snps_path <- file.path(context$paths$processed, "methyl_set_removed_snps.rds")
  prog$update(5, paste("Saving SNP cleaned-up methyl set under", methyl_set_removed_snps_path))
  saveRDS(methyl_set_removed_snp, methyl_set_removed_snps_path)

  prog$complete()
}
