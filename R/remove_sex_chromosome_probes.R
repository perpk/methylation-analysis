remove_sex_chromosome_probes <- function(context = NULL, methyl_set_filename = "methyl_set_removed_snps.rds") {
  print("remove sex chromosome probes")
  methyl_set <- readRDS(file.path(context$paths$processed, methyl_set_filename))
  autosomal_probes <- !(seqnames(methyl_set) %in% c("chrX", "chrY"))
  methyl_set_final <- methyl_set[autosomal_probes, ]
  print(paste("Probes after removing sex chromosomes:", nrow(methyl_set_final)))
  saveRDS(methyl_set_final, file.path(context$paths$processed, "methyl_set_filtered_chrom.rds"))
}
