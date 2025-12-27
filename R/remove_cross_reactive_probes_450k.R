remove_cross_reactive_probes_450k <- function(context = NULL, methyl_set_filename = "methyl_set_filtered_chrom.rds") {

  prog <- .create_progress_manager(4)

  prog$update(1, "Fetching set of Illumina 450k cross-reactive probes")

  ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

  cross_reactive_probes <- read.csv("https://github.com/sirselim/illumina450k_filtering/raw/master/48639-non-specific-probes-Illumina450k.csv", header = TRUE)

  cross_reactive <- cross_reactive_probes$TargetID

  prog$update(2, "Reading methyl set raw data")
  methyl_set_path <- file.path(context$paths$processed, methyl_set_filename)
  methyl_set <- readRDS(methyl_set_path)

  print(paste("Original probes in dataset:", nrow(getBeta(methyl_set))))

  prog$update(3, "Removing cross-reactive probes")
  methyl_set <- methyl_set[!rownames(methyl_set) %in% cross_reactive, ]

  print(paste("Probes after removing cross-reactive probes:", nrow(methyl_set)))

  methyl_set_cross_reactive_clean_path <- file.path(context$paths$processed, "methyl_set_removed_cross_reactive.rds")
  prog$update(4, paste("Saving cleaned-up methyl set under", methyl_set_cross_reactive_clean_path))
  saveRDS(methyl_set, methyl_set_cross_reactive_clean_path)

  prog$complete()

  rm(list = ls())
  gc(full = T)
}
