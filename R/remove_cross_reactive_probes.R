library(maxprobes)
remove_cross_reactive_probes <- function(
  context = NULL, 
  methyl_set = NULL, 
  methyl_set_filename = NULL
) {
  prog <- .create_progress_manager(3)

  print("remove cross reactive probes")

  prog$update(1, "Fetching set of cross-reactive probes")

  annlib <- NULL
  if (context$platform == "450k") {
    print("Platform is 450K")
    annlib <- IlluminaHumanMethylation450kanno.ilmn12.hg19
  } else if (context$platform == "EPIC") {
    print("Platform is EPIC")
    annlib <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19
  } else {
    stop("Platform must be specified as '450K' or 'EPIC'")
  }
  ann <- getAnnotation(annlib)

  # cross_reactive_probes <- read.csv("https://github.com/sirselim/illumina450k_filtering/raw/master/48639-non-specific-probes-Illumina450k.csv", header = TRUE)

  # cross_reactive <- cross_reactive_probes$TargetID

  prog$update(2, "Reading methyl set raw data")
  if (is.null(methyl_set)) {
    methyl_set <- readRDS(methyl_set_filename)
  }

  print(paste("Original probes in dataset:", nrow(getBeta(methyl_set))))

  prog$update(3, "Removing cross-reactive probes")
  # MsetEx_noXloci <- maxprobes::dropXreactiveLoci(MsetEx)
  methyl_set <- maxprobes::dropXreactiveLoci(methyl_set)
  # methyl_set <- methyl_set[!rownames(methyl_set) %in% cross_reactive, ]

  print(paste("Probes after removing cross-reactive probes:", nrow(methyl_set)))

  prog$complete()

  methyl_set_cross_reactive_clean_path <- file.path(context$paths$processed, "methyl_set_removed_cross_reactive.rds")
  methyl_set_cross_reactive_clean_container <- new("ResultsContainer", filename = methyl_set_cross_reactive_clean_path, object = methyl_set, future = NULL)   

  return(
    list(
      methyl_set_cross_reactive_clean_container = methyl_set_cross_reactive_clean_container
    )
  )
}
