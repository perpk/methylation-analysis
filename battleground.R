source('R/progress_mgr.R')

source('R/project_context.R')
context <- .load_methylation_project(base_dir = "/Volumes/Elements/methyl-pipe-out", project_id = "ppmi_20260103_120001", platform = "EPIC")

methyl_set_filepath <- file.path(context$paths$raw_data, "methyl_set.rds")
methyl_set <- readRDS(methyl_set_filepath)

library(minfi)
source('R/qc.R')

qc <- getQC(methyl_set)
plotQC(qc)

prog <- .create_progress_manager(1)
threshold <- .determine_qc_threshold(qc, "auto", context, prog)
bad_samples <- .identify_failed_samples(qc, threshold, prog)

rm(list = setdiff(ls(), c("targets", "bad_samples", "context")))

rg_set_filepath <- file.path(context$paths$raw_data, "rg_set.rds")
rg_set <- readRDS(rg_set_filepath)
rg_set_clean <- rg_set[, -bad_samples$indices, drop = FALSE]
targets_clean <- readRDS(file.path(context$paths$qc, "targets_clean.rds"))
rm(list = c("rg_set"))
gc(full=TRUE)
densityPlot(rg_set_clean, sampGroups = targets_clean$Sample_Group, main = "Beta Value Distribution")

rm(list = ls())
gc(full = TRUE)

cross_reactive_probes <- read.csv("https://github.com/sirselim/illumina450k_filtering/raw/master/48639-non-specific-probes-Illumina450k.csv", header = TRUE)

dim(cross_reactive_probes)

xloci <- maxprobes::xreactive_probes(array_type = "450K")
length(xloci)
class(xloci)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annlib <- IlluminaHumanMethylation450kanno.ilmn12.hg19
ann <- getAnnotation(annlib)
class(ann)
attr(ann, "annotation")
ann@package


m <- readRDS("/Volumes/Elements/methyl-pipe-out/ppmi_20260103_120001/qc/methyl_set_clean.rds")
print(annotation(m))

rm(m)
gc(full = TRUE)

context <- .load_methylation_project(base_dir = "/Volumes/Elements/methyl-pipe-out", project_id = "ppmi_20260103_120001", platform = "EPIC")
context$platform <- "EPIC"
context$platform

saveRDS(context, file.path(context$base_dir, "project_context.rds"))
