context <- .load_methylation_project("/Volumes/Elements/methyl-pipe-out", "GSE165081_20251230_222519")

rg_set <- readRDS(file.path(context$paths$processed, "rg_set_filtered.rds"))
targets <- readRDS(file.path(context$paths$processed, "targets_remove_mismatch.rds"))
densityPlot(rg_set, sampGroups = targets$Sample_Group, main = "Beta Value Distribution")

beta_matrix <- readRDS(file.path(context$paths$results, "beta_matrix_bmiq.rds"))

