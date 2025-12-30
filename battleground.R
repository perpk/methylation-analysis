source('R/project_context.R')

project_context <- .load_methylation_project("/Volumes/Elements/methyl-pipe-out", "GSE111629_20251226_102044")

targets <- readRDS(file.path(project_context$paths$qc, "targets_s_mismatch_cells.rds"))

# BCells, CD4T, CD8T,
