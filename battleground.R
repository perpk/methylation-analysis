project_to_load <- "GSE111629_20251226_102044"
project_location <- "/Volumes/Elements/methyl-pipe-out"
project_name <- "GSE111629"

source("R/project_context.R")
context <- .load_methylation_project(project_location, project_to_load, platform = "450k", cohorts = cohorts)

cells_only_dmp <- read.csv(file.path(context$paths$results, "CellsOnly_dmp_annotated_results.csv"), row.names = 1)
pcs_cells_dmp <- read.csv(file.path(context$paths$results, "PCplusCells_dmp_annotated_results.csv"), row.names = 1)
pcs_only_dmp <- read.csv(file.path(context$paths$results, "PCs_Only_dmp_annotated_results.csv"), row.names = 1)

cells_only_dmp_sign <- cells_only_dmp[abs(cells_only_dmp["logFC"]) > 0.2 & cells_only_dmp["adj.P.Val"] < 0.05, ]
pcs_cells_dmp_sign <- pcs_cells_dmp[abs(pcs_cells_dmp["logFC"]) > 0.2 & pcs_cells_dmp["adj.P.Val"] < 0.05, ]
pcs_only_dmp_sign <- pcs_only_dmp[abs(pcs_only_dmp["logFC"]) > 0.2 & pcs_only_dmp["adj.P.Val"] < 0.05, ]

dim(pcs_cells_dmp_sign)
dim(cells_only_dmp_sign)
dim(pcs_only_dmp_sign)

dfs <- list(pcs_cells_dmp_sign, cells_only_dmp_sign, pcs_only_dmp_sign)

library(dplyr)
i <- map(dfs, ~.$ProbeID) %>% reduce(intersect)
i
