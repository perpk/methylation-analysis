rm(list = ls())
gc(full = TRUE)

library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(FlowSorted.Blood.EPIC)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE111629_20260709_222156"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "targets.rds"))

source("R/intermediate_data_proxy.R")

source("R/extract_methyl_set.R")
res_extract_methyl_set <- intermediate_data_proxy(
  extract_methyl_set, project_context,
  targets = targets
)

source("R/cell_cnt_estimate.R")
res_cell_cnt_estimate <- intermediate_data_proxy(
  cell_cnt_estimate,
  project_context,
  rg_set = res_extract_methyl_set$rg_set_container@object,
  targets = targets
)

saveRDS(res_cell_cnt_estimate$targets_container@object, file.path(project_context$paths$results, "targets_s_mismatch_cells.rds"))

head(res_cell_cnt_estimate$targets_container@object)
