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

project_to_load <- "GSE145361_20260712_222124"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "targets.rds"))

beta_matrix <- readRDS(file.path(project_context$paths$results, "beta_matrix.rds"))

source("R/intermediate_data_proxy.R")

source("R/apply_BMIQ.R")
source("R/results_container.R")
bmiq_res <- intermediate_data_proxy(
  apply_BMIQ,
  project_context,
  beta_matrix = beta_matrix,
  beta_matrix_file = file.path(project_context$paths$results, "beta_matrix.rds"),
  plot = FALSE,
  save_bmiq = TRUE
)

is.null(bmiq_res)

saveRDS(bmiq_res$beta_bmiq_container@object, file = file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))
saveRDS(bmiq_res$m_bmiq_container@object, file = file.path(project_context$paths$results, "m_values_bmiq.rds"))
