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
library(tidyverse)

install.packages(c("tidyverse", "mirai"))

source("R/progress_mgr.R")
source("R/project_context.R")


project_name <- "GSE145361"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)
mode <- results_mode()$memory_only

source("R/intermediate_data_proxy.R")
source("R/results_container.R")
project_context <- create_methylation_project(project_name, project_location, platform = platform, cohorts = cohorts, mode = mode)

data_folder <- paste0(project_name, "_RAW")

targets <- read.metharray.sheet(data_folder, pattern = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")

gse <- getGEO(project_name, GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])
all_idat_files <- list.files(data_folder, pattern = "\\.idat\\.gz$", full.names = FALSE)
basenames <- unique(gsub("_(Grn|Red).idat.gz", "", all_idat_files))
targets <- data.frame(
    Sample_Name = gsub("_.*", "", basenames),
    Basename = file.path(data_folder, basenames),
    stringsAsFactors = FALSE
)

targets <- merge(targets, pdata, by.x = "Sample_Name", by.y = "row.names", all.x = TRUE)

library(tidyverse)
targets$Sample_Group <- factor(targets$`disease-state:ch1`)
targets <- targets %>% mutate(Sample_Group = case_when(Sample_Group == "Parkinson's disease" ~ "PD", Sample_Group == "Control" ~ "Control"))
targets$Sample_Group <- as.factor(targets$Sample_Group)


source("R/extract_methyl_set.R")
res <- extract_methyl_set(
  context = project_context,
  targets = targets
)

is.null(res$m_set_container@object)
is.null(res$rg_set_container@object)

source("R/qc.R")
res_qc <- qc(context = project_context,
  methyl_set = res$m_set_container@object,
  rg_set = res$rg_set_container@object,
  targets = targets
)

saveRDS(res_qc$failed_samples_results@object, file = file.path(project_context$paths$qc, "failed_samples.rds"))

saveRDS(res_qc$bisulfite_thresholds_results@object, file = file.path(project_context$paths$qc, "bisulfite_thresholds.rds"))
saveRDS(res_qc$qc_results@object, file = file.path(project_context$paths$qc, "qc_results.rds"))
saveRDS(res_qc$failed_samples_results@object, file = file.path(project_context$paths$qc, "failed_samples.rds"))


m_set_clean <- res$m_set_container@object
rg_set_clean <- res$rg_set_container@object

source("R/bg_correction_dye_bias_norm.R")
res_bg_corr <- bg_correction_dye_bias_norm(
  context = project_context,
  rg_set = rg_set_clean
)

is.null(res_bg_corr$rg_set_container@object)

source("R/probe_qc.R")
res_probe_qc <- probe_qc(
  context = project_context,
  rg_set = res_bg_corr$rg_set_container@object
)

removed_probes_container <- res_probe_qc$removed_probes_container@object
qc_summary_container <- res_probe_qc$qc_summary_container@object

saveRDS(removed_probes_container, file = res_probe_qc$removed_probes_container@filename)
saveRDS(qc_summary_container, file = res_probe_qc$qc_summary_container@filename)
