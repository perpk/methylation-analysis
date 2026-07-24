rm(list = ls())
gc(full = TRUE)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "ppmi_20260721_075730"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "EPIC"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "targets_after_cell_count_estimation.rds"))
targets_reduced <- readRDS(file.path(project_context$paths$processed, "targets_remove_mismatch.rds"))

head(rownames(targets))
head(rownames(targets_reduced))

rownames(targets_reduced) <- targets_reduced$Sample_Name
rownames(targets) <- targets$Sample_Name

targets_cells <- targets[rownames(targets_reduced) %in% rownames(targets), ]

targets_cells %>% colnames()

library(dplyr)
library(stringr)

head(targets_cells)


targets_cells$Sex <- targets_cells$SEX

idat_folder_loc <- "./ppmi/Project120_IDATS_n524final_toLONI_030718"
source("R/extract_scandate_from_idat.R")
scan_dates <- extract_scandate_from_idat(
    file_path = idat_folder_loc
)
scan_dates %>% head()
scan_dates <- scan_dates[!duplicated(scan_dates[c("SentrixID", "ScanDate")]), ]

targets_cells$Basename
targets_cells$Sentrix_ID <- sapply(targets_cells$Basename, function(x) {
    basename(x) %>% str_extract("\\d+_R\\d{2}C\\d{2}")
})

enriched_targets <- merge(
    x = targets_cells,
    y = scan_dates,
    by.x = "Sentrix_ID",
    by.y = "SentrixID",
    all.x = FALSE
)
enriched_targets %>% colnames()
saveRDS(
    enriched_targets,
    file.path(project_context$paths$results, "targets_s_mismatch_cells_scandate.rds")
)

m_values <- readRDS(file.path(project_context$paths$results, "m_values_bmiq_no_outliers.rds"))

m_values_samples <- colnames(m_values)

head(m_values_samples)

enriched_targets %>%
    rownames() %>%
    head()

rownames(enriched_targets) <- enriched_targets$Sentrix_ID

common <- intersect(rownames(enriched_targets), m_values_samples)
length(common)

harmonized_targets <- enriched_targets[rownames(enriched_targets) %in% common, ]
dim(harmonized_targets)

harmonized_m_values <- m_values[, colnames(m_values) %in% common]
dim(harmonized_m_values)
harmonized_targets %>% head()
saveRDS(
    harmonized_targets,
    file.path(project_context$paths$processed, "ppmi_harmonized_targets.rds")
)

saveRDS(
    harmonized_m_values,
    file.path(project_context$paths$processed, "ppmi_harmonized_m_values.rds")
)
