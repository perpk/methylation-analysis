rm(list = ls())
gc(full = TRUE)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE111629_20260709_222156"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$results, "targets_s_mismatch_cells.rds"))
targets_reduced <- readRDS(file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"))

head(rownames(targets))
head(rownames(targets_reduced))

rownames(targets_reduced) <- targets_reduced$Sample_Name
rownames(targets) <- targets$Sample_Name

targets_cells <- targets[rownames(targets) %in% rownames(targets_reduced), ]

library(dplyr)
library(stringr)

head(targets_cells)

targets_cells$Sex <- targets_cells$`gender:ch1`

targets_cells <- targets_cells %>%
    mutate(
        Age_Group =
            case_when(
                `age:ch1` >= 30 & `age:ch1` < 50 ~ "30-49",
                `age:ch1` >= 50 & `age:ch1` < 60 ~ "50-59",
                `age:ch1` >= 60 & `age:ch1` < 70 ~ "60-69",
                `age:ch1` >= 70 & `age:ch1` < 80 ~ "70-79",
                `age:ch1` >= 80 ~ "80+",
                TRUE ~ NA_character_
            )
    )

idat_folder_loc <- "./GSE111629_RAW"
source("R/extract_scandate_from_idat.R")
scan_dates <- extract_scandate_from_idat(
    file_path = idat_folder_loc
)
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

saveRDS(
    enriched_targets,
    file.path(project_context$paths$results, "targets_s_mismatch_cells_scandate.rds")
)

m_values <- readRDS(file.path(project_context$paths$results, "m_values_bmiq_no_outliers.rds"))

m_values_samples <- colnames(m_values)

head(m_values_samples)

rownames(enriched_targets) <- enriched_targets$Basename %>% str_extract("GSM\\d{7}_\\d{10}_R\\d{2}C\\d{2}")

common <- intersect(rownames(enriched_targets), m_values_samples)
length(common)


harmonized_targets <- enriched_targets[rownames(enriched_targets) %in% common, ]
dim(harmonized_targets)

harmonized_m_values <- m_values[, colnames(m_values) %in% common]
dim(harmonized_m_values)

saveRDS(
    harmonized_targets,
    file.path(project_context$paths$processed, "GSE111629_harmonized_targets.rds")
)

saveRDS(
    harmonized_m_values,
    file.path(project_context$paths$processed, "GSE111629_harmonized_m_values.rds")
)
