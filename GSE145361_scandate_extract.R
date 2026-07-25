rm(list = ls())
gc(full = TRUE)

source("R/progress_mgr.R")
source("R/project_context.R")

project_to_load <- "GSE145361_20260714_080452"
project_location <- "/root/workspace/methyl-pipe-out"
platform <- "450k"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)

targets <- readRDS(file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"))
targets_cells <- readRDS(file.path(project_context$paths$processed, "targets_s_mismatch_cells.rds"))

dim(targets)
dim(targets_cells)

head(rownames(targets))
head(rownames(targets_cells))
head(targets$Sample_Name)

library(dplyr)
library(stringr)

targets_cells %>%
    rownames() %>%
    str_detect("\\d{12}_R\\d{2}C\\d{2}") %>%
    sum() - dim(targets_cells)[1]
targets %>%
    pull(Sample_Name) %>%
    str_detect("GSM") %>%
    sum() - dim(targets)[1]

cols <- c("Sample_Name", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
data <- targets_cells %>% select(all_of(cols))

merged <- merge(
    x = targets,
    y = data,
    by.x = "Sample_Name",
    by.y = "Sample_Name",
    all.x = TRUE
)

merged <- merged %>% mutate(Sample_Group = case_when(is.na(Sample_Group) ~ "Control", TRUE ~ Sample_Group))
head(merged)

dim(merged)
View(merged)

merged[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")] %>%
    is.na() %>%
    sum()

merged$Sentrix_ID <- sapply(merged$Basename, function(x) {
    basename(x) %>% str_extract("\\d+_R\\d{2}C\\d{2}")
})
merged$Sentrix_ID %>%
    is.na() %>%
    sum()
merged$Sentrix_ID
scan_dates$SentrixID

idat_folder_loc <- "/workspace/methylation-analysis/GSE145361_RAW"
source("R/extract_scandate_from_idat.R")
scan_dates <- extract_scandate_from_idat(
    file_path = idat_folder_loc, format = "%d/%m/%Y"
)


# idat <- readIDAT("/workspace/methylation-analysis/GSE145361_RAW/GSM4315521_3998888024_R01C01_Grn.idat.gz")
# idat$RunInfo
# run_metadata <- idat$RunInfo
# scan_row <- which(run_metadata[, "BlockType"] == "Scan")[1]
# scan_date_string <- run_metadata[scan_row, "RunTime"]
# scan_date_string
# as.POSIXct(scan_date_string, format = "%d/%m/%Y") %>%
#     as.character() %>%
#     str_extract("^(\\d{4}-\\d{2})")



scan_dates$ScanDate %>%
    is.na() %>%
    sum()
scan_dates$SentrixID %>%
    is.na() %>%
    sum()

scan_dates <- scan_dates[!duplicated(scan_dates[c("SentrixID", "ScanDate")]), ]
dim(scan_dates)
dim(merged)
head(scan_dates)
enriched_targets <- merge(
    x = merged,
    y = scan_dates,
    by.x = "Sentrix_ID",
    by.y = "SentrixID",
    all.x = TRUE,
    all.y = FALSE
)

intersect(merged$Sentrix_ID, scan_dates$SentrixID) %>% length()

enriched_targets$ScanDate %>%
    is.na() %>%
    sum()
scan_dates$ScanDate %>%
    is.na() %>%
    sum()
dim(enriched_targets)
dim(merged)
View(targets)
View(enriched_targets)

m_values <- readRDS(file.path(project_context$paths$results, "m_values_bmiq.rds"))
dim(m_values)

m_values_samples <- colnames(m_values)

duplicates <- enriched_targets$Basename %>%
    str_extract("GSM\\d{7}_\\d{12}_R\\d{2}C\\d{2}") %>%
    duplicated()
View(enriched_targets[duplicates, ])

dim(enriched_targets[enriched_targets$Basename %>% str_extract("GSM\\d{7}_\\d{12}_R\\d{2}C\\d{2}"), ])

enriched_targets$Basename %>% str_extract("GSM\\d+_\\d+_R\\d{2}C\\d{2}")

rownames(enriched_targets) <- enriched_targets$Basename %>% str_extract("GSM\\d+_\\d+_R\\d{2}C\\d{2}")
enriched_targets %>%
    rownames() %>%
    head()

common <- intersect(rownames(enriched_targets), m_values_samples)
length(common)

harmonized_targets <- enriched_targets[rownames(enriched_targets) %in% common, ]
dim(harmonized_targets)

harmonized_m_values <- m_values[, colnames(m_values) %in% common]
dim(harmonized_m_values)

saveRDS(
    harmonized_targets,
    file.path(project_context$paths$processed, "GSE145361_harmonized_targets.rds")
)

saveRDS(
    harmonized_m_values,
    file.path(project_context$paths$processed, "GSE145361_harmonized_m_values.rds")
)
