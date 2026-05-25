rm(list = ls())
gc(full = TRUE)

targets <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260513_110353/processed/targets_remove_mismatch.rds")
targets_cells <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260415_170143/targets_s_mismatch_cells.rds")

dim(targets)
dim(targets_cells)

head(rownames(targets))
head(rownames(targets_cells))
head(targets$Sample_Name)

head(targets$Sample_Name)

library(dplyr)
library(stringr)

targets_cells %>% rownames() %>% str_detect("\\d{12}_R\\d{2}C\\d{2}") %>% sum() - dim(targets_cells)[1]
dim(targets_cells)
rownames(targets_cells) %>% head()

targets %>% pull(Sample_Name) %>% str_detect("\\d{4}_\\d{12}_R\\d{2}C\\d{2}") %>% sum() - dim(targets)[1]
dim(targets)

rownames(targets) <- targets %>% pull(Sample_Name) %>% str_extract("\\d{12}_R\\d{2}C\\d{2}")
head(rownames(targets))

colnames(targets_cells)

cols <- c("Sample_Name", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
data <- targets_cells %>% select(all_of(cols))

merged <- merge(
    x = targets,
    y = data,
    by.x = "Sample_Name",
    by.y = "Sample_Name",
    all.x = FALSE
)

library(terra)

merged[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]  %>% is.na() %>% sum()
cell_counts[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")] %>% is.na() %>% sum()

intersect(rownames(targets), rownames(targets_cells)) %>% length()

merged[is.na(merged[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]), ]


merged$Sentrix_ID <- sapply(merged$Basename, function(x) {
    basename(x) %>% str_extract("\\d{12}_R\\d{2}C\\d{2}")
})

idat_folder_loc <- "/Volumes/Elements/methylation-analysis/ppmi/Project120_IDATS_n524final_toLONI_030718"
source("R/extract_scandate_from_idat.R")
scan_dates <- extract_scandate_from_idat(
    file_path=idat_folder_loc
)
scan_dates <- scan_dates[!duplicated(scan_dates[c("SentrixID", "ScanDate")]), ]

enriched_targets <- merge(
    x = merged,
    y = scan_dates,
    by.x = "Sentrix_ID",
    by.y = "SentrixID",
    all.x = FALSE
)
enriched_targets %>% pull(ScanDate) %>% is.na() %>% sum()

m_values <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260513_110353/results/m_values_bmiq.rds")
dim(m_values)

m_values_samples <- colnames(m_values)
str(m_values_samples)

rownames(enriched_targets) <-  enriched_targets$Sentrix_ID

common <- intersect(rownames(enriched_targets), m_values_samples)

harmonized_targets <- enriched_targets[rownames(enriched_targets) %in% common,  ]
dim(harmonized_targets)

harmonized_m_values <- m_values[, colnames(m_values) %in% common]
dim(harmonized_m_values)

saveRDS(
    harmonized_targets, 
    "/Volumes/Elements/vastai/ppmi/ppmi_20260513_110353/processed/PPMI_harmonized_targets_local.rds"
)

saveRDS(
    harmonized_m_values, 
    "/Volumes/Elements/vastai/ppmi/ppmi_20260513_110353/processed/PPMI_harmonized_m_values_local.rds"
)

### ____________________

library(illuminaio)

idat <- readIDAT("/Volumes/Elements/methylation-analysis/ppmi/Project120_IDATS_n524final_toLONI_030718/200973410143_R06C01_Red.idat")