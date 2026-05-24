targets <- readRDS("/Volumes/Elements/vastai/gse111629/GSE111629_20260515_083253/processed/targets_remove_mismatch.rds")
targets_cells <- readRDS("/Volumes/Elements/vastai/gse111629/GSE111629_20260416_183131/targets_s_mismatch_cells.rds")

dim(targets)
dim(targets_cells)

head(rownames(targets))
head(rownames(targets_cells))
head(targets$Sample_Name)

library(dplyr)
library(stringr)

targets_cells %>% rownames() %>% str_detect("\\d{12}_R\\d{2}C\\d{2}") %>% sum() - dim(targets_cells)[1]targets %>% pull(Sample_Name) %>% str_detect("\\d{4}_\\d{12}_R\\d{2}C\\d{2}") %>% sum() - dim(targets)[1]
targets %>% pull(Sample_Name) %>% str_detect("GSM") %>% sum() - dim(targets)[1]

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


merged[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]  %>% is.na() %>% sum()

merged$Sentrix_ID <- sapply(merged$Basename, function(x) {
    basename(x) %>% str_extract("\\d{10}_R\\d{2}C\\d{2}")
})
merged$Sentrix_ID %>% is.na() %>% sum()

idat_folder_loc <- "/Volumes/Elements/methylation-analysis/GSE111629_RAW"
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
dim(enriched_targets)
dim(merged)

m_values <- readRDS("/Volumes/Elements/vastai/gse111629/GSE111629_20260515_083253/results/m_values_bmiq.rds")
dim(m_values)

m_values_samples <- colnames(m_values)

rownames(enriched_targets) <- enriched_targets$Basename %>% str_extract("GSM\\d{7}_\\d{10}_R\\d{2}C\\d{2}")
enriched_targets %>% rownames() %>% head()

common <- intersect(rownames(enriched_targets), m_values_samples)
length(common)

harmonized_targets <- enriched_targets[rownames(enriched_targets) %in% common,  ]
dim(harmonized_targets)

harmonized_m_values <- m_values[, colnames(m_values) %in% common]
dim(harmonized_m_values)

saveRDS(
    harmonized_targets,
    "/Volumes/Elements/vastai/gse111629/GSE111629_20260515_083253/processed/GSE111629_harmonized_targets.rds"
)

saveRDS(
    harmonized_m_values,
    "/Volumes/Elements/vastai/gse111629/GSE111629_20260515_083253/processed/GSE111629_harmonized_m_values.rds"
)
