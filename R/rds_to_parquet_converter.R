install.packages("arrow")

library(arrow)
library(dplyr)

rootDir <- "/Volumes/Elements/vastai"

projects <- list(
    list(path="gse111629/GSE111629_20260416_183131", harmonize="gse")
)

files <- list(
    list(m_values = "m_values_bmiq.rds", meta = "targets_remove_mismatch.rds")
)
m <- readRDS(file.path(rootDir, projects[[1]]$path, files[[1]]$meta))
meta <- harmonize_meta_m_values_gse(m_vals, m)
rownames(m)
head(gsub(m$Basename, pattern = "GSE111629_RAW/", replacement = ""))

harmonize_meta_m_values_ppmi <- function(m_values, meta) {
    rownames(meta) <- gsub(meta$Sample_Name, pattern = "^\\d{4}_", replacement = "")
    meta <- meta[match(colnames(m_values), rownames(meta)), ]
    return(meta)
}

harmonize_meta_m_values_gse <- function(m_values, meta) {
    rownames(meta) <- gsub(meta$Basename, pattern = "GSE111629_RAW/", replacement = "")
    meta <- meta[match(colnames(m_values), rownames(meta)), ]
    return(meta)
}

for (project in projects) {
    for (file in files) {
        m_vals <- readRDS(file.path(rootDir, project$path, file$m_values))
        meta <- readRDS(file.path(rootDir, project$path, file$meta))
        if (project$harmonize == "ppmi") {
            meta <- harmonize_meta_m_values_ppmi(m_vals, meta)
        }
        if (project$harmonize == "gse") {
            meta <- harmonize_meta_m_values_gse(m_vals, meta)
        }

        if (all(rownames(meta) == colnames(m_vals))) {
            combined <- cbind(meta, t(m_vals))
            write_parquet(combined, file.path(rootDir, project$path, "data.parquet"))
        } else {
            stop("Row names of meta do not match column names of m_values.")
        }
    }
}

head(rownames(meta))
head(colnames(m_vals_t))

dim(meta)
dim(m_vals)

head(colnames(m_vals))
head(colnames(meta))
View(meta)
head(gsub(meta$Sample_Name, pattern = paste0(meta$PATNO, "_"), replacement = ""))

meta %>% 
    mutate(Sample_Name = gsub(Sample_Name, pattern = "^\\d{4}_", replacement = "")) %>%
    head()
