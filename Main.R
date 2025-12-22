source('R/project_context.R')
source('R/qc.R')
source('R/progress_mgr.R')
source('R/normalize_extract_methylset.R')
source('R/principal_component_analysis.R')
source('R/plot_PCA.R')
source('R/outlier_analysis.R')
source('R/biological_gender_mismatch_analysis.R')
source('R/clean_data.R')

library(GEOquery)
library(tidyverse)
library(minfi)

project_context <- create_methylation_project("GSE111629", "/Volumes/Elements/methyl-pipe-out")

project_context <- .load_methylation_project("/Volumes/Elements/methyl-pipe-out", "GSE111629_20251211_103624")

targets <- read.metharray.sheet("GSE111629_RAW", pattern = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
gse <- getGEO('GSE111629', GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])
all_idat_files <- list.files("GSE111629_RAW", pattern = "\\.idat\\.gz$", full.names = FALSE)
basenames <- unique(gsub("_(Grn|Red).idat.gz", "", all_idat_files))
targets <- data.frame(
  Sample_Name = gsub("_.*", "", basenames),
  Basename = file.path("GSE111629_RAW", basenames),
  stringsAsFactors = FALSE
)
targets <- merge(targets, pdata, by.x = "Sample_Name", by.y = "row.names", all.x = TRUE)
targets$Sample_Group <- factor(targets$`disease state:ch1`)
targets <- targets %>% mutate(Sample_Group = case_when(Sample_Group == "Parkinson's disease (PD)" ~ "PD", Sample_Group == "PD-free control" ~ "Control"))
targets$Sample_Group <- as.factor(targets$Sample_Group)

rm(list=setdiff(ls(), c("project_context", "context", "targets")))
gc(full=T)

qc_res <- qc(targets, project_context)

methyl_set_res <- normalize_extract_methylset(project_context)

col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- "gender:ch1"
col_map[["Age"]] <- "age:ch1"
keys <- c("Sample_Group", "Gender", "Age")

source('R/principal_component_analysis.R')
pca_df <- principal_component_analysis(project_context, col_map, keys, 5)

pca_persisted_filename<-"pca_df_before_cleanup.rds"
pca_vars <- c("Sample_Group" = "By Diagnosis", "Gender" = "By Biological Gender", "Age" = "By Age")
convert_f <- function(df) {
  (transform(df, Age = as.numeric(Age)))
}
source('R/plot_PCA.R')
plot_PCA(
  context = project_context,
  pca_results_rds_filename = pca_persisted_filename,
  pca_vars = pca_vars,
  convert_fun = convert_f,
  continuously_scaled = c("Age"),
  "pca_plot_before_",
  pairplot_color_by="Age"
)

source('R/outlier_analysis.R')
outlier_analysis(context=project_context, pca=pca_df, targets_filename="targets_clean.rds", sample_metadata=c("Sample_Group", "gender:ch1", "age:ch1"))

source('R/biological_gender_mismatch_analysis.R')
mismatches <- biological_gender_mismatch_analysis(context=project_context, recorded_sex_col = "gender:ch1", targets_filename="targets_clean.rds")
mismatches

data_to_clean <- list(
  beta_matrix_data = list(
    filepath = file.path(project_context$paths$processed, "beta_matrix.rds"),
    clean_filename = file.path(project_context$paths$processed, "beta_matrix_s_mismatch.rds")
  ),
  m_values_data = list(
    filepath = file.path(project_context$paths$processed, "m_values.rds"),
    clean_filename = file.path(project_context$paths$processed, "m_values_s_mismatch.rds")
  ),
  rg_set_data = list(
    filepath = file.path(project_context$paths$qc, "rg_set_clean.rds"),
    clean_filename = file.path(project_context$paths$qc, "rg_set_s_mismatch.rds")
  ),
  methyl_set_data = list(
    filepath = file.path(project_context$paths$normalized, "methyl_set.rds"),
    clean_filename = file.path(project_context$paths$normalized, "methyl_set_s_mismatch.rds")
  )
)

if (!exists("mismatches")) {
  sex_mismatch_log <- file.path(project_context$paths$results, "sex_mismatch_log.csv")
  mismatches <- read.csv(sex_mismatch_log, row.names=1)
}

samples_to_remove <- rownames(mismatches)

source('R/clean_data.R')
clean_data(samples_to_remove = samples_to_remove, data_to_clean = data_to_clean)

targets_to_clean <-   list(
  targets_data = list(
    filepath = file.path(project_context$paths$qc, "targets_clean.rds"),
    clean_filename = file.path(project_context$paths$qc, "targets_s_mismatch.rds"),
    clean_rows = TRUE,
    by_col = "Sample_Name"
  )
)

rows_to_remove <- mismatches[['Sample_Name']]
clean_data(samples_to_remove = rows_to_remove, data_to_clean = targets_to_clean)
