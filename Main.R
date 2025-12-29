library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ggplot2)

load_project <- TRUE
project_to_load <- "GSE111629_20251226_102044"

source('R/progress_mgr.R')
source('R/project_context.R')
if (load_project == FALSE) {
  #### Create a new project
  project_context <- create_methylation_project("GSE111629", "/Volumes/Elements/methyl-pipe-out")
} else {
  #### Load an existing project
  project_context <- .load_methylation_project("/Volumes/Elements/methyl-pipe-out", project_to_load)
}

### Read array data from filesystem

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

### Read raw data and write to disk
source('R/extract_methyl_set.R')
extract_methyl_set(context = project_context, targets = targets)

### QC - Outlier detection and removal
source('R/qc.R')
qc(context = project_context, targets = targets, )

### Remove sex-mismatched samples
source('R/biological_gender_mismatch_analysis.R')
biological_gender_mismatch_analysis(context = project_context, recorded_sex_col = "gender:ch1")

### Remove cross-reactive probes, sex-chromosome related probes and single nucleotide polymorphisms (SNPs)
### Order matters, firstly SNPs must be removed then probes on XY chromosomes and thus keeping only those on autosomal and finally filtering of cross-reactive probes.
#### 1. Single Nucleotide Polymorphisms
source('R/remove_snp_450k.R')
remove_snp_450k(context = project_context)

#### 2. Sex-Chromosome Probe filtering
source('R/remove_sex_chromosome_probes.R')
remove_sex_chromosome_probes(context = project_context)

#### 3. Cross Reactive Probes
source('R/remove_cross_reactive_probes_450k.R')
remove_cross_reactive_probes_450k(context = project_context)

### Filter rg_set according to the remaining methyl_set probes after it has been cleared off from SNPs, X/Y-Chromosome- and cross-reactive probes.
source('R/filter_rg_set.R')
filter_rg_set(context = project_context)

### Perform background correction and dye-bias normalization on rg_set and extract new methyl_set & beta-matrix based on the filtered rg_set from previous step
source('R/bg_correction_dye_bias_norm.R')
bg_correction_dye_bias_norm(context = project_context)

### Beta-Mixture Quantile (BMIQ) Normalization for 450k
source('R/apply_BMIQ.R')
apply_BMIQ(context = project_context)


### Principal Component Analysis
source('R/principal_component_analysis.R')
col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- "gender:ch1"
col_map[["Age"]] <- "age:ch1"
keys <- c("Sample_Group", "Gender", "Age")
principal_component_analysis(project_context, col_maps = col_map, keys = keys, npc = 5)

#### Plot PCA
source('R/plot_PCA.R')
pca_vars <- c("Sample_Group" = "By Diagnosis", "Gender" = "By Biological Gender", "Age" = "By Age")
convert_f <- function(df) {
  (transform(df, Age = as.numeric(Age)))
}
plot_PCA(
  context = project_context,
  pca_results_rds_filename = "pca_df.rds",
  pca_vars = pca_vars,
  convert_fun = convert_f,
  continuously_scaled = c("Age"),
  "pca_plot_",
  pairplot_color_by="Age"
)

#### Outlier detection from PCA
source('R/outlier_analysis.R')
outlier_analysis(context=project_context, sample_metadata=c("Sample_Group", "gender:ch1", "age:ch1"))

#### PCA Outlier assessment (TODO Doesn't really work out)
source('R/pca_outlier_assessment.R')
pca_outlier_assessment(context = project_context)
metrics <- readRDS(file.path(project_context$paths$qc, "metrics_pca_outlier_assessment.rds"))

### Calculate Delta-Beta values
beta_matrix <- readRDS(file.path(project_context$paths$results, "beta_matrix_bmiq.rds"))
targets <- readRDS(file.path(project_context$paths$processed, "targets_remove_mismatch.rds"))
rownames(targets) <- targets$Basename %>% str_remove("GSE111629_RAW/")

pd_samples <- rownames(targets[targets$Sample_Group == 'PD',])
hc_samples <- rownames(targets[targets$Sample_Group == 'Control',])

colnames_pd <- which(colnames(beta_matrix) %in% pd_samples)
colnames_hc <- which(colnames(beta_matrix) %in% hc_samples)

mean_beta_pd <- rowMeans(beta_matrix[, colnames_pd], na.rm = TRUE)
mean_beta_hc <- rowMeans(beta_matrix[, colnames_hc], na.rm = TRUE)
delta_beta <- mean_beta_pd - mean_beta_hc

beta_means <- data.frame(
  mean_beta_pd,
  mean_beta_hc,
  delta_beta
)

write.csv(beta_means, file.path(project_context$paths$results, "beta_means.csv"))

rm(beta_matrix)
rm(targets)
gc(full = TRUE)

### Cell Count Estimation
source('R/cell_cnt_estimate.R')
cell_cnt_estimate(context = project_context)

### Plot Cell proportion per cohort and write results to files for each cell type
source('R/plot_cell_proportions.R')
plot_cell_proportions(context = project_context)
