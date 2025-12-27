library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

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









