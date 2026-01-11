library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

project_to_load <- "GSE111629_20251226_102044"
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
targets$Sample_Group <- factor(targets$`disease state:ch1`)
targets <- targets %>% mutate(Sample_Group = case_when(Sample_Group == "Parkinson's disease (PD)" ~ "PD", Sample_Group == "PD-free control" ~ "Control"))
targets$Sample_Group <- as.factor(targets$Sample_Group)

rm(list = setdiff(ls(), c("targets")))
gc(full = T)

source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = project_name)

cohorts <- list(
  PD_vs_Control = c("PD", "Control")
)

source("./pre_process_eda.R")
pre_process_eda(
  project_to_load = project_to_load,
  targets = targets,
  data_folder = data_folder,
  project_location = "/Volumes/Elements/methyl-pipe-out",
  var_mapping = var_mappings,
  cohorts = cohorts
)
