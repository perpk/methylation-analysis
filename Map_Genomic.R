## Here we map CpGs to genomic features. Ultimately, we will do DL based on engineered features.

peg1_m_values <- readRDS("/root/workspace/restored/GSE111629_harmonized_m_values_cleaned.rds")
ppmi_m_values <- readRDS("/root/workspace/restored/PPMI_harmonized_m_values_cleaned.rds")
sgpd_m_values <- readRDS("/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/processed/GSE145361_harmonized_m_values_cleaned_corrected.rds")
    
peg1_targets <- readRDS("/root/workspace/restored/peg1_harmonized_targets.rds")
ppmi_targets <- readRDS("/root/workspace/restored/ppmi_harmonized_targets.rds")
sgpd_targets <- readRDS("/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/processed/GSE145361_harmonized_targets.rds")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(tidyr)
library(arrow)

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

mapping <- data.frame(
  Probe_ID = rownames(anno),
  Gene = anno$UCSC_RefGene_Name,
  Region = anno$UCSC_RefGene_Group,
  stringsAsFactors = FALSE
)

mappingEPIC <- data.frame(
  Probe_ID = rownames(annoEPIC),
  Gene = annoEPIC$UCSC_RefGene_Name,
  Region = annoEPIC$UCSC_RefGene_Group,
  stringsAsFactors = FALSE
)

clean_mapping <- mapping %>%
  # Drop probes that do not map to any known gene (intergenic regions)
  filter(Gene != "") %>%
  
  # The magic function: splits the strings and multiplies the rows simultaneously
  # Because Gene and Region have identical counts of semicolons, they stay aligned
  separate_rows(Gene, Region, sep = ";") %>%
  
  # Remove duplicate mappings caused by multiple transcript variants
  distinct(Probe_ID, Gene, Region) %>%
  
  # Create the final unified feature name for the neural network
  # e.g., "SNCA_TSS200"
  mutate(Feature_Name = paste(Gene, Region, sep = "_")) %>%
  
  # Keep only what we need for cuDF
  select(Probe_ID, Feature_Name)
    
clean_mapping %>% head()
write_feather(clean_mapping, "/workspace/results/450_probe_to_feature_mapping.feather")
clean_mapping_epic %>% head()
write_feather(clean_mapping_epic, "/workspace/results/EPIC_probe_to_feature_mapping.feather")
