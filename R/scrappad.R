rm(list = ls())
gc(full = TRUE)

sgpd_targets <- readRDS("/Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/processed/GSE145361_harmonized_targets.rds")

library(tidyverse)
sgpd_targets %>% head()

sgpd_targets %>% colnames()

sgpd_targets$`characteristics_ch1` %>% head()
