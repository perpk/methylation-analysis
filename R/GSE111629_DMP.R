rm(list = ls())
gc(full = TRUE)

parent_path <- "/Volumes/saucepan/methylation-project"
sub_path <- "GSE111629_20260722_083339"

m_values <- readRDS(
    file.path(
        paste0(parent_path, "/", sub_path, "/results"),
        "GSE111629_combat_m_values.rds"
    )
)

targets <- readRDS(
    file.path(
        paste0(parent_path, "/", sub_path, "/processed"),
        "GSE111629_harmonized_targets.rds"
    )
)

library(ggplot2)
library(tidyverse)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


is.factor(targets$Sample_Group)
# targets$Sample_Group <- as.factor(targets$Sample_Group)

is.factor(targets$Sex)
targets$Sex <- as.factor(make.names(targets$Sex))

is.factor(targets$Age_Group)
targets$Age_Group <- as.factor(make.names(targets$Age_Group))

is.factor(targets$ScanDate)
targets$ScanDate <- as.factor(make.names(targets$ScanDate))


design <- model.matrix(
    ~ Sample_Group + Sex + Age_Group + ScanDate + CD8T + CD4T + Bcell + Gran + Mono + NK,
    data = targets
)

fit <- lmFit(m_values, design)
cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
design %>% head()
