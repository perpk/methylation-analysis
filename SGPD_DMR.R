rm(list = ls())
gc(full = TRUE)

library(tidyverse)

data_dir <- "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/"

m_values_bmiq_no_outliers <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_m_values.rds"))
targets <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_targets.rds"))

dmr <- function(m_values, targets, design) {
    library(DMRcate)

    annot <- cpg.annotate(
        datatype = "array",
        object = m_values,
        what = "M",
        arraytype = "450K", # "450K"
        analysis.type = "differential",
        design = design,
        coef = 2, # The column index in 'design' for DiagnosisPD
        fdr = 1 # FDR threshold for individual probes (tuning parameter)
    )

    dmrc_output <- dmrcate(
        annot,
        lambda = 1000, # Bandwidth in base pairs (1000 is standard for arrays)
        C = 2, # Statistical scaling factor (2 is recommended for arrays)
        min.cpgs = 3, # A region must have at least 3 CpGs to be considered a DMR
        pcutoff = 0.005 # FDR threshold for DMRs (tuning parameter)
    )
    results_ranges <- extractRanges(dmrc_output)
    results_df <- as.data.frame(results_ranges)
    results_df <- results_df[order(results_df$min_smoothed_fdr), ]

    (results_df)
}

targets$Slide <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets$Array <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])
targets$Sex <- targets$`gender:ch1`

design <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Array, data = targets)
dmr_results <- dmr(m_values_bmiq_no_outliers, targets, design)
