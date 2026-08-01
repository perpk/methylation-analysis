## DMR for PPMI, SGPD (GSE145361), PEG1 (GSE111629)
# DMR shall be done on non-combat normalized data since limma is used by DMRcate and pass batch information as covariate in the model.
# Also, the original pre_process_eda.R script had a flaw where after outlier removal the BMIQ normalized beta matrix was re-normalized,
# hence the part of outlier removal and normalization is performed here again to get the correct beta matrix and conclusively the
# correct M-values for DMR analysis.

# I need the methyl sets after removal of cross-reactive probes for each cohort and their targets with all cell information and
# scan dates
# The methyl set for SGPD (GSE145361) is missing though ¯\_(ツ)_/¯

# PPMI
## wget https://drive.usercontent.google.com/download?id=107zGQw0RH24365Z4WPPMdFATrc4Cr4dZ&export=download&authuser=0&confirm=t&uuid=af61a6c5-47b6-43a2-899c-ec2b50f5604f&at=AFYLz4M7uoIsJDIKyLenfcTXmpOC:1785576419974 -O ppmi_methyl_set_removed_cross_reactive.rds
## wget https://drive.usercontent.google.com/download?id=1OsaZ_Pqu0RXZ3H2VGsZoTRj87Io5rv1v&export=download&authuser=0&confirm=t&uuid=a1b4d8ac-bbe5-4b5d-88ed-d14f6ddb26ae&at=AFYLz4MqvQESEao_7nEDGsWxs2ad:1785576747081 -O ppmi_harmonized_targets.rds

# PEG1
## wget https://drive.usercontent.google.com/download?id=1vHhRtama9o2xSGde8bCmHWRq_spy93yc&export=download&authuser=0&confirm=t&uuid=959c754e-707a-4948-b864-ab5ca0342c4d&at=AFYLz4PM9Ezut4c1xYnKbD4CcQ1U:1785576540968 -O peg1_methyl_set_removed_cross_reactive.rds
## wget https://drive.usercontent.google.com/download?id=16jH7E6Io9tb3FqymrDfbsv-8_-53OAQs&export=download&authuser=0&confirm=t&uuid=25b11ad7-14d5-443a-858d-cec63d81f836&at=AFYLz4P2Pvmm3Hu8oRG25dI1UYDs:1785576678453 -O peg1_harmonized_targets.rds

# SGPD
## TODO

library(tidyverse)

targets_peg1 <- readRDS("~/Downloads/GSE111629_harmonized_targets.rds")

targets_peg1 %>% head()

targets_ppmi <- readRDS("~/Downloads/ppmi_harmonized_targets.rds")

targets_ppmi %>% head()

table(targets_ppmi$Sample_Group, targets_ppmi$Age_Group)

sum(is.na(targets_ppmi$Age_Group))


### DMRcate
targets$Sample_Group <- factor(targets$Sample_Group, levels = c("Control", "PD"))
targets$Sex <- factor(targets$Sex)
targets$Batch <- factor(targets$Batch) # e.g., Array slide or scan date


design <- model.matrix(~ Sample_Group + Age + Sex + CD8T + Mono + Batch, data = targets)

colnames(design)

# Probe-level annotation
ppmi_annot <- cpg.annotate(
    datatype = "array",
    object = ppmi_m_values,
    what = "M",
    arraytype = "EPIC", # "450K"
    analysis.type = "differential",
    design = design,
    coef = 2, # The column index in 'design' for DiagnosisPD
    fdr = 0.05 # FDR threshold for individual probes (tuning parameter)
)

# Find DMRs
dmrc_output <- dmrcate(
    ppmi_annot,
    lambda = 1000, # Bandwidth in base pairs (1000 is standard for arrays)
    C = 2, # Statistical scaling factor (2 is recommended for arrays)
    min.cpgs = 3 # A region must have at least 3 CpGs to be considered a DMR
)

results_ranges <- extractRanges(dmrc_output)

results_df <- as.data.frame(results_ranges)

results_df <- results_df[order(results_df$min_smoothed_fdr), ]

head(results_df)

# 1. Assign colors to your phenotypes
# We create a vector of colors corresponding to each sample's diagnosis
# Ensure this matches the order of samples in your original matrix and targets dataframe
pal <- c("blue", "red")
sample_colors <- pal[as.factor(targets$Diagnosis)]

# 2. Plot the top-ranked DMR (dmr = 1)
DMR.plot(
    ranges = results_ranges, # The GRanges object extracted from dmrcate()
    dmr = 1, # The index of the DMR to plot (1 = most significant)
    CpGs = myannotation, # The annotated object from the cpg.annotate() step
    phen.col = sample_colors, # The color assignments for your samples
    what = "Beta", # Plots biological proportions instead of M-values
    arraytype = "EPIC", # Specify your array type
    genome = "hg19"
) # The reference genome your data was aligned to
