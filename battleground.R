project_to_load <- "GSE111629_20251226_102044"
project_location <- "/Volumes/Elements/methyl-pipe-out"
project_name <- "GSE111629"

source("R/project_context.R")

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)
context <- .load_methylation_project(project_location, project_to_load, platform = "450k", cohorts = cohorts)

rg_set <- readRDS(file.path(context$paths$raw_data, "rg_set.rds"))

library(minfi)
library(minfiData)
library(IlluminaHumanMethylation450kmanifest)

# qc_data <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION I")
# qcReport(rg_set, pdf = "GSE111629_Bisulfite_Conversion_QC_Report.pdf")

# 1. Get addresses for both types
ctrls_i <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION I")
ctrls_ii <- getControlAddress(rg_set, controlType = "BISULFITE CONVERSION II")

# 2. Calculate average intensities
# BC I usually uses Green for the 'Converted' signal
bc1_scores <- colMeans(getGreen(rg_set)[ctrls_i, ])

# BC II usually uses Red for the 'Converted' signal
bc2_scores <- colMeans(getRed(rg_set)[ctrls_ii, ])

# 3. Identify failures (using your 2 SD logic)
fail_i <- names(bc1_scores[bc1_scores < (mean(bc1_scores) - 2 * sd(bc1_scores))])
fail_ii <- names(bc2_scores[bc2_scores < (mean(bc2_scores) - 2 * sd(bc2_scores))])

# 4. Combine into a unique list of problematic samples
all_failed_samples <- unique(c(fail_i, fail_ii))

print(all_failed_samples)
