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
## wget https://drive.usercontent.google.com/download?id=1PIaVMFefCEJxiW2dVaCFoMxRLFkBQCWE&export=download&authuser=0&confirm=t&uuid=93ebaf6d-e67e-4a09-a4ae-18571de2bcf1&at=AFYLz4OcFQMvDnTb0MgWknT8I61X:1786282868617 -O ppmi_pca_df_with_outliers.rds

# PEG1
## wget https://drive.usercontent.google.com/download?id=1vHhRtama9o2xSGde8bCmHWRq_spy93yc&export=download&authuser=0&confirm=t&uuid=959c754e-707a-4948-b864-ab5ca0342c4d&at=AFYLz4PM9Ezut4c1xYnKbD4CcQ1U:1785576540968 -O peg1_methyl_set_removed_cross_reactive.rds
## wget https://drive.usercontent.google.com/download?id=16jH7E6Io9tb3FqymrDfbsv-8_-53OAQs&export=download&authuser=0&confirm=t&uuid=25b11ad7-14d5-443a-858d-cec63d81f836&at=AFYLz4P2Pvmm3Hu8oRG25dI1UYDs:1785576678453 -O peg1_harmonized_targets.rds

# SGPD
## TODO

rm(list = ls())
gc(full = TRUE)

library(tidyverse)

dmr <- function(m_values, targets, design) {
    library(DMRcate)

    ppmi_annot <- cpg.annotate(
        datatype = "array",
        object = m_values,
        what = "M",
        arraytype = "EPICv1", # "450K"
        analysis.type = "differential",
        design = design,
        coef = 2, # The column index in 'design' for DiagnosisPD
        fdr = 1 # FDR threshold for individual probes (tuning parameter)
    )

    dmrc_output <- dmrcate(
        ppmi_annot,
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

qc_dmr <- function(m_values, design) {
    library(limma)

    # until Bcell λ=1.007
    # with Bcell and scandate at the end λ=1.055

    fit <- lmFit(m_values, design)
    fit <- eBayes(fit)
    stats <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    p_values <- stats$P.Value
    chisq <- qchisq(1 - p_values, 1)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

    print(paste("Genomic Inflation Factor (Lambda):", round(lambda, 3)))

    expected_p <- ppoints(length(p_values))
    observed_p <- sort(p_values)

    png(file.path("/root/workspace/data/", "PPMI_QQ_plot.png"), width = 800, height = 800)

    plot(-log10(expected_p), -log10(observed_p),
        main = paste("QQ Plot of PD Methylation\nLambda =", round(lambda, 3)),
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        pch = 16, cex = 0.5, col = "darkblue"
    )

    abline(0, 1, col = "red", lwd = 2)
    dev.off()
}

plot_dmr <- function(targets, results_ranges, ppmi_annot) {
    # 1. Assign colors to your phenotypes
    # We create a vector of colors corresponding to each sample's diagnosis
    # Ensure this matches the order of samples in your original matrix and targets dataframe
    pal <- c("blue", "red")
    sample_colors <- pal[as.factor(targets$Sample_Group)]

    # 2. Plot the top-ranked DMR (dmr = 1)
    DMR.plot(
        ranges = results_ranges, # The GRanges object extracted from dmrcate()
        dmr = 1, # The index of the DMR to plot (1 = most significant)
        CpGs = ppmi_annot, # The annotated object from the cpg.annotate() step
        phen.col = sample_colors, # The color assignments for your samples
        what = "Beta", # Plots biological proportions instead of M-values
        arraytype = "EPICv1", # Specify your array type
        genome = "hg19"
    ) # The reference genome your data was aligned to
}

.get_probe_design_vector <- function(probe_names, platform) {
    # Returns: 1 for Type I probes, 2 for Type II probes
    anno <- getAnnotation(platform)

    # Match probes to annotation (handle missing probes)
    matched_probes <- intersect(probe_names, rownames(anno))

    if (length(matched_probes) < length(probe_names)) {
        warning(paste(
            length(probe_names) - length(matched_probes),
            "probes not found in annotation"
        ))
    }

    # Get probe types
    probe_types <- rep(NA, length(probe_names))
    names(probe_types) <- probe_names

    # Type I = 1, Type II = 2 (BMIQ convention)
    probe_types[matched_probes] <- ifelse(anno[matched_probes, "Type"] == "I", 1, 2)

    # Check for NAs
    if (any(is.na(probe_types))) {
        warning("Some probes could not be classified. Using alternative method...")

        # Alternative: Use probe name patterns
        # Type I: typically start with cg and have specific Infinium I design
        # Type II: typically start with ch or some cg with Infinium II design
        # This is less reliable but works as fallback
        na_probes <- probe_names[is.na(probe_types)]
        probe_types[na_probes] <- ifelse(grepl("^ch", na_probes), 2, 1)
    }

    return(probe_types)
}

## Raw m-values from methylset (after removal of cross-reactive probes) yield with the design matrix:
### ~ Sample_Group + Age_Group + Sex + CD8T + CD4T + NK + Bcell a λ = 1.007 and analysis yields 20 DMRs
### The ComBat dataset yields around 0.7 and is considered genetically deflated. The same set is also thought to be over-normalized (BMIQ accidentally ran twice).

targets_ppmi <- readRDS("/root/workspace/data/ppmi_harmonized_targets.rds")
methyl_set_ppmi <- readRDS("/root/workspace/data/ppmi_methyl_set_removed_cross_reactive.rds")

library(minfi)
m_values_ppmi <- getM(methyl_set_ppmi)
beta_values_ppmi <- getBeta(methyl_set_ppmi)

# Load the pca df with the outliers.
ppmi_pca_df_with_outliers <- readRDS("/root/workspace/data/ppmi_pca_df_with_outliers.rds")

outlier_samples <- rownames(ppmi_pca_df_with_outliers)[ppmi_pca_df_with_outliers$Is_Outlier]
beta_matrix_no_outliers <- beta_values_ppmi[, !colnames(beta_values_ppmi) %in% outlier_samples]

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
platform <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19

design.v <- .get_probe_design_vector(rownames(beta_matrix_no_outliers), platform)
beta_bmiq <- matrix(NA, nrow = nrow(beta_matrix_no_outliers), ncol = ncol(beta_matrix_no_outliers))

library(wateRmelon)
for (i in 1:ncol(beta_matrix_no_outliers)) {
    print(paste("Processing sample", i, "of", ncol(beta_matrix_no_outliers)))
    beta.v <- beta_matrix_no_outliers[, i]

    tryCatch(
        {
            bmiq_result <- BMIQ(
                beta.v = beta.v,
                design.v = design.v,
                nfit = 10000,
                plots = FALSE,
                pri = TRUE
            )
            print(bmiq_result$nbeta[1:10])

            # Store the normalized values
            beta_bmiq[, i] <- bmiq_result$nbeta
        },
        error = function(e) {
            error_message <- paste("BMIQ failed for sample", colnames(beta_matrix_no_outliers)[i], ":", e$message)
            print(error_message)
            warning(error_message)
            beta_bmiq[, i] <- beta.v # Keep original values if BMIQ fails
        }
    )
}

library(lumi)
rownames(beta_bmiq) <- rownames(beta_matrix_no_outliers)
colnames(beta_bmiq) <- colnames(beta_matrix_no_outliers)
m_values_bmiq_no_outliers <- beta2m(beta_bmiq)

### Batch analysis
library(tidyverse)

targets_ppmi %>% head()

targets_ppmi$Slide <- targets_ppmi$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets_ppmi$Array <- targets_ppmi$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])

targets_ppmi %>% head()

targets_original_ppmi <- readRDS("/root/workspace/data/targets_original_ppmi.rds")
targets_original_ppmi %>% head()

targets_test <- merge(
    targets_ppmi,
    targets_original_ppmi %>% select(Sample_Name, ENROLL_AGE),
    by = "Sample_Name",
)
targets_test %>% head()

m_values_bmiq_no_outliers <- readRDS("/root/workspace/methyl-pipe-out/ppmi_m_values_bmiq_no_outliers.rds")
# design <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate + Array, data = targets_ppmi)
design <- model.matrix(~ Sample_Group + ENROLL_AGE.y + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate + Array, data = targets_test)
qc_dmr(m_values_bmiq_no_outliers, design)

### remove batch effects via limma

table(targets_test$Array, targets_test$ScanDate)

dd <- model.matrix(~Sample_Group, data = targets_test)
ppmi_covariates <- as.matrix(targets_test[, c("CD8T", "CD4T", "Bcell", "Mono", "NK", "Neu")])
sex_numeric <- as.numeric(as.factor(targets_test$Sex)) - 1
age_numeric <- as.numeric(targets_test$ENROLL_AGE.y)
ppmi_covariates <- cbind(ppmi_covariates, Sex = sex_numeric)
ppmi_covariates <- cbind(ppmi_covariates, Age = age_numeric)

m_matrix_clean <- removeBatchEffect(
    x = m_values_bmiq_no_outliers,
    batch = targets_test$ScanDate,
    batch2 = targets_test$Array,
    covariates = ppmi_covariates,
    design = dd
)

saveRDS(
    m_matrix_clean,
    file.path("/root/workspace/data", "PPMI_harmonized_m_values_cleaned.rds")
)

library(arrow)

targets_test %>% dim()
m_matrix_clean %>% dim()

rownames(targets_test) <- colnames(m_matrix_clean)

all(rownames(targets_test) %in% colnames(m_matrix_clean))

all(colnames(m_matrix_clean) %in% rownames(targets_test))

combined <- cbind(targets_test, t(m_matrix_clean))
dim(combined)
write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path("/root/workspace/methyl-pipe-out", "ppmi_data_test.parquet"))

library(lumi)
beta_matrix_clean <- m2beta(m_matrix_clean)

saveRDS(
    beta_matrix_clean,
    file.path("/root/workspace/methyl-pipe-out", "ppmi_harmonized_beta_values_cleaned.rds")
)

saveRDS(
    m_values_bmiq_no_outliers,
    file.path("/root/workspace/methyl-pipe-out", "ppmi_m_values_bmiq_no_outliers.rds")
)

### end batch investigation
m_values_bmiq_no_outliers <- beta2m(beta_bmiq)

dmr_results <- dmr(m_values_bmiq_no_outliers, targets_ppmi, design)

dmr_results %>% dim()

plot_dmr(targets_ppmi, results_ranges = dmr_results, ppmi_annot = cpg.annotate(
    datatype = "array",
    object = m_values_bmiq_no_outliers,
    what = "M",
    arraytype = "EPICv1", # "450K"
    analysis.type = "differential",
    design = design,
    coef = 2, # The column index in 'design' for DiagnosisPD
    fdr = 1 # FDR threshold for individual probes (tuning parameter)
))

head(dmr_results)

write.csv(dmr_results, file = "/root/workspace/data/ppmi_dmr_results.csv", row.names = TRUE)

ppmi_pca_df_with_outliers %>% head()
targets_ppmi %>% head()

targets_test <- targets_ppmi

targets_test %>%
    pull(Array) %>%
    unique()


targets_test$PC1 <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "PC1"]
targets_test$PC2 <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "PC2"]
targets_test$PC3 <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "PC3"]
targets_test$PC4 <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "PC4"]
targets_test$PC5 <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "PC5"]
targets_test$Age <- ppmi_pca_df_with_outliers[match(rownames(targets_test), rownames(ppmi_pca_df_with_outliers)), "Age"]

design_test <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate + Array, data = targets_test)
qc_dmr(m_values_bmiq_no_outliers, design_test)

library(sva)

mod <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate, data = targets_test)
mod0 <- model.matrix(~ Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate, data = targets_test)

n.sv <- num.sv(m_values_bmiq_no_outliers, mod, method = "leek")
svobj <- sva(m_values_bmiq_no_outliers, mod, mod0, n.sv = n.sv)

design_sva <- cbind(mod, svobj$sv)
colnames(design_sva) <- make.names(colnames(design_sva))
