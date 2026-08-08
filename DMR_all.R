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
## wget https://drive.usercontent.google.com/download?id=1yZBVUPqHNHEHsEihT6ThKkFiOQHWIsFL&export=download&authuser=0&confirm=t&uuid=1adbefc9-f93a-4183-aaf8-999764f88169&at=AFYLz4OlqCDuH8dLNdbQFZcINEJ_:1785682138681 -O peg1_pca_df_with_outliers.rds

# SGPD
## TODO see SGPD_Batch_Investigation.R

rm(list = ls())
gc(full = TRUE)

library(tidyverse)
library(minfi)

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

    png(file.path("/root/workspace/data/", "GSE111629_QQ_plot.png"), width = 800, height = 800)

    plot(-log10(expected_p), -log10(observed_p),
        main = paste("QQ Plot of PD Methylation\nLambda =", round(lambda, 3)),
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        pch = 16, cex = 0.5, col = "darkblue"
    )

    abline(0, 1, col = "red", lwd = 2)
    dev.off()
}

plot_dmr <- function(targets, results_ranges, annot) {
    # 1. Assign colors to your phenotypes
    CpGs <- annot # The annotated object from the cpg.annotate() step
    # Ensure this matches the order of samples in your original matrix and targets dataframe
    pal <- c("blue", "red")
    sample_colors <- pal[as.factor(targets$Sample_Group)]

    # 2. Plot the top-ranked DMR (dmr = 1)
    DMR.plot(
        ranges = results_ranges, # The GRanges object extracted from dmrcate()
        dmr = 1, # The index of the DMR to plot (1 = most significant)
        CpGs = annot, # The annotated object from the cpg.annotate() step
        phen.col = sample_colors, # The color assignments for your samples
        what = "Beta", # Plots biological proportions instead of M-values
        arraytype = "450K", # Specify your array type
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

targets_peg1 <- readRDS("/root/workspace/data/peg1_harmonized_targets.rds")
methyl_set_peg1 <- readRDS("/root/workspace/data/peg1_methyl_set_removed_cross_reactive.rds")

m_values_peg1 <- getM(methyl_set_peg1)
beta_values_peg1 <- getBeta(methyl_set_peg1)

# Load the pca df with the outliers.
peg1_pca_df_with_outliers <- readRDS("/root/workspace/data/peg1_pca_df_with_outliers.rds")

outlier_samples <- rownames(peg1_pca_df_with_outliers)[peg1_pca_df_with_outliers$Is_Outlier]
beta_matrix_no_outliers <- beta_values_peg1[, !colnames(beta_values_peg1) %in% outlier_samples]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
platform <- IlluminaHumanMethylation450kanno.ilmn12.hg19

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

# design <- model.matrix(~ Sample_Group + Age_Group + Sex + CD8T + CD4T + NK + Bcell, data = targets_peg1)

m_values_bmiq_no_outliers <- beta2m(beta_bmiq)

design <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets_test)

dmr_results <- dmr(m_values_bmiq_no_outliers, targets_peg1, design)

dmr_results %>% dim()

head(dmr_results)
View(dmr_results)

write.csv(dmr_results, file = "peg1_dmr_results.csv", row.names = TRUE)

plot_dmr(targets_peg1, dmr_results, getAnnotation(platform))

peg1_pca_df_with_outliers %>% head()
targets_peg1 %>% head()

targets_test <- targets_peg1

targets_test %>%
    pull(Array) %>%
    unique()


targets_test$PC1 <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "PC1"]
targets_test$PC2 <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "PC2"]
targets_test$PC3 <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "PC3"]
targets_test$PC4 <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "PC4"]
targets_test$PC5 <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "PC5"]
targets_test$Age <- peg1_pca_df_with_outliers[match(rownames(targets_test), rownames(peg1_pca_df_with_outliers)), "Age"]

targets_test$Slide <- targets_test$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets_test$Array <- targets_test$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])

targets_test$Slide


design_test <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array, data = targets_test)
design_pca <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets_test)
qc_dmr(m_values_bmiq_no_outliers, design_pca)

# Remove batch effects via limma
dd <- model.matrix(~Sample_Group, data = targets_test)
cell_types <- as.matrix(targets_test[, c("CD8T", "CD4T", "Bcell", "Mono", "NK", "Gran", "PC2", "PC4")])
m_matrix_clean <- removeBatchEffect(
    x = m_values_bmiq_no_outliers,
    batch = targets_test$Array,
    covariates = cell_types,
    design = dd
)

saveRDS(
    m_matrix_clean,
    file.path("/root/workspace/data", "GSE111629_harmonized_m_values_cleaned.rds")
)

library(arrow)

all(rownames(targets_test) %in% colnames(m_matrix_clean))

all(colnames(m_matrix_clean) %in% rownames(targets_test))

combined <- cbind(targets_test, t(m_matrix_clean))
dim(combined)
write_parquet(combined, write_statistics = FALSE, use_dictionary = FALSE, file.path("/root/workspace/data", "GSE111629_data_test.parquet"))

library(lumi)
beta_matrix_clean <- m2beta(m_matrix_clean)

saveRDS(
    beta_matrix_clean,
    file.path("/root/workspace/data", "GSE111629_harmonized_beta_values_cleaned.rds")
)



# 1.052 all PCs + Age + Sex

pc_pvalues <- sapply(1:5, function(i) {
    pc_col <- paste0("PC", i)

    # Run a t-test comparing the PC scores of PD vs Control
    test_result <- t.test(targets_test[[pc_col]] ~ targets_test$Sample_Group)

    return(test_result$p.value)
})

png(file.path("/root/workspace/data/", "GSE111629_PC_boxplots.png"), width = 1200, height = 800)
# Set up a 2x3 grid for plotting multiple charts at once
par(mfrow = c(2, 3))

# Plot a boxplot for each of the first 5 PCs
for (i in 1:5) {
    pc_col <- paste0("PC", i)
    boxplot(targets_test[[pc_col]] ~ targets_test$Sample_Group,
        main = paste(pc_col, "vs Sample Group"),
        ylab = "PC Score",
        xlab = "",
        col = c("lightblue", "lightcoral")
    )
}
dev.off()

# Reset the plotting grid
par(mfrow = c(1, 1))

names(pc_pvalues) <- paste0("PC", 1:5)

print(pc_pvalues)

table(targets_test$Sample_Group, targets_test$Array)
table(targets_test$Sample_Group, targets_test$ScanDate)


## SVA if necessary (λ doesn't improve with covariates)
library(sva)

mod <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate, data = targets_test)
mod0 <- model.matrix(~ Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate, data = targets_test)

n.sv <- num.sv(m_values_bmiq_no_outliers, mod, method = "leek")
svobj <- sva(m_values_bmiq_no_outliers, mod, mod0, n.sv = n.sv)

design_sva <- cbind(mod, svobj$sv)
colnames(design_sva) <- make.names(colnames(design_sva))
