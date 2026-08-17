rm(list = ls())
gc(full = TRUE)

library(tidyverse)
library(minfi)

m_values <- readRDS(file.path("/root/workspace/methyl-pipe-out", "peg1_m_values_bmiq_no_outliers.rds"))
peg1_pca_df_with_outliers <- readRDS("/root/workspace/data/peg1_pca_df_with_outliers.rds")
targets <- readRDS("/root/workspace/data/peg1_harmonized_targets.rds")

targets$PC2 <- peg1_pca_df_with_outliers[match(rownames(targets), rownames(peg1_pca_df_with_outliers)), "PC2"]
targets$PC4 <- peg1_pca_df_with_outliers[match(rownames(targets), rownames(peg1_pca_df_with_outliers)), "PC4"]
targets$Age <- peg1_pca_df_with_outliers[match(rownames(targets), rownames(peg1_pca_df_with_outliers)), "Age"]

targets$Slide <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets$Array <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])


library(sva)

mod <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets)
mod0 <- model.matrix(~ Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets)

sva_results <- sva(m_values, mod, mod0)

sv_matrix <- sva_results$sv

colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))

sv_matrix %>% head()

saveRDS(sv_matrix, file.path("/root/workspace/methyl-pipe-out", "sv_matrix_peg1.rds"))

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


design <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets)
design_sva_full <- cbind(design, sv_matrix)

library(limma)
fit <- lmFit(m_values, design_sva_full)
fit <- eBayes(fit)
stats <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
p_values <- stats$P.Value
chisq <- qchisq(1 - p_values, 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

print(paste("Genomic Inflation Factor (Lambda):", round(lambda, 3)))

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
dmr_results <- dmr(m_values, targets, design_sva_full)

dmr_results %>%
    as.data.frame() %>%
    View()

overlaps_prev_run <- c(
    "HKR1",
    "GDF7",
    "CTBP1-AS2",
    "CTBP1",
    "ARHGEF4",
    "HCAR1",
    "H2AFJ",
    "LSP1",
    "HOXA-AS3", "HOXA3", "RP1-170O19.22", "HOXA5",
    "MEIOB",
    "ARHGEF10",
    "MTTP", "TRMT10A",
    "ZNF43",
    "ENKUR",
    "TMEM184C",
    "ID2",
    "GTPBP3",
    "ANK1",
    "PRRT1",
    "ZBTB22",
    "MCF2L",
    "PHACTR2",
    "NOTCH4",
    "SLFN12",
    "DDR1",
    "RP11-328M4.3",
    "RP11-328M4.2",
    "METTL22",
    "SULT1C2",
    "MFSD6L",
    "C22orf34"
)

# Extract and split comma-separated genes from overlapping.genes column
genes_from_dmr <- dmr_results$overlapping.genes %>%
    strsplit(",") %>%
    unlist() %>%
    trimws()

# Find intersection with previous run
common_genes <- intersect(genes_from_dmr, overlaps_prev_run)
print(common_genes)
