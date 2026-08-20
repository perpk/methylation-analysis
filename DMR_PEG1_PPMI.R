rm(list = ls())
gc(full = TRUE)

library(tidyverse)

m_values_combined <- readRDS("/workspace/results/consolidated_m_values.rds")
targets_combined <- readRDS("/workspace/results/consolidated_targets.rds")

samples_to_exclude <- targets_combined %>%
    filter(Cohort == "SGPD") %>%
    rownames()

m_values_peg1_ppmi <- m_values_combined[, !colnames(m_values_combined) %in% samples_to_exclude]
targets_peg1_ppmi <- targets_combined[!rownames(targets_combined) %in% samples_to_exclude, ]

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

    png(file.path("/workspace/results/GSE145361_QQ_plot.png"), width = 800, height = 800)

    plot(-log10(expected_p), -log10(observed_p),
        main = paste("QQ Plot of PD Methylation\nLambda =", round(lambda, 3)),
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        pch = 16, cex = 0.5, col = "darkblue"
    )

    abline(0, 1, col = "red", lwd = 2)
    dev.off()
}



design <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Cohort, data = targets_peg1_ppmi)

qc_dmr(m_values_peg1_ppmi, design)

mod <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Cohort, data = targets_peg1_ppmi)
mod0 <- model.matrix(~ Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Cohort, data = targets_peg1_ppmi)

library(sva)
sva_results <- sva(m_values_peg1_ppmi, mod, mod0)
sva_results$n.sv
sv_matrix <- sva_results$sv
colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))

design_sva_full <- cbind(design, sv_matrix)

qc_dmr(m_values_peg1_ppmi, design_sva_full)

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

    results <- list(
        results_ranges = results_ranges,
        results_df = results_df
    )
    return(results)
}

dmr_results <- dmr(m_values_peg1_ppmi, targets_peg1_ppmi, design_sva_full)
