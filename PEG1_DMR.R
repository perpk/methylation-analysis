rm(list = ls())
gc(full = TRUE)

targets_peg1 <- readRDS("/root/workspace/data/peg1_harmonized_targets.rds")
peg1_pca_df_with_outliers <- readRDS("/root/workspace/data/peg1_pca_df_with_outliers.rds")
peg1_m_values <- readRDS(file.path("/root/workspace/methyl-pipe-out", "peg1_m_values_bmiq_no_outliers.rds"))

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

targets_peg1$PC2 <- peg1_pca_df_with_outliers[match(rownames(targets_peg1), rownames(peg1_pca_df_with_outliers)), "PC2"]
targets_peg1$PC4 <- peg1_pca_df_with_outliers[match(rownames(targets_peg1), rownames(peg1_pca_df_with_outliers)), "PC4"]

targets_peg1$Slide <- targets_peg1$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets_peg1$Array <- targets_peg1$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])

targets_peg1 %>% head()

targets_peg1$Age <- targets_peg1$`age:ch1`

design_pca <- model.matrix(~ Sample_Group + Age + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + ScanDate + Array + PC2 + PC4, data = targets_peg1)
dmr_results <- dmr(peg1_m_values, targets_peg1, design_pca)

dmr_results$results_ranges %>%
    as.data.frame() %>%
    View()

dmr_df <- as.data.frame(dmr_results$results_ranges)
write.csv(dmr_df, file = "/root/workspace/data/peg1_dmr_results.csv", row.names = FALSE)
saveRDS(dmr_results$results_ranges, file = "/root/workspace/data/GSE111629_DMR_results_ranges.rds")

pal <- c("blue", "red")

sample_colors <- pal[as.factor(targets_peg1$Sample_Group)]
names(sample_colors) <- rownames(targets_peg1)

library(lumi)
beta <- m2beta(peg1_m_values)

pd_samples <- targets_peg1[targets_peg1$Sample_Group == "PD", ]
control_samples <- targets_peg1[targets_peg1$Sample_Group == "Control", ]

pd_samples$Sample_Name <- rownames(pd_samples)
control_samples$Sample_Name <- rownames(control_samples)

pd_means <- rowMeans(beta[, pd_samples$Sample_Name], na.rm = TRUE)
control_means <- rowMeans(beta[, control_samples$Sample_Name], na.rm = TRUE)

beta_aggregated <- cbind(PD = pd_means, Control = control_means)

dmr_ranges <- dmr_results$results_ranges

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_gr <- GRanges(
    seqnames = ann$chr,
    ranges = IRanges(start = ann$pos, end = ann$pos),
    strand = ann$strand,
    probe_id = rownames(ann)
)

top_dmr <- dmr_ranges[1]
overlaps <- findOverlaps(ann_gr, top_dmr)
dmr_probes <- ann_gr$probe_id[queryHits(overlaps)]

library(DMRcate)
library(tidyverse)
calculate_effect_size_all_dmrs <- function(dmr_ranges, beta, targets) {
    dmr_ranges_df <- as.data.frame(dmr_ranges)

    dmr_effect_sizes <- lapply(seq_along(dmr_ranges), function(index) {
        valid_dmr_probes <- dmr_probes(dmr_ranges[index], beta)

        if (length(valid_dmr_probes) == 0) {
            delta_beta <- NA_real_
        } else {
            dmr_beta_subset <- beta[valid_dmr_probes, , drop = FALSE]

            pd_samples <- rownames(targets[targets$Sample_Group == "PD", , drop = FALSE])
            control_samples <- rownames(targets[targets$Sample_Group == "Control", , drop = FALSE])

            mean_pd <- mean(dmr_beta_subset[, pd_samples, drop = FALSE], na.rm = TRUE)
            mean_control <- mean(dmr_beta_subset[, control_samples, drop = FALSE], na.rm = TRUE)
            delta_beta <- mean_pd - mean_control
        }

        data.frame(
            dmr_index = index,
            deltabeta = delta_beta,
            FDR = dmr_ranges_df$min_smoothed_fdr[index],
            genes = dmr_ranges_df$overlapping.genes[index],
            chromosome = as.character(seqnames(dmr_ranges[index])),
            stringsAsFactors = FALSE
        )
    })

    bind_rows(dmr_effect_sizes)
}

dmr_effect_sizes_df <- calculate_effect_size_all_dmrs(dmr_ranges, beta, targets_peg1)
