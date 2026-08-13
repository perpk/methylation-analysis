ppmi_m_values_bmiq_no_outliers <- readRDS(file.path("/root/workspace/methyl-pipe-out", "ppmi_m_values_bmiq_no_outliers.rds"))
targets_ppmi <- readRDS("/root/workspace/data/ppmi_harmonized_targets.rds")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

dmr <- function(m_values, targets, design) {
    library(DMRcate)

    annot <- cpg.annotate(
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
targets_original_ppmi <- readRDS("/root/workspace/data/targets_original_ppmi.rds")

targets_test <- merge(
    targets_ppmi,
    targets_original_ppmi %>% select(Sample_Name, ENROLL_AGE),
    by = "Sample_Name",
)
targets_test %>% head()

design <- model.matrix(~ Sample_Group + ENROLL_AGE.y + Sex + CD8T + CD4T + Bcell + Mono + NK + Neu + ScanDate + Array, data = targets_test)

dmr_results <- dmr(ppmi_m_values_bmiq_no_outliers, targets_test, design)
dmr_results$results_ranges %>%
    as.data.frame() %>%
    View()

saveRDS(dmr_results$results_ranges, "/root/workspace/data/PPMI_DMR_results_ranges.rds")
saveRDS(dmr_results$results_df, "/root/workspace/data/PPMI_DMR_results_df.rds")
write.csv(dmr_results$results_df, "/root/workspace/data/PPMI_DMR_results_df.csv", row.names = FALSE)

library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann_gr <- GRanges(
    seqnames = ann$chr,
    ranges = IRanges(start = ann$pos, end = ann$pos),
    strand = ann$strand,
    probe_id = rownames(ann)
)

top_dmr <- dmr_results$results_ranges[1]
overlaps <- findOverlaps(ann_gr, top_dmr)
top_dmr_probes <- ann_gr$probe_id[queryHits(overlaps)]

get_dmr_probes <- function(dmr_range, beta) {
    overlaps <- findOverlaps(ann_gr, dmr_range)
    probe_ids <- ann_gr$probe_id[queryHits(overlaps)]
    intersect(rownames(beta), probe_ids)
}

library(lumi)
beta <- m2beta(ppmi_m_values_bmiq_no_outliers)


library(DMRcate)
library(tidyverse)
calculate_effect_size_all_dmrs <- function(dmr_ranges, beta, targets) {
    dmr_ranges_df <- as.data.frame(dmr_ranges)

    dmr_effect_sizes <- lapply(seq_along(dmr_ranges), function(index) {
        valid_dmr_probes <- get_dmr_probes(dmr_ranges[index], beta)

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

dmr_effect_sizes <- calculate_effect_size_all_dmrs(dmr_results$results_ranges, beta, targets_ppmi)
dmr_effect_sizes %>% head()

write.csv(dmr_effect_sizes, "/root/workspace/data/PPMI_DMR_effect_sizes.csv", row.names = FALSE)

pd_samples <- targets_ppmi[targets_ppmi$Sample_Group == "PD", ]
control_samples <- targets_ppmi[targets_ppmi$Sample_Group == "Control", ]

pd_samples$Sample_Name <- rownames(pd_samples)
control_samples$Sample_Name <- rownames(control_samples)

pd_means <- rowMeans(beta[, pd_samples$Sample_Name], na.rm = TRUE)
control_means <- rowMeans(beta[, control_samples$Sample_Name], na.rm = TRUE)

beta_aggregated <- cbind(PD = pd_means, Control = control_means)


# 14, 15, 7, 1; 16, 10, 5, 11

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR14.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 14,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR15.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 15,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR7.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 7,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR1.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 1,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

# 14, 15, 7, 1; 16, 10, 5, 11

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR16.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 16,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR10.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 10,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR5.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 5,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()

png(file = file.path("/root/workspace/data", "PPMI_Top_DMR11.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_results$results_ranges,
    dmr = 11,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "EPICv1",
    genome = "hg19"
)

dev.off()
