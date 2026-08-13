rm(list = ls())
gc(full = TRUE)

library(tidyverse)

data_dir <- "/root/workspace/methyl-pipe-out/GSE145361_20260804_162813/"

m_values_bmiq_no_outliers <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_m_values.rds"))
targets <- readRDS(file.path(data_dir, "processed/GSE145361_harmonized_targets.rds"))
dmr_results_ranges <- readRDS(file.path(data_dir, "results/GSE145361_DMR_results_ranges.rds"))

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

targets$Slide <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[1])

targets$Array <- targets$Sentrix_ID %>%
    strsplit(split = "_") %>%
    sapply(function(x) x[2])
targets$Sex <- targets$`gender:ch1`

design <- model.matrix(~ Sample_Group + Sex + CD8T + CD4T + Bcell + Mono + NK + Gran + Array, data = targets)
dmr_results <- dmr(m_values_bmiq_no_outliers, targets, design)

write.csv(dmr_results$results_df, file.path(data_dir, "results/GSE145361_DMR_results.csv"), row.names = FALSE)
saveRDS(dmr_results$results_ranges, file.path(data_dir, "results/GSE145361_DMR_results_ranges.rds"))

pal <- c("blue", "red")

sample_colors <- pal[as.factor(targets$Sample_Group)]
names(sample_colors) <- rownames(targets)

dmr_ranges <- dmr_results$results_ranges

library(lumi)
beta <- m2beta(m_values_bmiq_no_outliers)

pd_samples <- targets[targets$Sample_Group == "PD", ]
control_samples <- targets[targets$Sample_Group == "Control", ]

pd_samples$Sample_Name <- rownames(pd_samples)
control_samples$Sample_Name <- rownames(control_samples)

pd_means <- rowMeans(beta[, pd_samples$Sample_Name], na.rm = TRUE)
control_means <- rowMeans(beta[, control_samples$Sample_Name], na.rm = TRUE)

beta_aggregated <- cbind(PD = pd_means, Control = control_means)

dmr_ranges <- dmr_results_ranges
png(file = file.path(data_dir, "results/SGPD_Top_DMR.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 1,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR2.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 2,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR3.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 3,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR4.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 4,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR5.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 5,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR6.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 6,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR17.png"), width = 10, height = 8, units = "in", res = 300)
library(DMRcate)
DMR.plot(
    ranges = dmr_ranges,
    dmr = 17,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR11.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 11,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR9.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 9,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR24.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 24,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR8.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 8,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR21.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 21,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

png(file = file.path(data_dir, "results/SGPD_Top_DMR18.png"), width = 10, height = 8, units = "in", res = 300)

DMR.plot(
    ranges = dmr_ranges,
    dmr = 18,
    CpGs = beta_aggregated,
    phen.col = c("PD" = "#1f77b4", "Control" = "#d62728"),
    what = "Beta",
    arraytype = "450K",
    genome = "hg19"
)

dev.off()

library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Map probes to genomic position once, reused for every DMR boxplot below
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_gr <- GRanges(
    seqnames = ann$chr,
    ranges = IRanges(start = ann$pos, end = ann$pos),
    strand = ann$strand,
    probe_id = rownames(ann)
)

plot_dmr_boxplot <- function(dmr_index, dmr_ranges, beta, targets, ann_gr, out_dir) {
    overlaps <- findOverlaps(ann_gr, dmr_ranges[dmr_index])
    # Some annotated probes may have been filtered out of beta during QC
    dmr_cpg_probes <- intersect(ann_gr$probe_id[queryHits(overlaps)], rownames(beta))

    plot_data <- beta[dmr_cpg_probes, , drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("probe_id") %>%
        pivot_longer(cols = -probe_id, names_to = "Sample_Name", values_to = "Beta") %>%
        left_join(rownames_to_column(targets["Sample_Group"], "Sample_Name"), by = "Sample_Name")

    dmr_boxplot <- ggplot(plot_data, aes(x = Sample_Group, y = Beta, fill = Sample_Group)) +
        geom_boxplot(outlier.alpha = 0.4) +
        geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
        scale_fill_manual(values = c("PD" = "#1f77b4", "Control" = "#d62728")) +
        labs(
            title = paste("DMR", dmr_index, "Beta Value Distribution"),
            x = "Sample Group",
            y = "Beta Value"
        ) +
        theme_minimal() +
        theme(legend.position = "none")

    png(file = file.path(out_dir, paste0("SGPD_DMR", dmr_index, "_Boxplot.png")), width = 6, height = 6, units = "in", res = 300)
    print(dmr_boxplot)
    dev.off()
}

dmr_boxplot_indices <- c(5, 17, 1, 11, 24, 8, 21, 18, 6)
for (dmr_index in dmr_boxplot_indices) {
    plot_dmr_boxplot(dmr_index, dmr_ranges, beta, targets, ann_gr, file.path(data_dir, "results"))
}

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

dmr_effect_sizes <- calculate_effect_size_all_dmrs(dmr_ranges, beta, targets)

dmr_effect_sizes <- dmr_effect_sizes[order(-abs(dmr_effect_sizes$deltabeta)), ]
dmr_effect_sizes %>% head()

write.csv(dmr_effect_sizes, file.path(data_dir, "results/SGPD_DMR_Effect_Sizes.csv"), row.names = FALSE)

## ==== No...
library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # Or EPIC, depending on the cohort

# 1. Get the coordinates for your top DMR
top_dmr <- dmr_ranges[1]

# 2. Load the standard Illumina annotation to get the coordinates of all probes
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_gr <- GRanges(
    seqnames = ann$chr,
    ranges = IRanges(start = ann$pos, end = ann$pos),
    strand = ann$strand,
    probe_id = rownames(ann)
)

# 3. Find which probes overlap with your top DMR
overlaps <- findOverlaps(ann_gr, top_dmr)
dmr_probes <- ann_gr$probe_id[queryHits(overlaps)]

library(reshape2)
library(dplyr)

# 1. Subset your Beta matrix to only include the probes in DMR #1
dmr_available_probes <- intersect(rownames(beta), dmr_probes)
dmr_beta_matrix <- beta[dmr_available_probes, ]

dmr_beta_matrix %>% dim()

write.csv(dmr_beta_matrix, file.path(data_dir, "results/SGPD_Top_DMR_Beta_Matrix.csv"), row.names = TRUE)

# 2. Melt the matrix into a long dataframe
plot_data <- melt(dmr_beta_matrix)
colnames(plot_data) <- c("Probe", "Sample_Name", "Beta")

# 3. Merge in your target data so you have the Sample_Group labels
# (Assuming targets_test has a column called 'Sample' matching the matrix column names)
plot_data <- merge(plot_data, targets[, c("Sample_Name", "Sample_Group")], by.x = "Sample_Name", by.y = "row.names")

# 4. Merge in the genomic positions of the probes so the x-axis is ordered correctly
probe_positions <- data.frame(Probe = ann_gr$probe_id, Position = start(ann_gr))
plot_data <- merge(plot_data, probe_positions, by = "Probe")

library(ggplot2)

mplot <- ggplot(plot_data, aes(x = Position, y = Beta, color = Sample_Group, fill = Sample_Group)) +
    # Add the individual sample points (made slightly transparent)
    geom_jitter(alpha = 0.4, width = 10, size = 1.5) +

    # Add the smoothed trend line (the biological average)
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +

    # Force the Y-axis to strict biological percentages
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +

    # Apply a soft pastel palette for the groups
    scale_color_brewer(palette = "Pastel1") +
    scale_fill_brewer(palette = "Pastel1") +

    # Clean up the labels
    labs(
        title = "Methylation Profile of Top DMR",
        x = paste("Genomic Position (Chromosome", seqnames(top_dmr), ")"),
        y = "Methylation (Beta Value)"
    ) +

    # Apply a clean, light-themed grid layout
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "#e5e5e5"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom"
    )


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(patchwork)

# 1. Pull the gene coordinates from the UCSC hg19 database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb)

# 2. Find which genes overlap your specific DMR region
overlapping_genes <- subsetByOverlaps(genes_gr, top_dmr)

# 3. Convert the genomic ranges object into a dataframe for ggplot
gene_df <- as.data.frame(overlapping_genes)

# 4. Map the Entrez IDs to readable Gene Symbols (e.g., "SNCA")
gene_df$symbol <- mapIds(org.Hs.eg.db,
    keys = gene_df$gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)

library(ggplot2)
library(patchwork)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# ==========================================
# STEP 1: CALCULATE SHARED BOUNDARIES
# ==========================================
# We calculate the absolute min and max positions from your methylation data first,
# so we can force both plots to adhere to this exact window.
plot_xmin <- min(plot_data$Position) - 50
plot_xmax <- max(plot_data$Position) + 50


# ==========================================
# STEP 2: BUILD THE METHYLATION PLOT
# ==========================================
meth_plot <- ggplot(plot_data, aes(x = Position, y = Beta, color = Sample_Group, fill = Sample_Group)) +
    geom_jitter(alpha = 0.4, width = 10, size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    scale_color_brewer(palette = "Pastel1") +
    scale_fill_brewer(palette = "Pastel1") +
    # INJECTION: Force the x-axis to the shared boundaries here
    coord_cartesian(xlim = c(plot_xmin, plot_xmax)) +
    labs(
        title = "Methylation Profile of Top DMR",
        y = "Methylation (Beta Value)"
    ) +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "#e5e5e5"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )


# ==========================================
# STEP 3: PREP THE GENE ANNOTATION DATA
# ==========================================
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb)
overlapping_genes <- subsetByOverlaps(genes_gr, top_dmr)
gene_df <- as.data.frame(overlapping_genes)

gene_df$symbol <- mapIds(org.Hs.eg.db,
    keys = gene_df$gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)

# INJECTION: Fix the missing gene names (NAs) immediately after mapping IDs
gene_df$symbol <- ifelse(is.na(gene_df$symbol),
    paste0("GeneID:", gene_df$gene_id),
    gene_df$symbol
)


# ==========================================
# STEP 4: BUILD THE ANNOTATION TRACK
# ==========================================
anno_plot <- ggplot(gene_df) +
    geom_segment(aes(x = start, xend = end, y = 1, yend = 1),
        linewidth = 3, color = "#6c757d", lineend = "butt"
    ) +
    # Adjust y position and vjust so the text sits comfortably inside the panel
    geom_text(aes(x = (start + end) / 2, y = 1.0, label = symbol),
        size = 4, fontface = "bold", vjust = -1
    ) +
    # Expand the y-axis upper limit to give the text plenty of headroom
    coord_cartesian(xlim = c(plot_xmin, plot_xmax), ylim = c(0.5, 1.8)) +
    labs(x = paste("Genomic Position (Chromosome", seqnames(top_dmr), ")")) +
    theme_minimal() +
    theme(
        panel.grid.major.x = element_line(color = "#e5e5e5"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
    )


# ==========================================
# STEP 5: STITCH AND SAVE
# ==========================================
final_combined_plot <- meth_plot / anno_plot + plot_layout(heights = c(4, 1))
final_combined_plot

ggsave(file.path(data_dir, "results/SGPD_Custom_Annotated_DMR_Fixed.pdf"),
    plot = final_combined_plot,
    width = 10,
    height = 7
)
