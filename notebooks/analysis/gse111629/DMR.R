rm(list = ls())
gc(full = TRUE)
install.packages("interp")

library(DMRcate)

rootDir <- "/Volumes/Elements/vastai/gse111629/GSE111629_20260416_183131"

targets <- readRDS(file.path(rootDir, "targets_s_mismatch_cells.rds"))

targets$Sex <- targets$"gender:ch1"

library(dplyr)
targets <- targets %>%
    mutate(
        Age_Group =
            case_when(
                targets$"age:ch1" >= 30 & targets$"age:ch1" < 50 ~ "30-49",
                targets$"age:ch1" >= 50 & targets$"age:ch1" < 60 ~ "50-59",
                targets$"age:ch1" >= 60 & targets$"age:ch1" < 70 ~ "60-69",
                targets$"age:ch1" >= 70 & targets$"age:ch1" < 80 ~ "70-79",
                targets$"age:ch1" >= 80 ~ "80+",
            )
    )

m_vals <- readRDS(file.path(rootDir, "m_values_bmiq.rds"))

beta_vals <- readRDS(file.path(rootDir, "beta_matrix_bmiq.rds"))

pheno_colors <- c("Control" = "forestgreen", "PD" = "magenta") 
sample_colors <- pheno_colors[targets$Sample_Group]

targets <- targets[match(colnames(m_vals), rownames(targets)), ]

mod <- model.matrix(~ Sample_Group + Age_Group + Sex + CD4T + Bcell, data = targets)

annotated <- cpg.annotate(
    object = m_vals,
    datatype = "array",
    what = "M",
    analysis.type = "differential",
    design = mod,
    coef = 2,
    fdr = 0.05,
    arraytype = "450K"
)

dmrcoutput <- dmrcate(annotated, lambda = 1000, C = 2)
BiocManager::install("DMRcatedata")
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results_ranges_df <- as.data.frame(results.ranges)
dim(results_ranges_df)
results_ranges_df <- results_ranges_df[results_ranges_df$no.cpgs >= 3, ]
dim(results_ranges_df)

write.csv(results_ranges_df, file.path(rootDir, "DMR_results_hg19_GSE111629.csv"), row.names = FALSE)
View(results_ranges_df)

results_ranges_df <- read.csv(file.path(rootDir, "DMR_results_hg19_GSE111629.csv"))

library(ggplot2)

chrom_cpg_counts <- results_ranges_df %>%
    mutate(seqnames = as.character(seqnames)) %>%
    group_by(seqnames) %>%
    summarise(total_cpgs = sum(no.cpgs, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_cpgs)) %>%
    mutate(seqnames = factor(seqnames, levels = seqnames))

ggplot(chrom_cpg_counts, aes(x = seqnames, y = total_cpgs)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = total_cpgs), vjust = -0.3, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_minimal() +
    labs(title = "Total CpG Counts per Chromosome", x = "Chromosome", y = "Total CpG Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

chr_of_interest <- "chr6"
dmr_indices <- which(results_ranges_df$seqnames == chr_of_interest)
dmr_to_plot <- dmr_indices[1] 
dmr_chr <- results_ranges_df$seqnames[dmr_to_plot]
dmr_start <- results_ranges_df$start[dmr_to_plot]
dmr_end <- results_ranges_df$end[dmr_to_plot]

DMR.plot(
    ranges = results.ranges, 
    dmr = dmr_to_plot,                    
    CpGs = beta_vals, 
    what = "Beta",              
    arraytype = "450K",         
    phen.col = sample_colors, 
    genome = "hg19"
)

gr_probes <- makeGRangesFromDataFrame(beta_means_filtered,
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "pos",
    end.field = "pos"
)

gr_regions <- makeGRangesFromDataFrame(results_ranges_df,
    keep.extra.columns = TRUE
)

overlaps <- findOverlaps(gr_probes, gr_regions)

annotated_results <- cbind(
    beta_means_filtered[queryHits(overlaps), ],
    results_ranges_df[subjectHits(overlaps), ]
)

# Re-create after filtering according to Delta-beta > 7%
gr_regions <- makeGRangesFromDataFrame(results_ranges_df[subjectHits(overlaps), ],
    keep.extra.columns = TRUE
)

BiocManager::install("Gviz")
BiocManager::install("GenomicRanges")


itrack <- IdeogramTrack(genome = "hg19", chromosome = "chr7")

gtrack <- GenomeAxisTrack()

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

grtrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
    chromosome = "chr7",
    name = "Gene Model",
    stacking = "dense"
)

atrack <- AnnotationTrack(gr_regions, name = "CpGs", fill = "red")

beta_vals <- readRDS(file.path(rootDir, "results/beta_matrix_bmiq.rds"))
beta_vals <- beta_vals[subjectHits(overlaps), ]

cohorts <- as.factor(targets$Sample_Group)

dtrack <- DataTrack(
    range = gr_regions,
    data = t(beta_vals),
    name = "Beta",
    groups = cohorts,
    type = "boxplot",
    legend = TRUE
)

plotTracks(list(itrack, gtrack, grtrack, atrack, dtrack),
    from = 151442351,
    to = 151442967,
    chromosome = "chr7",
    background.panel = "#f6f6f6", # Light grey background
    col.title = "black",
    col.axis = "black"
)

gr_regions
