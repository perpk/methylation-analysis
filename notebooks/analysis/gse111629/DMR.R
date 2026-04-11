rm(list = ls())
gc(full = TRUE)
install.packages("interp")

library(DMRcate)

rootDir <- "/Volumes/Elements/vastai/gse111629"

targets <- readRDS(file.path(rootDir, "qc/targets_s_mismatch_cells.rds"))

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

m_vals <- readRDS(file.path(rootDir, "results/m_values_bmiq.rds"))

head(m_vals)

targets <- targets[match(colnames(m_vals), rownames(targets)), ]

dim(targets)
dim(m_vals)
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

write.csv(results_ranges_df, file.path(rootDir, "local/DMR_results_hg19_GSE111629.csv"), row.names = FALSE)
View(results_ranges_df)

beta_means_filtered <- read.csv(file.path(rootDir, "results/beta_means_annotated_filtered.csv"))
colnames(beta_means_filtered)

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
