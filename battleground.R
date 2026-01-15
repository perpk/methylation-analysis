rm(list = ls())
gc(full = TRUE)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(karyoploteR)
library(GenomicRanges)

ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

dmp_ann_res <- read.csv("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/dmp_annotated_results.csv", row.names = 2)
dmp <- dmp_ann_res[abs(dmp_ann_res$logFC) > 0.3 & dmp_ann_res$adj.P.Val < 0.05, ]
beta_means <- read.csv("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/beta_means.csv", row.names = 1)
targets <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/qc/targets_s_mismatch_cells_pcs.rds")

dmp <- merge(
  dmp,
  beta_means,
  by = 0
)

as.factor(dmp_ann_res$chr)
dim(dmp_ann_res[dmp_ann_res$chr %in% c("chrX", "chrY"), ])
xy_chrom_probes <- rownames(dmp_ann_res[dmp_ann_res$chr %in% c("chrX", "chrY"), ])

beta_matrix_noob <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/beta_matrix_noob.rds")

str(beta_matrix_noob)
head(rownames(beta_matrix_noob))

dim(beta_matrix_noob[rownames(beta_matrix_noob) %in% xy_chrom_probes, ])

methyl_set_norm <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/processed/methyl_set_norm.rds")
gr_set <- ratioConvert(methyl_set_norm, type = "Illumina")
gr_set <- mapToGenome(gr_set)

is_sex_chr <- seqnames(gr_set) %in% c("chrX", "chrY")
sum(is_sex_chr)

gr_set <- gr_set[!is_sex_chr, ]
cat("Probes after sex chromosome removal:", nrow(gr_set), "\n")
beta_matrix_noob <- getBeta(gr_set)

dim(beta_matrix_noob[rownames(beta_matrix_noob) %in% xy_chrom_probes, ])

m_values_bmiq <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/m_values_bmiq.rds")
beta_matrix_bmiq <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/beta_matrix_bmiq.rds")

library(tidyverse)
dim(beta_matrix_bmiq[rownames(beta_matrix_bmiq) %in% xy_chrom_probes, ])
dim(m_values_bmiq[rownames(m_values_bmiq) %in% xy_chrom_probes, ])
rownames(dmp) <- dmp$Row.names
dmp <- dmp[, -1]

dmp$start <- dmp$pos
dmp$end <- dmp$pos
dmp$type <- ifelse(dmp$logFC > 0, "hyper", "hypo")
dmp$name <- rownames(dmp)

dmp_vis <- data.frame(
  chr = dmp$chr,
  start = dmp$start,
  stop = dmp$end,
  type = dmp$type,
  name = dmp$name
)

dmp_gr <- makeGRangesFromDataFrame(
  dmp_vis,
  keep.extra.columns = TRUE # Keep 'type' and 'name' columns
)

hyper_dmp <- dmp_gr[dmp_gr$type == "hyper"]
hypo_dmp <- dmp_gr[dmp_gr$type == "hypo"]

kp <- plotKaryotype(
  genome = "hg19",
  plot.type = 2,
  main = "DMP Distribution Across Chromosomes",
  cex = 1.0
)

kpPoints(
  kp,
  data = hyper_dmp,
  y = 0.7,
  pch = 19,
  col = "red",
  cex = 0.5
)
kpPoints(
  kp,
  data = hypo_dmp,
  y = 0.3,
  pch = 19,
  col = "blue",
  cex = 0.5
)

legend(x = "bottomright", fill = c("red", "blue"), legend = c("Hypermethylated", "Hypomethylated"))

# this one contains sex-chromosome probes
m_set_before <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/processed/methyl_set_removed_snps.rds")
m_set_after <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/processed/methyl_set_filtered_chrom.rds")

probes_before <- !(as.character(seqnames(m_set_before)) %in% c("chrX", "chrY"))
num_autosomal_before <- sum(probes_before)
num_sex_before <- sum(!probes_before)
num_sex_before

probes_after <- !(as.character(seqnames(m_set_after)) %in% c("chrX", "chrY"))
num_autosomal_after <- sum(probes_after)
num_sex_after <- sum(!probes_after)
num_sex_after

head(rownames(m_set_before))

m_set_before_rows <- rownames(m_set_before)
m_set_after_rows <- rownames(m_set_after)

class(m_set_after_rows)

dmp[dmp$chr == "chrX", ]

"cg07660381" %in% m_set_after_rows
"cg07660381" %in% m_set_before_rows

rg_set_filtered <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/processed/rg_set_filtered.rds")

probe_id <- "cg07660381"

cat("In rgSet:", probe_id %in% featureNames(rg_set_filtered), "\n")

beta_noob <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/results/beta_matrix_noob.rds")
beta_noob_rows <- rownames(beta_noob)

probe_id %in% beta_noob_rows

methyl_set_norm <- readRDS("/Volumes/Elements/methyl-pipe-out/GSE111629_20251226_102044/processed/methyl_set_norm.rds")

probe_id %in% rownames(methyl_set_norm)

library(tidyverse)

m_set_norm <- mapToGenome(methyl_set_norm)

norm_probes <- !(as.character(seqnames(m_set_norm)) %in% c("chrX", "chrY"))
num_autosomal_before <- sum(norm_probes)
num_sex_before <- sum(!norm_probes)

num_autosomal_before
num_sex_before

if(probe_id %in% featureNames(rg_set_filtered)) {
  cat("Probe IS in rg_set according to featureNames\n")
  # Get its information
  probe_info <- getProbeInfo(rg_set_filtered, probe_id)
  print(probe_info)
} else {
  cat("Probe is NOT in rg_set featureNames\n")
  
  # Check if it's hidden somewhere else
  cat("Checking all slots...\n")
  slot_names <- slotNames(rg_set_filtered)
  for(slot in slot_names) {
    obj <- slot(rg_set_filtered, slot)
    if(is.matrix(obj) || is.data.frame(obj)) {
      if(probe_id %in% rownames(obj)) {
        cat("Found in slot:", slot, "\n")
      }
    }
  }
}

cat("\nRGChannelSet structure:\n")
str(rg_set_filtered, max.level = 2

# Check the discrepancy


cat("=== Probe Name Investigation ===\n")

# 1. Check different ways to get probe names
cat("Length of @NAMES slot:", length(rg_set_filtered@NAMES), "\n")
cat("featureNames() length:", length(featureNames(rg_set_filtered)), "\n")
cat("nrow() of rg_set:", nrow(rg_set_filtered), "\n")

# 2. Are they the same?
names_from_slot <- rg_set_filtered@NAMES
names_from_feature <- featureNames(rg_set_filtered)

cat("\nAre they identical?", identical(names_from_slot, names_from_feature), "\n")

# 3. Look for our probe in both
probe_id <- "cg07660381"
cat("\nLooking for probe", probe_id, ":\n")
cat("In @NAMES slot?", probe_id %in% names_from_slot, "\n")
cat("In featureNames?", probe_id %in% names_from_feature, "\n")

# 4. If it's in @NAMES but not featureNames, show surrounding probes
if(probe_id %in% names_from_slot && !probe_id %in% names_from_feature) {
  idx <- which(names_from_slot == probe_id)
  cat("\nProbe found at index", idx, "in @NAMES slot\n")
  cat("Surrounding probes in @NAMES:\n")
  print(names_from_slot[max(1, idx-2):min(length(names_from_slot), idx+2)])
}
