rm(list = ls())
gc(full = TRUE)

library(minfi)

rg_set <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260406/local/rg_set.rds")

counts_epic <- estimateCellCounts2(
    rgSet = rg_set,
    compositeCellType = "Blood",
    referencePlatform = "IlluminaHumanMethylationEPIC"
)
