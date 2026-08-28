sgpd_idat_path <- "/Volumes/saucepan/methylation-project/data/GSE145361_RAW/GSE145361_RAW"
random_file_to_dissect <- "GSM4315405_3998888019_R01C01_Grn.idat.gz"

library(illuminaio)
library(dplyr)

idat_data <- readIDAT(paste0(sgpd_idat_path, "/", random_file_to_dissect))

class(idat_data)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
targets <- read.metharray.sheet(sgpd_idat_path, pattern = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")

library(GEOquery)
gse <- getGEO("GSE145361", GSEMatrix = TRUE, getGPL = FALSE)

gse[[1]]
pdata <- pData(gse[[1]])

pdata %>% head()

# # 1. Install required packages if you don't have them
# # BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# # install.packages("arrow")
# install.packages("arrow")

# library(minfi)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(arrow)

# # 2. Extract the official Bioconductor annotation object
# ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# # 3. Convert the complex DataFrame to a flat base R data.frame
# manifest_df <- as.data.frame(ann)

# # 4. FIX THE OFFSET ISSUE:
# # Move the rownames (cg probes) into an explicit column before export
# manifest_df$IlmnID <- rownames(manifest_df)

# # 5. Subset only the columns the PyTorch Graph needs to save memory
# manifest_clean <- manifest_df[, c("IlmnID", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]

# # 6. Rename columns to strictly match the Python pipeline's expectations
# colnames(manifest_clean) <- c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")

# # 7. Clean up empty Bioconductor strings to prevent Python parsing errors
# manifest_clean$UCSC_RefGene_Name[manifest_clean$UCSC_RefGene_Name == ""] <- "nan"
# manifest_clean$UCSC_RefGene_Group[manifest_clean$UCSC_RefGene_Group == ""] <- "nan"

# # 8. Export directly to Parquet (No row.names, no CSV shifting!)
# write_parquet(manifest_clean, "infinium450k_manifest.parquet")

# cat("Manifest successfully exported to infinium450k_manifest.parquet\n")
