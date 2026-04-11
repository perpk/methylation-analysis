rootDir <- "/Volumes/Elements/vastai/gse111629"

beta_means <- read.csv(file.path(rootDir, "results/beta_means.csv"))

head(beta_means)

dim(beta_means[abs(beta_means$delta_beta) > 0.07, ])

library(minfi)

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

df_annotated <- merge(beta_means, ann[, c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island", "DMR", "Regulatory_Feature_Name", "Regulatory_Feature_Group")],
    by.x = "X", by.y = "Name"
)

beta_means_filtered <- as.data.frame(df_annotated[abs(df_annotated$delta_beta) > 0.07, ])

write.csv(beta_means_filtered, file.path(rootDir, "results/beta_means_annotated_filtered.csv"), row.names = FALSE)

View(beta_means_filtered)

colnames(ann)

rm(beta_means)
gc(full = TRUE)
