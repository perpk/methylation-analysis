rm(list = ls())
gc(full = TRUE)

library(arrow)
library(dplyr)

rootDir <- "/Volumes/Elements/vastai/combat"
gse111629_data <- read_parquet(file.path(rootDir, "GSE111629_data.parquet"))
ppmi_data <- read_parquet(file.path(rootDir, "PPMI_data.parquet"))
gse145361_data <- read_parquet(file.path(rootDir, "GSE145361_data.parquet"))

dim(ppmi_data[, startsWith(colnames(ppmi_data), "cg")])

ppmi_data_n <- ppmi_data[, startsWith(colnames(ppmi_data), "cg")]
ppmi_data_n$Sample_Group <- ppmi_data$Sample_Group
ppmi_data_n$SEX <- ppmi_data$SEX
ppmi_data_n$Age_Group <- ppmi_data$Age_Group
ppmi_data_n$Sample_Name <- ppmi_data$Sample_Name
ppmi_data_n$Cohort <- "PPMI"

gse111629_data <- gse111629_data %>%
  mutate(
    Age_Group =
      case_when(
        as.numeric(`age:ch1`) >= 30 & as.numeric(`age:ch1`) < 50 ~ "30-49",
        as.numeric(`age:ch1`) >= 50 & as.numeric(`age:ch1`) < 60 ~ "50-59",
        as.numeric(`age:ch1`) >= 60 & as.numeric(`age:ch1`) < 70 ~ "60-69",
        as.numeric(`age:ch1`) >= 70 & as.numeric(`age:ch1`) < 80 ~ "70-79",
        as.numeric(`age:ch1`) >= 80 ~ "80+",
      )
  )

gse111629_data_n <- gse111629_data[, startsWith(colnames(gse111629_data), "cg")]
gse111629_data_n$Sample_Group <- gse111629_data$Sample_Group
gse111629_data_n$SEX <- gse111629_data[["gender:ch1"]]
gse111629_data_n$Age_Group <- gse111629_data$Age_Group
gse111629_data_n$Sample_Name <- gse111629_data$Sample_Name
gse111629_data_n$Cohort <- "GSE111629"

ppmi_probes <- colnames(ppmi_data_n)[startsWith(colnames(ppmi_data_n), "cg")]
gse111629_probes <- colnames(gse111629_data_n)[startsWith(colnames(gse111629_data_n), "cg")]
common_probes <- intersect(ppmi_probes, gse111629_probes)

ppmi_data_common <- ppmi_data_n[, c(common_probes, "Sample_Group", "SEX", "Age_Group", "Sample_Name", "Cohort")]
gse111629_data_common <- gse111629_data_n[, c(common_probes, "Sample_Group", "SEX", "Age_Group", "Sample_Name", "Cohort")]

gse145361_data %>% dim()

gse145361_data %>% colnames() %>% head(50)
gse145361_data$hyb_protocol %>% table()

gse145361_data_n <- gse145361_data[, startsWith(colnames(gse145361_data), "cg")]
gse145361_data_n$Sample_Group <- gse145361_data$Sample_Group
gse145361_data_n$SEX <- gse145361_data$`gender:ch1`
gse145361_data_n$Age_Group <- NA
gse145361_data_n$Sample_Name <- gse145361_data$Sample_Name
gse145361_data_n$Cohort <- "GSE145361"

dim(ppmi_data_common)
dim(gse111629_data_common)

gse145361_reduced <- gse145361_data_n[, colnames(gse145361_data_n) %in% colnames(gse111629_data_common)]
dim(gse145361_reduced)

gse145361_data_common <- gse145361_data_n[, colnames(gse145361_data_n) %in% colnames(gse111629_data_common)]

dim(gse145361_data_common)

merged_data <- bind_rows(ppmi_data_common, gse111629_data_common, gse145361_data_common) %>%
  as.data.frame()

dim(merged_data)

rm(list = setdiff(ls(), c("merged_data")))
gc(full = TRUE)

m_values <- merged_data[, startsWith(colnames(merged_data), "cg")]
pca <- prcomp(m_values, center = TRUE, scale. = FALSE)

npc <- 2
pca_df <- data.frame(matrix(NA, nrow = nrow(m_values), ncol = npc))
dim(pca_df)
rownames(pca_df) <- rownames(m_values)
colnames(pca_df) <- sapply(1:npc, function(pca_df) {
    (paste0("PC", pca_df))
})
for (index in 1:npc) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}

pca_df$Sample_Group <- merged_data$Sample_Group
pca_df$SEX <- merged_data$SEX
pca_df$Age_Group <- merged_data$Age_Group
pca_df$Sample_Name <- merged_data$Sample_Name
pca_df$Cohort <- merged_data$Cohort

head(pca_df)

library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point() +
  labs(title = "PCA of Methylation Data Colored by Cohort") +
  theme_minimal()
