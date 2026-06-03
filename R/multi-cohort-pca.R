rm(list = ls())
gc(full = TRUE)

library(arrow)
library(dplyr)

rootDir <- "/workspace"
gse111629_data <- read_parquet(file.path(rootDir, "GSE111629_data.parquet"))
ppmi_data <- read_parquet(file.path(rootDir, "PPMI_data.parquet"))
gse145361_data <- read_parquet(file.path(rootDir, "GSE145361_data.parquet"))

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

gse145361_data_n <- gse145361_data[, startsWith(colnames(gse145361_data), "cg")]
gse145361_data_n$Sample_Group <- gse145361_data$Sample_Group
gse145361_data_n$SEX <- gse145361_data$`gender:ch1`
gse145361_data_n$Age_Group <- NA
gse145361_data_n$Sample_Name <- gse145361_data$Sample_Name
gse145361_data_n$Cohort <- "GSE145361"

gse145361_reduced <- gse145361_data_n[, colnames(gse145361_data_n) %in% colnames(gse111629_data_common)]
dim(gse145361_reduced)

gse145361_data_common <- gse145361_data_n[, colnames(gse145361_data_n) %in% colnames(gse111629_data_common)]

gse111629_data_common <- gse111629_data_common[, colnames(gse111629_data_common) %in% colnames(gse145361_data_common)]
dim(gse111629_data_common)

ppmi_data_common <- ppmi_data_common[, colnames(ppmi_data_common) %in% colnames(gse111629_data_common)]
dim(ppmi_data_common)

merged_data <- bind_rows(ppmi_data_common, gse145361_data_common, gse111629_data_common) %>%
  as.data.frame()
merged_data %>%
  colnames() %>%
  startsWith("cg") %>%
  merged_data[.] %>%
  is.na() %>%
  sum()

rm(list = setdiff(ls(), c("merged_data")))
gc(full = TRUE)

saveRDS(merged_data, "/workspace/merged_methylation_data.rds")
# PCA rows must be samples, probes columns
m_values <- merged_data[, startsWith(colnames(merged_data), "cg")]
pca <- prcomp(m_values, center = TRUE, scale. = FALSE)
dim(pca$x)

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
pcaplot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point() +
  labs(title = "PCA of Methylation Data Colored by Cohort") +
  theme_minimal()

ggsave("/workspace/PCA_cohorts.png",
  pcaplot,
  dpi = 300,
  limitsize = FALSE
)

saveRDS(pca_df, "/workspace/multi_cohort_pca_df.rds")

# Reference based ComBat

batch <- merged_data$Cohort
batch %>% table()

mod <- model.matrix(~ as.factor(merged_data$Sample_Group), data = merged_data)

ref_id <- "GSE145361"

colnames(m_values) %>% head()

library(sva)
harmonized_m_values <- ComBat(
  dat = t(m_values),
  batch = batch,
  mod = mod,
  ref.batch = ref_id
)

saveRDS(harmonized_m_values, "/workspace/harmonized_m_values.rds")

harmonized_m_values %>% head()
pca <- prcomp(t(harmonized_m_values), center = TRUE, scale. = FALSE)
dim(m_values)
pca_df_harmonized_cohorts <- data.frame(matrix(NA, nrow = nrow(t(harmonized_m_values)), ncol = npc))
dim(pca_df_harmonized_cohorts)
rownames(pca_df_harmonized_cohorts) <- colnames(harmonized_m_values)
colnames(pca_df_harmonized_cohorts) <- sapply(1:npc, function(index) {
  (paste0("PC", index))
})
for (index in 1:npc) {
  pca_df_harmonized_cohorts[[paste0("PC", index)]] <- pca$x[, index]
}
pca_df_harmonized_cohorts %>% dim()
dim(merged_data)
pca_df_harmonized_cohorts$Sample_Group <- merged_data$Sample_Group
pca_df_harmonized_cohorts$SEX <- merged_data$SEX
pca_df_harmonized_cohorts$Age_Group <- merged_data$Age_Group
pca_df_harmonized_cohorts$Sample_Name <- merged_data$Sample_Name
pca_df_harmonized_cohorts$Cohort <- merged_data$Cohort

head(pca_df_harmonized_cohorts)

saveRDS(pca_df_harmonized_cohorts, "/workspace/multi_cohort_pca_df_harmonized.rds")

library(ggplot2)
pcaplot <- ggplot(pca_df_harmonized_cohorts, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point() +
  labs(title = "PCA of Methylation Data Colored by Cohort - Harmonized via Combat w. GSE145361 as batch id") +
  theme_minimal()

ggsave("/workspace/PCA_cohorts_harmonized.png",
  pcaplot,
  dpi = 300,
  limitsize = FALSE
)

