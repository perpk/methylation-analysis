rm(list = ls())
gc(full = TRUE)

path <- "/Volumes/Elements/vastai/ppmi/ppmi_20260403_164555"
methyl_set_norm <- readRDS(file.path(paste0(path, "/processed"), "methyl_set_norm.rds"))
targets <- readRDS(file.path(paste0(path, "/qc"), "targets_clean.rds"))

library(minfi)
methyl_set_genomic <- mapToGenome(methyl_set_norm)
predicted_sex <- getSex(methyl_set_genomic, cutoff = -2)
sex_check <- data.frame(
    Sample_Name = targets$Sample_Name,
    Recorded_Sex = targets[["SEX"]],
    Predicted_Sex = predicted_sex$predictedSex,
    X_chr_median = predicted_sex$xMed,
    Y_chr_median = predicted_sex$yMed
)
sex_check <- sex_check %>% mutate(
    Predicted_Sex = case_when(Predicted_Sex == "M" ~ "Male", Predicted_Sex == "F" ~ "Female")
)
mismatches <- sex_check[sex_check$Recorded_Sex != sex_check$Predicted_Sex, ]
print(mismatches)

sex_check
