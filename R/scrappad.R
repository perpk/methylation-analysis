rm(list = ls())
gc(full = TRUE)

library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(FlowSorted.Blood.EPIC)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(tidyverse)

source("R/progress_mgr.R")
source("R/project_context.R")
source("R/intermediate_data_proxy.R")
source("R/results_container.R")

# wget https://drive.usercontent.google.com/download?id=1D9PyaHhZNBB2-DW80ZaAgcIlomfgl6lB&export=download&authuser=0&confirm=t&uuid=9608423f-4d31-4d7f-aa4f-8ee7efe1e455&at=AFYLz4Pmx-cmnG7juuQM_dl3Y9BE:1785519150140 -O targets_s_mismatch_cells_scandate.rds
# wget https://drive.usercontent.google.com/download?id=1045WmyC8psedUMwgmtEy1o-_YOZJGyBk&export=download&authuser=0&confirm=t&uuid=bef8834a-7bc2-4675-8dca-ad6eea651173&at=AFYLz4PDSABhzUFn-lDr-54S7O5D:1785519212713 -O beta_matrix.rds
# wget https://drive.usercontent.google.com/download?id=1LTlKTziuwnYND3ULSr4TtYex9GcVTwQO&export=download&authuser=0&confirm=t&uuid=216af4b7-b769-42ad-976a-d4ace364b2a1&at=AFYLz4MM5Gf62YfvUs1gYoxixHWB:1785519395943 -O beta_matrix_bmiq.rds
# wget https://drive.usercontent.google.com/download?id=1pKEXWGd-WtXHLQc8oHGSSKQpm_kV4D4E&export=download&authuser=0&confirm=t&uuid=89cca26a-db7d-4454-abb2-4a7edb022dfc&at=AFYLz4PiJ0lNZ9F0sLFpaPhXT8mU:1785519465129 -O m_values_bmiq.rds
# wget https://drive.usercontent.google.com/download?id=1PIaVMFefCEJxiW2dVaCFoMxRLFkBQCWE&export=download&authuser=0&confirm=t&uuid=535c8b79-860e-41c9-b722-b1828a528727&at=AFYLz4ORcayXZgTrw9A2li1-sDGs:1785519568230 -O pca_df_with_outliers.rds



# methyl_set_removed_cross_reactive <- readRDS("/root/data/methyl_set_removed_cross_reactive.rds")
# targets <- readRDS("/root/data/targets_s_mismatch_cells_scandate.rds")

# beta_matrix <- getBeta(methyl_set_removed_cross_reactive)
# m_matrix <- getM(methyl_set_removed_cross_reactive)

# targets %>% head()

# design <- model.matrix(
#   ~ 0 + Sample_Group + SEX,
#   data = targets
# )

# is.factor(targets$Age_Group)
# targets$Age_Group <- as.factor(targets$Age_Group)

# fit <- lmFit(m_matrix, design)
# cont.matrix <- makeContrasts(Parkinsons_vs_Control = Sample_GroupPD - Sample_GroupControl, levels = design)
# fit_contrasts <- contrasts.fit(fit, cont.matrix)
# fit2 <- eBayes(fit_contrasts)
# results <- topTable(fit2, coef = "Parkinsons_vs_Control", number = Inf, adjust.method = "BH")

# table(results$adj.P.Val < 0.05)

project_to_load <- "GSE145361_20260804_162813"
project_location <- "/root/workspace/methyl-pipe-out"

cohorts <- list(
    PD_vs_Control = c("PD", "Control")
)

source("R/project_context.R")
project_context <- .load_methylation_project(project_location, project_to_load, platform = "450k", cohorts = cohorts)


source("R/outlier_remove_redo_BMIQ.R")
outlier_removed_bmiq_res <- outlier_remove_redo_BMIQ(
  context = project_context,
  pca = NULL,
  beta_matrix = NULL,
  pca_filename = file.path(project_context$paths$results, "pca_df.rds"),
  beta_matrix_filename = file.path(project_context$paths$results, "beta_matrix.rds")
)
