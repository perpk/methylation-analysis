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

beta_matrix <- readRDS("/root/data/beta_matrix.rds")
# beta_matrix_bmiq <- readRDS("/root/data/beta_matrix_bmiq.rds")
# m_values_bmiq <- readRDS("/root/data/m_values_bmiq.rds")

targets <- readRDS("/root/data/targets_s_mismatch_cells_scandate.rds")

pca_df_with_outliers <- readRDS("/root/data/pca_df_with_outliers.rds")

outlier_samples <- rownames(pca_df_with_outliers)[pca_df_with_outliers$Is_Outlier]
beta_matrix_no_outliers <- beta_matrix[, !colnames(beta_matrix) %in% outlier_samples]

beta_matrix %>% dim()
beta_matrix_no_outliers %>% dim()

beta_bmiq <- matrix(NA, nrow = nrow(beta_matrix_no_outliers), ncol = ncol(beta_matrix_no_outliers))
rownames(beta_bmiq) <- rownames(beta_matrix_no_outliers)
colnames(beta_bmiq) <- colnames(beta_matrix_no_outliers)

.get_probe_design_vector <- function(probe_names, platform = IlluminaHumanMethylationEPICanno.ilm10b4.hg19) {
  # Returns: 1 for Type I probes, 2 for Type II probes
  anno <- getAnnotation(platform)

  # Match probes to annotation (handle missing probes)
  matched_probes <- intersect(probe_names, rownames(anno))

  if (length(matched_probes) < length(probe_names)) {
    warning(paste(
      length(probe_names) - length(matched_probes),
      "probes not found in annotation"
    ))
  }

  # Get probe types
  probe_types <- rep(NA, length(probe_names))
  names(probe_types) <- probe_names

  # Type I = 1, Type II = 2 (BMIQ convention)
  probe_types[matched_probes] <- ifelse(anno[matched_probes, "Type"] == "I", 1, 2)

  # Check for NAs
  if (any(is.na(probe_types))) {
    warning("Some probes could not be classified. Using alternative method...")

    # Alternative: Use probe name patterns
    # Type I: typically start with cg and have specific Infinium I design
    # Type II: typically start with ch or some cg with Infinium II design
    # This is less reliable but works as fallback
    na_probes <- probe_names[is.na(probe_types)]
    probe_types[na_probes] <- ifelse(grepl("^ch", na_probes), 2, 1)
  }

  return(probe_types)
}

design.v <- .get_probe_design_vector(rownames(beta_matrix_no_outliers))

print(paste("Type I probes:", sum(design.v == 1)))
print(paste("Type II probes:", sum(design.v == 2)))

beta_bmiq <- matrix(NA, nrow = nrow(beta_matrix_no_outliers), ncol = ncol(beta_matrix_no_outliers))
rownames(beta_bmiq) <- rownames(beta_matrix_no_outliers)
colnames(beta_bmiq) <- colnames(beta_matrix_no_outliers)

library(wateRmelon)
for (i in 1:ncol(beta_matrix_no_outliers)) {
  print(paste("Processing sample", i, "of", ncol(beta_matrix_no_outliers)))

  # Extract beta values for this sample
  beta.v <- beta_matrix_no_outliers[, i]

  # Apply BMIQ to this sample
  tryCatch(
    {
      bmiq_result <- BMIQ(
        beta.v = beta.v,
        design.v = design.v,
        nfit = 10000,
        plots = FALSE,
        pri = TRUE
      )
      print(bmiq_result$nbeta[1:10])

      # Store the normalized values
      beta_bmiq[, i] <- bmiq_result$nbeta
    },
    error = function(e) {
      warning(paste("BMIQ failed for sample", colnames(beta_matrix_no_outliers)[i], ":", e$message))
      beta_bmiq[, i] <- beta.v # Keep original values if BMIQ fails
    }
  )
}

saveRDS(beta_bmiq, "/root/data/beta_matrix_bmiq_new.rds")

m_bmiq <- log2(beta_bmiq / (1 - beta_bmiq))
m_bmiq[is.infinite(m_bmiq)] <- NA

m_raw <- log2(beta_matrix_no_outliers / (1 - beta_matrix_no_outliers))

targets %>% dim()
m_bmiq %>% dim()

targets <- targets %>% filter(rownames(targets) %in% m_bmiq %>% colnames())

m_bmiq %>%
  colnames() %>%
  head()

targets %>%
  pull(Sample_Name) %>%
  head()

targets_harmonized <- targets[targets$Sentrix_ID %in% colnames(m_bmiq), ]

targets_harmonized %>% head()
targets %>% head()

design <- model.matrix(
  ~ Sample_Group + Age_Group + Sex,
  data = targets_harmonized
)
fit <- lmFit(m_raw, design)
fit2 <- eBayes(fit)
results <- topTable(fit2, coef = "Sample_GroupPD", number = Inf, adjust.method = "BH")

table(results$adj.P.Val < 0.05)
