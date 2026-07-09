library(GEOquery)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library(FlowSorted.Blood.EPIC)
library(ggplot2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

source("R/intermediate_data_proxy.R")

pre_process_eda <- function(
  project_name = NULL,
  project_to_load = NULL,
  targets = NULL,
  data_folder = NULL,
  project_location = NULL,
  var_mapping = NULL,
  platform = NULL,
  qc_threshold = 10.5,
  cohorts = NULL,
  mode = results_mode()$memory_only
) {
  if (is.null(platform)) {
    stop("Platform must be specified as '450K' or 'EPIC'")
  }
  source("R/progress_mgr.R")

  source("R/project_context.R")

  print(paste("mode is:", mode))
  if (is.null(project_to_load)) {
    # Create a new project
    project_context <- create_methylation_project(project_name, project_location, platform = platform, cohorts = cohorts, mode = mode)
  } else {
    # Load an existing project
    project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
  }

  print(paste("Saving initial targets to", file.path(project_context$paths$processed, "targets.rds")))
  saveRDS(targets, file.path(project_context$paths$processed, "targets.rds"))

  ### Read raw data and write to disk
  #### Returns a list with:
  #### - rg_set_container: A ResultsContainer object containing the raw RG set
  #### - m_set_container: A ResultsContainer object containing the raw methylation set
  source("R/extract_methyl_set.R")
  res_extract_methyl_set <- intermediate_data_proxy(
    extract_methyl_set, project_context,
    targets = targets
  )

  ### Sample QC - Outlier detection and removal: Here samples are removed based on the median methylated and unmethylated signal intensities. Samples with a median methylated or unmethylated signal intensity below the specified threshold (default: 10.5) are flagged as outliers and removed from the dataset.
  #### It is crucial to do this step alongside cleaning up the data from problematic samples in order to not have outliers skewing the normalization and thus the downstream analyses. It is also important to do this step before the biological gender mismatch analysis since the outliers can
  #### Returns a list with ResultsContainer objects containing the filtered RG set and methylation set after sample QC:
  #### - qc_results: A ResultsContainer object containing the filtered RG set and methylation set after sample QC
  #### - rg_set_results: A ResultsContainer object containing the filtered RG set
  #### - targets_results: A ResultsContainer object containing the filtered targets data frame
  #### - methyl_set_results: A ResultsContainer object containing the filtered methylation set
  #### - bisulfite_thresholds_results: A ResultsContainer object containing the bisulfite conversion control thresholds
  source("R/qc.R")
  res_qc <- intermediate_data_proxy(
    qc,
    project_context,
    targets = targets,
    rg_set = res_extract_methyl_set$rg_set_container@object,
    methyl_set = res_extract_methyl_set$m_set_container@object,
    rg_set_filename = res_extract_methyl_set$rg_set_container@filename,
    methyl_set_filename = res_extract_methyl_set$m_set_container@filename,
    qc_threshold = qc_threshold
  )

  if (project_context$mode == results_mode()$memory_only) {
    print(paste("Saving targets after sample QC to", file.path(project_context$paths$processed, "targets_after_qc.rds")))
    saveRDS(res_qc$targets_results@object, file.path(project_context$paths$processed, "targets_after_qc.rds"))
  }

  ### Perform background correction and dye-bias normalization on rg_set and extract new methyl_set & beta-matrix based on the filtered rg_set from previous step
  ### Here, preprocessNoob is used and by doing so on the rgset, the methyl_set emerges.
  #### Returns a list with ResultsContainer objects containing the normalized methylation set and RG set:
  #### - methyl_set_container: A ResultsContainer object containing the normalized methylation set
  #### - rg_set_container: A ResultsContainer object containing the normalized RG set
  source("R/bg_correction_dye_bias_norm.R")
  res_bg_corr <- intermediate_data_proxy(
    bg_correction_dye_bias_norm,
    project_context,
    rg_set = res_qc$rg_set_results@object,
    rg_set_filename = res_qc$rg_set_results@filename
  )

  # Probe QC - Detection p-value based probe filtering
  # This is performed on the rg_set and is done after normalization since the probes are necessary to provide an unbiased normalization
  # Respectively the probes are removed from methyl_set as well and a new, filtered version gets persisted to the filesystem.
  # The new filename for the filtered methyl_set is "methyl_set_probe_qc.rds"
  ## Returns a list with ResultsContainer objects containing the candidate probes to remove and a summary of the probe QC results:
  ## - removed_probes_container: A ResultsContainer object containing the candidate probes to remove based on detection p-values
  ## - qc_summary_container: A ResultsContainer object containing a summary of the probe QC results
  source("R/probe_qc.R")
  res_probe_qc <- intermediate_data_proxy(
    probe_qc,
    project_context,
    rg_set = res_bg_corr$rg_set_container@object,
    rg_set_filename = res_bg_corr$rg_set_container@filename
  )

  ### Remove sex-mismatched samples
  ### This operation is performed on the methyl_set. Also, the mismatched probes are removed from the rgset as well.
  #### Returns a list with ResultsContainer objects containing the filtered methylation set and RG set after removing sex-mismatched samples:
  #### - methyl_set_container: A ResultsContainer object containing the filtered methylation set after removing sex-mismatched samples
  #### - rg_set_container: A ResultsContainer object containing the filtered RG set after removing sex-mismatched samples
  #### - targets_container: A ResultsContainer object containing the filtered targets data frame after removing sex-mismatched samples
  source("R/biological_gender_mismatch_analysis.R")
  res_bio_gender_mismatch <- intermediate_data_proxy(
    biological_gender_mismatch_analysis,
    project_context,
    recorded_sex_col = var_mapping$gender_var,
    methyl_set = res_bg_corr$methyl_set_container@object,
    rg_set = res_bg_corr$rg_set_container@object,
    targets = res_qc$targets_results@object,
  )

  if (project_context$mode == results_mode()$memory_only) {
    print(paste("Saving targets after biological gender mismatch analysis to", file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds")))
    saveRDS(res_bio_gender_mismatch$targets_container@object, file.path(project_context$paths$processed, "targets_after_bio_gender_mismatch.rds"))
  }

  removed_pdp <- res_probe_qc$removed_probes_container@object
  if (is.null(removed_pdp)) {
    removed_pdp <- readRDS(res_probe_qc$removed_probes_container@filename)
  }
  rg_set <- res_bio_gender_mismatch$rg_set_container@object
  if (is.null(rg_set)) {
    rg_set <- readRDS(res_bio_gender_mismatch$rg_set_container@filename)
  }
  m_set <- res_bio_gender_mismatch$methyl_set_container@object
  if (is.null(m_set)) {
    m_set <- readRDS(res_bio_gender_mismatch$methyl_set_container@filename)
  }

  rg_set <- rg_set[!rownames(rg_set) %in% removed_pdp[["probe_id"]], ]
  m_set <- m_set[!rownames(m_set) %in% rownames(removed_pdp), ]

  if (project_context$mode == results_mode()$disk_only ||
    project_context$mode == results_mode()$disk_and_memory) {
    saveRDS(rg_set, file.path(project_context$paths$processed, "rg_set_remove_probe_qc.rds"))
    saveRDS(m_set, file.path(project_context$paths$processed, "methyl_set_remove_probe_qc.rds"))
  }

  # Remove cross-reactive probes, sex-chromosome related probes and single nucleotide polymorphisms (SNPs)
  # Order matters, firstly SNPs must be removed then probes on XY chromosomes and thus keeping only those on autosomal and finally filtering of cross-reactive probes.
  # All of the following, enumerated operations are performed on the methyl_set.
  ## 1. Single Nucleotide Polymorphisms
  source("R/remove_snp.R")
  res_snp <- intermediate_data_proxy(
    remove_snp,
    project_context,
    methyl_set = m_set,
    methyl_set_file =
      file.path(
        project_context$paths$processed,
        "methyl_set_remove_probe_qc.rds"
      )
  )

  ### 2. Sex-Chromosome Probe filtering
  source("R/remove_sex_chromosome_probes.R")
  res_sex_chromosome <- intermediate_data_proxy(
    remove_sex_chromosome_probes,
    project_context,
    methyl_set = res_snp$methyl_set_removed_snps_container@object,
    methyl_set_filename = res_snp$methyl_set_removed_snps_container@filename
  )

  #### 3. Cross Reactive Probes
  source("R/remove_cross_reactive_probes.R")
  res_cross_reactive <- intermediate_data_proxy(
    remove_cross_reactive_probes,
    project_context,
    methyl_set = res_sex_chromosome$mismatch_container@object,
    methyl_set_filename = res_sex_chromosome$mismatch_container@filename
  )

  methyl_set_final <- res_cross_reactive$methyl_set_cross_reactive_clean_container@object
  if (is.null(methyl_set_final)) {
    methyl_set_final <- readRDS(res_cross_reactive$methyl_set_cross_reactive_clean_container@filename)
  }

  ### The filtering of the dataset is complete and beta & m-values are now extracted from the methyl_set and written to disk as "beta_matrix.rds" and "m_matrix.rds" respectively.
  beta_matrix <- getBeta(methyl_set_final)
  m_matrix <- getM(methyl_set_final)

  beta_matrix_filepath <- NULL
  m_matrix_filepath <- NULL
  if (project_context$mode == results_mode()$disk_only ||
    project_context$mode == results_mode()$disk_and_memory) {
    beta_matrix_filepath <- file.path(project_context$paths$results, "beta_matrix.rds")
    m_matrix_filepath <- file.path(project_context$paths$results, "m_matrix.rds")
    saveRDS(beta_matrix, beta_matrix_filepath)
    saveRDS(m_matrix, m_matrix_filepath)
  }

  beta_matrix_container <- new("ResultsContainer", filename = beta_matrix_filepath, object = beta_matrix, future = NULL)
  m_matrix_container <- new("ResultsContainer", filename = m_matrix_filepath, object = m_matrix, future = NULL)

  # Beta-Mixture Quantile (BMIQ) Normalization
  source("R/apply_BMIQ.R")
  bmiq_res <- intermediate_data_proxy(
    apply_BMIQ,
    project_context,
    beta_matrix = beta_matrix_container@object,
    beta_matrix_file = beta_matrix_container@filename,
    plot = TRUE
  )

  # Principal Component Analysis
  source("R/principal_component_analysis.R")
  col_map <- list()
  col_map[["Sample_Group"]] <- "Sample_Group"
  col_map[["Gender"]] <- var_mapping$gender_var
  col_map[["Age"]] <- var_mapping$age_var
  keys <- c("Sample_Group", "Gender", "Age")

  res_pca <- intermediate_data_proxy(
    principal_component_analysis,
    project_context,
    m_values = bmiq_res$m_bmiq_container@object,
    m_values_filename = bmiq_res$m_bmiq_container@filename,
    targets = res_bio_gender_mismatch$targets_container@object,
    targets_filename = res_bio_gender_mismatch$targets_container@filename,
    col_maps = col_map,
    keys = keys,
    npc = 5
  )

  #### Plot PCA
  source("R/plot_PCA.R")
  pca_vars <- c("Sample_Group" = "By Diagnosis", "Gender" = "By Biological Gender", "Age" = "By Age")
  convert_f <- function(df) {
    (transform(df, Age = as.numeric(Age)))
  }
  plot_PCA(
    context = project_context,
    pca_df = res_pca$pca_df_container@object,
    pca_results_rds_filename = res_pca$pca_df_container@filename,
    pca_vars = pca_vars,
    convert_fun = convert_f,
    continuously_scaled = c("Age"),
    "pca_plot_",
    pairplot_color_by = "Age"
  )

  ## Outlier detection from PCA
  source("R/outlier_analysis.R")

  res_outlier <- intermediate_data_proxy(
    outlier_analysis,
    project_context,
    pca = res_pca$pca_df_container@object,
    pca_filename = res_pca$pca_df_container@filename,
    targets = res_bio_gender_mismatch$targets_container@object,
    targets_filename = res_bio_gender_mismatch$targets_container@filename,
    beta_bmiq = bmiq_res$beta_bmiq_container@object,
    beta_bmiq_filename = bmiq_res$beta_bmiq_container@filename
  )

  source("R/outlier_remove_redo_BMIQ.R")
  outlier_removed_bmiq_res <- outlier_remove_redo_BMIQ(
    context = project_context,
    pca = res_outlier$pca_outliers_container@object,
    beta_matrix = bmiq_res$beta_bmiq_container@object,
    pca_filename = res_outlier$pca_outliers_container@filename,
    beta_matrix_filename = bmiq_res$beta_bmiq_container@filename
  )

  source("R/apply_BMIQ.R")
  plot_BMIQ(
    project_context,
    beta_before = bmiq_res$beta_bmiq_container@object,
    beta_after = outlier_removed_bmiq_res$beta_bmiq_container@object,
    beta_before_filename = bmiq_res$beta_bmiq_container@filename,
    beta_after_filename = outlier_removed_bmiq_res$beta_bmiq_container@filename
  )

  ## Calculate Delta-Beta values
  print("Calculate Delta-Beta values")

  beta_matrix <- outlier_removed_bmiq_res$beta_bmiq_container@object
  if (is.null(beta_matrix)) {
    beta_matrix <- readRDS(outlier_removed_bmiq_res$beta_bmiq_container@filename)
  }
  targets <- res_bio_gender_mismatch$targets_container@object
  if (is.null(targets)) {
    targets <- readRDS(res_bio_gender_mismatch$targets_container@filename)
  }
  library(dplyr)
  library(stringr)

  rownames(targets) <- targets$Basename %>% str_remove(paste0(data_folder, "/"))

  pd_samples <- rownames(targets[targets$Sample_Group == "PD", ])
  hc_samples <- rownames(targets[targets$Sample_Group == "Control", ])

  colnames_pd <- which(colnames(beta_matrix) %in% pd_samples)
  colnames_hc <- which(colnames(beta_matrix) %in% hc_samples)

  mean_beta_pd <- rowMeans(beta_matrix[, colnames_pd], na.rm = TRUE)
  mean_beta_hc <- rowMeans(beta_matrix[, colnames_hc], na.rm = TRUE)
  delta_beta <- mean_beta_pd - mean_beta_hc

  beta_means <- data.frame(
    mean_beta_pd,
    mean_beta_hc,
    delta_beta
  )

  write.csv(beta_means, file.path(project_context$paths$results, "beta_means.csv"))

  # Cell Count Estimation
  source("R/cell_cnt_estimate.R")
  res_cell_cnt_estimate <- intermediate_data_proxy(
    cell_cnt_estimate,
    project_context,
    rg_set = res_extract_methyl_set$rg_set_container@object,
    targets = targets
  )

  # Plot Cell proportion per cohort and write results to files for each cell type
  source("R/plot_cell_proportions.R")
  plot_cell_proportions(
    context = project_context,
    targets = res_cell_cnt_estimate$targets_container@object,
    targets_file = res_cell_cnt_estimate$targets_container@filename,
    cell_types = getCellTypesForPlatform(project_context$platform)
  )

  if (project_context$mode == results_mode()$disk_and_memory) {
    library(mirai)

    all_res <- list(
      res_extract_methyl_set = res_extract_methyl_set,
      res_qc = res_qc,
      res_bg_corr = res_bg_corr,
      res_probe_qc = res_probe_qc,
      res_bio_gender_mismatch = res_bio_gender_mismatch,
      res_snp = res_snp,
      res_sex_chromosome = res_sex_chromosome,
      res_cross_reactive = res_cross_reactive,
      bmiq_res = bmiq_res,
      res_pca = res_pca,
      res_outlier = res_outlier,
      res_cell_cnt_estimate = res_cell_cnt_estimate,
    )
    for (name in names(all_res)) {
      future_obj <- all_res[[name]]
      if (is(future_obj, "list")) {
        for (item in future_obj) {
          if (is(item, "ResultsContainer") && !is.null(item@future)) {
            cat(paste("Waiting for future of", item@filename, "to complete...\n"))
            value(item@future)
          }
        }
      } else if (is(future_obj, "ResultsContainer") && !is.null(future_obj@future)) {
        cat(paste("Waiting for future of", future_obj@filename, "to complete...\n"))
        value(future_obj@future)
      }
    }
  }
}
