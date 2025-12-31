pca_outlier_assessment <- function(context = NULL) {
  prog <- .create_progress_manager(7)

  raw_rg_set_filepath <- file.path(context$paths$raw_data, "rg_set.rds")
  prog$update(1, paste("Reading raw rg_set from", raw_rg_set_filepath))
  rg_set <- readRDS(raw_rg_set_filepath)

  outlier_filepath <- file.path(context$paths$results, "outlier_log.csv")
  prog$update(2, paste("Reading outlier log from ", outlier_filepath))
  outliers <- read.csv(outlier_filepath, row.names = 1)

  if (nrow(outliers) == 0) {
    print("No outliers available")
    prog$complete()
    return()
  }

  prog$update(3, "Calculating detection p-values")
  detP <- detectionP(rg_set)
  failure_rates <- colMeans(detP > 0.01)

  prog$update(4, "Calculating bisulfite conversion efficiency")
  bs <- bscon(rg_set)

  prog$update(5, "Consolidating metrics")
  metrics <- data.frame(
    Sample = colnames(rg_set),
    DetectionFailure = failure_rates,
    BisulfiteConversion = bs,
    IsOutlier = colnames(rg_set) %in% rownames(outliers),
    Detection_raw_p = NA_real_,
    Bisulfite_raw_p = NA_real_,
    Detection_fdr = NA_real_,
    Bisulfite_fdr = NA_real_,
    Detection_bonf = NA_real_,
    Bisulfite_bonf = NA_real_,
    stringsAsFactors = FALSE
  )

  prog$update(6, "Running statistical tests")
  det_test <- t.test(DetectionFailure ~ IsOutlier, data = metrics)
  bs_test <- t.test(BisulfiteConversion ~ IsOutlier, data = metrics)

  metrics$Detection_raw_p <- det_test$p.value
  metrics$Bisulfite_raw_p <- bs_test$p.value

  all_pvals <- c(det_test$p.value, bs_test$p.value)
  fdr_adjusted <- p.adjust(all_pvals, method = "fdr")
  bonf_adjusted <- p.adjust(all_pvals, method = "bonferroni")

  metrics$Detection_fdr <- fdr_adjusted[1]
  metrics$Bisulfite_fdr <- fdr_adjusted[2]

  metrics$Detection_bonf <- bonf_adjusted[1]
  metrics$Bisulfite_bonf <- bonf_adjusted[2]

  metrics$TechnicalFailure_FDR <- with(metrics, IsOutlier & (
    (Detection_fdr < 0.05 & DetectionFailure > 0.05) |
      (Bisulfite_fdr < 0.05 &
         BisulfiteConversion < 90)
  ))

  metrics$FailedTest <- NA_character_
  for(i in which(metrics$IsOutlier)) {
    failed_tests <- c()
    if(metrics$Detection_fdr[i] < 0.05 & metrics$DetectionFailure[i] > 0.05) {
      failed_tests <- c(failed_tests, "Detection")
    }
    if(metrics$Bisulfite_fdr[i] < 0.05 & metrics$BisulfiteConversion[i] < 90) {
      failed_tests <- c(failed_tests, "Bisulfite")
    }
    metrics$FailedTest[i] <- paste(failed_tests, collapse = ", ")
  }

  metrics_filepath <- file.path(context$paths$qc, "metrics_pca_outlier_assessment.rds")
  prog$update(7, paste("Saving metrics under", metrics_filepath))
  saveRDS(metrics, metrics_filepath)

  prog$complete()
  rm(list = ls())
  gc(full = T)
}
