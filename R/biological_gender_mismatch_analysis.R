biological_gender_mismatch_analysis <- function(context=NULL,
                                                methyl_set_filename = "methyl_set_clean.rds",
                                                rg_set_filename = "rg_set_clean.rds",
                                                targets_filename="targets_clean.rds",
                                                targets_sample_col="Sample_Name",
                                                recorded_sex_col=NULL) {
  prog <- .create_progress_manager(5)

  methyl_set_file <- file.path(context$paths$qc, methyl_set_filename)
  rg_set_file <- file.path(context$paths$qc, rg_set_filename)
  prog$update(1, paste("Reading Methylset from filesystem", methyl_set_file))
  methyl_set <- readRDS(methyl_set_file)
  rg_set <- readRDS(rg_set_file)

  prog$update(2, "Map Methylation Data to the Genome")
  methyl_set_genomic <- mapToGenome(methyl_set)

  prog$update(3, "Predict Sex")
  predicted_sex <- getSex(methyl_set_genomic, cutoff = -2)

  targets_filepath <- file.path(context$paths$qc, targets_filename)
  targets <- readRDS(targets_filepath)

  prog$update(4, "Map Methylation Data to the Genome")
  sex_check <- data.frame(
    Sample_Name = targets$Sample_Name,
    Recorded_Sex = targets[[recorded_sex_col]],
    Predicted_Sex = predicted_sex$predictedSex,
    X_chr_median = predicted_sex$xMed,
    Y_chr_median = predicted_sex$yMed
  )

  prog$update(5, "Fetch mismatches and log them persistently")
  sex_check <- sex_check %>% mutate(
    Predicted_Sex = case_when(Predicted_Sex == "M" ~ "Male", Predicted_Sex == "F" ~ "Female")
  )

  mismatches <- sex_check[sex_check$Recorded_Sex != sex_check$Predicted_Sex, ]
  print(mismatches)
  sex_mismatch_log <- file.path(context$paths$results, "sex_mismatch_log.csv")
  print(paste("Writing mismatch results to file", sex_mismatch_log))
  write.csv(mismatches, sex_mismatch_log)

  mismatch_plot <- ggplot(sex_check, aes(x = X_chr_median, y = Y_chr_median,
                        color = Recorded_Sex,
                        shape = Predicted_Sex)) +
    geom_point(aes(shape = Predicted_Sex), size = 3, alpha = 0.5) +  scale_color_manual(values = c("Female" = "red", "Male" = "blue")) +
    labs(x = "X Chromosome Median Intensity",
         y = "Y Chromosome Median Intensity",
         title = "Sex Prediction Quality Control",
         color = "Recorded Sex",
         shape = "Predicted Sex") +
    theme_minimal() +
    theme(legend.position = "bottom")

  mismatch_plot_file <- file.path(context$paths$plots, "sex_mismatch_plot.png")
  ggsave(filename = mismatch_plot_file,
         plot = mismatch_plot,
         width = 8,
         height = 6,
         dpi = 300)
  print(paste("Saved:", mismatch_plot_file))
  print(paste("Found", nrow(mismatches), "sex mismatches:"))
  print(mismatches)

  s <- colnames(methyl_set)
  indices <- which(s %in% rownames(mismatches))
  methyl_set_remove_mismatch <- methyl_set[, -indices]
  rg_set_remove_mismatch <- rg_set[, -indices]

  s <- targets[[targets_sample_col]]
  indices <- which(s %in% mismatches[[targets_sample_col]])
  targets_remove_mismatch <- targets[-indices, ]

  saveRDS(methyl_set_remove_mismatch, file.path(context$paths$processed, "methyl_set_remove_mismatch.rds"))
  saveRDS(rg_set_remove_mismatch, file.path(context$paths$processed, "rg_set_remove_mismatch.rds"))
  saveRDS(targets_remove_mismatch, file.path(context$paths$processed, "targets_remove_mismatch.rds"))

  prog$complete()
}
