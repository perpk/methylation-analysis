source("R/progress_mgr.R")
source("R/project_context.R")
source("R/intermediate_data_proxy.R")
source("R/results_container.R")

projects_location <- "/Volumes/saucepan/methylation-project"
ppmi_project <- "ppmi_20260721_075730"


ppmi_targets <- readRDS(file.path(paste0(projects_location, "/", ppmi_project, "/processed"), "targets_after_cell_count_estimation.rds"))
gse111629_targets <- readRDS(file.path(paste0(projects_location, "/", "GSE111629_20260722_083339", "/processed"), "targets_after_cell_count_estimation.rds"))
gse145361_targets <- readRDS(file.path(paste0(projects_location, "/", "GSE145361_20260714_080452", "/processed"), "targets.rds"))

library(tidyverse)

ppmi_targets %>% head()

# --- Harmonize Sex / Sample_Group / Age_Group across cohorts for summary barplots ---

gse111629_targets <- gse111629_targets %>%
    mutate(
        Age_Group = case_when(
            as.numeric(`age:ch1`) >= 30 & as.numeric(`age:ch1`) < 50 ~ "30-49",
            as.numeric(`age:ch1`) >= 50 & as.numeric(`age:ch1`) < 60 ~ "50-59",
            as.numeric(`age:ch1`) >= 60 & as.numeric(`age:ch1`) < 70 ~ "60-69",
            as.numeric(`age:ch1`) >= 70 & as.numeric(`age:ch1`) < 80 ~ "70-79",
            as.numeric(`age:ch1`) >= 80 ~ "80+"
        )
    )

gse145361_targets$Age_Group <- NA_character_

age_group_levels <- c("30-49", "50-59", "60-69", "70-79", "80+")

cohort_summary <- bind_rows(
    ppmi_targets %>% transmute(Cohort = "PPMI", Sex = SEX, Sample_Group = Sample_Group, Age_Group = Age_Group),
    gse111629_targets %>% transmute(Cohort = "GSE111629", Sex = `gender:ch1`, Sample_Group = Sample_Group, Age_Group = Age_Group),
    gse145361_targets %>% transmute(Cohort = "GSE145361", Sex = `gender:ch1`, Sample_Group = Sample_Group, Age_Group = Age_Group)
) %>%
    mutate(
        Cohort = factor(Cohort, levels = c("PPMI", "GSE111629", "GSE145361")),
        Age_Group = factor(Age_Group, levels = age_group_levels)
    )

# 1. Total number of samples per cohort
sample_counts <- cohort_summary %>% count(Cohort, name = "n_samples")

total_samples_barplot <- ggplot(sample_counts, aes(x = Cohort, y = n_samples, fill = Cohort)) +
    geom_col() +
    geom_text(aes(label = n_samples), vjust = -0.4) +
    labs(title = "Total Number of Samples by Cohort", x = "Cohort", y = "Number of Samples") +
    theme_minimal()

ggsave(paste0(projects_location, "/cohort_total_samples_barplot.png"), total_samples_barplot, width = 7, height = 5, dpi = 300)

# 2. Distribution of Age_Group per cohort
age_group_barplot <- ggplot(
    cohort_summary %>% filter(!is.na(Age_Group)),
    aes(x = Age_Group, fill = Cohort)
) +
    geom_bar(position = "dodge") +
    labs(title = "Distribution of Age Group by Cohort", x = "Age Group", y = "Number of Samples") +
    theme_minimal()
ggsave(paste0(projects_location, "/cohort_age_group_barplot.png"), age_group_barplot, width = 7, height = 5, dpi = 300)

# 3. Distribution of Sample_Group per cohort
sample_group_barplot <- ggplot(cohort_summary, aes(x = Sample_Group, fill = Cohort)) +
    geom_bar(position = "dodge") +
    labs(title = "Distribution of Sample Group by Cohort", x = "Sample Group", y = "Number of Samples") +
    theme_minimal()

ggsave(paste0(projects_location, "/cohort_sample_group_barplot.png"), sample_group_barplot, width = 7, height = 5, dpi = 300)

# 4. Distribution of Sex per cohort
sex_barplot <- ggplot(cohort_summary, aes(x = Sex, fill = Cohort)) +
    geom_bar(position = "dodge") +
    labs(title = "Distribution of Sex by Cohort", x = "Sex", y = "Number of Samples") +
    theme_minimal()

ggsave(paste0(projects_location, "/cohort_sex_barplot.png"), sex_barplot, width = 7, height = 5, dpi = 300)

# ========================================================================================================================================================================

### Sample QC Summary Statistics

rm(list = ls())
gc(full = TRUE)

methyl_set_before_probe_qc <- readRDS("/Volumes/saucepan/methylation-project/ppmi_20260721_075730/processed/methyl_set_remove_mismatch.rds")
methyl_set_after_probe_qc <- readRDS("/Volumes/saucepan/methylation-project/ppmi_20260721_075730/processed/methyl_set_remove_probe_qc.rds")

dim(methyl_set_before_probe_qc)
dim(methyl_set_after_probe_qc)

rg_set_before_probe_qc <- readRDS("/Volumes/saucepan/methylation-project/ppmi_20260721_075730/processed/rg_set_remove_mismatch.rds")
rg_set_after_probe_qc <- readRDS("/Volumes/saucepan/methylation-project/ppmi_20260721_075730/processed/rg_set_remove_probe_qc.rds")

dim(rg_set_before_probe_qc)
dim(rg_set_after_probe_qc)
1052641 - 1043999
866836 - 858194

rg_set_after_probe_qc %>%
    rownames() %>%
    head()

methyl_set_after_probe_qc %>%
    rownames() %>%
    head()

library(minfi)
ann <- getAnnotation(rg_set_after_probe_qc)

probe_det_p <- readRDS("/Volumes/saucepan/methylation-project/ppmi_20260721_075730/qc/removed_probes_detection_p.rds")

probe_det_p %>%
    head()

# ========================================================================================================================================================================
# Sample QC Summary Statistics
library(tidyverse)

rm(list = ls())
gc(full = TRUE)

projects_location <- "/Volumes/saucepan/methylation-project"
ppmi_project <- "ppmi_20260721_075730"

ppmi_targets <- readRDS(file.path(paste0(projects_location, "/", ppmi_project, "/processed"), "targets_after_cell_count_estimation.rds"))

qc_results <- readRDS(file.path(paste0(projects_location, "/", ppmi_project, "/qc"), "qc_metrics.rds"))

qc_results %>%
    dim()

failed_samples_gse111629 <- readRDS(file.path(paste0(projects_location, "/", "GSE111629_20260722_083339", "/qc"), "failed_samples.rds"))
failed_samples_gse145361 <- readRDS(file.path(paste0(projects_location, "/", "GSE145361_20260714_080452", "/qc"), "failed_samples.rds"))
failed_samples_ppmi <- readRDS(file.path(paste0(projects_location, "/", ppmi_project, "/qc"), "failed_samples.rds"))

failed_samples_ppmi$Beadcount_Failure == FALSE %>% count()
failed_samples_ppmi %>% dim()

failed_samples_gse111629$Beadcount_Failure %>% table()
failed_samples_gse111629 %>% dim()
failed_samples_gse145361$Beadcount_Failure %>% table()
failed_samples_gse145361 %>% dim()
# --- Failed samples per cohort by failure category (Intensity_Failure, Bisulfite_Failure) ---

failed_samples_summary <- bind_rows(
    failed_samples_ppmi %>% mutate(Cohort = "PPMI"),
    failed_samples_gse111629 %>% mutate(Cohort = "GSE111629"),
    failed_samples_gse145361 %>% mutate(Cohort = "GSE145361")
) %>%
    select(Cohort, Intensity_Failure, Bisulfite_Failure, Beadcount_Failure) %>%
    pivot_longer(
        cols = c(Intensity_Failure, Bisulfite_Failure, Beadcount_Failure),
        names_to = "Failure_Category",
        values_to = "Failed"
    ) %>%
    mutate(
        Failure_Category = factor(
            Failure_Category,
            levels = c("Intensity_Failure", "Bisulfite_Failure", "Beadcount_Failure")
        )
    ) %>%
    filter(Failed) %>%
    count(Cohort, Failure_Category, name = "n_failed") %>%
    complete(Cohort, Failure_Category, fill = list(n_failed = 0)) %>%
    mutate(Cohort = factor(Cohort, levels = c("PPMI", "GSE111629", "GSE145361")))

failed_samples_summary$Failure_Category == "Beadcount_Failure"

failed_samples_barplot <- ggplot(failed_samples_summary, aes(x = Cohort, y = n_failed, fill = Failure_Category)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n_failed), position = position_dodge(width = 0.9), vjust = -0.4) +
    labs(
        title = "Failed Samples by Cohort and Failure Category",
        x = "Cohort", y = "Number of Failed Samples", fill = "Failure Category"
    ) +
    theme_minimal()

ggsave(paste0(projects_location, "/cohort_failed_samples_barplot.png"), failed_samples_barplot, width = 7, height = 5, dpi = 300)

View(failed_samples_summary)
