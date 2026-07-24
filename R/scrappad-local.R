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
