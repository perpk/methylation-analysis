library(dplyr)
library(minfi)
library(dplyr)

data_folder <- "ppmi"

demographics <- read.csv("./ppmi/Demographics_24Nov2025.csv")

participant_status <- read.csv("./ppmi/Participant_Status_24Nov2025.csv")


participant_status <- participant_status %>% filter(!ENROLL_STATUS %in% c("Excluded", "Declined", "Screen_failed"))
dim(participant_status)

patient_meta <- merge(
  x = demographics[, c("PATNO", "SEX")],
  y = participant_status[, c("ENROLL_STATUS", "PATNO", "COHORT", "ENROLL_AGE", "ENRLPINK1", "ENRLPRKN", "ENRLSRDC", "ENRLNORM", "ENRLOTHGV", "ENRLHPSM", "ENRLLRRK2", "ENRLRBD", "ENRLSNCA", "ENRLGBA")],
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = FALSE
)
dim(patient_meta)


ppmi_meth_120_txt <- read.delim("./ppmi/Project120_IDATS_n524final_toLONI_030718/PPMI_Meth_n524_for_LONI030718.txt", header = T, sep = "\t")

ppmi_meth_120_meta <- merge(
  x = patient_meta,
  y = ppmi_meth_120_txt,
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = FALSE
)
dim(ppmi_meth_120_meta)

ppmi_meth_120_meta$Basename <- paste0(ppmi_meth_120_meta$Sentrix.ID, "_", ppmi_meth_120_meta$Sentrix.Position)
ppmi_meth_120_meta$Sample_Name <- paste0(ppmi_meth_120_meta$PATNO, "_", ppmi_meth_120_meta$Sentrix.ID, "_", ppmi_meth_120_meta$Sentrix.Position)
ppmi_meth_120_meta$Sample_Group <- ppmi_meth_120_meta$COHORT
ppmi_meth_120_meta <- ppmi_meth_120_meta %>% dplyr::rename(Sentrix_Position = Sentrix.Position, Sentrix_ID = Sentrix.ID)

# 0 -> F; 1 -> M
ppmi_meth_120_meta$SEX <- ifelse(ppmi_meth_120_meta$SEX == 0, "Female", "Male")

# 1 -> PD; 2 -> Control; 3 -> SWEDD; 4 -> Prodromal
ppmi_meth_120_meta <- ppmi_meth_120_meta %>%
  mutate(Sample_Group = case_when(
    Sample_Group == 1 ~ "PD",
    Sample_Group == 2 ~ "Control",
    Sample_Group == 3 ~ "SWEDD",
    Sample_Group == 4 ~ "Prodromal"
  ))

ppmi_meth_120_meta <- ppmi_meth_120_meta %>%
  mutate(
    Age_Group =
      case_when(
        ENROLL_AGE >= 30 & ENROLL_AGE < 50 ~ "30-49",
        ENROLL_AGE >= 50 & ENROLL_AGE < 60 ~ "50-59",
        ENROLL_AGE >= 60 & ENROLL_AGE < 70 ~ "60-69",
        ENROLL_AGE >= 70 & ENROLL_AGE < 80 ~ "70-79",
        ENROLL_AGE >= 80 ~ "80+",
      )
  )

write.csv(ppmi_meth_120_meta, "./ppmi/Project120_IDATS_n524final_toLONI_030718/sample_sheet.csv", row.names = F)

rm(list = setdiff(ls(), c("data_folder")))
gc(full = TRUE)

project_name <- "ppmi"
data_folder <- "ppmi"

targets <- read.metharray.sheet(data_folder, pattern = "sample_sheet")

source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = project_name)

source("./pre_process_eda.R")
pre_process_eda(
  project_name = project_name,
  targets = targets,
  data_folder = data_folder,
  project_location = "/Volumes/Elements/methyl-pipe-out",
  var_mapping = var_mappings,
  platform = "EPIC",
  qc_threshold = "auto"
)

project_to_load <- "ppmi_20260103_120001"

source("./pre_process_eda.R")
pre_process_eda(
  project_to_load = project_to_load,
  targets = targets,
  data_folder = data_folder,
  project_location = "/Volumes/Elements/methyl-pipe-out",
  var_mapping = var_mappings,
  platform = "EPIC",
  qc_threshold = "auto"
)
