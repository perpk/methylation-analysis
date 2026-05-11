library(dplyr)
library(minfi)
library(dplyr)

data_folder <- "ppmi"

demographics <- read.csv("./ppmi/Demographics_24Nov2025.csv")

participant_status <- read.csv("./ppmi/Participant_Status_24Nov2025.csv")

participant_status <- participant_status %>% filter(!ENROLL_STATUS %in% c("Screen failed", "Withdrew", "Withdrew Deceased", "Baseline Withdraw", "Excluded", "Declined"))

patient_meta <- merge(
  x = demographics[, c("PATNO", "SEX")],
  y = participant_status[, c("ENROLL_STATUS", "PATNO", "COHORT", "ENROLL_AGE", "ENRLPINK1", "ENRLPRKN", "ENRLSRDC", "ENRLNORM", "ENRLOTHGV", "ENRLHPSM", "ENRLLRRK2", "ENRLRBD", "ENRLSNCA", "ENRLGBA")],
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = TRUE
)

# family_history <- read.csv("./ppmi/condition/Family_History_24Nov2025.csv")
# multi_patno <- family_history %>% group_by(PATNO) %>% summarise(count = n()) %>% filter(count > 1)
# single_patno <- family_history %>% group_by(PATNO) %>% summarise(count = n()) %>% filter(count == 1)

# trans_rows <- family_history[family_history$PATNO %in% multi_patno$PATNO, ] %>% filter(EVENT_ID == "TRANS")
# single_rows <- family_history[family_history$PATNO %in% single_patno$PATNO, ]

# clean_fam_hist <- rbind(trans_rows, single_rows)

# patient_meta <- merge(
#   x = patient_meta,
#   y = clean_fam_hist[, c("PATNO", "ANYFAMPD")],
#   by.x = "PATNO",
#   by.y = "PATNO",
#   all.x = TRUE
# )

# patient_meta <- patient_meta %>% filter(!PATNO %in% c("141081","140041"))

ppmi_meth_120_txt <- read.delim("./ppmi/Project120_IDATS_n524final_toLONI_030718/PPMI_Meth_n524_for_LONI030718.txt", header = T, sep = "\t")
View(ppmi_meth_120_txt)
ppmi_meth_120_meta <- merge(
  x = patient_meta,
  y = ppmi_meth_120_txt,
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = FALSE
)

colnames(ppmi_meth_120_meta)
View(ppmi_meth_120_meta)


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
# saveRDS(targets, file = file.path(data_folder, "targets.rds"))
# targets <- readRDS(file.path(data_folder, "targets.rds"))
source("./meta_vars_mapping.R")
var_mappings <- meta_vars_mapping(dataset = project_name)

cohorts <- list(
  PD_vs_Control = c("PD", "Control"),
  SWEDD_vs_Control = c("SWEDD", "Control")
)
project_to_load <- "ppmi"

head(targets)

targets_pd_lrrk2 <- targets %>% filter(Sample_Group == "PD" & ENRLLRRK2 == 1)

nrow(targets[targets$ANYFAMPD == 1, ])

source("./pre_process_eda.R")
pre_process_eda(
  project_to_load = project_to_load,
  targets = targets,
  data_folder = data_folder,
  project_location = "/root/workspace/methyl-pipe-out",
  var_mapping = var_mappings,
  platform = "EPIC",
  qc_threshold = "auto",
  cohorts = cohorts
)

rootDir <- "/Users/kpax/Documents/study/temp"
library(dplyr)

beta_matrix <- readRDS(file.path(rootDir, "beta_matrix_bmiq.rds"))
targets <- readRDS(file.path(rootDir, "targets_remove_mismatch.rds"))
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

write.csv(beta_means, file.path(rootDir, "beta_means.csv"))


###

library(limma)

idat_file <- "./ppmi/Project120_IDATS_n524final_toLONI_030718/200973410121_R01C01_Grn.idat"
manifest_file <- "./ppmi/Project120_IDATS_n524final_toLONI_030718/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
idat <- read.idat(c(idat_file), bgxfile = manifest_file, dateinfo=TRUE, skip=7)

rgset <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260415_170143/rg_set_ppmi.rds")
annotation_data <- pData(rgset)
View(annotation_data)

library(illuminaio)
data <- readIDAT(idat_file)
print(data$RunInfo)


run_metadata <- data$RunInfo
scan_row <- which(run_metadata[, "BlockType"] == "Scan")[1]
scan_date_string <- run_metadata[scan_row, "RunTime"]
scan_date <- as.POSIXct(scan_date_string, format="%m/%d/%Y")
print(scan_date)
data$Unknowns$MostlyA
scan_date <- run_metadata[run_metadata[, "BlockType"] == "Scan", "RunTime"][1]
scan_date
