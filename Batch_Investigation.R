demographics <- read.csv("./ppmi/Demographics_24Nov2025.csv")

participant_status <- read.csv("./ppmi/Participant_Status_24Nov2025.csv")

participant_status <- participant_status %>% filter(!ENROLL_STATUS %in% c("Screen failed", "Withdrew", "Withdrew Deceased", "Baseline Withdraw", "Excluded", "Declined", "Withdraw Deceased"))
dim(participant_status)

patient_meta <- merge(
  x = demographics[, c("PATNO", "SEX")],
  y = participant_status[, c("ENROLL_STATUS", "PATNO", "COHORT", "ENROLL_AGE", "ENRLPINK1", "ENRLPRKN", "ENRLSRDC", "ENRLNORM", "ENRLOTHGV", "ENRLHPSM", "ENRLLRRK2", "ENRLRBD", "ENRLSNCA", "ENRLGBA")],
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = FALSE
)

View(patient_meta)
dim(patient_meta)

ppmi_meth_120_txt <- read.delim("./ppmi/Project120_IDATS_n524final_toLONI_030718/PPMI_Meth_n524_for_LONI030718.txt", header = T, sep = "\t")
dim(ppmi_meth_120_txt)

ppmi_meth_120_meta <- merge(
  x = patient_meta,
  y = ppmi_meth_120_txt,
  by.x = "PATNO",
  by.y = "PATNO",
  all.x = FALSE
)
View(ppmi_meth_120_meta)

ppmi_meth_120_meta$Basename <- paste0(ppmi_meth_120_meta$Sentrix.ID, "_", ppmi_meth_120_meta$Sentrix.Position)
ppmi_meth_120_meta$Sample_Name <- paste0(ppmi_meth_120_meta$PATNO, "_", ppmi_meth_120_meta$Sentrix.ID, "_", ppmi_meth_120_meta$Sentrix.Position)
ppmi_meth_120_meta$Sample_Group <- ppmi_meth_120_meta$COHORT
ppmi_meth_120_meta <- ppmi_meth_120_meta %>% dplyr::rename(Sentrix_Position = Sentrix.Position, Sentrix_ID = Sentrix.ID)

View(ppmi_meth_120_meta)
dim(ppmi_meth_120_meta)

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
dim(ppmi_meth_120_meta)

ppmi_scan_dates <- read.csv("./ppmi/ppmi_scan_dates.csv", header = TRUE)

dim(ppmi_scan_dates)
View(ppmi_scan_dates)

summary(duplicated(ppmi_scan_dates[c("SentrixID","ScanDate")]))
ppmi_scan_dates[c("SentrixID","ScanDate")]

duplicate_indices <- duplicated(ppmi_scan_dates[c("SentrixID","ScanDate")])

ppmi_scan_dates_harm <- ppmi_scan_dates[!duplicate_indices, ]

ppmi_meth_120_meta_scandates <- merge(
  x = ppmi_meth_120_meta,
  y = ppmi_scan_dates_harm,
  by.x = "Basename",
  by.y = "SentrixID",
  all = FALSE
)
dim(ppmi_meth_120_meta_scandates)
View(ppmi_meth_120_meta_scandates)

beta_matrix <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260415_170143/beta_matrix_bmiq.rds")
targets <- readRDS("/Volumes/Elements/vastai/ppmi/ppmi_20260415_170143/targets_remove_mismatch.rds")

pca <- prcomp(t(beta_matrix))

head(targets)
targets$Sentrix_ID <- paste0(targets$Slide,"_", targets$Array)

dim(targets)
targets_harm <- merge(
  x = targets,
  y = ppmi_scan_dates_harm,
  by.x = "Sentrix_ID",
  by.y = "SentrixID",
  all = FALSE
)
dim(targets_harm)

head(targets_harm)

pca_df <- data.frame(matrix(NA, nrow = ncol(beta_matrix), ncol = 10))
rownames(pca_df) <- colnames(beta_matrix)
colnames(pca_df) <- sapply(1:10, function(i) paste0("PC", i))
for (index in 1:10) {
    pca_df[[paste0("PC", index)]] <- pca$x[, index]
}

library(ggplot2)

targets_harm$Sentrix_ID <- paste0(targets_harm$Slide,"_", targets_harm$Array)

pca_df_with_dates <- merge(
  x = pca_df,
  y = targets_harm,
  by.x = "row.names",
  by.y = "Sentrix_ID",
  all = FALSE
)

head(pca_df_with_dates)

ggplot(pca_df_with_dates, aes(x = PC1, y = PC2, color = ScanDate)) +
  geom_point() +
  labs(title = "PCA of Methylation Data Colored by Scan Date") +
  theme_minimal()

ggplot(pca_df_with_dates, aes(x = PC1, y = PC3, color = ScanDate)) +
  geom_point() +
  labs(title = "PCA of Methylation Data Colored by Scan Date") +
  theme_minimal()

head(targets_harm)

library(stringr)
library(patchwork)

pca_df_with_dates$ScanDate <- pca_df_with_dates$ScanDate %>% str_extract("^(\\d{4}-\\d{2})")

color_by <- "ScanDate"
pca_pairplot <- function(pca_df_with_dates, color_by = NULL) {
  n_pcs <- 6
  pc_columns <- grep("^PC\\d+$", colnames(pca_df_with_dates), value = TRUE)
  pc_columns <- pc_columns[1:min(n_pcs, length(pc_columns))]
  combos <- expand.grid(x = pc_columns, y = pc_columns)
  plots <- list()
  plot_idx <- 1

  for (i in 1:nrow(combos)) {
      x_var <- as.character(combos$x[i])
      y_var <- as.character(combos$y[i])

      if (x_var == y_var) {
        # Diagonal: Density plot
        p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]])) +
          geom_density(fill = "steelblue", alpha = 0.5) +
          labs(x = x_var, y = "Density") +
          theme_minimal() +
          theme(axis.text = element_text(size = 6))
      } else {
        # Off-diagonal: Scatter plot
        if (!is.null(color_by) && color_by %in% colnames(pca_df_with_dates)) {
          p <- ggplot(pca_df_with_dates, aes(
            x = .data[[x_var]], y = .data[[y_var]],
            color = as.factor(.data[[color_by]])
          )) +
            geom_point(size = 0.5, alpha = 0.5) +
            scale_color_viridis_d(name = color_by)
        } else {
          p <- ggplot(pca_df_with_dates, aes(x = .data[[x_var]], y = .data[[y_var]])) +
            geom_point(size = 0.5, alpha = 0.5, color = "steelblue")
        }
        p <- p +
          labs(x = x_var, y = y_var) +
          theme_minimal() +
          theme(
            axis.text = element_text(size = 6),
            legend.position = "right"
          )
      }

      plots[[plot_idx]] <- p
      plot_idx <- plot_idx + 1
  }

  pplot <- wrap_plots(plots, ncol = length(pc_columns)) +
      plot_annotation(title = paste("PCA Pairplot - First", length(pc_columns), "PCs"))

  ggsave(paste0("./ppmi/ppmi_scan_date_pca_pairplot_", color_by, ".png"), pplot, width = 3 * length(pc_columns), height = 3 * length(pc_columns), dpi = 300)

