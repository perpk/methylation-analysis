rm(list = ls())
gc(full = TRUE)

csv_paths <- c(
    "/Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/processed/SGPD_DMR_Effect_Sizes_corrected.csv"
)

# Load and combine all CSV files
data_list <- lapply(csv_paths, read.csv)

# Add dataset identifier and standardize column names
shape_map <- c("PPMI" = 16, "PEG1" = 17, "SGPD" = 18)
for (i in seq_along(csv_paths)) {
    data_list[[i]]$dataset <- basename(dirname(csv_paths[i]))
    # Extract cohort prefix from CSV filename
    filename <- basename(csv_paths[i])
    cohort <- ifelse(grepl("PPMI", filename), "PPMI",
                 ifelse(grepl("PEG1", filename), "PEG1",
                    ifelse(grepl("SGPD", filename), "SGPD", NA)))
    data_list[[i]]$cohort <- cohort
    data_list[[i]]$shape <- shape_map[cohort]
    # Standardize column names to lowercase
    names(data_list[[i]]) <- tolower(names(data_list[[i]]))
}

# Ensure all dataframes have the same columns before combining
all_cols <- Reduce(union, lapply(data_list, names))
data_list <- lapply(data_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) {
        df[[col]] <- NA
    }
    df[, all_cols]
})

data_combined <- do.call(rbind, data_list)

# Prepare data for plotting
data_combined$color <- ifelse(data_combined$deltabeta < 0, "blue", "red")
data_combined$size <- abs(data_combined$deltabeta)
plot_dir <- file.path(getwd(), "documentation", "stats", "plots")

# Create scatterplot per cohort
cohorts <- unique(data_combined$cohort)
for (cohort in cohorts) {
    cohort_data <- data_combined[data_combined$cohort == cohort, ]
    
    png(file.path(plot_dir, paste0("chromosome_effect_sizes_", cohort, ".png")), width = 14, height = 8, units = "in", res = 300)
    par(mar = c(5, 5, 4, 2))
    
    plot(as.numeric(factor(cohort_data$chromosome)),
        seq_len(nrow(cohort_data)),
        col = cohort_data$color,
        cex = 1 + 4 * (cohort_data$size / max(cohort_data$size)),
        pch = 16,
        xlab = "Chromosome",
        ylab = "Region Index",
        main = paste("DMR Effect Sizes Distribution by Chromosome -", cohort),
        xaxt = "n"
    )
    
    axis(1,
        at = 1:length(unique(cohort_data$chromosome)),
        labels = sort(unique(cohort_data$chromosome))
    )
    
    legend("topright",
        legend = c("Δβ < 0", "Δβ > 0"),
        col = c("blue", "red"),
        pch = 16,
        title = "Effect Direction"
    )
    
    dev.off()
}

# Summary statistics by chromosome
summary_stats <- aggregate(deltabeta ~ chromosome,
    data = data_combined,
    FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x))
)
print(summary_stats)
