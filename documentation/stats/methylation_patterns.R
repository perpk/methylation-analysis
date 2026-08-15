# /Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv
# /Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv
# /Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/results/SGPD_DMR_Effect_Sizes.csv

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
    library(readr)
    library(ggplot2)
})

csv_paths <- c(
    "/Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/processed/SGPD_DMR_Effect_Sizes_corrected.csv"
)

# Saves plot under: <project-root>/documentation/stats/plots
plot_dir <- file.path(getwd(), "documentation", "stats", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

all_counts <- data.frame(
    Cohort = character(),
    Sign = character(),
    Count = numeric(),
    stringsAsFactors = FALSE
)

for (path in csv_paths) {
    if (!file.exists(path)) {
        warning(sprintf("File not found: %s", path))
        next
    }

    cohort <- strsplit(basename(path), "_")[[1]][1]
    if (is.na(cohort) || cohort == "") {
        warning(sprintf("Could not parse cohort from filename: %s", path))
        next
    }

    dat <- read_csv(path, show_col_types = FALSE)

    # fixed column name as requested
    if (!("deltabeta" %in% names(dat))) {
        warning(sprintf("Column 'deltabeta' not found in: %s", path))
        next
    }

    delta <- suppressWarnings(as.numeric(dat[["deltabeta"]]))
    # delta <- delta[!is.na(delta)]

    pos_n <- sum(delta > 0, na.rm = TRUE)
    neg_n <- sum(delta < 0, na.rm = TRUE)

    all_counts <- rbind(
        all_counts,
        data.frame(Cohort = cohort, Sign = "Hypomethylated", Count = neg_n, stringsAsFactors = FALSE),
        data.frame(Cohort = cohort, Sign = "Hypermethylated", Count = pos_n, stringsAsFactors = FALSE)
    )
}

# Remove any malformed rows that can create NA bars
all_counts <- all_counts[complete.cases(all_counts[, c("Cohort", "Sign", "Count")]), , drop = FALSE]
all_counts <- all_counts[all_counts$Cohort != "" & all_counts$Sign != "", , drop = FALSE]

all_counts$Sign <- factor(all_counts$Sign, levels = c("Hypomethylated", "Hypermethylated"))

p <- ggplot(all_counts, aes(x = Cohort, y = Count, fill = Sign)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
        aes(label = Count),
        position = position_dodge(width = 0.75),
        vjust = -0.25,
        size = 3.5
    ) +
    scale_fill_manual(
        values = c("Hypomethylated" = "#D55E00", "Hypermethylated" = "#009E73"),
        na.translate = FALSE
    ) +
    labs(
        title = "Delta-beta sign counts by cohort",
        x = "Cohort",
        y = "Number of DMRs",
        fill = NULL
    ) +
    theme_minimal(base_size = 12)

print(p)

ggsave(
    filename = file.path(plot_dir, "cohorts_delta_beta_sign_counts.png"),
    plot = p,
    width = 9,
    height = 5,
    dpi = 300
)

## === Percentages

# /Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv
# /Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv
# /Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/results/SGPD_DMR_Effect_Sizes.csv

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
    library(readr)
    library(ggplot2)
    library(dplyr)
})

csv_paths <- c(
    "/Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv",
    "/Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/processed/SGPD_DMR_Effect_Sizes_corrected.csv"
)

# Create pie charts for percentages
all_pie_data <- data.frame(
    Cohort = character(),
    Sign = character(),
    Count = numeric(),
    stringsAsFactors = FALSE
)

for (path in csv_paths) {
    if (!file.exists(path)) {
        warning(sprintf("File not found: %s", path))
        next
    }

    cohort <- strsplit(basename(path), "_")[[1]][1]
    if (is.na(cohort) || cohort == "") {
        warning(sprintf("Could not parse cohort from filename: %s", path))
        next
    }

    dat <- read_csv(path, show_col_types = FALSE)

    if (!("deltabeta" %in% names(dat))) {
        warning(sprintf("Column 'deltabeta' not found in: %s", path))
        next
    }

    delta <- suppressWarnings(as.numeric(dat[["deltabeta"]]))

    pos_n <- sum(delta > 0, na.rm = TRUE)
    neg_n <- sum(delta < 0, na.rm = TRUE)

    all_pie_data <- rbind(
        all_pie_data,
        data.frame(Cohort = cohort, Sign = "Hypomethylated", Count = neg_n, stringsAsFactors = FALSE),
        data.frame(Cohort = cohort, Sign = "Hypermethylated", Count = pos_n, stringsAsFactors = FALSE)
    )
}

all_pie_data <- all_pie_data[complete.cases(all_pie_data[, c("Cohort", "Sign", "Count")]), , drop = FALSE]
all_pie_data <- all_pie_data[all_pie_data$Cohort != "" & all_pie_data$Sign != "", , drop = FALSE]
all_pie_data$Sign <- factor(all_pie_data$Sign, levels = c("Hypomethylated", "Hypermethylated"))

# Create pie charts
n_cohorts <- length(unique(all_pie_data$Cohort))
png(
    filename = file.path(plot_dir, "cohorts_delta_beta_percentages.png"),
    width = 4 * n_cohorts,
    height = 5,
    units = "in",
    res = 300
)

par(mfrow = c(1, n_cohorts), mar = c(2, 2, 3, 2))

cohorts <- unique(all_pie_data$Cohort)
for (cohort in cohorts) {
    subset_data <- all_pie_data[all_pie_data$Cohort == cohort, ]
    counts <- subset_data$Count
    labels <- subset_data$Sign

    total <- sum(counts)
    percentages <- (counts / total) * 100

    pie(counts,
        labels = sprintf("%s\n(%.1f%%)", labels, percentages),
        main = cohort,
        col = c("Hypomethylated" = "#D55E00", "Hypermethylated" = "#009E73")[as.character(labels)],
        cex = 1
    )
}

dev.off()

print("Pie charts saved.")
