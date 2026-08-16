# ============================================================
# Cohort x Chromosome DeltaBeta Heatmap
# ============================================================
# Reads DMR effect-size CSVs from multiple cohorts, extracts
# chromosome + deltabeta, aggregates per chromosome per cohort,
# and plots a heatmap (cohorts vs chromosomes, fill = deltabeta)
# using a viridis color gradient.
# ============================================================

# ---- 0. Packages -------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tidyr", "stringr", "ggplot2", "viridis", "forcats")
missing_pkgs  <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(forcats)

# ---- 1. Inputs -----------------------------------------------------
csv_paths <- c(
  "/Volumes/saucepan/methylation-project/ppmi_20260721_075730/results/PPMI_DMR_effect_sizes.csv",
  "/Volumes/saucepan/methylation-project/GSE111629_20260722_083339/results/PEG1_DMR_Effect_Sizes.csv",
  "/Volumes/saucepan/methylation-project/GSE145361_20260804_162813/GSE145361_20260804_162813/processed/SGPD_DMR_Effect_Sizes_corrected.csv"
)

# How to collapse multiple DMRs per chromosome into one value per cell.
# Change to "median", "sum", "max", or "min" if preferred.
agg_fun <- "mean"

# ---- 2. Helper: flexible column finder ------------------------------
# Finds a column whose name matches a pattern, ignoring case/underscores.
find_col <- function(df, patterns) {
  clean_names <- names(df) |> str_to_lower() |> str_remove_all("[^a-z0-9]")
  for (p in patterns) {
    p_clean <- str_remove_all(str_to_lower(p), "[^a-z0-9]")
    hit <- which(clean_names == p_clean)
    if (length(hit) > 0) return(names(df)[hit[1]])
    hit <- which(str_detect(clean_names, p_clean))
    if (length(hit) > 0) return(names(df)[hit[1]])
  }
  NA_character_
}

# ---- 3. Read + standardize each cohort file --------------------------
read_cohort_file <- function(path) {
  cohort <- str_split_fixed(basename(path), "_", 2)[, 1]

  df <- read_csv(path, show_col_types = FALSE)

  chr_col   <- find_col(df, c("chromosome", "chr", "seqnames"))
  delta_col <- find_col(df, c("deltabeta", "delta_beta", "dbeta", "meandiff", "effectsize"))

  if (is.na(chr_col))   stop(sprintf("Could not find a chromosome column in: %s", path))
  if (is.na(delta_col)) stop(sprintf("Could not find a deltabeta column in: %s", path))

  df |>
    transmute(
      cohort      = cohort,
      chromosome  = as.character(.data[[chr_col]]),
      deltabeta   = as.numeric(.data[[delta_col]])
    )
}

all_data <- bind_rows(lapply(csv_paths, read_cohort_file))

# ---- 4. Normalize chromosome labels & order ---------------------------
normalize_chr <- function(x) {
  x <- str_trim(x)
  x <- str_remove(x, regex("^chr", ignore_case = TRUE))
  x <- str_to_upper(x)
  paste0("chr", x)
}

chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

all_data <- all_data |>
  filter(!is.na(chromosome), !is.na(deltabeta)) |>
  mutate(
    chromosome = normalize_chr(chromosome),
    chromosome = factor(chromosome, levels = chr_levels)
  ) |>
  filter(!is.na(chromosome))  # drops anything not matching chr1-22/X/Y/M

# ---- 5. Aggregate: one value per cohort x chromosome -------------------
agg_data <- all_data |>
  group_by(cohort, chromosome) |>
  summarise(
    deltabeta = match.fun(agg_fun)(deltabeta, na.rm = TRUE),
    n_dmrs    = n(),
    .groups = "drop"
  ) |>
  complete(cohort, chromosome)  # fill missing combos with NA so grid is complete

# ---- 6. Plot heatmap ----------------------------------------------------
p <- ggplot(agg_data, aes(x = chromosome, y = fct_rev(cohort), fill = deltabeta)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(
    option    = "viridis",
    name      = paste0(str_to_title(agg_fun), " \u0394\u03B2"),
    na.value  = "grey90"
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Cohort \u00d7 Chromosome Δβ Association",
    x     = "Chromosome",
    y     = "Cohort"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

print(p)

# ---- 7. Save output -------------------------------------------------
out_dir <- "documentation/stats/plots"
if (!dir.exists(out_dir)) dir.create(out_dir)

ggsave(file.path(out_dir, "cohort_chromosome_deltabeta_heatmap.png"),
       plot = p, width = 10, height = 4.5, dpi = 300)

# also save the underlying aggregated table for reference
write_csv(agg_data, file.path(out_dir, "cohort_chromosome_deltabeta_summary.csv"))

message("Done. Heatmap saved to: ", file.path(out_dir, "cohort_chromosome_deltabeta_heatmap.png"))