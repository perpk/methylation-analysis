source('R/project_context.R')

project_context <- .load_methylation_project("/Volumes/Elements/methyl-pipe-out", "GSE111629_20251211_103624")

beta_val_filepath <- file.path(project_context$paths$processed, "beta_matrix.rds")
beta_matrix <- readRDS(beta_val_filepath)

beta_matrix_reduced <- beta_matrix[1:100,]

targets_filepath <- file.path(project_context$paths$qc, "targets_clean.rds")
targets <- readRDS(targets_filepath)

col_map <- list()
col_map[["Sample_Group"]] <- "Sample_Group"
col_map[["Gender"]] <- "gender:ch1"
col_map[["Age"]] <- "age:ch1"
keys <- c("Sample_Group", "Gender", "Age")

npc <- 5

pca <- prcomp(t(beta_matrix_reduced))

pca_df <- data.frame(matrix(NA, nrow = ncol(beta_matrix[1:100,]), ncol = npc))
colnames(pca_df) <- sapply(1:npc, function(pca_df) {
  (paste0("PC", pca_df))
})
for (index in 1:npc) {
  pca_df[[paste0("PC", index)]] <- pca$x[,index]
}

for (k in keys) {
  pca_df[[k]] <- targets[[col_map[[k]]]]
}

pca_vars <- c("Gender" = "By Biological Gender")

ggplot(pca_df, aes(x = PC1, y = PC2,
                   color = as.factor(.data[["Gender"]]))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "test",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Gender") +
  theme_minimal() +
  theme(legend.position = "right")


###


methyl_set_file <- file.path(project_context$paths$normalized, "methyl_set.rds")
methyl_set <- readRDS(methyl_set_file)

###

d_1 <- data.frame(c(10:100))
d_2 <- data.frame(c(10:100))

samples_to_remove <- c(10, 30, 60, 80)

datasets <- list(methyl_set = d_1, beta_matrix = d_2)
head(datasets$methyl_set)

results <- lapply(datasets, function(df) {
  df[-samples_to_remove, , drop = F]
})

###

bm <- readRDS(file.path(project_context$paths$processed, "beta_matrix.rds"))
sex_mismatch_log <- file.path(project_context$paths$results, "sex_mismatch_log.csv")
mismatches <- read.csv(sex_mismatch_log, row.names=1)

samples_to_remove <- rownames(mismatches)
s <- colnames(bm)
indices <- which(s %in% samples_to_remove)
indices
clean <- bm[,-indices]
dim(clean)
dim(bm)


##

t <- readRDS(file.path(project_context$paths$qc, "targets_clean.rds"))
sex_mismatch_log <- file.path(project_context$paths$results, "sex_mismatch_log.csv")
mismatches <- read.csv(sex_mismatch_log, row.names=1)

samples_to_remove <- mismatches[["Sample_Name"]]

s <- t[["Sample_Name"]]

indices <- which(s %in% samples_to_remove)
indices
samples_to_remove
s
mismatches
t
head(t)

t_1 <- t[-indices,]
dim(t_1)
