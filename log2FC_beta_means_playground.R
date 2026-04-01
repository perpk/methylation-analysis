beta_means <- read.csv(file.path(context$paths$results, "beta_means.csv"))
beta_means$ProbeID <- beta_means$X
beta_means <- beta_means[-1, ]

dmp_pc1_pc3 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC1_PC3.csv")
dmp_pc1_3 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC1-3.csv")
dmp_pc1_2 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC1-2.csv")
dmp_pc1 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC1.csv")
dmp_pc2_pc3 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC2_PC3.csv")
dmp_pc2 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC2.csv")
dmp_pc3 <- read_csv("/Volumes/Elements/GSE111629_temp/dmp_annotated_results_PC3.csv")

dmp_pc1_pc3 <- dmp_pc1_pc3 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc1_3 <- dmp_pc1_3 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc1_2 <- dmp_pc1_2 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc1 <- dmp_pc1 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc2_pc3 <- dmp_pc2_pc3 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc2 <- dmp_pc2 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)
dmp_pc3 <- dmp_pc3 %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 0.2)

dmp_pc1_pc3 <- merge(dmp_pc1_pc3, beta_means, by = "ProbeID")
dmp_pc1_3 <- merge(dmp_pc1_3, beta_means, by = "ProbeID")
dmp_pc1_2 <- merge(dmp_pc1_2, beta_means, by = "ProbeID")
dmp_pc1 <- merge(dmp_pc1, beta_means, by = "ProbeID")
dmp_pc2_pc3 <- merge(dmp_pc2_pc3, beta_means, by = "ProbeID")
dmp_pc2 <- merge(dmp_pc2, beta_means, by = "ProbeID")
dmp_pc3 <- merge(dmp_pc3, beta_means, by = "ProbeID")
dmp_pc1_pc3 <- dmp_pc1_pc3 %>% arrange(desc(delta_beta))
dmp_pc1_3 <- dmp_pc1_3 %>% arrange(desc(delta_beta))
dmp_pc1_2 <- dmp_pc1_2 %>% arrange(desc(delta_beta))
dmp_pc1 <- dmp_pc1 %>% arrange(desc(delta_beta))
dmp_pc2_pc3 <- dmp_pc2_pc3 %>% arrange(desc(delta_beta))
dmp_pc2 <- dmp_pc2 %>% arrange(desc(delta_beta))
dmp_pc3 <- dmp_pc3 %>% arrange(desc(delta_beta))

write_csv(dmp_pc1_pc3, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc1_pc3_filtered_deltaBeta.csv"))
write_csv(dmp_pc1_3, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc1_3_filtered_deltaBeta.csv"))
write_csv(dmp_pc1_2, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc1_2_filtered_deltaBeta.csv"))
write_csv(dmp_pc1, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc1_filtered_deltaBeta.csv"))
write_csv(dmp_pc2_pc3, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc2_pc3_filtered_deltaBeta.csv"))
write_csv(dmp_pc2, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc2_filtered_deltaBeta.csv"))
write_csv(dmp_pc3, file.path("/Volumes/Elements/GSE111629_temp/", "dmp_pc3_filtered_deltaBeta.csv"))
