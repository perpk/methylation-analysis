run_combat <- function(m_values, meta_df, batch_colname, mod_colnames) {
    library(sva)
    batch <- meta_df[[batch_colname]]
    mod_matrix <- model.matrix(~ as.factor(Sample_Group), data=meta_df)
    combat_m_values <- ComBat(dat=m_values, batch=batch, mod=mod_matrix)
    return(combat_m_values)
}