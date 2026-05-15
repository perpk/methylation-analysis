run_combat <- function(m_values, meta_df, batch_colname, mod_colnames) {
    library(sva)
    batch <- as.factor(meta_df[[batch_colname]])
    mod_matrix <- model.matrix(~ ., data=meta_df[mod_colnames])
    combat_m_values <- ComBat(dat=m_values, batch=batch, mod=mod_matrix)
    return(combat_m_values)
}