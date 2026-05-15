enrich_meta_with_batch <- function(meta_df, batch_df, x_colname, y_colname) {
  enriched_meta <- merge(
    x = meta_df,
    y = batch_df,
    by.x = x_colname,
    by.y = y_colname,
    all.x = FALSE
  )
  return (enriched_meta)
}