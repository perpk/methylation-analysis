library(rlist)

meta_vars_mapping <- function(dataset) {
  mapping <- list()
  if (dataset == "GSE165081") {
    mapping <- list.append(
      mapping,
      gender_var = "Sex",
      age_var = "Age"
    )
  } else if (dataset == "GSE111629") {
    mapping <- list.append(
      mapping,
      gender_var = "gender:ch1",
      age_var = "age:ch1"
    )
  }
  return (mapping)
}
