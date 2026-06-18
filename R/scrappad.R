setClass(
  "ResultsContainer",
  slots = list(
    filename = "character",
    object = "ANY",
    future = "ANY"
  )
)

library(mirai)

result_one <- new("ResultsContainer", filename = "result_one.rds", object = data.frame(x = 1:10, y = rnorm(10)), future = NULL)
result_one@future = mirai({
    Sys.sleep(5)
    o <- data.frame(x = 1:10, y = rnorm(10))
    saveRDS(o, file = f)
}, o = result_one@object, f = result_one@filename)
unresolved(result_one@future)
result_one@future[]

test <- list(
  result_one = 1
)
test$result_one

