intermediate_data_proxy <- function(pipeline_function, auto_clean = FALSE, project_context, ...) {
    library(mirai)
    source("R/project_context.R")

    results <- pipeline_function(project_context, ...)

    ## In case of disk-only, the memory must be cleaned up,
    ### since presumably the user does not want to keep the
    ### results in memory. In case of disk-and-memory, the
    ### results are kept in memory, but also saved to disk
    ### asynchronously. In case of memory-only, the results
    ### are kept in memory and not saved to disk.
    if (project_context$mode == results_mode()$disk_only) {
        print("Disk only mode: Saving results to disk (sync).")
        updated_results <- list()

        updated_results <- lapply(results, function(result) {
            saveRDS(result@object, file = result@filename)
            if (auto_clean) {
                rm(list = deparse(substitute(result@object)))
            }
            if (is(result, "ResultsContainer")) {
                return(new("ResultsContainer",
                    filename = result@filename,
                    object = result@object,
                    future = NULL
                ))
            } else {
                return(result)
            }
        })
        gc(full = TRUE)
        return(updated_results)
    }
    if (project_context$mode == results_mode()$disk_and_memory) {
        print("Disk and memory mode: Saving results to disk (async) and keeping them in memory.")
        updated_results <- lapply(results, function(result) {
            if (is(result, "ResultsContainer")) {
                future <- mirai(
                    {
                        cat(paste("Saving", f, "to disk asynchronously.\n"))
                        saveRDS(o, file = f)
                        return(f)
                    },
                    o = result@object,
                    f = result@filename
                )
                return(new("ResultsContainer",
                    filename = result@filename,
                    object = result@object,
                    future = future
                ))
            } else {
                return(result)
            }
        })
        return(updated_results)
    }
    if (project_context$mode == results_mode()$memory_only) {
        print("Memory only mode: Keeping results in memory.")
        return(results)
    }
}
