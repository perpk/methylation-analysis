
intermediate_data_proxy <- function(pipeline_function, project_context, ...) {
    library(mirai)

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
        for (result in results) {
            if (is(result, "ResultsContainer")) {
                cat(paste("Saving", result@filename, "to disk synchronously.\n"))
                saveRDS(result@object, file = result@filename)
                updated_result <- new("ResultsContainer", 
                    filename = result@filename, 
                    object = NULL,
                    future = NULL)
                rm(result@object)
                gc(full = TRUE)
                updated_results <- c(updated_results, list(updated_result))
            }
        }
        return(updated_results)
    }
    if (project_context$mode == results_mode()$disk_and_memory) {
        print("Disk and memory mode: Saving results to disk (async) and keeping them in memory.")
        updated_results <- list()
        for (result in results) {
            if (is(result, "ResultsContainer")) {
                _future <- mirai({
                    cat(paste("Saving", _f, "to disk asynchronously.\n"))
                    saveRDS(_o, file = _f)
                    return(_f)
                }, _o = result@object, _f = result@filename)
                updated_result <- new("ResultsContainer", 
                    filename = result@filename, 
                    object = result@object, 
                    future = _future)
                updated_results <- c(updated_results, list(updated_result))
            }
        }
        return(updated_results)
    }
    if (project_context$mode == results_mode()$memory_only) {
        print("Memory only mode: Keeping results in memory.")
        return(results)
    }
}