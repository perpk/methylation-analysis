#' Progress manager for methylation pipeline
#'
#' @param total_steps Total number of steps
#' @param enabled Whether to show progress
#' @return A progress manager object
#'
#' @keywords internal
.create_progress_manager <- function(total_steps, enabled = TRUE) {
  if (!enabled) {
    return(list(
      update = function(step, msg) invisible(NULL),
      complete = function() invisible(NULL)
    ))
  }

  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent :current/:total (:elapsed)",
    total = total_steps,
    clear = FALSE,
    width = 70,
    show_after = 0
  )

  list(
    update = function(step_num, step_name, details = "") {
      pb$tick(tokens = list(step = step_name))
      if (details != "") {
        cat(paste0("\n  → ", details))
      }
    },
    complete = function() {
      cat("\n\n✨ Pipeline completed successfully!\n")
    }
  )
}
