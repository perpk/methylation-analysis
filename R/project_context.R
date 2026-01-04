.load_methylation_project <- function(base_dir, project_id, platform = NULL) {
  if (is.null(platform)) {
    stop("Platform must be specified when loading a project (Either 450k or EPIC).")
  }
  filepath <- paste(base_dir, "/", project_id, sep = "")
  if (!file.exists(filepath)) {
    stop("Project not found", filepath)
  }
  project_path <- file.path(filepath, "project_context.rds")

  project_context <- readRDS(project_path)
  if (is.null(project_context$platform)) {
    project_context$platform <- platform
  }
  (project_context)
}

create_methylation_project <- function(project_name, output_dir, keep_intermediates = T, platform = NULL) {
  if (is.null(platform)) {
    stop("Platform must be specified when creating a project (Either 450k or EPIC).")
  }
  prog <- .create_progress_manager(1)
  prog$update(1, "Creating Project Context...", "Preparing file and folder structures")
  project_id <- paste0(project_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  cat("Creating methylation project:", project_id, "\n")

  base_dir <- file.path(output_dir, project_id)

  dirs <- c(
    "raw_data",
    "qc",
    "normalized",
    "processed",
    "results",
    "plots",
    "logs",
    "temp"
  )

  paths <- sapply(dirs, function(d) {
    path <- file.path(base_dir, d)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
  })

  context <- list(
    project_id = project_id,
    project_name = project_name,
    base_dir = base_dir,
    paths = as.list(paths),
    platform = platform
  )

  class(context) <- "methylation_project"
  saveRDS(context, file.path(base_dir, "project_context.rds"))
  prog$complete()

  cat("Created methylation project:", project_id, "\n")
  cat("Location:", base_dir, "\n")

  return(context)
}

#' Get or set current project context
#' @keywords internal
.project_context <- local({
  current_context <- NULL

  list(
    get = function() {
      if (is.null(current_context)) {
        stop("No active project context. Use set_project_context() first.")
      }
      current_context
    },
    set = function(context) {
      if (!inherits(context, "methylation_project")) {
        stop("Context must be a methylation_project object")
      }
      current_context <<- context
      invisible(context)
    },
    clear = function() {
      current_context <<- NULL
    }
  )
})
