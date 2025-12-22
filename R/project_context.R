.load_methylation_project <- function(base_dir, project_id) {
  filepath <- paste(base_dir, "/", project_id, sep="")
  if (!file.exists(filepath)) {
    stop("Project not found", filepath)
  }
  project_path <- file.path(filepath, "project_context.rds")

  project_context <- readRDS(project_path)
}

create_methylation_project <- function(project_name, output_dir, keep_intermediates=T) {

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
    config = list(
      keep_intermediates = keep_intermediates,
      created = Sys.time()
    ),
    state = list(
      current_step = "initialized",
      steps_completed = character(0),
      objects_in_memory = list(),
      files_on_disk = list()
    )
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
