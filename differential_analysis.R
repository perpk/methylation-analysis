differential_analysis <- function(
  project_name = NULL,
  project_to_load = NULL,
  project_location = NULL,
  platform = NULL,
  cohorts = NULL,
  design_formula
) {
    if (is.null(platform)) {
        stop("Platform must be specified as '450K' or 'EPIC'")
    }
    source("R/progress_mgr.R")
    source("R/project_context.R")
    if (is.null(project_to_load)) {
        #### Create a new project
        project_context <- create_methylation_project(project_name, project_location, platform = platform, cohorts = cohorts)
    } else {
        #### Load an existing project
        project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
    }

    if (is.null(project_context$design_formula)) {
        project_context$design_formula <- design_formula
        print(paste("Design Formula is ", design_formula))
    }

    source("R/differential_probe_analysis.R")
    differential_probe_analysis(context = project_context)
}
