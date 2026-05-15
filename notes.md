> project_location = "/root/workspace/methyl-pipe-out"
> platform = "EPIC"
> cohorts <- list(
  PD_vs_Control = c("PD", "Control"),
  SWEDD_vs_Control = c("SWEDD", "Control")
)
> project_context <- .load_methylation_project(project_location, project_to_load, platform = platform, cohorts = cohorts)
