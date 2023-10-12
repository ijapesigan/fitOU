#' Data Analysis - `data/bivariate_ou.rda`
#'
DataProcessBivariateOU <- function() {
  rproj <- rprojroot::is_rstudio_project
  data_dir <- rproj$find_file(
    "data"
  )
  dir.create(
    path = data_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  bivariate_ou <- read.csv(
    rproj$find_file(
      ".setup",
      "data-raw",
      "bivariate-ou.txt"
    )
  )
  save(
    bivariate_ou,
    file = file.path(
      data_dir,
      "bivariate_ou.rda"
    ),
    compress = "xz"
  )
}
DataProcessBivariateOU()
rm(DataProcessBivariateOU)
