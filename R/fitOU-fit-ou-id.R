#' Fit the Ornstein–Uhlenbeck Model for Each Individual
#'
#' This is a wrapper function that makes fitting the Ornstein–Uhlenbeck model
#' for each individual
#' convenient using the `dynr` package.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' 
#' @inheritParams FitOU
#' @param ncores Positive integer.
#'   Number of cores to use.
#' @inherit FitOU references details
#'
#' @examples
#' \dontrun{
#' FitOUID(
#'   data = bivariate_ou,
#'   observed = c("y1", "y2"),
#'   id = "id",
#'   time = "time",
#'   verbose = FALSE
#' )
#' }
#' @family Fit Ornstein–Uhlenbeck Model Functions
#' @keywords fitOU fit
#' @export
FitOUID <- function(data,
                    observed,
                    id,
                    time,
                    mu0 = NULL,
                    sigma0 = NULL,
                    mu_start = NULL,
                    phi_start = NULL,
                    sigma_start = NULL,
                    theta_start = NULL,
                    ...,
                    ncores = NULL) {
  stopifnot(
    length(unique(data[ , id])) > 1
  )
  data <- as.data.frame(data)
  data <- split(x = data, f = data[ , id])
  if (is.null(ncores)) {
    cl <- NULL
  } else {
    if (ncores > 1) {
      cl <- parallel::makeCluster(ncores)
    } else {
      cl <- NULL
    }
  }
  return(
    pbapply::pblapply(
      X = data,
      FUN = FitOU,
      observed = observed,
      id = id,
      time = time,
      mu0 = mu0,
      sigma0 = sigma0,
      mu_start = mu_start,
      phi_start = phi_start,
      sigma_start = sigma_start,
      theta_start = theta_start,
      ...,
      cl = cl
    )
  )
}
