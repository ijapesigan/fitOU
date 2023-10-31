#' Fit the Ornstein–Uhlenbeck Model
#'
#' This is a wrapper function that makes fitting the Ornstein–Uhlenbeck model
#' convenient using the `dynr` package.
#'
#' @details The measurement model is given by
#' \deqn{
#'   \mathbf{y}_{i, t}
#'   =
#'   \boldsymbol{\nu}
#'   +
#'   \boldsymbol{\Lambda}
#'   \boldsymbol{\eta}_{i, t}
#'   +
#'   \boldsymbol{\varepsilon}_{i, t}
#'   \quad
#'   \mathrm{with}
#'   \quad
#'   \boldsymbol{\varepsilon}_{i, t}
#'   \sim
#'   \mathcal{N}
#'   \left(
#'   \mathbf{0},
#'   \boldsymbol{\Theta}
#'   \right)
#' }
#'
#' where \eqn{\mathbf{y}_{i, t}}, \eqn{\boldsymbol{\eta}_{i, t}},
#' and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#' are random variables and \eqn{\boldsymbol{\nu}},
#' \eqn{\boldsymbol{\Lambda}},
#' and \eqn{\boldsymbol{\Theta}} are model parameters.
#' \eqn{\mathbf{y}_{i, t}} is a vector of observed random variables
#' at time \eqn{t} and individual \eqn{i},
#' \eqn{\boldsymbol{\eta}_{i, t}} is a vector of latent random variables
#' at time \eqn{t} and individual \eqn{i},
#' and \eqn{\boldsymbol{\varepsilon}_{i, t}}
#' is a vector of random measurement errors
#' at time \eqn{t} and individual \eqn{i},
#' while \eqn{\boldsymbol{\nu}} is a vector of intercept,
#' \eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
#' and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
#' \eqn{\boldsymbol{\varepsilon}}.
#'
#' The dynamic structure is given by
#'
#' \deqn{
#'   \mathrm{d} \boldsymbol{\eta}_{i, t}
#'   =
#'   \boldsymbol{\Phi}
#'   \left(
#'   \boldsymbol{\mu}
#'   -
#'   \boldsymbol{\eta}_{i, t}
#'   \right)
#'   \mathrm{d}t
#'   +
#'   \boldsymbol{\Sigma}^{\frac{1}{2}}
#'   \mathrm{d}
#'   \mathbf{W}_{i, t}
#' }
#'
#' where \eqn{\boldsymbol{\mu}} is the long-term mean or equilibrium level,
#' \eqn{\boldsymbol{\Phi}} is the rate of mean reversion,
#' determining how quickly the variable returns to its mean,
#' \eqn{\boldsymbol{\Sigma}} is the matrix of volatility
#' or randomness in the process, and \eqn{\mathrm{d}\boldsymbol{W}}
#' is a Wiener process or Brownian motion,
#' which represents random fluctuations.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param data Data frame.
#'   A data frame object of data for potentially
#'   multiple subjects that contain
#'   a column of subject ID numbers
#'   (i.e., an ID variable),
#'   a column indicating subject-specific measurement occasions
#'   (i.e., a TIME variable),
#'   at least one column of observed values.
#' @param observed Character vector.
#'   A vector of character strings
#'   of the names of the observed variables in the data.
#' @param id Character string.
#'   A character string of the name of the ID variable in the data.
#' @param time Character string.
#'   A character string of the name of the TIME variable in the data.
#' @param mu0 Numeric vector.
#'   Mean of initial latent variable values
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
#'   If `mu0 = NULL`, a vector of zeros is used.
#' @param sigma0 Numeric matrix.
#'   Covariance matrix
#'   of initial latent variable values
#'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
#'   If `sigma0 = NULL`, an identity matrix is used.
#' @param mu_start Numeric vector.
#'   Starting values of the `mu` vector,
#'   that is,
#'   the long-term mean or equilibrium level.
#'   If `mu_start = NULL`,
#'   a vector means of the observed variables is used.
#' @param phi_start Numeric matrx.
#'   Starting values of the `phi` matrix,
#'   that is,
#'   the rate of mean reversion,
#'   determining how quickly the variable returns to its mean.
#'   If `phi_start = NULL`, a matrix of zeros is used.
#' @param sigma_start Numeric matrx.
#'   Starting values of the `sigma` matrix,
#'   that is,
#'   the matrix of volatility
#'   or randomness in the process.
#'   If `sigma_start = NULL`, an identity matrix is used.
#' @param theta_start Numeric matrx.
#'   Starting values of the `theta` matrix,
#'   that is,
#'   the measurement error covariance matrix
#'   (\eqn{\boldsymbol{\Theta}}).
#'   If `theta_start = NULL`, an identity matrix is used.
#' @param sigma_diag Logical.
#'   If `sigma_diag = TRUE`,
#'   estimate only the diagonals of \eqn{\boldsymbol{\Sigma}}.
#' @param center Logical.
#'   If `center = TRUE`, mean center by `id`.
#' @param scale Logical.
#'   If `scale = TRUE`, standardize by `id`.
#' @param ub Numeric vector.
#'   Optional.
#'   The upper bounds for \eqn{\boldsymbol{\Phi}}.
#' @param lb Numeric vector.
#'   Optional.
#'   The lower bounds for \eqn{\boldsymbol{\Phi}}.
#' @param ... Additional arguments to pass to [dynr::dynr.cook()].
#'
#' @references
#'   Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
#'   Continuous-time dynamic models:
#'   Connections to structural equation models and other discrete-time models.
#'   In R. H. Hoyle (Ed.), Handbook of structural equation modeling (2nd ed.).
#'   The Guilford Press.
#'
#'   Ou, L., Hunter, M. D., & Chow, S.-M. (2019).
#'   What's for dynr:
#'   A package for linear and nonlinear dynamic modeling in R.
#'   *The R Journal*, *11*(1), 91.
#'   \doi{10.32614/rj-2019-012}
#'
#'   Uhlenbeck, G. E., & Ornstein, L. S. (1930).
#'   On the theory of the brownian motion.
#'   *Physical Review*, *36*(5), 823–841.
#'   \doi{doi.org/10.1103/physrev.36.823}
#'
#' @examples
#' \dontrun{
#' FitOU(
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
FitOU <- function(data,
                  observed,
                  id,
                  time,
                  mu0 = NULL,
                  sigma0 = NULL,
                  mu_start = NULL,
                  phi_start = NULL,
                  sigma_start = NULL,
                  theta_start = NULL,
                  sigma_diag = FALSE,
                  center = FALSE,
                  scale = FALSE,
                  lb = NULL,
                  ub = NULL,
                  ...) {
  y_names <- observed
  k <- length(y_names)
  eta_names <- paste0("eta_", seq_len(k))
  iden <- diag(k)
  null_vec <- rep(x = 0, times = k)
  if (center | scale) {
    if (scale) {
      foo <- function(x) {
        (
          x - mean(
            x = x,
            na.rm = TRUE
          )
        ) / sd(
          x = x,
          na.rm = TRUE
        )
      }
    } else {
      foo <- function(x) {
        x - mean(
          x = x,
          na.rm = TRUE
        )
      }
    }
    for (i in seq_len(k)) {
      data[, y_names[i]] <- stats::ave(
        x = data[, y_names[i]],
        data[, id],
        FUN = foo
      )
    }
  }
  # data
  dynr_data <- dynr::dynr.data(
    dataframe = data,
    id = id,
    time = time,
    observed = y_names
  )
  # initial
  if (is.null(mu0)) {
    mu0 <- null_vec
  }
  if (is.null(sigma0)) {
    sigma0 <- iden
  }
  row_idx <- rep(x = 1:k, each = k)
  col_idx <- rep(x = 1:k, times = k)
  dynr_initial <- dynr::prep.initial(
    values.inistate = mu0,
    params.inistate = paste0("mu0_", seq_len(k)),
    values.inicov = sigma0,
    params.inicov = matrix(
      data = paste0(
        "sigma0_",
        pmin(row_idx, col_idx),
        pmax(row_idx, col_idx)
      ),
      nrow = k
    )
  )
  # measurement
  dynr_measurement <- dynr::prep.measurement(
    values.load = iden,
    params.load = matrix(data = "fixed", nrow = k, ncol = k),
    state.names = eta_names,
    obs.names = y_names
  )
  # dynamics
  formula <- character(k)
  mu_names <- character(k)
  phi_names <- matrix(
    data = "0",
    nrow = k,
    ncol = k
  )
  for (i in seq_len(k)) {
    terms <- character(k)
    for (j in 1:k) {
      terms[j] <- paste0("phi_", i, j, " * (mu_", j, " - eta_", j, ")")
      phi_names[i, j] <- paste0("phi_", i, j)
    }
    formula[i] <- paste0("eta_", i, " ~ ", paste(terms, collapse = " + "))
    mu_names[i] <- paste0("mu_", i)
  }
  formula_list <- list()
  for (i in formula) {
    formula_obj <- stats::as.formula(i)
    formula_list <- c(formula_list, list(formula_obj))
  }
  dim(phi_names) <- NULL
  if (is.null(mu_start)) {
    mu_start <- colMeans(
      data[, observed, drop = FALSE],
      na.rm = TRUE
    )
  }
  names(mu_start) <- mu_names
  if (is.null(phi_start)) {
    phi_start <- rep(x = 0.1, times = k * k)
  } else {
    dim(phi_start) <- NULL
  }
  names(phi_start) <- phi_names
  dynr_dynamics <- dynr::prep.formulaDynamics(
    formula = formula_list,
    startval = c(
      mu_start,
      phi_start
    ),
    isContinuousTime = TRUE
  )
  # noise
  if (is.null(sigma_start)) {
    sigma_start <- iden
  }
  if (is.null(theta_start)) {
    theta_start <- iden
  }
  if (sigma_diag) {
    sigma_params <- matrix(
      data = "fixed",
      nrow = k,
      ncol = k
    )
    diag(sigma_params) <- paste0(
      "sigma_",
      1:k,
      1:k
    )
    sigma_start_diag <- diag(sigma_start)
    sigma_start <- matrix(
      data = 0,
      nrow = k,
      ncol = k
    )
    diag(sigma_start) <- sigma_start_diag
  } else {
    sigma_params <- matrix(
      data = paste0(
        "sigma_",
        pmin(row_idx, col_idx),
        pmax(row_idx, col_idx)
      ),
      nrow = k
    )
  }
  theta <- matrix(data = "fixed", nrow = k, ncol = k)
  diag(theta) <- paste0("theta_", seq_len(k), seq_len(k))
  dynr_noise <- dynr::prep.noise(
    values.latent = sigma_start,
    params.latent = sigma_params,
    values.observed = theta_start,
    params.observed = theta
  )
  # model
  model <- dynr::dynr.model(
    data = dynr_data,
    initial = dynr_initial,
    measurement = dynr_measurement,
    dynamics = dynr_dynamics,
    noise = dynr_noise,
    outfile = paste0(tempfile(), ".c")
  )
  if (!is.null(lb)) {
    model$lb[phi_names] <- lb
  }
  if (!is.null(ub)) {
    model$ub[phi_names] <- ub
  }
  # fit
  fit <- dynr::dynr.cook(
    dynrModel = model,
    ...
  )
  return(
    list(
      model = model,
      fit = fit
    )
  )
}
