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
#' @param model Output of [fitOU::ModelOU()].
#' @param retry Positive integer.
#'   Maximum number of reruns.
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
#' data <- DataOU(
#'   data = bivariate_ou,
#'   observed = c("y1", "y2"),
#'   id = "id",
#'   time = "time"
#' )
#' model <- ModelOU(
#'   data = data,
#'   observed = c("y1", "y2"),
#'   id = "id",
#'   time = "time"
#' )
#' FitOU(
#'   model = model,
#'   verbose = FALSE
#' )
#' }
#' @family Fit Ornstein–Uhlenbeck Model Functions
#' @keywords fitOU fit
#' @export
FitOU <- function(model,
                  retry = NULL,
                  ...) {
  fit <- dynr::dynr.cook(
    dynrModel = model,
    ...
  )
  if (is.null(retry)) {
    return(
      fit
    )
  } else {
    exitflag <- fit$exitflag
    if (exitflag %in% c(1, 2, 3, 4)) {
      return(fit)
    }
    # loop if exitflag is not 1, 2, 3, 4
    i <- 0
    while (!(exitflag %in% c(1, 2, 3, 4))) {
      model$xstart <- fit$transformed.parameters + stats::runif(
        n = length(model$xstart),
        min = -.02,
        max = .02
      )
      fit <- dynr::dynr.cook(
        dynrModel = model,
        ...
      )
      exitflag <- fit$exitflag
      i <- i + 1
      if (i == retry) {
        warning(
          paste0(
            "Maximum number of retries reached.\n"
          )
        )
        return(fit)
      }
    }
  }
}
