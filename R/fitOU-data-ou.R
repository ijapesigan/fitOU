.Subset <- function(data,
                    id,
                    time,
                    observed) {
  data <- as.matrix(data)
  data <- data[order(data[, time]), , drop = FALSE]
  data <- data[order(data[, id]), , drop = FALSE]
  data <- cbind(
    id = data[, id],
    time = data[, time],
    data[
      ,
      observed,
      drop = FALSE
    ]
  )
  ids <- unique(data[, "id"])
  data <- lapply(
    X = ids,
    FUN = function(id) {
      data <- data[
        which(data[, "id"] == id), ,
        drop = FALSE
      ]
      data[is.nan(data)] <- NA
      return(data)
    }
  )
  return(data)
}

.InsertNA <- function(data) {
  times <- unique(
    unlist(
      lapply(
        X = data,
        FUN = function(i) {
          return(
            i[, "time"]
          )
        }
      )
    )
  )
  return(
    lapply(
      X = data,
      FUN = function(i) {
        timei <- i[, "time"]
        timej <- times[!(times %in% timei)]
        if (length(timej > 0)) {
          observed <- matrix(
            data = NA,
            nrow = length(timej),
            ncol = ncol(i) - 2
          )
          varnames <- colnames(i)
          varnames <- varnames[!(varnames %in% c("id", "time"))]
          colnames(observed) <- varnames
          data <- rbind(
            i,
            cbind(
              id = unique(i[, "id"]),
              times = timej,
              observed
            )
          )
          data <- data[order(data[, "time"]), , drop = FALSE]
          return(
            data
          )
        } else {
          return(i)
        }
      }
    )
  )
}

.InitialNA <- function(data) {
  data <- lapply(
    X = data,
    FUN = function(i) {
      has_na <- any(
        is.na(
          i[1, , drop = FALSE]
        )
      )
      if (has_na & nrow(i) <= 1) {
        return(NULL)
      }
      while (has_na) {
        i <- i[-1, , drop = FALSE]
        if (nrow(i) <= 1) {
          return(NULL)
        }
        has_na <- any(
          is.na(
            i[1, ]
          )
        )
        if (has_na & nrow(i) == 1) {
          return(NULL)
        }
      }
      return(i)
    }
  )
  data[sapply(data, is.null)] <- NULL
  return(data)
}

.InitialNAN <- function(data) {
  data <- lapply(
    X = data,
    FUN = function(i) {
      has_nan <- any(
        is.nan(
          i[1, ]
        )
      )
      if (has_nan & nrow(i) <= 1) {
        return(NULL)
      }
      while (has_nan) {
        i <- i[-1, , drop = FALSE]
        if (nrow(i) <= 1) {
          return(NULL)
        }
        has_nan <- any(
          is.nan(
            i[1, ]
          )
        )
        if (has_nan & nrow(i) == 1) {
          return(NULL)
        }
      }
      return(i)
    }
  )
  data[sapply(data, is.null)] <- NULL
  return(data)
}


.ScaleID <- function(data,
                     center = TRUE,
                     scale = TRUE,
                     skip = NULL) {
  return(
    lapply(
      X = data,
      FUN = function(i,
                     center,
                     scale,
                     skip) {
        varnames <- colnames(i)
        if (is.null(skip)) {
          scale_data <- i[
            ,
            varnames[
              !(varnames %in% c("id", "time"))
            ]
          ]
        } else {
          stopifnot(
            skip %in% varnames
          )
          skip_data <- i[, skip, drop = FALSE]
          scale_data <- i[
            ,
            varnames[
              !(varnames %in% c("id", "time", skip))
            ]
          ]
        }
        if (center) {
          scaled_data <- apply(
            X = scale_data,
            MARGIN = 2,
            FUN = function(y) {
              y - mean(y, na.rm = TRUE)
            }
          )
        }
        if (scale) {
          # prevent NaN
          scaled_data <- apply(
            X = scale_data,
            MARGIN = 2,
            FUN = function(y) {
              (
                y - mean(y, na.rm = TRUE)
              ) / sd(
                y,
                na.rm = TRUE
              )^as.logical(
                sd(
                  y,
                  na.rm = TRUE
                )
              )
            }
          )
        }
        if (is.null(skip)) {
          data <-  cbind(
            id = i[, "id"],
            time = i[, "time"],
            scaled_data
          )
        } else {
          data <-  cbind(
            id = i[, "id"],
            time = i[, "time"],
            scaled_data,
            skip_data
          )
        }
        return(data)
      },
      center = center,
      scale = scale,
      skip = skip
    )
  )
}

#' Preliminary Data Preparation
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param insert_na Logical.
#'   Insert `NA` to `observed` variables for existing `time` points.
#' @param center Logical.
#'   If `center = TRUE`, mean center by `id`.
#' @param scale Logical.
#'   If `scale = TRUE`, standardize by `id`.
#' @param skip Character vector.
#'   A vector of character strings
#'   of the names of the observed variables to skip centering/scaling.
#' @param initial_na Logical.
#'   Iteratively remove rows where any observed variable
#'   for the first time point has `NA`.
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
#'
#' @examples
#' data <- DataOU(
#'   data = bivariate_ou,
#'   observed = c("y1", "y2"),
#'   id = "id",
#'   time = "time"
#' )
#' summary(data)
#'
#' @family Fit Ornstein–Uhlenbeck Model Functions
#' @keywords fitOU fit
#' @export
DataOU <- function(data,
                   observed,
                   id,
                   time,
                   insert_na = FALSE,
                   center = FALSE,
                   scale = FALSE,
                   skip = NULL,
                   initial_na = TRUE) {
  data <- .Subset(
    data = data,
    id = id,
    time = time,
    observed = observed
  )
  if (insert_na) {
    data <- .InsertNA(
      data = data
    )
  }
  if (any(c(center, scale))) {
    data <- .ScaleID(
      data = data,
      center = center,
      scale = scale,
      skip = skip
    )
  }
  if (initial_na) {
    data <- .InitialNA(
      data = data
    )
    data <- .InitialNAN(
      data = data
    )
  }
  data <- do.call(
    what = "rbind",
    args = data
  )
  data <- data[
    ,
    c(
      "id",
      "time",
      observed
    ),
    drop = FALSE
  ]
  return(
    as.data.frame(
      data
    )
  )
}
