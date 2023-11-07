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
      return(
        data[
          which(data[, "id"] == id), ,
          drop = FALSE
        ]
      )
    }
  )
  return(
    data
  )
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
            i[1, ] # i[1, , drop = FALSE]
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
          i[1, ] # i[1, , drop = FALSE]
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
                     scale = TRUE) {
  return(
    lapply(
      X = data,
      FUN = function(i,
                     center,
                     scale) {
        varnames <- colnames(i)
        varnames <- varnames[!(varnames %in% c("id", "time"))]
        data <- scale(
          x = i[, varnames, drop = FALSE],
          center = center,
          scale = scale
        )
        colnames(data) <- varnames
        return(
          cbind(
            id = i[, "id"],
            time = i[, "time"],
            data
          )
        )
      },
      center = center,
      scale = scale
    )
  )
}

#' Preliminary Data Preparation
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param insert_na Logical.
#'   Insert `NA` to `observed` variables for existing `time` points.
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
#' @param center Logical.
#'   If `center = TRUE`, mean center by `id`.
#' @param scale Logical.
#'   If `scale = TRUE`, standardize by `id`.
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
                   initial_na = TRUE,
                   center = FALSE,
                   scale = FALSE) {
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
  if (initial_na) {
    data <- .InitialNA(
      data = data
    )
    data <- .InitialNAN(
      data = data
    )
  }
  if (center | scale) {
    data <- .ScaleID(
      data = data,
      center = center,
      scale = scale
    )
  }
  return(
    as.data.frame(
      do.call(
        what = "rbind",
        args = data
      )
    )
  )
}
