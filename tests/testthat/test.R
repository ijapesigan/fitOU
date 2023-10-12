## ---- test
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste(text, "eigen"),
      {
        testthat::expect_true(
          TRUE
        )
      }
    )
  },
  text = "test"
)
