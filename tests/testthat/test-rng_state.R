## fit_lppls() seeds its own random starts in modes F2 and MPL, so a calibration
## is reproducible. It must still leave the caller's random number stream as it
## found it.

## A cheap calibration with random restarts (num_searches > 1), so the internal
## seeding is exercised.
small_fit <- function(mode = "F2", ...) {
  t     <- 1:200
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8)
  fit_lppls(log_price = log_p, fh = 10, hold_out = 0,
            lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
            tc_init = 205, m_init = 0.5, o_init = 8, num_searches = 3,
            mode = mode, ...)
}


test_that("an F2 or MPL fit leaves the caller's random number stream untouched", {
  for (mode in c("F2", "MPL")) {
    set.seed(42)
    expected <- stats::runif(3)          # the caller's next draws without a fit
    set.seed(42)
    small_fit(mode)
    expect_identical(stats::runif(3), expected)
  }
})


test_that("the trace replay leaves the stream untouched too", {
  set.seed(7)
  expected <- stats::runif(1)
  set.seed(7)
  small_fit("F2", tp = c(1, 0, 0))
  expect_identical(stats::runif(1), expected)
})


test_that("a fit made with no seed in place leaves none behind", {
  saved <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit(if (!is.null(saved)) assign(".Random.seed", saved, envir = globalenv()))
  if (!is.null(saved)) rm(".Random.seed", envir = globalenv())

  small_fit("F2")
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})


test_that("the internal seeding still makes F2 calibrations reproducible", {
  set.seed(1)
  a <- small_fit("F2")
  set.seed(999)
  b <- small_fit("F2")
  expect_identical(a$fit[[1]], b$fit[[1]])
  expect_identical(a$fit[[2]], b$fit[[2]])
})


test_that("F1 still draws its restarts from the caller's stream", {
  set.seed(42)
  before <- .Random.seed
  small_fit("F1")
  expect_false(identical(.Random.seed, before))   # it consumed random numbers
})
