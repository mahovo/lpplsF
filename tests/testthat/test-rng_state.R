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


## ---- Functions with a `seed` argument -----------------------------------------
## They seed the generator only when given a seed, and then restore the caller's
## stream. Without one they draw from it, and a calibration in between must not
## reset it: before fit_lppls() restored the stream, every noise level or
## bootstrap replicate after the first repeated the same draw.

sweep_signal <- function() {
  eval_lppls(1:200, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
             tc = 221, m = 0.6, omega = 8)
}
sweep_args <- list(fh = 10, tc_init = 205, m_init = 0.5, o_init = 8,
                   num_searches = 1)

## A noise generator that keeps a copy of every draw it makes.
recording_noise <- function() {
  draws <- list()
  list(noise = function(n) {
         z <- stats::rnorm(n)
         draws[[length(draws) + 1L]] <<- z
         z
       },
       draws = function() draws)
}

test_that("lppls_noise_sweep with a seed leaves the caller's stream untouched", {
  set.seed(42)
  expected <- stats::runif(3)
  set.seed(42)
  lppls_noise_sweep(sweep_signal(), sd_seq = c(0.005, 0.01),
                    fit_args = sweep_args, seed = 3630)
  expect_identical(stats::runif(3), expected)
})


test_that("with a seed, every noise level sees the same innovations", {
  rec <- recording_noise()
  lppls_noise_sweep(sweep_signal(), sd_seq = c(0.005, 0.01), noise = rec$noise,
                    fit_args = sweep_args, seed = 3630)
  expect_identical(rec$draws()[[1]], rec$draws()[[2]])
})


test_that("without a seed, the noise levels draw independently", {
  rec <- recording_noise()
  set.seed(42)
  lppls_noise_sweep(sweep_signal(), sd_seq = c(0.005, 0.01, 0.02),
                    noise = rec$noise, fit_args = sweep_args)
  d <- rec$draws()
  expect_length(d, 3L)
  expect_false(identical(d[[1]], d[[2]]))
  expect_false(identical(d[[2]], d[[3]]))
})


## A rolling calibration for lr_trend_test(): single regime, eight windows.
trend_rolling <- function() {
  set.seed(3630)
  N  <- 300
  lp <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                   tc = N + 1, m = 0.6, omega = 8) + stats::rnorm(N) * 0.01
  roll <- lppls_rolling(lp, t1_seq = seq(1, by = 20, length.out = 8),
                        hold_out = 0, fh = 30,
                        fit_args = list(lower = c(0.1, 6, -0.03, 0.8),
                                        upper = c(0.9, 13, -1e-14, 1e6),
                                        m_init = 0.5, o_init = 8,
                                        num_searches = 1))
  list(roll = roll, log_p = lp)
}

test_that("lr_trend_test with a seed leaves the caller's stream untouched", {
  skip_on_cran()
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_if_not_installed("sandwich")
  x <- trend_rolling()
  set.seed(42)
  expected <- stats::runif(3)
  set.seed(42)
  lr_trend_test(x$roll, log_price = x$log_p, B = 3, seed = 1)
  expect_identical(stats::runif(3), expected)
})


test_that("without a seed, lr_trend_test's bootstrap replicates differ", {
  skip_on_cran()
  skip_if_not_installed("tsgarch")
  skip_if_not_installed("xts")
  skip_if_not_installed("sandwich")
  x <- trend_rolling()
  set.seed(42)
  res <- lr_trend_test(x$roll, log_price = x$log_p, B = 4)
  expect_true(res$valid)
  expect_length(unique(res$null), length(res$null))
})
