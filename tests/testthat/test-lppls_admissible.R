# Tests for lppls_admissible(), and for the out-of-range flag print.lppls_fit()
# raises when the reported estimate itself fails a condition.

lower <- c(0.1, 6, -1e14, 0.8)
upper <- c(0.9, 13, -1e-14, 1e6)

## A cheap, well-behaved F2 calibration, as in test-print_lppls_fit.R. The
## warnings it may raise about out-of-range searches are not under test here.
small_fit <- function(...) {
  t     <- 1:200
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8)
  suppressWarnings(fit_lppls(log_price = log_p, fh = 10, hold_out = 0,
                             tc_init = 205, m_init = 0.5, o_init = 8,
                             num_searches = 2, mode = "F2", ...))
}

## Printed output joined into one whitespace-normalised string, for phrases that
## wrapping may split.
flat_print <- function(x) gsub("\\s+", " ", paste(capture.output(print(x)), collapse = " "))


test_that("each condition is checked on its own", {
  est <- data.frame(
    m     = c(0.5,   0.05,  0.5,   0.5,  0.5,   0.5),
    omega = c(9,     9,     14,    9,    9,     9),
    B     = c(-0.01, -0.01, -0.01, 0.01, -0.01, NA),
    D     = c(1,     1,     1,     1,    0.5,   1)
  )
  r <- lppls_admissible(est, lower = lower, upper = upper)
  expect_identical(r$m_ok,       c(TRUE, FALSE, TRUE,  TRUE,  TRUE,  TRUE))
  expect_identical(r$omega_ok,   c(TRUE, TRUE,  FALSE, TRUE,  TRUE,  TRUE))
  expect_identical(r$B_ok,       c(TRUE, TRUE,  TRUE,  FALSE, TRUE,  FALSE))
  expect_identical(r$D_ok,       c(TRUE, TRUE,  TRUE,  TRUE,  FALSE, TRUE))
  expect_identical(r$admissible, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))

  ## The estimates come back unchanged, with the checks added.
  expect_identical(r[names(est)], est)
})


test_that("the conditions are those of fit_lppls(): B strictly below upper[3]", {
  est <- data.frame(m = c(0.1, 0.9), omega = c(6, 13), B = c(-1e-14, -1),
                    D = c(0.8, 0.8))
  r <- lppls_admissible(est, lower = lower, upper = upper)
  expect_identical(r$m_ok,     c(TRUE, TRUE))    # bounds inclusive
  expect_identical(r$omega_ok, c(TRUE, TRUE))
  expect_identical(r$B_ok,     c(FALSE, TRUE))   # B = upper[3] is out of range
  expect_identical(r$D_ok,     c(TRUE, TRUE))    # D = lower[4] is in range
})


test_that("one estimate can be given as a list, checked against the default bounds", {
  r <- lppls_admissible(list(tc = 210, m = 0.5, omega = 9, B = c(B = -0.01), D = 1))
  expect_identical(nrow(r), 1L)
  expect_true(r$admissible)
  expect_false(lppls_admissible(list(m = 0.5, omega = 9, B = -0.01, D = 0.5))$admissible)
})


test_that("missing columns and malformed bounds are errors", {
  expect_error(lppls_admissible(data.frame(m = 0.5, omega = 9, B = -0.01)),
               "has no D")
  expect_error(lppls_admissible(data.frame(m = 0.5, omega = 9, B = -0.01, D = 1),
                                lower = c(0.1, 6)),
               "length 4")
  expect_error(lppls_admissible(1:3), "lppls_fit")
})


test_that("a fit is checked against the bounds it was run with", {
  f <- small_fit(lower = lower, upper = upper)

  r <- lppls_admissible(f)
  expect_identical(nrow(r), 1L)
  expect_equal(r$tc, f$fit[[1]]$tc)
  expect_true(r$admissible)

  ## Every per-tc fit, including any the B filter set aside.
  all_fits <- if (length(f$fit) >= 3L) f$fit[[3]] else f$fit[[2]]
  expect_identical(nrow(lppls_admissible(f, which = "fits")), nrow(all_fits))

  ## Other bounds can be applied after the fact, one side at a time.
  strict <- lppls_admissible(f, lower = c(lower[1:3], 1e3))
  expect_false(strict$D_ok)
  expect_false(strict$admissible)

  ## A fit that does not record its bounds needs them supplied.
  f$fit_args <- NULL
  expect_error(lppls_admissible(f), "does not record its bounds")
  expect_true(lppls_admissible(f, lower = lower, upper = upper)$admissible)
})


test_that("a rolling calibration is checked window by window", {
  roll <- structure(list(
    table = data.frame(T1 = 1:3, m = c(0.5, 0.5, 0.95), omega = 9,
                       B = c(-0.01, 0.01, -0.01), D = 1),
    fit_args = list(lower = lower, upper = upper)),
    class = "lppls_rolling")
  r <- lppls_admissible(roll)
  expect_identical(r$T1, 1:3)
  expect_identical(r$admissible, c(TRUE, FALSE, FALSE))

  ## Bounds lppls_rolling() was not given are fit_lppls()'s defaults.
  roll$fit_args <- list()
  expect_identical(lppls_admissible(roll)$admissible, c(TRUE, FALSE, FALSE))
})


test_that("print() flags an estimate whose damping is out of range", {
  ## D is reported, not filtered, so a floor no fit can reach leaves the
  ## estimate where it was, but now out of range.
  f <- small_fit(lower = c(lower[1:3], 1e3), upper = upper)
  expect_match(flat_print(f),
               "! the estimate is out of range: D = [0-9.]+ \\(needs D >= 1000\\)")

  s <- summary(f)
  expect_false(s$admissible[["D_ok"]])
  expect_false(s$admissible[["admissible"]])
})


test_that("print() flags an estimate whose B is out of range", {
  ## Every B here is near -0.025, so none is below -1: no fit passes the B
  ## filter, fit_lppls() falls back to the unfiltered best fit, and print()
  ## says the estimate is out of range.
  f <- small_fit(lower = lower, upper = c(upper[1:2], -1, upper[4]))
  expect_match(flat_print(f), "! the estimate is out of range: B = \\S+ \\(needs B < -1\\)")
  expect_false(summary(f)$admissible[["B_ok"]])
})


test_that("an estimate in range raises no flag", {
  f <- small_fit(lower = lower, upper = upper)
  expect_false(grepl("out of range:", flat_print(f)))
  expect_true(summary(f)$admissible[["admissible"]])
})
