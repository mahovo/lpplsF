# Tests for SSE (Sum of Squared Errors) function

test_that("SSE returns a single numeric value", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  # Generate synthetic data
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                 par$tc, par$m, par$omega)

  result <- SSE(par, log_p, t)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_false(is.na(result))
})

test_that("SSE is zero when model perfectly fits data", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  # Generate data from the exact model
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                 par$tc, par$m, par$omega)

  result <- SSE(par, log_p, t)

  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("SSE is positive for noisy data", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  # Generate data with noise
  set.seed(123)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                 par$tc, par$m, par$omega) + rnorm(100, 0, 0.01)

  result <- SSE(par, log_p, t)

  expect_gt(result, 0)
})

test_that("SSE validates input lengths", {
  t <- 1:100
  log_p <- 1:50  # Wrong length
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  expect_error(SSE(par, log_p, t), "'log_p' and 't' must have the same length")
})

test_that("SSE increases with larger noise", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  base <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                par$tc, par$m, par$omega)

  set.seed(123)
  log_p_small_noise <- base + rnorm(100, 0, 0.001)
  set.seed(456)
  log_p_large_noise <- base + rnorm(100, 0, 0.1)

  sse_small <- SSE(par, log_p_small_noise, t)
  sse_large <- SSE(par, log_p_large_noise, t)

  # SSE should be larger with more noise (on average)
  expect_gt(sse_large, sse_small)
})

test_that("SSE is symmetric in residuals", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  base <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                par$tc, par$m, par$omega)

  # Add positive and negative deviations
  log_p_positive <- base + 0.01
  log_p_negative <- base - 0.01

  sse_pos <- SSE(par, log_p_positive, t)
  sse_neg <- SSE(par, log_p_negative, t)

  # Both should produce the same SSE (since errors are squared)
  expect_equal(sse_pos, sse_neg, tolerance = 1e-10)
})

test_that("SSE handles NA values gracefully", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                 par$tc, par$m, par$omega)
  log_p[50] <- NA  # Introduce an NA

  result <- SSE(par, log_p, t)

  # Should still return a numeric value (NA removed)
  expect_type(result, "double")
  expect_false(is.na(result))
})

# ---------------------------------------------------------------------------
# SSE — threading of higher-order parameters (T1/T2/omega2/omega3) via `par`
# ---------------------------------------------------------------------------

test_that("SSE forwards T1 and omega2 for mode 2 (not silently defaulted)", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5, omega = 9)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega)

  expect_false(isTRUE(all.equal(
    SSE(c(par, list(T1 = 300, omega2 = 5)), log_p, t, mode = 2),
    SSE(c(par, list(T1 = 500, omega2 = 5)), log_p, t, mode = 2)
  )))                                                        # T1 threaded
  expect_false(isTRUE(all.equal(
    SSE(c(par, list(T1 = 300, omega2 = 0)), log_p, t, mode = 2),
    SSE(c(par, list(T1 = 300, omega2 = 5)), log_p, t, mode = 2)
  )))                                                        # omega2 threaded
})

test_that("SSE forwards T2 for mode 3", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5, omega = 9,
              T1 = 300, omega2 = 5, omega3 = 2)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega)
  expect_false(isTRUE(all.equal(
    SSE(c(par, list(T2 = 1000)), log_p, t, mode = 3),
    SSE(c(par, list(T2 = 1990)), log_p, t, mode = 3)
  )))
})

test_that("mode 0 ignores higher-order params (behaviour unchanged)", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5, omega = 9)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega)
  expect_equal(SSE(c(par, list(T1 = 300)), log_p, t, mode = 0),
               SSE(par, log_p, t, mode = 0))
})

# ---------------------------------------------------------------------------
# SSE — validation
# ---------------------------------------------------------------------------

test_that("SSE errors when par omits a core parameter", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5)  # no omega
  expect_error(SSE(par, log_p, t), "missing required parameter")
})

test_that("SSE errors on a non-list par, non-numeric input, and bad mode", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5, omega = 9)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega)
  expect_error(SSE("nope", log_p, t),                 "'par' must be a named list")
  expect_error(SSE(par, as.character(log_p), t),      "must be numeric")
  expect_error(SSE(par, log_p, t, mode = 5),          "'mode' must be 0, 1, 2 or 3")
})

test_that("SSE warns (and defaults) when a higher-order mode's params are missing", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, tc = 150, m = 0.5, omega = 9)
  log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega)
  expect_warning(SSE(par, log_p, t, mode = 2), "eval_lppls\\(\\) default")   # T1/omega2 absent
})

# ---------------------------------------------------------------------------
# contour_plot_sse — validation and single-warning behaviour
# ---------------------------------------------------------------------------

test_that("contour_plot_sse validates vars (valid names, no fix-and-scan overlap)", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(tc = 150, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001)
  expect_error(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "nope"),
                     lower = c(0.1, 6), upper = c(0.9, 13), cp = FALSE),
    "must name eval_lppls parameters")
  expect_error(
    contour_plot_sse(log_p, t, c(par, list(m = 0.5)), vars = list(x = "m", y = "omega"),
                     lower = c(0.1, 6), upper = c(0.9, 13), cp = FALSE),
    "also fixed in 'par'")
})

test_that("contour_plot_sse validates the lower/upper bounds", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(tc = 150, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001)
  expect_error(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "omega"),
                     lower = c(0.9, 6), upper = c(0.1, 13), cp = FALSE),   # lower >= upper (x)
    "strictly less")
  expect_error(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "omega"),
                     lower = 0.1, upper = c(0.9, 13), cp = FALSE),          # length-1 lower
    "length 2")
})

test_that("contour_plot_sse errors when par + vars omit a core parameter", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(tc = 150, A = 4, B = -0.015, C1 = 0.001)                     # missing C2
  expect_error(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "omega"),
                     lower = c(0.1, 6), upper = c(0.9, 13), cp = FALSE),
    "missing from 'par'")
})

test_that("contour_plot_sse warns once (not per grid cell) for a missing higher-order param", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(tc = 150, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001)         # mode 2, no T1/omega2
  w <- capture_warnings(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "omega"),
                     lower = c(0.1, 6), upper = c(0.9, 13), cp = FALSE, mode = 2)
  )
  expect_length(w, 1L)                                       # once, not 101 * 101
  expect_match(w, "eval_lppls\\(\\) default")
})

test_that("contour_plot_sse runs mode 2 silently when T1/omega2 are supplied", {
  t <- 1:100; log_p <- rnorm(100)
  par <- list(tc = 150, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001, T1 = 300, omega2 = 5)
  expect_silent(
    contour_plot_sse(log_p, t, par, vars = list(x = "m", y = "omega"),
                     lower = c(0.1, 6), upper = c(0.9, 13), cp = FALSE, mode = 2)
  )
})
