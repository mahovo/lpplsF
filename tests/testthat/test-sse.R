# Tests for SSE (Sum of Squared Errors) function

test_that("SSE returns a single numeric value", {
  t <- 1:100
  par <- list(A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
              tc = 150, m = 0.5, omega = 9)

  # Generate synthetic data
  log_p <- LPPLS(t, par$A, par$B, par$C1, par$C2,
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
  log_p <- LPPLS(t, par$A, par$B, par$C1, par$C2,
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
  log_p <- LPPLS(t, par$A, par$B, par$C1, par$C2,
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

  base <- LPPLS(t, par$A, par$B, par$C1, par$C2,
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

  base <- LPPLS(t, par$A, par$B, par$C1, par$C2,
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

  log_p <- LPPLS(t, par$A, par$B, par$C1, par$C2,
                 par$tc, par$m, par$omega)
  log_p[50] <- NA  # Introduce an NA

  result <- SSE(par, log_p, t)

  # Should still return a numeric value (NA removed)
  expect_type(result, "double")
  expect_false(is.na(result))
})
