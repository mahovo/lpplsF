# Tests for internal helper functions

# We need to test internal functions by calling them with :::
# or by testing them indirectly through exported functions

test_that("calculate_damping returns expected values", {
  # Access internal function
  calc_damp <- lpplsF:::calculate_damping

  # Test with known values
  m <- 0.5
  B <- -0.015
  omega <- 9
  C1 <- 0.001
  C2 <- 0.002

  result <- calc_damp(m, B, omega, C1, C2)

  # Manual calculation
  expected <- m * abs(B) / (omega * sqrt(C1^2 + C2^2))

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("calculate_damping handles edge cases", {
  calc_damp <- lpplsF:::calculate_damping

  # Very small oscillation amplitude
  result <- calc_damp(m = 0.5, B = -0.015, omega = 9, C1 = 1e-10, C2 = 1e-10)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("lppls_f helper computes (tc-t)^m correctly", {
  lppls_f <- lpplsF:::lppls_f

  t <- c(10, 20, 30)
  tc <- 100
  m <- 0.5

  result <- lppls_f(t, tc, m)
  expected <- (tc - t)^m

  expect_equal(result, expected)
})

test_that("lppls_g helper computes cosine term correctly", {
  lppls_g <- lpplsF:::lppls_g

  t <- 50
  tc <- 100
  m <- 0.5
  omega <- 9

  result <- lppls_g(t, tc, m, omega)
  expected <- (tc - t)^m * cos(omega * log(tc - t))

  expect_equal(result, expected)
})

test_that("lppls_h helper computes sine term correctly", {
  lppls_h <- lpplsF:::lppls_h

  t <- 50
  tc <- 100
  m <- 0.5
  omega <- 9

  result <- lppls_h(t, tc, m, omega)
  expected <- (tc - t)^m * sin(omega * log(tc - t))

  expect_equal(result, expected)
})

test_that("sum functions aggregate correctly", {
  sum_f <- lpplsF:::sum_f
  sum_g <- lpplsF:::sum_g
  sum_h <- lpplsF:::sum_h

  t <- 1:10
  tc <- 50
  m <- 0.5
  omega <- 9

  # sum_f
  result_f <- sum_f(t, tc, m)
  expected_f <- sum((tc - t)^m)
  expect_equal(result_f, expected_f)

  # sum_g
  result_g <- sum_g(t, tc, m, omega)
  expected_g <- sum((tc - t)^m * cos(omega * log(tc - t)))
  expect_equal(result_g, expected_g)

  # sum_h
  result_h <- sum_h(t, tc, m, omega)
  expected_h <- sum((tc - t)^m * sin(omega * log(tc - t)))
  expect_equal(result_h, expected_h)
})

test_that("beta_calculator produces correct linear coefficients", {
  skip_if_not_installed("symengine")

  # Create the calculator
  create_beta_calc <- lpplsF:::create_beta_calculator
  beta_calc <- create_beta_calc()

  # Generate known LPPLS data
  t <- 1:100
  true_params <- list(A = 4, B = -0.015, C1 = 0.002, C2 = 0.001)
  tc <- 150
  m <- 0.5
  omega <- 9

  log_p <- eval_lppls(t, true_params$A, true_params$B, true_params$C1,
                 true_params$C2, tc, m, omega)

  # Calculate beta given known nonlinear params
  beta <- beta_calc(log_p, t, tc, m, omega)

  # Should recover true linear parameters
  expect_equal(beta["A"], true_params$A, tolerance = 1e-6)
  expect_equal(beta["B"], true_params$B, tolerance = 1e-6)
  expect_equal(beta["C1"], true_params$C1, tolerance = 1e-6)
  expect_equal(beta["C2"], true_params$C2, tolerance = 1e-6)
})

test_that("validate_params correctly identifies valid parameters", {
  validate <- lpplsF:::validate_params

  lower <- c(0.1, 6, -1e14, 0.8)
  upper <- c(0.9, 13, -1e-14, 1e6)

  # Valid parameters
  valid_params <- list(m = 0.5, omega = 9, B = -0.01, D = 1.0)
  result <- validate(valid_params, lower, upper)
  expect_true(result$valid)
  expect_true(result$B_valid)
  expect_true(result$D_valid)
  expect_true(result$m_valid)
  expect_true(result$omega_valid)
})

test_that("validate_params correctly identifies invalid parameters", {
  validate <- lpplsF:::validate_params

  lower <- c(0.1, 6, -1e14, 0.8)
  upper <- c(0.9, 13, -1e-14, 1e6)

  # Invalid B (positive)
  invalid_B <- list(m = 0.5, omega = 9, B = 0.01, D = 1.0)
  result <- validate(invalid_B, lower, upper)
  expect_false(result$valid)
  expect_false(result$B_valid)

  # Invalid D (too low)
  invalid_D <- list(m = 0.5, omega = 9, B = -0.01, D = 0.5)
  result <- validate(invalid_D, lower, upper)
  expect_false(result$valid)
  expect_false(result$D_valid)

  # Invalid m (out of bounds)
  invalid_m <- list(m = 0.05, omega = 9, B = -0.01, D = 1.0)
  result <- validate(invalid_m, lower, upper)
  expect_false(result$valid)
  expect_false(result$m_valid)

  # Invalid omega (out of bounds)
  invalid_omega <- list(m = 0.5, omega = 15, B = -0.01, D = 1.0)
  result <- validate(invalid_omega, lower, upper)
  expect_false(result$valid)
  expect_false(result$omega_valid)
})
