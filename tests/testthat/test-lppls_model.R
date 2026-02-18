# Tests for LPPLS model function
# Following test-driven development approach

test_that("LPPLS function returns correct dimensions", {
  t <- 1:100
  result <- LPPLS(t, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
                  tc = 150, m = 0.5, omega = 9)

  expect_length(result, 100)
  expect_type(result, "double")
})

test_that("LPPLS function handles scalar input", {
  result <- LPPLS(t = 50, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
                  tc = 150, m = 0.5, omega = 9)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_false(is.na(result))
})

test_that("LPPLS function produces NaN for t >= tc", {
  t <- 1:200
  result <- LPPLS(t, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
                  tc = 150, m = 0.5, omega = 9)

  # Values before tc should be valid
  expect_true(all(!is.nan(result[1:149])))

  # Values at and after tc will produce NaN due to (tc-t)^m with tc-t <= 0
  expect_warning(
    LPPLS(150:200, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
          tc = 150, m = 0.5, omega = 9),
    "Some time values are >= tc"
  )
})

test_that("LPPLS function validates input types", {
  expect_error(
    LPPLS(t = "not numeric", A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
          tc = 150, m = 0.5, omega = 9),
    "'t' must be a numeric vector"
  )
})

test_that("LPPLS mode 0 (Filimonov) produces expected behavior", {
  # Test that the function follows the expected formula
  t <- 50
  A <- 4
  B <- -0.015
  C1 <- 0.001
  C2 <- 0.002
  tc <- 150
  m <- 0.5
  omega <- 9

  result <- LPPLS(t, A, B, C1, C2, tc, m, omega, mode = 0)

  # Manual calculation
  tau <- tc - t
  d <- omega * log(tau)
  expected <- A + tau^m * (B + C1 * cos(d) + C2 * sin(d))

  expect_equal(result, expected)
})

test_that("LPPLS mode parameter is validated", {
  expect_error(
    LPPLS(1:10, A = 4, B = -0.015, C1 = 0.001, C2 = 0.001,
          tc = 150, m = 0.5, omega = 9, mode = 5),
    "LPPLS mode must be 0, 1, 2 or 3"
  )
})

test_that("LPPLS produces monotonically changing values approaching tc", {
  # For a bubble, prices should generally increase approaching tc
  t <- 1:100
  result <- LPPLS(t, A = 4, B = -0.015, C1 = 0.0001, C2 = 0.0001,
                  tc = 150, m = 0.5, omega = 9)

  # The trend should be increasing (ignoring oscillations)
  # Check that the moving average is increasing
  window <- 10
  ma <- stats::filter(result, rep(1/window, window), sides = 1)
  ma <- ma[!is.na(ma)]

  # Most differences should be positive for an upward trend
  diffs <- diff(ma)
  expect_true(sum(diffs > 0) > length(diffs) * 0.5)
})
