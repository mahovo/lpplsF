# Tests for lr_trend_test() (P-11) and the fit_args slot on lppls_rolling().
#
# lr_trend_test() runs a parametric AR(1)-GARCH(1,1) bootstrap of the chi^2_np bias
# trend, which internally re-rolls the calibration B times. To keep the suite fast we
# mock fit_lppls with a deterministic stub (as in test-lppls_lagrange.R), so the LPPLS
# optimisation is instant; the AR-GARCH fit and the bootstrap machinery run for real.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

## Replace lpplsF:::fit_lppls with `stub` for the calling test_that() block.
local_mock_fit_lppls <- function(stub, envir = parent.frame()) {
  ns <- asNamespace("lpplsF")
  original <- ns$fit_lppls
  locked <- bindingIsLocked("fit_lppls", ns)
  if (locked) unlockBinding("fit_lppls", ns)
  assign("fit_lppls", stub, envir = ns)
  do.call("on.exit", list(substitute({
    assign("fit_lppls", original, envir = ns)
    if (locked) lockBinding("fit_lppls", ns)
  }), add = TRUE), envir = envir)
}

## Generic namespace mock: replace lpplsF:::`name` with `value` for the block,
## restored on exit. Used to control .fit_ar_garch() in the guard tests.
local_mock_ns <- function(name, value, envir = parent.frame()) {
  ns <- asNamespace("lpplsF")
  original <- get(name, envir = ns)
  locked <- bindingIsLocked(name, ns)
  if (locked) unlockBinding(name, ns)
  assign(name, value, envir = ns)
  do.call("on.exit", list(substitute({
    if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
    assign(name, original, envir = ns)
    if (locked) lockBinding(name, ns)
  }), add = TRUE), envir = envir)
}

## Stub returning fixed LPPLS parameters `theta` as the best fit for every window.
## tc is set window-locally (n_total + 1) so that, for a series whose true tc = N+1,
## the fitted curve reproduces the true LPPLS and the residuals equal the noise.
make_theta_stub <- function(theta) {
  function(log_price, fh = 30, hold_out = 0, tc_init = NULL, mode = "F2", ...) {
    th <- theta; th$tc <- length(log_price) + 1
    list(fit = list(c(th, list(D = 0, value = 0))))
  }
}

theta <- list(A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001, tc = 501, m = 0.6, omega = 8)
syn_fit_args <- list(lower = c(0.1, 6, -0.03, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
                     m_init = 0.2, o_init = 12, num_searches = 1, mode = "F2", factr = 1e7)

## Single-regime series (well specified: residuals of the true fit = noise).
make_single <- function(N = 500, seed = 1) {
  set.seed(seed)
  eval_lppls(1:N, theta$A, theta$B, theta$C1, theta$C2, theta$tc, theta$m, theta$omega) +
    rnorm(N) * 0.05
}
build_rolling <- function(y) {
  lppls_rolling(y, t1_seq = seq(1, 201, by = 10), hold_out = 50, fh = 30,
                fit_args = syn_fit_args)
}

# ---------------------------------------------------------------------------
# fit_args slot on lppls_rolling
# ---------------------------------------------------------------------------

test_that("lppls_rolling stores the fit_args used to build it", {
  local_mock_fit_lppls(make_theta_stub(theta))
  roll <- build_rolling(make_single())
  expect_true("fit_args" %in% names(roll))
  expect_equal(roll$fit_args$mode, "F2")
  expect_equal(roll$fit_args$lower, syn_fit_args$lower)
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("lr_trend_test rejects a non-lppls_rolling object", {
  expect_error(lr_trend_test(list(1, 2), log_price = rnorm(50)),
               "lppls_rolling")
})

test_that("lr_trend_test rejects log_price of the wrong length", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  roll <- build_rolling(make_single())
  expect_error(lr_trend_test(roll, log_price = rnorm(10), B = 4),
               "log_price")
})

# ---------------------------------------------------------------------------
# Return structure and ranges
# ---------------------------------------------------------------------------

test_that("lr_trend_test returns an lr_trend_test object with the expected fields", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  y <- make_single()
  res <- lr_trend_test(build_rolling(y), log_price = y, B = 10, seed = 7)

  expect_s3_class(res, "lr_trend_test")
  for (f in c("beta1", "lambda", "p_ols", "p_hac", "p_boot", "z", "null",
              "ar1", "garch", "persist", "valid", "B"))
    expect_true(f %in% names(res), info = f)
  expect_equal(unname(res$lambda), unname(-res$beta1))
})

test_that("lr_trend_test p-values and z are in range; null has the right length", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  y <- make_single()
  res <- lr_trend_test(build_rolling(y), log_price = y, B = 12, seed = 7)

  expect_gte(res$p_hac, 0);  expect_lte(res$p_hac, 1)
  expect_gte(res$p_boot, 0); expect_lte(res$p_boot, 1)
  expect_gte(res$z, 0)
  expect_length(res$null, res$B)
  expect_true(all(is.finite(res$null)))
})

# ---------------------------------------------------------------------------
# The alpha + beta < 1 validity guard
# ---------------------------------------------------------------------------

test_that("valid = TRUE and the bootstrap runs when the fit is stationary (persist < 1)", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  local_mock_ns(".fit_ar_garch",                      # stationary fit: alpha + beta = 0.9
                function(r) list(ar1 = 0.3, om0 = 1e-4, alpha1 = 0.10, beta1 = 0.80))
  y <- make_single()
  res <- lr_trend_test(build_rolling(y), log_price = y, B = 8, seed = 1)
  expect_lt(res$persist, 1)
  expect_true(res$valid)
  expect_length(res$null, res$B)
})

test_that("valid = FALSE (warns, skips bootstrap) when the fit is IGARCH (persist >= 1)", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  local_mock_ns(".fit_ar_garch",                      # IGARCH: alpha + beta = 1.05
                function(r) list(ar1 = 0.1, om0 = 1e-4, alpha1 = 0.15, beta1 = 0.90))
  y <- make_single()
  expect_warning(
    res <- lr_trend_test(build_rolling(y), log_price = y, B = 8, seed = 1),
    "persist|misspecif|IGARCH"
  )
  expect_false(res$valid)
  expect_true(is.na(res$p_boot))
})

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

test_that("lr_trend_test is reproducible given a seed", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  y <- make_single()
  roll <- build_rolling(y)
  a <- lr_trend_test(roll, log_price = y, B = 10, seed = 42)
  b <- lr_trend_test(roll, log_price = y, B = 10, seed = 42)
  expect_equal(a$p_boot, b$p_boot)
  expect_equal(a$null, b$null)
})

test_that("print.lr_trend_test runs without error", {
  skip_if_not_installed("tsgarch"); skip_if_not_installed("sandwich")
  local_mock_fit_lppls(make_theta_stub(theta))
  y <- make_single()
  res <- lr_trend_test(build_rolling(y), log_price = y, B = 8, seed = 1)
  expect_output(print(res), "bias-trend test")
})
