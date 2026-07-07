# Tests for lppls_lagrange(), which post-processes an lppls_rolling() object.
#
# lppls_lagrange() itself does no fitting: it takes a rolling calibration and
# derives chi^2_np / chi^2_lambda / T1*. The rolling objects used below are built
# with a mocked fit_lppls (deterministic, instant), except for the final
# skip_on_cran() integration test which uses the real fitter end-to-end.
#
# Window-slicing, tc_init, hold_out/fh forwarding and the rolling return
# structure are covered in test-lppls_rolling.R; here we focus on the Lagrange
# math, the return structure, and the chi^2 plots.

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

## Generate clean LPPLS bubble data
make_bubble_data <- function(n = 300, tc = 320, noise_sd = 0.005, seed = 42) {
  set.seed(seed)
  t <- seq_len(n)
  log_p <- eval_lppls(
    t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
    tc = tc, m = 0.6, omega = 8
  ) + rnorm(n) * noise_sd
  list(t = t, log_p = log_p)
}

## Replace lpplsF::fit_lppls with a deterministic stub for the calling
## test_that() block.
local_mock_fit_lppls <- function(stub, envir = parent.frame()) {
  ns <- asNamespace("lpplsF")
  original <- ns$fit_lppls
  was_locked <- bindingIsLocked("fit_lppls", ns)
  if (was_locked) unlockBinding("fit_lppls", ns)
  assign("fit_lppls", stub, envir = ns)
  do.call(
    "on.exit",
    list(
      substitute({
        assign("fit_lppls", original, envir = ns)
        if (was_locked) lockBinding("fit_lppls", ns)
      }),
      add = TRUE
    ),
    envir = envir
  )
}

## Fixed best fit: window mean as A, zero oscillation amplitudes -> deterministic
## SSE per window.
make_fixed_stub <- function() {
  function(log_price, fh = 30, hold_out = 0, tc_init = NULL,
           mode = "F2", ...) {
    n_total <- length(log_price)
    n_cal <- n_total - hold_out
    A_hat <- mean(log_price[seq_len(n_cal)])
    list(fit = list(list(
      tc = if (!is.null(tc_init)) tc_init else (n_cal + floor(fh / 2)),
      m = 0.5, omega = 8, A = A_hat, B = -0.01, C1 = 0, C2 = 0, D = 1,
      value = 0
    )))
  }
}

## Build a rolling calibration from clean bubble data using a supplied stub.
build_rolling <- function(stub = make_fixed_stub(), n = 120,
                          t1_seq = c(1, 10, 20, 30, 40), hold_out = 5,
                          fh = 20, envir = parent.frame()) {
  data <- make_bubble_data(n = n)
  local_mock_fit_lppls(stub, envir = envir)
  lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    fit_args = list(num_searches = 1, mode = "F2")
  )
}

## Find the first layer in a ggplot whose geom inherits from `class_name`.
find_layer <- function(p, class_name) {
  idx <- which(vapply(
    p$layers,
    function(l) inherits(l$geom, class_name),
    logical(1)
  ))
  if (length(idx) == 0L) {
    stop(sprintf("No layer of geom class '%s' found", class_name))
  }
  p$layers[[idx[1]]]
}

# ----------------------------------------------------------------------------
# Input validation
# ----------------------------------------------------------------------------

test_that("lppls_lagrange rejects a non-lppls_rolling object", {
  expect_error(lppls_lagrange(list()),
               "must be an object of class 'lppls_rolling'")
  expect_error(lppls_lagrange(1:10),
               "must be an object of class 'lppls_rolling'")
  expect_error(lppls_lagrange("abc"),
               "must be an object of class 'lppls_rolling'")
})

test_that("lppls_lagrange errors when n_cal - p is non-positive", {
  ## Window at T1 = 44 (hold_out = 0) has n_cal = 7, so p = 7 -> denom 0.
  local_mock_fit_lppls(make_fixed_stub())
  roll <- lppls_rolling(
    log_price = rnorm(50),
    t1_seq = c(1, 44),
    hold_out = 0,
    fh = 10,
    fit_args = list(num_searches = 1)
  )
  expect_error(lppls_lagrange(roll, p = 7), "Non-positive denominator")
})

# ----------------------------------------------------------------------------
# Return value structure
# ----------------------------------------------------------------------------

test_that("lppls_lagrange returns an lppls_lagrange-classed list with all fields", {
  t1_seq <- c(1, 10, 20, 30)
  roll <- build_rolling(t1_seq = t1_seq, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  expect_s3_class(res, "lppls_lagrange")
  expect_named(
    res,
    c("table", "rolling", "trend", "intercept", "slope", "lambda",
      "T2", "p", "t1_star"),
    ignore.order = TRUE
  )

  expect_s3_class(res$table, "data.frame")
  expect_named(res$table, c("T1", "n", "chi_sq_np", "chi_sq_lambda"))
  expect_equal(nrow(res$table), length(t1_seq))

  expect_s3_class(res$rolling, "lppls_rolling")
  expect_identical(res$rolling, roll)
  expect_s3_class(res$trend, "lm")

  expect_length(res$intercept, 1L)
  expect_length(res$slope, 1L)
  expect_length(res$lambda, 1L)
  expect_length(res$T2, 1L)
  expect_length(res$p, 1L)
  expect_length(res$t1_star, 1L)
})

test_that("T2 and table$T1 / table$n are carried through from the rolling object", {
  t1_seq <- c(1, 11, 21, 31)
  hold_out <- 7
  roll <- build_rolling(t1_seq = t1_seq, hold_out = hold_out, n = 100)
  res <- lppls_lagrange(roll, p = 0)

  expect_equal(res$T2, roll$T2)
  expect_equal(res$T2, 100 - hold_out)
  expect_equal(res$table$T1, t1_seq)
  expect_equal(res$table$n, roll$table$n_cal)
})

# ----------------------------------------------------------------------------
# Mathematical correctness: chi^2_np, trend, lambda, chi^2_lambda, t1*
# ----------------------------------------------------------------------------

test_that("chi_sq_np[k] = sse_cal[k] / (n_cal[k] - p)", {
  t1_seq <- c(1, 10, 25, 40)
  p_val <- 7
  roll <- build_rolling(t1_seq = t1_seq, hold_out = 8, fh = 20, n = 120)
  res <- lppls_lagrange(roll, p = p_val)

  for (k in seq_along(t1_seq)) {
    expect_equal(
      res$table$chi_sq_np[k],
      roll$table$sse_cal[k] / (roll$table$n_cal[k] - p_val),
      tolerance = 1e-12
    )
  }
})

test_that("p = 0 reproduces the thesis form (no parameter correction)", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 100)
  res0 <- lppls_lagrange(roll, p = 0)
  res7 <- lppls_lagrange(roll, p = 7)

  ## Same fits + same windows -> same SSEs; chi_np differs only by denom.
  for (k in seq_len(nrow(res0$table))) {
    n_cal <- res0$table$n[k]
    expect_equal(
      res0$table$chi_sq_np[k] * n_cal,
      res7$table$chi_sq_np[k] * (n_cal - 7),
      tolerance = 1e-12
    )
  }
})

test_that("intercept and slope equal coefficients of lm(chi_sq_np ~ T1)", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30, 40), n = 120)
  res <- lppls_lagrange(roll, p = 0)

  manual <- stats::lm(res$table$chi_sq_np ~ res$table$T1)
  expect_equal(res$intercept, unname(stats::coef(manual)[[1]]),
               tolerance = 1e-12)
  expect_equal(res$slope, unname(stats::coef(manual)[[2]]),
               tolerance = 1e-12)
})

test_that("lambda equals -slope", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30, 40), n = 120)
  res <- lppls_lagrange(roll, p = 0)
  expect_equal(res$lambda, -res$slope, tolerance = 1e-14)
})

test_that("chi_sq_lambda = chi_sq_np - lambda * (T2 - T1)", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30, 40), n = 120)
  res <- lppls_lagrange(roll, p = 0)

  expected <- res$table$chi_sq_np - res$lambda * (res$T2 - res$table$T1)
  expect_equal(res$table$chi_sq_lambda, expected, tolerance = 1e-12)
})

test_that("t1_star = t1_seq[which.min(chi_sq_lambda)]", {
  t1_seq <- c(1, 10, 20, 30, 40)
  roll <- build_rolling(t1_seq = t1_seq, n = 120)
  res <- lppls_lagrange(roll, p = 0)

  expect_equal(res$t1_star, t1_seq[which.min(res$table$chi_sq_lambda)])
  expect_true(res$t1_star %in% t1_seq)
})

test_that("chi_sq_lambda equals OLS residuals plus (beta_0 + beta_1 * T2)", {
  ## Derived identity: lambda regularization subtracts the OLS trend up to a
  ## constant shift.
  roll <- build_rolling(t1_seq = seq(1, 60, by = 6), hold_out = 10, n = 150)
  res <- lppls_lagrange(roll, p = 0)

  resid <- unname(stats::residuals(res$trend))
  expected <- resid + res$intercept + res$slope * res$T2
  expect_equal(res$table$chi_sq_lambda, expected, tolerance = 1e-10)

  ## Consequently t1_star is also the T1 of the minimum residual:
  expect_equal(res$t1_star, res$table$T1[which.min(resid)])
})

test_that("the minimizer is invariant to vertical shifts of chi_sq_lambda", {
  ## chi^2_lambda(T1) = chi^2_np(T1) + beta1 * T2 - beta1 * T1; the beta1 * T2
  ## term is constant in T1 and so does not move the argmin.
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30, 40), n = 120)
  res <- lppls_lagrange(roll, p = 0)

  alt <- res$table$chi_sq_np - res$slope * res$table$T1
  expect_equal(res$table$T1[which.min(alt)], res$t1_star)
})

test_that("lambda has the expected sign for a decreasing chi_sq_np trend", {
  ## Build a stub whose per-window SSE makes chi^2_np (= SSE / n_cal at p = 0)
  ## fall as T1 grows: offset A from the window mean by sqrt(n_cal), so the sum
  ## of squared offsets is n_cal^2 and chi^2_np = n_cal, which decreases in T1.
  stub_decreasing <- function(log_price, fh = 30, hold_out = 0,
                              tc_init = NULL, ...) {
    n_total <- length(log_price)
    n_cal <- n_total - hold_out
    A_hat <- mean(log_price[seq_len(n_cal)]) - sqrt(n_cal)
    list(fit = list(list(
      tc = if (!is.null(tc_init)) tc_init else (n_cal + floor(fh / 2)),
      m = 0.5, omega = 8, A = A_hat, B = 0, C1 = 0, C2 = 0, D = 1, value = 0
    )))
  }

  roll <- build_rolling(stub = stub_decreasing, t1_seq = seq(1, 50, by = 10),
                        n = 120)
  res <- lppls_lagrange(roll, p = 0)

  expect_true(res$slope < 0)
  expect_true(res$lambda > 0)
})

# ----------------------------------------------------------------------------
# print method
# ----------------------------------------------------------------------------

test_that("print.lppls_lagrange returns its argument invisibly", {
  roll <- build_rolling(t1_seq = c(1, 10, 20), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  out <- capture.output(returned <- print(res))
  expect_identical(returned, res)
  expect_true(any(grepl("LPPLS Lagrange regularization", out)))
  expect_true(any(grepl("Optimal T1\\*", out)))
  expect_true(any(grepl("lambda\\*", out)))
})

# ----------------------------------------------------------------------------
# Plot methods: chi_sq_np_plot / chi_sq_lambda_plot
# ----------------------------------------------------------------------------

test_that("chi_sq_np_plot and chi_sq_lambda_plot reject non-lppls_lagrange input", {
  expect_error(chi_sq_np_plot(list()),
               "must be an object of class 'lppls_lagrange'")
  expect_error(chi_sq_lambda_plot(list()),
               "must be an object of class 'lppls_lagrange'")
  expect_error(chi_sq_np_plot(1:10),
               "must be an object of class 'lppls_lagrange'")
})

test_that("chi_sq_np_plot returns a ggplot mapping T1 -> chi_sq_np over x$table", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_np_plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, res$table)
  expect_equal(rlang::as_label(p$mapping$x), "T1")
  expect_equal(rlang::as_label(p$mapping$y), "chi_sq_np")
})

test_that("chi_sq_np_plot's geom_abline has the fitted intercept and slope", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_np_plot(res)
  abline_layer <- find_layer(p, "GeomAbline")
  expect_equal(abline_layer$data$intercept, res$intercept, tolerance = 1e-14)
  expect_equal(abline_layer$data$slope, res$slope, tolerance = 1e-14)
  expect_equal(abline_layer$aes_params$colour, "red")
})

test_that("chi_sq_np_plot has exactly one geom_point and one geom_abline", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_np_plot(res)
  geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geom_classes == "GeomPoint"), 1L)
  expect_equal(sum(geom_classes == "GeomAbline"), 1L)
})

test_that("chi_sq_lambda_plot returns a ggplot mapping T1 -> chi_sq_lambda", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_lambda_plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, res$table)
  expect_equal(rlang::as_label(p$mapping$x), "T1")
  expect_equal(rlang::as_label(p$mapping$y), "chi_sq_lambda")
})

test_that("chi_sq_lambda_plot's geom_vline xintercept equals t1_star", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_lambda_plot(res)
  vline_layer <- find_layer(p, "GeomVline")
  expect_equal(vline_layer$data$xintercept, res$t1_star)
  expect_equal(vline_layer$aes_params$colour, "green3")
})

test_that("chi_sq_lambda_plot has exactly one geom_point and one geom_vline", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  p <- chi_sq_lambda_plot(res)
  geom_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geom_classes == "GeomPoint"), 1L)
  expect_equal(sum(geom_classes == "GeomVline"), 1L)
})

test_that("chi^2 plots build without errors or warnings", {
  roll <- build_rolling(t1_seq = c(1, 10, 20, 30), hold_out = 5, n = 80)
  res <- lppls_lagrange(roll, p = 0)

  expect_silent(ggplot2::ggplot_build(chi_sq_np_plot(res)))
  expect_silent(ggplot2::ggplot_build(chi_sq_lambda_plot(res)))
})

# ----------------------------------------------------------------------------
# Integration test: real fit_lppls, end-to-end rolling -> lagrange
# ----------------------------------------------------------------------------

test_that("end-to-end run with the real fit_lppls is internally consistent", {
  skip_on_cran()

  data <- make_bubble_data(n = 220, tc = 240)
  t1_seq <- seq(1, 60, by = 15)   # 5 windows
  hold_out <- 20
  fh <- 30
  p_val <- 7

  roll <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    fit_args = list(
      lower = c(0.1, 6, -0.03, 0.8),
      upper = c(0.9, 13, -1e-14, 1e6),
      m_init = 0.5, o_init = 8,
      num_searches = 1, mode = "F2"
    )
  )
  res <- lppls_lagrange(roll, p = p_val)

  expect_s3_class(res, "lppls_lagrange")
  expect_equal(res$T2, length(data$log_p) - hold_out)
  expect_equal(res$table$T1, t1_seq)

  ## chi_sq_np recomputed from the rolling object's stored fits
  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    window <- data$log_p[T1:length(data$log_p)]
    n_cal <- length(window) - hold_out
    expect_equal(res$table$n[k], n_cal)

    sse_k <- SSE(roll$fits[[k]], window[seq_len(n_cal)], seq_len(n_cal),
                 mode = 0)
    expect_equal(res$table$chi_sq_np[k], sse_k / (n_cal - p_val),
                 tolerance = 1e-10)
  }

  ## Trend coefficients match a manual lm()
  manual <- stats::lm(res$table$chi_sq_np ~ res$table$T1)
  expect_equal(res$intercept, unname(stats::coef(manual)[[1]]),
               tolerance = 1e-10)
  expect_equal(res$slope, unname(stats::coef(manual)[[2]]),
               tolerance = 1e-10)
  expect_equal(res$lambda, -res$slope, tolerance = 1e-14)

  ## chi_sq_lambda matches its definition
  expected_chi_lambda <- res$table$chi_sq_np -
    res$lambda * (res$T2 - res$table$T1)
  expect_equal(res$table$chi_sq_lambda, expected_chi_lambda,
               tolerance = 1e-10)

  ## t1_star is the argmin and lives inside t1_seq
  expect_equal(res$t1_star, t1_seq[which.min(res$table$chi_sq_lambda)])

  ## Plots build cleanly on real output
  expect_silent(ggplot2::ggplot_build(chi_sq_np_plot(res)))
  expect_silent(ggplot2::ggplot_build(chi_sq_lambda_plot(res)))
})
