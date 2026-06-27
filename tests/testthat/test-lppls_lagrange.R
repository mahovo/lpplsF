# Tests for lppls_lagrange()

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

## Replace lpplsF::fit_lppls with a deterministic stub for the duration
## of the calling test_that() block. Returns the deterministic best fit
## that lppls_lagrange() will see for each window.
##
## The stub records every call in .calls so tests can inspect what arguments
## lppls_lagrange forwarded.
local_mock_fit_lppls <- function(stub, .calls = NULL,
                                 envir = parent.frame()) {
  ns <- asNamespace("lpplsF")
  original <- ns$fit_lppls
  was_locked <- bindingIsLocked("fit_lppls", ns)

  if (was_locked) unlockBinding("fit_lppls", ns)

  wrapped <- function(...) {
    args <- list(...)
    if (!is.null(.calls)) {
      .calls$args[[length(.calls$args) + 1L]] <- args
    }
    stub(...)
  }

  assign("fit_lppls", wrapped, envir = ns)

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

## A simple stub that returns a fixed best fit using the window mean as A
## and zero oscillation amplitudes. This makes SSE deterministic.
make_fixed_stub <- function() {
  function(log_price, fh = 30, hold_out = 0, tc_init = NULL,
           mode = "F2", ...) {
    n_total <- length(log_price)
    n_cal <- n_total - hold_out
    A_hat <- mean(log_price[seq_len(n_cal)])
    list(
      fit = list(
        list(
          tc = if (!is.null(tc_init)) tc_init else (n_cal + floor(fh / 2)),
          m = 0.5,
          omega = 8,
          A = A_hat,
          B = -0.01,
          C1 = 0,
          C2 = 0,
          D = 1,
          value = 0
        )
      )
    )
  }
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

test_that("lppls_lagrange rejects non-numeric log_price", {
  expect_error(
    lppls_lagrange(log_price = "abc", t1_seq = c(1, 2)),
    "'log_price' must be a numeric vector"
  )
})

test_that("lppls_lagrange rejects t1_seq that is not numeric or too short", {
  log_p <- rnorm(50)

  expect_error(
    lppls_lagrange(log_price = log_p, t1_seq = "abc"),
    "'t1_seq' must be a numeric vector of length >= 2"
  )

  expect_error(
    lppls_lagrange(log_price = log_p, t1_seq = c(5)),
    "'t1_seq' must be a numeric vector of length >= 2"
  )
})

test_that("lppls_lagrange rejects non-positive or non-integer t1_seq", {
  log_p <- rnorm(50)

  expect_error(
    lppls_lagrange(log_price = log_p, t1_seq = c(0, 5)),
    "'t1_seq' must contain positive integers"
  )

  expect_error(
    lppls_lagrange(log_price = log_p, t1_seq = c(-1, 5)),
    "'t1_seq' must contain positive integers"
  )

  expect_error(
    lppls_lagrange(log_price = log_p, t1_seq = c(1.5, 5)),
    "'t1_seq' must contain positive integers"
  )
})

test_that("lppls_lagrange rejects t1_seq values >= T2", {
  log_p <- rnorm(50)
  ## T2 = 50 - 10 = 40, so 45 must fail
  expect_error(
    lppls_lagrange(
      log_price = log_p,
      t1_seq = c(5, 45),
      hold_out = 10
    ),
    "must be < T2"
  )

  ## T2 = 50 - 0 = 50, so 50 (==T2) must fail (boundary)
  expect_error(
    lppls_lagrange(
      log_price = log_p,
      t1_seq = c(5, 50),
      hold_out = 0
    ),
    "must be < T2"
  )
})

test_that("lppls_lagrange errors when n_cal - p is non-positive", {
  log_p <- rnorm(50)
  ## A window starting at T1=44 with hold_out=0 has n_cal=7, so p=7
  ## gives denom=0 and must error.
  expect_error(
    suppressWarnings(lppls_lagrange(
      log_price = log_p,
      t1_seq = c(1, 44),
      hold_out = 0,
      p = 7,
      fit_args = list(num_searches = 1)
    )),
    "Non-positive denominator"
  )
})

# ----------------------------------------------------------------------------
# Reserved fit_args handling
# ----------------------------------------------------------------------------

test_that("lppls_lagrange warns about reserved fit_args and strips them", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  expect_warning(
    lppls_lagrange(
      log_price = data$log_p,
      t1_seq = c(1, 10),
      hold_out = 5,
      fh = 20,
      p = 0,
      fit_args = list(log_price = "bad", fh = 999, hold_out = 999,
                      tc_init = 999, fb = TRUE,
                      num_searches = 1)
    ),
    "controlled by lppls_lagrange"
  )

  ## The actual fit_lppls calls should have used lppls_lagrange's values,
  ## not the bogus ones from fit_args.
  expect_true(length(calls$args) >= 1L)
  first <- calls$args[[1L]]
  expect_equal(first$fh, 20)
  expect_equal(first$hold_out, 5)
  expect_false(identical(first$log_price, "bad"))
})

test_that("lppls_lagrange defaults mode to 'F2' when not provided", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  for (a in calls$args) expect_equal(a$mode, "F2")
})

test_that("lppls_lagrange honours user-supplied mode in fit_args", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1, mode = "MPL")
  )

  for (a in calls$args) expect_equal(a$mode, "MPL")
})

# ----------------------------------------------------------------------------
# Per-window argument forwarding
# ----------------------------------------------------------------------------

test_that("lppls_lagrange slices log_price[T1:n_full] for each window", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  t1_seq <- c(1, 5, 15)
  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_length(calls$args, length(t1_seq))
  for (k in seq_along(t1_seq)) {
    expect_equal(
      calls$args[[k]]$log_price,
      data$log_p[t1_seq[k]:length(data$log_p)]
    )
  }
})

test_that("default tc_init per window is n_cal + floor(fh / 2)", {
  data <- make_bubble_data(n = 80)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  t1_seq <- c(1, 10, 20)
  hold_out <- 5
  fh <- 21
  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  for (k in seq_along(t1_seq)) {
    window_len <- length(data$log_p) - t1_seq[k] + 1
    n_cal <- window_len - hold_out
    expect_equal(calls$args[[k]]$tc_init, n_cal + floor(fh / 2))
  }
})

test_that("explicit tc_init overrides the per-window default", {
  data <- make_bubble_data(n = 80)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 5, 10),
    hold_out = 5,
    fh = 30,
    tc_init = 1234,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  for (a in calls$args) expect_equal(a$tc_init, 1234)
})

test_that("fit_lppls receives hold_out and fh as supplied to lppls_lagrange", {
  data <- make_bubble_data(n = 80)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 5, 10),
    hold_out = 7,
    fh = 11,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  for (a in calls$args) {
    expect_equal(a$hold_out, 7)
    expect_equal(a$fh, 11)
  }
})

# ----------------------------------------------------------------------------
# Return value structure
# ----------------------------------------------------------------------------

test_that("lppls_lagrange returns an lppls_lagrange-classed list with all fields", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_s3_class(res, "lppls_lagrange")
  expect_named(
    res,
    c("table", "fits", "trend", "intercept", "slope", "lambda",
      "T2", "p", "t1_star"),
    ignore.order = TRUE
  )

  expect_s3_class(res$table, "data.frame")
  expect_named(res$table, c("T1", "n", "chi_sq_np", "chi_sq_lambda"))
  expect_equal(nrow(res$table), length(t1_seq))

  expect_length(res$fits, length(t1_seq))
  expect_s3_class(res$trend, "lm")

  expect_length(res$intercept, 1L)
  expect_length(res$slope, 1L)
  expect_length(res$lambda, 1L)
  expect_length(res$T2, 1L)
  expect_length(res$p, 1L)
  expect_length(res$t1_star, 1L)
})

test_that("T2 equals length(log_price) - hold_out", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  for (hold_out in c(0, 5, 10)) {
    res <- lppls_lagrange(
      log_price = data$log_p,
      t1_seq = c(1, 10, 20),
      hold_out = hold_out,
      fh = 20,
      p = 0,
      fit_args = list(num_searches = 1)
    )
    expect_equal(res$T2, length(data$log_p) - hold_out)
  }
})

test_that("table$n[k] equals (n_full - T1[k] + 1) - hold_out", {
  data <- make_bubble_data(n = 100)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 11, 21, 31)
  hold_out <- 7
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 15,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  for (k in seq_along(t1_seq)) {
    window_len <- length(data$log_p) - t1_seq[k] + 1
    expect_equal(res$table$n[k], window_len - hold_out)
  }
})

test_that("table$T1 matches the supplied t1_seq exactly", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(2, 7, 13, 19)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 15,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_equal(res$table$T1, t1_seq)
})

# ----------------------------------------------------------------------------
# Mathematical correctness: chi^2_np, trend, lambda, chi^2_lambda, t1*
# ----------------------------------------------------------------------------

test_that("chi_sq_np[k] = SSE(fit_k, window[1:n_cal], 1:n_cal, mode=0) / (n_cal - p)", {
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 25, 40)
  hold_out <- 8
  p_val <- 7
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 20,
    p = p_val,
    fit_args = list(num_searches = 1)
  )

  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    window <- data$log_p[T1:length(data$log_p)]
    n_cal <- length(window) - hold_out
    sse_k <- SSE(res$fits[[k]], window[seq_len(n_cal)], seq_len(n_cal),
                 mode = 0)
    expect_equal(res$table$chi_sq_np[k], sse_k / (n_cal - p_val),
                 tolerance = 1e-12)
  }
})

test_that("p = 0 reproduces the thesis form (no parameter correction)", {
  data <- make_bubble_data(n = 100)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30)
  hold_out <- 5
  res0 <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )
  res7 <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 20,
    p = 7,
    fit_args = list(num_searches = 1)
  )

  ## Same fits + same windows → same SSEs; chi_np differs only by denom.
  for (k in seq_along(t1_seq)) {
    n_cal <- res0$table$n[k]
    expect_equal(
      res0$table$chi_sq_np[k] * n_cal,
      res7$table$chi_sq_np[k] * (n_cal - 7),
      tolerance = 1e-12
    )
  }
})

test_that("intercept and slope equal coefficients of lm(chi_sq_np ~ T1)", {
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30, 40)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  manual <- stats::lm(res$table$chi_sq_np ~ res$table$T1)
  expect_equal(res$intercept, unname(stats::coef(manual)[[1]]),
               tolerance = 1e-12)
  expect_equal(res$slope, unname(stats::coef(manual)[[2]]),
               tolerance = 1e-12)
})

test_that("lambda equals -slope", {
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30, 40),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_equal(res$lambda, -res$slope, tolerance = 1e-14)
})

test_that("chi_sq_lambda = chi_sq_np - lambda * (T2 - T1)", {
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30, 40)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expected <- res$table$chi_sq_np - res$lambda * (res$T2 - res$table$T1)
  expect_equal(res$table$chi_sq_lambda, expected, tolerance = 1e-12)
})

test_that("t1_star = t1_seq[which.min(chi_sq_lambda)]", {
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30, 40)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_equal(res$t1_star, t1_seq[which.min(res$table$chi_sq_lambda)])
  expect_true(res$t1_star %in% t1_seq)
})

test_that("chi_sq_lambda equals OLS residuals plus (beta_0 + beta_1 * T2)", {
  ## This is a derived identity: lambda regularization is equivalent to
  ## subtracting the OLS trend up to a constant shift. Verifies the math.
  data <- make_bubble_data(n = 150)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- seq(1, 60, by = 6)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 10,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  resid <- unname(stats::residuals(res$trend))
  expected <- resid + res$intercept + res$slope * res$T2
  expect_equal(res$table$chi_sq_lambda, expected, tolerance = 1e-10)

  ## Consequently t1_star is also the T1 of the minimum residual:
  expect_equal(res$t1_star, t1_seq[which.min(resid)])
})

test_that("the minimizer is invariant to vertical shifts of chi_sq_lambda", {
  ## chi^2_lambda(T1) = chi^2_np(T1) + beta1 * T2 - beta1 * T1
  ## The beta1 * T2 term is constant in T1 and so does not move the argmin.
  data <- make_bubble_data(n = 120)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30, 40)
  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  ## Form with constant `beta1 * T2` term dropped:
  ##   chi_lambda - beta1*T2 = chi_np - beta1*T1
  alt <- res$table$chi_sq_np - res$slope * res$table$T1
  expect_equal(t1_seq[which.min(alt)], res$t1_star)
})

test_that("lambda has the expected sign for a decreasing chi_sq_np trend", {
  ## If chi^2_np decreases in T1 (slope < 0), then lambda* = -beta1 > 0.
  ## Build a stub that returns SSE proportional to (n_cal)^2, so chi^2_np
  ## (= SSE / n_cal for p=0) grows linearly with n_cal; since n_cal
  ## decreases as T1 increases, chi^2_np decreases with T1.
  data <- make_bubble_data(n = 120)

  stub_decreasing <- function(log_price, fh = 30, hold_out = 0,
                              tc_init = NULL, ...) {
    n_total <- length(log_price)
    n_cal <- n_total - hold_out
    ## Choose B such that residuals scale with n_cal, yielding
    ## chi^2_np = SSE / n_cal increasing in n_cal => decreasing in T1.
    ##
    ## Simpler: return a fit whose A is offset from window mean by an
    ## amount proportional to sqrt(n_cal). Then sum of squared offsets
    ## = n_cal * (sqrt(n_cal))^2 = n_cal^2, and chi^2_np = n_cal.
    A_hat <- mean(log_price[seq_len(n_cal)]) - sqrt(n_cal)
    list(fit = list(list(
      tc = if (!is.null(tc_init)) tc_init else (n_cal + floor(fh / 2)),
      m = 0.5, omega = 8, A = A_hat,
      B = 0, C1 = 0, C2 = 0, D = 1, value = 0
    )))
  }

  local_mock_fit_lppls(stub_decreasing)

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = seq(1, 50, by = 10),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  ## chi^2_np ≈ n_cal which decreases as T1 grows → slope < 0, lambda > 0.
  expect_true(res$slope < 0)
  expect_true(res$lambda > 0)
})

# ----------------------------------------------------------------------------
# print method
# ----------------------------------------------------------------------------

test_that("print.lppls_lagrange returns its argument invisibly", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  out <- capture.output(returned <- print(res))
  expect_identical(returned, res)
  expect_true(any(grepl("LPPLS Lagrange regularization", out)))
  expect_true(any(grepl("Optimal T1\\*", out)))
  expect_true(any(grepl("lambda\\*", out)))
})

# ----------------------------------------------------------------------------
# Plot methods
# ----------------------------------------------------------------------------

test_that("chi_sq_np_plot and chi_sq_lambda_plot reject non-lppls_lagrange input", {
  expect_error(chi_sq_np_plot(list()),
               "must be an object of class 'lppls_lagrange'")
  expect_error(chi_sq_lambda_plot(list()),
               "must be an object of class 'lppls_lagrange'")
  expect_error(chi_sq_np_plot(1:10),
               "must be an object of class 'lppls_lagrange'")
})

test_that("chi_sq_np_plot returns a ggplot whose data is x$table", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_np_plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, res$table)
})

test_that("chi_sq_np_plot maps x to T1 and y to chi_sq_np", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_np_plot(res)
  expect_equal(rlang::as_label(p$mapping$x), "T1")
  expect_equal(rlang::as_label(p$mapping$y), "chi_sq_np")
})

test_that("chi_sq_np_plot's geom_abline has the fitted intercept and slope", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_np_plot(res)
  abline_layer <- find_layer(p, "GeomAbline")
  expect_equal(abline_layer$data$intercept, res$intercept, tolerance = 1e-14)
  expect_equal(abline_layer$data$slope, res$slope, tolerance = 1e-14)
  expect_equal(abline_layer$aes_params$colour, "red")
})

test_that("chi_sq_np_plot has exactly one geom_point and one geom_abline", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_np_plot(res)
  geom_classes <- vapply(p$layers,
                         function(l) class(l$geom)[1],
                         character(1))
  expect_equal(sum(geom_classes == "GeomPoint"), 1L)
  expect_equal(sum(geom_classes == "GeomAbline"), 1L)
})

test_that("chi_sq_lambda_plot returns a ggplot whose data is x$table", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_lambda_plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, res$table)
})

test_that("chi_sq_lambda_plot maps x to T1 and y to chi_sq_lambda", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_lambda_plot(res)
  expect_equal(rlang::as_label(p$mapping$x), "T1")
  expect_equal(rlang::as_label(p$mapping$y), "chi_sq_lambda")
})

test_that("chi_sq_lambda_plot's geom_vline xintercept equals t1_star", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_lambda_plot(res)
  vline_layer <- find_layer(p, "GeomVline")
  expect_equal(vline_layer$data$xintercept, res$t1_star)
  expect_equal(vline_layer$aes_params$colour, "green3")
})

test_that("chi_sq_lambda_plot has exactly one geom_point and one geom_vline", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  p <- chi_sq_lambda_plot(res)
  geom_classes <- vapply(p$layers,
                         function(l) class(l$geom)[1],
                         character(1))
  expect_equal(sum(geom_classes == "GeomPoint"), 1L)
  expect_equal(sum(geom_classes == "GeomVline"), 1L)
})

test_that("plots build without errors or warnings", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    p = 0,
    fit_args = list(num_searches = 1)
  )

  expect_silent(ggplot2::ggplot_build(chi_sq_np_plot(res)))
  expect_silent(ggplot2::ggplot_build(chi_sq_lambda_plot(res)))
})

# ----------------------------------------------------------------------------
# Integration test: real fit_lppls
# ----------------------------------------------------------------------------

test_that("end-to-end run with the real fit_lppls is internally consistent", {
  skip_on_cran()

  data <- make_bubble_data(n = 220, tc = 240)
  t1_seq <- seq(1, 60, by = 15)   # 5 windows
  hold_out <- 20
  fh <- 30
  p_val <- 7

  res <- lppls_lagrange(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    p = p_val,
    fit_args = list(
      lower = c(0.1, 6, -0.03, 0.8),
      upper = c(0.9, 13, -1e-14, 1e6),
      m_init = 0.5, o_init = 8,
      num_searches = 1, mode = "F2"
    )
  )

  expect_s3_class(res, "lppls_lagrange")
  expect_equal(res$T2, length(data$log_p) - hold_out)
  expect_length(res$fits, length(t1_seq))
  expect_equal(res$table$T1, t1_seq)

  ## chi_sq_np values recomputed from stored fits
  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    window <- data$log_p[T1:length(data$log_p)]
    n_cal <- length(window) - hold_out
    expect_equal(res$table$n[k], n_cal)

    sse_k <- SSE(res$fits[[k]], window[seq_len(n_cal)], seq_len(n_cal),
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
  expect_equal(res$t1_star,
               t1_seq[which.min(res$table$chi_sq_lambda)])

  ## Plots build cleanly on real output
  expect_silent(ggplot2::ggplot_build(chi_sq_np_plot(res)))
  expect_silent(ggplot2::ggplot_build(chi_sq_lambda_plot(res)))
})
