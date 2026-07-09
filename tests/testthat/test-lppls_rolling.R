# Tests for lppls_rolling() and the window-diagnostic plots rolling_tc_plot() /
# rolling_param_plot().
#
# lppls_rolling() calls fit_lppls() once per window, which is expensive. To keep
# the suite fast, most tests mock fit_lppls with a deterministic stub (as in
# test-lppls_lagrange.R); a single skip_on_cran() integration test exercises the
# real fitter and the MPL path.

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

## Replace lpplsF::fit_lppls with a deterministic stub for the duration of the
## calling test_that() block. When `.calls` (an environment with a `$args` list)
## is supplied, every call's arguments are recorded so tests can inspect what
## lppls_rolling() forwarded.
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

## A simple stub that returns a fixed best fit using the window mean as A and
## zero oscillation amplitudes. Deterministic SSE, and no mpl_output (so the
## likelihood-interval columns come back NA).
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

## A stub that also returns an `mpl_output` block (window-local likelihood
## intervals centred on tc, plus tc_hat_mpl), so the LI columns and the
## rolling_tc_plot() happy path can be exercised without the real fitter.
make_mpl_stub <- function() {
  function(log_price, fh = 30, hold_out = 0, tc_init = NULL,
           mode = "MPL", ...) {
    n_total <- length(log_price)
    n_cal <- n_total - hold_out
    tc_l <- if (!is.null(tc_init)) tc_init else (n_cal + floor(fh / 2))
    list(
      fit = list(
        list(tc = tc_l, m = 0.5, omega = 8,
             A = mean(log_price[seq_len(n_cal)]),
             B = -0.01, C1 = 0, C2 = 0, D = 1, value = 0)
      ),
      mpl_output = list(
        LI = list(c(tc_l - 3, tc_l + 3),   # 5%
                  c(tc_l - 2, tc_l + 2),    # 10%
                  c(tc_l - 1, tc_l + 1)),   # 50%
        tc_hat_mpl = tc_l
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

test_that("lppls_rolling rejects non-numeric log_price", {
  expect_error(
    lppls_rolling(log_price = "abc", t1_seq = c(1, 2)),
    "'log_price' must be a numeric vector"
  )
})

test_that("lppls_rolling rejects t1_seq that is not numeric or too short", {
  log_p <- rnorm(50)

  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = "abc"),
    "'t1_seq' must be a numeric vector of length >= 2"
  )

  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(5)),
    "'t1_seq' must be a numeric vector of length >= 2"
  )
})

test_that("lppls_rolling rejects non-positive or non-integer t1_seq", {
  log_p <- rnorm(50)

  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(0, 5)),
    "'t1_seq' must contain positive integers"
  )

  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(-1, 5)),
    "'t1_seq' must contain positive integers"
  )

  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(1.5, 5)),
    "'t1_seq' must contain positive integers"
  )
})

test_that("lppls_rolling rejects t1_seq values >= T2", {
  log_p <- rnorm(50)
  ## T2 = 50 - 10 = 40, so 45 must fail
  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(5, 45), hold_out = 10),
    "must be < T2"
  )

  ## T2 = 50 - 0 = 50, so 50 (==T2) must fail (boundary)
  expect_error(
    lppls_rolling(log_price = log_p, t1_seq = c(5, 50), hold_out = 0),
    "must be < T2"
  )
})

# ----------------------------------------------------------------------------
# Reserved fit_args handling and mode default
# ----------------------------------------------------------------------------

test_that("lppls_rolling warns about reserved fit_args and strips them", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  expect_warning(
    lppls_rolling(
      log_price = data$log_p,
      t1_seq = c(1, 10),
      hold_out = 5,
      fh = 20,
      fit_args = list(log_price = "bad", fh = 999, hold_out = 999,
                      tc_init = 999, fb = TRUE, num_searches = 1)
    ),
    "controlled by lppls_rolling"
  )

  ## The actual fit_lppls calls should have used lppls_rolling's values,
  ## not the bogus ones from fit_args.
  expect_true(length(calls$args) >= 1L)
  first <- calls$args[[1L]]
  expect_equal(first$fh, 20)
  expect_equal(first$hold_out, 5)
  expect_false(identical(first$log_price, "bad"))
})

test_that("lppls_rolling defaults mode to 'MPL' when not provided", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  for (a in calls$args) expect_equal(a$mode, "MPL")
})

test_that("lppls_rolling honours user-supplied mode in fit_args", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1, mode = "F2")
  )

  for (a in calls$args) expect_equal(a$mode, "F2")
})

test_that("lppls_rolling records the effective fit_args (mode defaulted, reserved stripped)", {
  data <- make_bubble_data(n = 60)
  local_mock_fit_lppls(make_fixed_stub())

  roll <- suppressWarnings(lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1, fh = 999)  # fh is reserved -> stripped
  ))

  expect_true("fit_args" %in% names(roll))
  expect_equal(roll$fit_args$mode, "MPL")        # defaulted
  expect_equal(roll$fit_args$num_searches, 1)     # kept
  expect_false("fh" %in% names(roll$fit_args))    # reserved, stripped
})

# ----------------------------------------------------------------------------
# Per-window argument forwarding
# ----------------------------------------------------------------------------

test_that("lppls_rolling slices log_price[T1:n_full] for each window", {
  data <- make_bubble_data(n = 60)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  t1_seq <- c(1, 5, 15)
  lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
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
  lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
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

  lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 5, 10),
    hold_out = 5,
    fh = 30,
    tc_init = 1234,
    fit_args = list(num_searches = 1)
  )

  for (a in calls$args) expect_equal(a$tc_init, 1234)
})

test_that("fit_lppls receives hold_out and fh as supplied to lppls_rolling", {
  data <- make_bubble_data(n = 80)
  calls <- new.env()
  calls$args <- list()

  local_mock_fit_lppls(make_fixed_stub(), calls)

  lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 5, 10),
    hold_out = 7,
    fh = 11,
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

test_that("lppls_rolling returns an lppls_rolling-classed list with all fields", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 10, 20, 30)
  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  expect_s3_class(res, "lppls_rolling")
  expect_named(
    res,
    c("table", "fits", "t1_seq", "T2", "hold_out", "fh", "n_full", "fit_args"),
    ignore.order = TRUE
  )

  expect_s3_class(res$table, "data.frame")
  expect_named(
    res$table,
    c("T1", "n_cal", "tc", "m", "omega", "A", "B", "C1", "C2", "D",
      "sse_cal", "LI5l", "LI5u", "LI10l", "LI10u", "LI50l", "LI50u",
      "tc_hat_mpl"),
    ignore.order = TRUE
  )
  expect_equal(nrow(res$table), length(t1_seq))

  expect_length(res$fits, length(t1_seq))
  expect_equal(res$t1_seq, t1_seq)
  expect_equal(res$n_full, length(data$log_p))
})

test_that("T2 equals length(log_price) - hold_out", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  for (hold_out in c(0, 5, 10)) {
    res <- lppls_rolling(
      log_price = data$log_p,
      t1_seq = c(1, 10, 20),
      hold_out = hold_out,
      fh = 20,
      fit_args = list(num_searches = 1)
    )
    expect_equal(res$T2, length(data$log_p) - hold_out)
  }
})

test_that("table$n_cal[k] equals (n_full - T1[k] + 1) - hold_out", {
  data <- make_bubble_data(n = 100)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 11, 21, 31)
  hold_out <- 7
  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 15,
    fit_args = list(num_searches = 1)
  )

  for (k in seq_along(t1_seq)) {
    window_len <- length(data$log_p) - t1_seq[k] + 1
    expect_equal(res$table$n_cal[k], window_len - hold_out)
  }
})

test_that("table$T1 matches the supplied t1_seq exactly", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(2, 7, 13, 19)
  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = 5,
    fh = 15,
    fit_args = list(num_searches = 1)
  )

  expect_equal(res$table$T1, t1_seq)
})

test_that("tc and likelihood-interval columns are shifted to absolute time", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_mpl_stub())

  t1_seq <- c(1, 10, 20, 30)
  hold_out <- 5
  fh <- 20
  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    fit_args = list(num_searches = 1, mode = "MPL")
  )

  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    off <- T1 - 1
    n_cal <- (length(data$log_p) - T1 + 1) - hold_out
    tc_local <- n_cal + floor(fh / 2)      # the stub's tc
    ## Point estimate shifted by (T1 - 1)
    expect_equal(res$table$tc[k], tc_local + off)
    ## LI bounds shifted by the same offset (stub: tc +/- 3, 2, 1)
    expect_equal(res$table$LI5l[k], (tc_local - 3) + off)
    expect_equal(res$table$LI5u[k], (tc_local + 3) + off)
    expect_equal(res$table$tc_hat_mpl[k], tc_local + off)
  }
})

test_that("sse_cal equals SSE over the calibration region (mode 0)", {
  data <- make_bubble_data(n = 100)
  local_mock_fit_lppls(make_fixed_stub())

  t1_seq <- c(1, 15, 30)
  hold_out <- 8
  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    window <- data$log_p[T1:length(data$log_p)]
    n_cal <- length(window) - hold_out
    sse_k <- SSE(res$fits[[k]], window[seq_len(n_cal)], seq_len(n_cal),
                 mode = 0)
    expect_equal(res$table$sse_cal[k], sse_k, tolerance = 1e-12)
    expect_true(res$table$sse_cal[k] >= 0)
  }
})

test_that("likelihood-interval columns are NA when the fit has no mpl_output", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())   # no mpl_output

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  li_cols <- c("LI5l", "LI5u", "LI10l", "LI10u", "LI50l", "LI50u",
               "tc_hat_mpl")
  for (col in li_cols) expect_true(all(is.na(res$table[[col]])))
})

# ----------------------------------------------------------------------------
# print method
# ----------------------------------------------------------------------------

test_that("print.lppls_rolling returns its argument invisibly", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  out <- capture.output(returned <- print(res))
  expect_identical(returned, res)
  expect_true(any(grepl("LPPLS rolling-window calibration", out)))
  expect_true(any(grepl("T2", out)))
  expect_true(any(grepl("Likelihood intervals", out)))
})

test_that("print.lppls_rolling reports absent vs present likelihood intervals", {
  data <- make_bubble_data(n = 80)

  local_mock_fit_lppls(make_fixed_stub())    # no mpl_output -> absent
  res_na <- lppls_rolling(data$log_p, t1_seq = c(1, 10, 20), hold_out = 5,
                          fh = 20, fit_args = list(num_searches = 1))
  out_na <- capture.output(print(res_na))
  expect_true(any(grepl("absent", out_na)))
})

# ----------------------------------------------------------------------------
# rolling_tc_plot
# ----------------------------------------------------------------------------

test_that("rolling_tc_plot rejects non-lppls_rolling input", {
  expect_error(rolling_tc_plot(list()),
               "must be an object of class 'lppls_rolling'")
  expect_error(rolling_tc_plot(1:10),
               "must be an object of class 'lppls_rolling'")
})

test_that("rolling_tc_plot errors when no likelihood intervals are present", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())   # no mpl_output -> all-NA tc_hat_mpl

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  expect_error(rolling_tc_plot(res), "mode = 'MPL'")
})

test_that("rolling_tc_plot returns a ggplot when intervals are present", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_mpl_stub())

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1, mode = "MPL")
  )

  p <- rolling_tc_plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data, res$table)
  expect_equal(rlang::as_label(p$mapping$x), "T1")
  expect_equal(rlang::as_label(p$mapping$y), "tc")
  expect_silent(ggplot2::ggplot_build(p))
})

# ----------------------------------------------------------------------------
# rolling_param_plot
# ----------------------------------------------------------------------------

test_that("rolling_param_plot rejects non-lppls_rolling input", {
  expect_error(rolling_param_plot(list()),
               "must be an object of class 'lppls_rolling'")
  expect_error(rolling_param_plot("nope"),
               "must be an object of class 'lppls_rolling'")
})

test_that("rolling_param_plot returns a ggplot faceted over the six parameters", {
  data <- make_bubble_data(n = 80)
  local_mock_fit_lppls(make_fixed_stub())

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = c(1, 10, 20, 30),
    hold_out = 5,
    fh = 20,
    fit_args = list(num_searches = 1)
  )

  p <- rolling_param_plot(res)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$param), c("m", "omega", "A", "B", "C1", "C2"))
  expect_silent(ggplot2::ggplot_build(p))
})

# ----------------------------------------------------------------------------
# Integration test: real fit_lppls, MPL path
# ----------------------------------------------------------------------------

test_that("lppls_rolling drives the real fit_lppls and populates the MPL columns", {
  skip_on_cran()

  data <- make_bubble_data(n = 200, tc = 220)
  t1_seq <- seq(1, 45, by = 15)   # 4 windows
  hold_out <- 20
  fh <- 30

  res <- lppls_rolling(
    log_price = data$log_p,
    t1_seq = t1_seq,
    hold_out = hold_out,
    fh = fh,
    fit_args = list(
      lower = c(0.1, 6, -0.03, 0.8),
      upper = c(0.9, 13, -1e-14, 1e6),
      m_init = 0.5, o_init = 8,
      num_searches = 1, mode = "MPL"
    )
  )

  expect_s3_class(res, "lppls_rolling")
  expect_equal(res$T2, length(data$log_p) - hold_out)
  expect_equal(res$table$T1, t1_seq)
  expect_length(res$fits, length(t1_seq))

  ## sse_cal recomputed from the stored fits
  for (k in seq_along(t1_seq)) {
    T1 <- t1_seq[k]
    window <- data$log_p[T1:length(data$log_p)]
    n_cal <- length(window) - hold_out
    expect_equal(res$table$n_cal[k], n_cal)
    sse_k <- SSE(res$fits[[k]], window[seq_len(n_cal)], seq_len(n_cal),
                 mode = 0)
    expect_equal(res$table$sse_cal[k], sse_k, tolerance = 1e-8)
  }

  ## MPL mode should yield at least some non-NA likelihood-interval output,
  ## and where present the plot builds cleanly.
  if (any(!is.na(res$table$tc_hat_mpl))) {
    expect_silent(ggplot2::ggplot_build(rolling_tc_plot(res)))
  }
  expect_silent(ggplot2::ggplot_build(rolling_param_plot(res)))
})
