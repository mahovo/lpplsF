## Diagnostics reported by print.lppls_fit(). The checks run against
## .fit_diagnostics(), which returns the quantities as values; print() is tested
## only for the lines and flags it derives from them.

## A cheap, well-behaved F2 calibration: interior optimum, few tc candidates.
clean_fit <- function(mode = "F2", ...) {
  t     <- 1:200
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8)
  fit_lppls(log_price = log_p, fh = 10, hold_out = 0,
            lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
            tc_init = 205, m_init = 0.5, o_init = 8, num_searches = 2,
            mode = mode, ...)
}

## A fit forced onto its m bound by squeezing the box around a value the data
## do not want, so the boundary checks have something to find.
pinned_fit <- function() {
  t     <- 1:200
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8)
  fit_lppls(log_price = log_p, fh = 10, hold_out = 0,
            lower = c(0.80, 6, -1e14, 0.8), upper = c(0.85, 13, -1e-14, 1e6),
            tc_init = 205, m_init = 0.82, o_init = 8, num_searches = 2,
            mode = "F2")
}


test_that("fit_lppls returns a classed object carrying its own settings", {
  f <- clean_fit()
  expect_s3_class(f, "lppls_fit")
  expect_true(is.list(f))              # still a plain list underneath
  expect_named(f$fit_args,
               c("mode", "n", "fh", "hold_out", "lower", "upper",
                 "num_searches", "mpl_cutoff"))
  expect_identical(f$fit_args$mode, "F2")
  expect_identical(f$fit_args$n, 200)
  expect_identical(f$fit_args$fh, 10)
  expect_identical(f$fit_args$upper[1], 0.9)
  ## $-access, which existing code relies on, is unaffected by the class.
  expect_true(is.numeric(f$fit[[1]]$tc))
})


test_that("print dispatches and reports the estimate", {
  f   <- clean_fit()
  out <- capture.output(print(f))
  expect_match(out[1], "LPPLS calibration \\(mode = F2\\)")
  expect_true(any(grepl("Window:.*n = 200.*fh = 10", out)))
  expect_true(any(grepl("Estimate:.*tc = .*m = .*omega = ", out)))
  expect_true(any(grepl("Filters:", out)))

  ## Returns x invisibly, checked without letting the output leak into the run.
  capture.output(v <- withVisible(print(f)))
  expect_false(v$visible)
  expect_identical(v$value, f)
})


test_that("an interior optimum is reported as interior and raises no boundary flag", {
  f <- clean_fit()
  d <- .fit_diagnostics(f)
  expect_length(d$pinned, 0)
  expect_equal(d$boundary_counts[["m"]], 0)

  out <- capture.output(print(f))
  expect_true(any(grepl("Bounds:.*estimate interior", out)))
  expect_false(any(grepl("pinned at an optimisation bound", out)))
})


test_that("a boundary solution is detected at the estimate and across the fits", {
  f <- pinned_fit()
  d <- .fit_diagnostics(f)

  expect_true("m" %in% d$pinned)
  ## Every per-tc conditional fit is pinned too, not just the point estimate.
  expect_equal(d$boundary_counts[["m"]], d$n_fits)

  out <- capture.output(print(f))
  expect_true(any(grepl("Bounds:.*estimate m at a bound", out)))
  expect_true(any(grepl("! m pinned at an optimisation bound", out)))
})


test_that("boundary counts respect the tolerance and both bounds", {
  tbl   <- data.frame(m = c(0.1, 0.5, 0.9), omega = c(6, 9.5, 12.9))
  lower <- c(0.1, 6, -1e14, 0.8)
  upper <- c(0.9, 13, -1e-14, 1e6)

  cnt <- .boundary_counts(tbl, lower, upper)
  expect_equal(cnt[["m"]], 2)       # 0.1 and 0.9 sit on the two bounds
  expect_equal(cnt[["omega"]], 1)   # only the one sitting on 6

  ## Widening the tolerance to 2% of the omega span (0.14) pulls in 12.9, which
  ## is 0.1 below its bound, but still not the interior 9.5.
  expect_equal(.boundary_counts(tbl, lower, upper, tol = 0.02)[["omega"]], 2)
})


test_that("basin split is measured as a gap, so a smooth profile is not flagged", {
  f <- clean_fit()
  d <- .fit_diagnostics(f)

  expect_equal(d$best_sse, min(f$fit[[2]]$value))
  expect_gte(d$basin_gap, 0)
  expect_lte(d$basin_gap, 1)
  ## A conditional SSE profile rises away from the optimal tc, so plenty of fits
  ## are far above the best. That must not by itself count as a basin split.
  expect_gt(sum(f$fit[[2]]$value > 1.10 * d$best_sse), 0)
  expect_false(d$split)
  expect_false(any(grepl("distinct inferior basin", capture.output(print(f)))))
})


test_that("two separated clusters are split at the gap and flagged", {
  f <- clean_fit()
  n <- nrow(f$fit[[2]])
  ## A good cluster near 1 and a distinctly worse one near 10.
  f$fit[[2]]$value <- c(1, 1.01, 1.02, rep(c(10, 10.01), length.out = n - 3))

  d <- .fit_diagnostics(f)
  expect_true(d$split)
  expect_gt(d$basin_gap, 0.5)
  expect_equal(d$n_inferior, n - 3)
  expect_equal(d$split_at, c(1.02, 10))
  expect_equal(d$best_sse, 1)

  out <- capture.output(print(f))
  expect_true(any(grepl("largest gap .* splits off", out)))
  expect_true(any(grepl("! .* sit in a distinct inferior basin", out)))
})


test_that("the split threshold is the documented argument", {
  f <- clean_fit()
  n <- nrow(f$fit[[2]])
  f$fit[[2]]$value <- c(1, 1.01, 1.02, rep(c(10, 10.01), length.out = n - 3))

  ## The gap here is almost the whole spread, so pin the assertions to it: the
  ## split must switch off exactly when gap_tol rises above the measured gap.
  g <- .fit_diagnostics(f)$basin_gap
  expect_gt(g, 0.9)
  expect_true(.fit_diagnostics(f, gap_tol = g - 0.01)$split)
  expect_false(.fit_diagnostics(f, gap_tol = g + 0.001)$split)
})


test_that("basin statistics survive degenerate SSE columns", {
  f <- clean_fit()

  ## All identical: no spread, so no gap can be computed, but nothing errors.
  f$fit[[2]]$value <- rep(2, nrow(f$fit[[2]]))
  d <- expect_silent(.fit_diagnostics(f))
  expect_equal(d$best_sse, 2)
  expect_null(d$basin_gap)
  expect_silent(capture.output(print(f)))

  ## Near-zero SSE from a noiseless series must not blow the statistic up: it is
  ## a share of the spread, so it stays within [0, 1].
  f$fit[[2]]$value <- c(1e-16, seq(0.001, 0.01, length.out = nrow(f$fit[[2]]) - 1))
  d <- .fit_diagnostics(f)
  expect_gte(d$basin_gap, 0)
  expect_lte(d$basin_gap, 1)
})


test_that("MPL diagnostics report coverage, the estimate and the intervals", {
  f <- clean_fit(mode = "MPL", mpl_cutoff = c(0.05, 0.5))
  d <- .fit_diagnostics(f)

  expect_equal(d$mpl$total, f$fit_args$fh)
  expect_lte(d$mpl$defined, d$mpl$total)
  expect_equal(d$mpl$mple, f$mpl_output$tc_hat_mpl)
  expect_equal(d$mpl$cutoff, c(0.05, 0.5))

  out <- capture.output(print(f))
  expect_true(any(grepl("MPL:.*defined at .*/.* tc", out)))
  ## Intervals are labelled with the cutoffs actually used, not a fixed set.
  expect_true(any(grepl("LI\\(0.05\\)", out)))
  expect_true(any(grepl("LI\\(0.50\\)", out)))
  expect_false(any(grepl("LI\\(0.10\\)", out)))
})


test_that("a truncated MPL curve is flagged: edge maximum and interval on the MPLE", {
  f  <- clean_fit(mode = "MPL", mpl_cutoff = c(0.05))
  fh <- f$fit_args$fh

  ## A curve rising to its last defined point, NA beyond: exactly the shape a
  ## det(X'X - H) <= 0 truncation produces.
  f$mpl_output$LL <- c(seq_len(fh - 3), rep(NA_real_, 3))
  f$mpl_output$R  <- f$mpl_output$LL - max(f$mpl_output$LL, na.rm = TRUE)
  mple            <- f$fit_args$n + (fh - 3)
  f$mpl_output$tc_hat_mpl <- mple
  f$mpl_output$LI <- list(c(f$fit_args$n + 1, mple))

  d <- .fit_diagnostics(f)
  expect_equal(d$mpl$defined, fh - 3)
  expect_identical(d$mpl$at_edge, "upper")
  expect_true(d$mpl$li_touches_mple)

  out <- capture.output(print(f))
  expect_true(any(grepl("MPL undefined at 3 of", out)))
  expect_true(any(grepl("MPLE is the upper edge", out)))
  expect_true(any(grepl("bounded by the MPLE itself", out)))
})


test_that("an interior MPL maximum with a straddling interval raises no MPL flag", {
  f  <- clean_fit(mode = "MPL", mpl_cutoff = c(0.05))
  fh <- f$fit_args$fh
  n  <- f$fit_args$n

  peak <- 5L
  f$mpl_output$LL <- -abs(seq_len(fh) - peak)          # maximum in the interior
  f$mpl_output$R  <- f$mpl_output$LL - max(f$mpl_output$LL)
  f$mpl_output$tc_hat_mpl <- n + peak
  f$mpl_output$LI <- list(c(n + peak - 2, n + peak + 2))

  d <- .fit_diagnostics(f)
  expect_true(is.na(d$mpl$at_edge))
  expect_false(d$mpl$li_touches_mple)

  out <- capture.output(print(f))
  expect_false(any(grepl("edge of the defined region", out)))
  expect_false(any(grepl("bounded by the MPLE itself", out)))
  expect_false(any(grepl("MPL undefined", out)))
})


test_that("non-MPL fits print without an MPL section", {
  out <- capture.output(print(clean_fit()))
  expect_false(any(grepl("MPL:", out)))
  expect_false(any(grepl("^\\s+LI", out)))
})


test_that("F1 fits describe their rows as restarts, not conditional fits", {
  t     <- 1:200
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8)
  f <- fit_lppls(log_price = log_p, fh = 10, hold_out = 0,
                 tc_init = 205, m_init = 0.5, o_init = 8, num_searches = 3,
                 mode = "F1")

  expect_s3_class(f, "lppls_fit")
  expect_identical(.fit_diagnostics(f)$mode, "F1")
  out <- capture.output(print(f))
  expect_true(any(grepl("mode = F1", out)))
  expect_false(any(grepl("conditional fits", out)))
})


test_that("objects predating fit_args degrade instead of failing", {
  f <- clean_fit()
  f$fit_args <- NULL                    # what an older cached .rds looks like

  d <- expect_silent(.fit_diagnostics(f))
  expect_true(is.na(d$mode))
  expect_null(d$pinned)
  expect_null(d$boundary_counts)
  ## Everything not needing the recorded bounds still works.
  expect_false(is.null(d$best_sse))
  expect_false(is.null(d$basin_gap))

  out <- capture.output(print(f))
  expect_true(any(grepl("Bounds:.*not recorded", out)))
  expect_true(any(grepl("Basins:", out)))
})
