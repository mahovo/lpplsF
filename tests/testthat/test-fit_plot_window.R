## fit_plot()'s window markers: T1 bounds the fitted curve, T2 marks the end,
## and `dates` maps every x position onto a calendar axis.

make_fit <- function() {
  list(A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
       tc = 1001, m = 0.6, omega = 8)
}
make_series <- function(n = 900) {
  set.seed(3630)
  eval_lppls(seq_len(n), 8, -0.025, 0.0015, 0.001, 1001, 0.6, 8) +
    stats::rnorm(n) * 0.05
}

## the fitted curve is layer 2; pull its x range back off the built plot
fit_curve_x <- function(p) {
  d <- ggplot2::ggplot_build(p)$data[[2]]
  range(d$x, na.rm = TRUE)
}

test_that("without T1 the fitted curve starts at the first plotted point", {
  lp <- make_series()
  p  <- fit_plot(make_fit(), seq_along(lp), lp)
  expect_equal(fit_curve_x(p)[1], 1)
})

test_that("T1 starts the fitted curve at T1, not at the first observation", {
  lp <- make_series()
  p  <- fit_plot(make_fit(), seq_along(lp), lp, T1 = 400)
  expect_equal(fit_curve_x(p)[1], 400)
})

test_that("T1 does not truncate the observed series", {
  lp <- make_series()
  p  <- fit_plot(make_fit(), seq_along(lp), lp, T1 = 400)
  obs <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(min(obs$x, na.rm = TRUE), 1)
  expect_equal(max(obs$x, na.rm = TRUE), length(lp))
})

test_that("T1 and T2 both appear as markers, with tc_hat", {
  lp <- make_series()
  p  <- fit_plot(make_fit(), seq_along(lp), lp, T1 = 400, T2 = 800)
  b  <- ggplot2::ggplot_build(p)
  vlines <- b$data[[3]]
  expect_equal(nrow(vlines), 3L)
  expect_setequal(round(sort(vlines$xintercept)), c(400, 800, 1001))
})

test_that("T1 works on a calendar axis and shifts nothing numerically", {
  lp <- make_series()
  d  <- as.Date("2010-01-01") + seq_along(lp) - 1
  p_i <- fit_plot(make_fit(), seq_along(lp), lp, T1 = 400)
  p_d <- fit_plot(make_fit(), seq_along(lp), lp, T1 = 400, dates = d)
  expect_equal(p_i$labels$x, "Time Index")
  expect_equal(p_d$labels$x, "Date")
  ## same fitted values, just relabelled on the x axis
  expect_equal(ggplot2::ggplot_build(p_i)$data[[2]]$y,
               ggplot2::ggplot_build(p_d)$data[[2]]$y)
  ## the T1 marker lands on the right calendar day (d[400], the 400th observation)
  vl        <- ggplot2::ggplot_build(p_d)$data[[3]]
  marked    <- as.Date(vl$xintercept, origin = "1970-01-01")
  expect_true(d[400] %in% marked)
})

test_that("a rolling fit shifted to the global clock is drawn correctly", {
  ## the case the guide needs: fit stored in window-local time, replayed globally
  lp   <- make_series(2000)
  T1   <- 601
  f    <- make_fit()
  f$tc <- 1001 - (T1 - 1)              # window-local tc
  g    <- f; g$tc <- f$tc + (T1 - 1)   # back onto the global clock
  p    <- fit_plot(g, seq_along(lp), lp, T1 = T1, T2 = 900)
  expect_equal(fit_curve_x(p)[1], T1)
  expect_lt(fit_curve_x(p)[2], 1001)   # curve stops below tc
})
