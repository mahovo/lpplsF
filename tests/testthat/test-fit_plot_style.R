## fit_plot()'s two styles. The "simulated" default must stay exactly as it was,
## because the thesis figures are generated through it.

est <- list(A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
            tc = 1001, m = 0.6, omega = 8)
series <- function(n = 900) {
  set.seed(3630)
  eval_lppls(seq_len(n), 8, -0.025, 0.0015, 0.001, 1001, 0.6, 8) +
    stats::rnorm(n) * 0.05
}
layers_of <- function(p) ggplot2::ggplot_build(p)$data
limits_of <- function(p) ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$get_limits()
line_colour <- function(p, i) unique(layers_of(p)[[i]]$colour)

test_that("the default style is simulated, and names the series 'sim data'", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800)
  expect_true("sim data" %in% limits_of(p))
  expect_false("prediction" %in% limits_of(p))
  expect_equal(line_colour(p, 1), "royalblue1")
  expect_equal(line_colour(p, 2), "red")
})

test_that("the simulated style is unchanged by the style argument existing", {
  lp <- series()
  a <- fit_plot(est, seq_along(lp), lp, T2 = 800)
  b <- fit_plot(est, seq_along(lp), lp, T2 = 800, style = "simulated")
  da <- layers_of(a); db <- layers_of(b)
  expect_equal(length(da), length(db))
  for (i in seq_along(da)) expect_equal(da[[i]], db[[i]])
})

test_that("the empirical style greys the data and splits the curve at T2", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800, style = "empirical")
  expect_equal(limits_of(p), c("data", "fit", "prediction"))
  expect_equal(line_colour(p, 1), "darkgray")
  expect_equal(line_colour(p, 2), "royalblue")
  expect_equal(line_colour(p, 3), "red")
})

test_that("the split meets at T2, with no gap and no overlap", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800, style = "empirical")
  d  <- layers_of(p)
  expect_equal(max(d[[2]]$x), 800)
  expect_equal(min(d[[3]]$x), 800)
  expect_equal(max(d[[3]]$x), floor(est$tc) - 1)
})

test_that("empirical keeps the window markers out of the legend", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800, T1 = 200, style = "empirical")
  expect_false(any(c("T1", "T2", "tc_hat") %in% limits_of(p)))
  expect_equal(length(layers_of(p)), 6L)   # 3 lines + 3 rules
})

test_that("simulated keeps the markers in the legend", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800, T1 = 200)
  expect_true(all(c("T1", "T2", "tc_hat") %in% limits_of(p)))
})

test_that("with no T2 there is nothing to split, so one line is drawn", {
  lp <- series()
  p  <- fit_plot(est, seq_along(lp), lp, style = "empirical")
  expect_false("prediction" %in% limits_of(p))
  expect_equal(line_colour(p, 2), "royalblue")
})

test_that("a tc inside the calibration window is not split either", {
  lp <- series()
  p <- fit_plot(est, seq_along(lp), lp, T2 = 1000, style = "empirical")
  expect_false("prediction" %in% limits_of(p))
})

test_that("data_label overrides the legend entry in either style", {
  lp <- series()
  p1 <- fit_plot(est, seq_along(lp), lp, data_label = "S&P 500")
  p2 <- fit_plot(est, seq_along(lp), lp, style = "empirical",
                 data_label = "S&P 500")
  expect_true("S&P 500" %in% limits_of(p1))
  expect_true("S&P 500" %in% limits_of(p2))
  expect_false("sim data" %in% limits_of(p1))
})

test_that("empirical works on a calendar axis", {
  lp <- series()
  d  <- as.Date("2010-01-01") + seq_along(lp) - 1
  p  <- fit_plot(est, seq_along(lp), lp, T2 = 800, style = "empirical", dates = d)
  expect_equal(p$labels$x, "Date")
  expect_equal(line_colour(p, 1), "darkgray")
  expect_equal(as.Date(max(layers_of(p)[[2]]$x), origin = "1970-01-01"), d[800])
})

test_that("an unknown style is rejected", {
  lp <- series(100)
  expect_error(fit_plot(est, seq_along(lp), lp, style = "pretty"), "should be one of")
})

test_that("fit_lppls passes style through to the attached plot", {
  set.seed(123)
  t  <- 1:200
  lp <- eval_lppls(t, 4, -0.015, 0.002, 0.001, 250, 0.5, 9) + stats::rnorm(200, 0, 0.01)
  f  <- suppressWarnings(fit_lppls(log_price = lp, fh = 20, hold_out = 20,
                                   tc_init = 210, num_searches = 1, mode = "F2",
                                   fp = TRUE, style = "empirical"))
  expect_true("data" %in% limits_of(f$fit_plot))
  expect_false("sim data" %in% limits_of(f$fit_plot))
})
