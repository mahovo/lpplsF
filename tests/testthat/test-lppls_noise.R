## lppls_noise_sweep() and its two plots.

make_signal <- function(N = 400) {
  eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
             tc = N + 1, m = 0.6, omega = 8)
}
quick_args <- list(fh = 20, hold_out = 20, num_searches = 1, mode = "F2")
mpl_args   <- list(fh = 20, hold_out = 20, num_searches = 1,
                   mode = "MPL", mpl_cutoff = c(0.05, 0.5))

test_that("a sweep returns one row per noise level, with the estimate columns", {
  sw <- lppls_noise_sweep(make_signal(), c(0.01, 0.05),
                          fit_args = quick_args, seed = 3630)
  expect_s3_class(sw, "lppls_noise_sweep")
  expect_s3_class(sw, "data.frame")
  expect_equal(nrow(sw), 2L)
  expect_equal(sw$sd, c(0.01, 0.05))
  expect_true(all(c("tc", "m", "omega", "A", "B", "C1", "C2", "D", "value")
                  %in% names(sw)))
})

test_that("likelihood-interval columns are named for their cutoff percentage", {
  sw <- suppressWarnings(
    lppls_noise_sweep(make_signal(), c(0.01, 0.02),
                      fit_args = mpl_args, seed = 3630))
  ## cutoffs 0.05 and 0.50 -> LI5l/LI5u and LI50l/LI50u, matching lppls_rolling()
  expect_true(all(c("LI5l", "LI5u", "LI50l", "LI50u") %in% names(sw)))
  expect_false(any(is.na(sw$tc_hat_mpl)))
})

test_that("a seed makes the sweep reproducible", {
  a <- lppls_noise_sweep(make_signal(), c(0.02, 0.04),
                         fit_args = quick_args, seed = 11)
  b <- lppls_noise_sweep(make_signal(), c(0.02, 0.04),
                         fit_args = quick_args, seed = 11)
  expect_equal(a$tc, b$tc)
  expect_equal(a$omega, b$omega)
})

test_that("the noise generator is pluggable", {
  heavy <- function(n) stats::rt(n, df = 3) / sqrt(3)   # unit-scale heavy tails
  sw <- lppls_noise_sweep(make_signal(), c(0.02),
                          noise = heavy, fit_args = quick_args, seed = 7)
  expect_equal(nrow(sw), 1L)
  expect_true(is.finite(sw$tc))
  ## a different generator on the same seed must give a different series
  g <- lppls_noise_sweep(make_signal(), c(0.02),
                         fit_args = quick_args, seed = 7)
  expect_false(isTRUE(all.equal(sw$value, g$value)))
})

test_that("bad input is rejected with a specific message", {
  sig <- make_signal(100)
  expect_error(lppls_noise_sweep("a", 0.01), "'signal' must be a numeric")
  expect_error(lppls_noise_sweep(sig, numeric(0)), "non-empty numeric")
  expect_error(lppls_noise_sweep(sig, c(-1, 1)), "must not contain negative")
  expect_error(lppls_noise_sweep(sig, 0.01, noise = "rnorm"), "must be a function")
  expect_error(lppls_noise_sweep(sig, 0.01, fit_args = "x"), "must be a list")
})

test_that("noise_param_plot draws one panel per requested parameter", {
  sw <- lppls_noise_sweep(make_signal(), c(0.01, 0.03, 0.05),
                          fit_args = quick_args, seed = 3630)
  p <- noise_param_plot(sw)
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(ggplot2::ggplot_build(p)$plot$data$param), 6L)

  p2 <- noise_param_plot(sw, params = c("m", "omega"))
  expect_equal(nlevels(ggplot2::ggplot_build(p2)$plot$data$param), 2L)
})

test_that("noise_tc_plot needs likelihood intervals and says so", {
  sw <- lppls_noise_sweep(make_signal(), c(0.01, 0.03),
                          fit_args = quick_args, seed = 3630)   # F2: no LI
  expect_error(noise_tc_plot(sw), "No likelihood-interval columns")
})

test_that("noise_tc_plot draws a band per likelihood interval", {
  sw <- suppressWarnings(
    lppls_noise_sweep(make_signal(), c(0.01, 0.02, 0.03),
                      fit_args = mpl_args, seed = 3630))
  p <- noise_tc_plot(sw, true_tc = 401)
  expect_s3_class(p, "ggplot")
  ## two cutoffs -> two bands in the ribbon layer's data
  ribbon <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(ribbon$group)), 2L)
})

test_that("the plots accept a plain data.frame saved before the class existed", {
  sw <- suppressWarnings(
    lppls_noise_sweep(make_signal(), c(0.01, 0.02),
                      fit_args = mpl_args, seed = 3630))
  bare <- as.data.frame(unclass(sw))            # drop the class, keep the columns
  expect_false(inherits(bare, "lppls_noise_sweep"))
  expect_s3_class(noise_param_plot(bare), "ggplot")
  expect_s3_class(noise_tc_plot(bare), "ggplot")
})

test_that("a table missing required columns is rejected by name", {
  expect_error(noise_param_plot(data.frame(x = 1)), "missing the column")
  expect_error(noise_tc_plot(data.frame(sd = 1)), "missing the column")
})
