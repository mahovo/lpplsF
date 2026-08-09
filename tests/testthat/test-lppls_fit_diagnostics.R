## Tier 2 (likelihood-interval integrity), Tier 3 (curvature at the estimate)
## and Tier 4 (summary / plot) diagnostics for an "lppls_fit".

mpl_fit <- function(...) {
  t     <- 1:200
  set.seed(4)
  log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                      tc = 221, m = 0.6, omega = 8) + stats::rnorm(200) * 0.01
  fit_lppls(log_price = log_p, fh = 20, hold_out = 0,
            lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
            tc_init = 210, m_init = 0.5, o_init = 8, num_searches = 2,
            mode = "MPL", mpl_cutoff = c(0.05, 0.5), ...)
}


test_that("the fit carries its input series", {
  f <- mpl_fit()
  expect_length(f$log_price, 200)
  expect_identical(f$log_price[seq_len(f$fit_args$n)], f$log_price)  # hold_out = 0
})


## ---- Tier 2: are the intervals closed by the likelihood or by the grid? -----

test_that("a contiguous interval is reported as contiguous and raises no flag", {
  f  <- mpl_fit()
  n  <- f$fit_args$n
  ## A single-peaked curve: every point above a cutoff forms one run.
  peak <- 10L
  f$mpl_output$LL <- -abs(seq_len(20) - peak)
  f$mpl_output$R  <- f$mpl_output$LL - max(f$mpl_output$LL)
  f$mpl_output$tc_hat_mpl <- n + peak
  f$mpl_output$LI <- list(c(n + peak - 2, n + peak + 2),   # |R| < 3.00
                          c(n + peak - 0, n + peak + 0))

  k <- .fit_diagnostics(f)$mpl$li_check
  expect_true(all(vapply(k, function(z) isTRUE(z$contiguous), logical(1))))
  expect_true(all(vapply(k, function(z) z$bridged == 0L, logical(1))))
  expect_false(any(grepl("range\\(\\) bridges", capture.output(print(f)))))
})


test_that("an interval bridging a hole is detected and counted", {
  f <- mpl_fit()
  n <- f$fit_args$n
  ## Two separated runs above log(0.05) = -3.00, with three points below the
  ## cutoff in between: range() spans them, so the interval is not what it seems.
  R <- c(0, -1, -1, -8, -8, -8, -1, -1, rep(-20, 12))
  f$mpl_output$R  <- R
  f$mpl_output$LL <- R
  f$mpl_output$tc_hat_mpl <- n + 1L
  f$mpl_output$LI <- list(c(n + 1, n + 8), c(n + 1, n + 3))

  k <- .fit_diagnostics(f)$mpl$li_check
  expect_equal(k[[1]]$n_in, 5)        # positions 1,2,3,7,8
  expect_equal(k[[1]]$span, 8)
  expect_false(k[[1]]$contiguous)
  expect_equal(k[[1]]$bridged, 3)

  out <- capture.output(print(f))
  expect_true(any(grepl("LI\\(0.05\\) spans 8 grid points but only 5", out)))
  expect_true(any(grepl("bridges 3", out)))
})


test_that("li_check is absent when the cutoffs were not recorded", {
  f <- mpl_fit()
  f$fit_args$mpl_cutoff <- NULL
  expect_null(.fit_diagnostics(f)$mpl$li_check)
  expect_silent(capture.output(print(f)))
})


## ---- Tier 3: curvature at the estimate --------------------------------------

test_that("lppls_curvature factorises det(X'X - H) through the Schur complement", {
  f <- mpl_fit()
  cv <- lppls_curvature(f)

  expect_s3_class(cv, "lppls_curvature")
  ## det(full) = det(linear block) * det(profile Hessian) is the identity the
  ## whole diagnostic rests on.
  expect_equal(cv$det_full, cv$det_linear * cv$det_profile, tolerance = 1e-6)
  expect_gt(cv$det_linear, 0)          # linear block is a crossproduct
  expect_equal(dim(cv$hessian_profile), c(2L, 2L))
  expect_length(cv$eigen_profile, 2L)
  expect_equal(prod(cv$eigen_profile), cv$det_profile, tolerance = 1e-6)
  expect_true(cv$classification %in%
                c("positive definite", "saddle", "negative definite"))
})


test_that("the profile Hessian matches a numerical Hessian of the profiled SSE", {
  f  <- mpl_fit()
  cv <- lppls_curvature(f)
  e  <- f$fit[[1]]
  n  <- f$fit_args$n
  t  <- seq_len(n)
  lp <- f$log_price[seq_len(n)]
  bc <- create_beta_calculator()

  ## g(m, omega) = 1/2 * SSE, profiled over the linear parameters.
  g <- function(m, om) {
    b <- bc(lp, t, e$tc, m, om)
    fitted <- eval_lppls(t, b["A"], b["B"], b["C1"], b["C2"], e$tc, m, om, mode = 0)
    0.5 * sum((lp - fitted)^2, na.rm = TRUE)
  }
  hm <- 1e-4; ho <- 1e-3
  gmm <- (g(e$m + hm, e$omega) - 2 * g(e$m, e$omega) + g(e$m - hm, e$omega)) / hm^2
  goo <- (g(e$m, e$omega + ho) - 2 * g(e$m, e$omega) + g(e$m, e$omega - ho)) / ho^2

  expect_equal(cv$hessian_profile[1, 1], gmm, tolerance = 1e-3)
  expect_equal(cv$hessian_profile[2, 2], goo, tolerance = 1e-3)
  expect_equal(cv$hessian_profile[1, 2], cv$hessian_profile[2, 1], tolerance = 1e-8)
})


test_that("lppls_curvature accepts a supplied series and validates its input", {
  f  <- mpl_fit()
  lp <- f$log_price
  ref <- lppls_curvature(f)

  f$log_price <- NULL                                  # an older fit
  expect_error(lppls_curvature(f), "does not carry its input series")
  expect_equal(lppls_curvature(f, log_price = lp)$det_full, ref$det_full)
  expect_error(lppls_curvature(f, log_price = lp[1:10]), "shorter than")
  expect_error(lppls_curvature(list(a = 1)), "lppls_fit")
})


test_that("curvature prints its classification", {
  out <- capture.output(print(lppls_curvature(mpl_fit())))
  expect_match(out[1], "LPPLS curvature at the estimate")
  expect_true(any(grepl("det\\(X'X - H\\)", out)))
  expect_true(any(grepl("profile Hessian is", out)))
})


## ---- Tier 4: summary and plot ------------------------------------------------

test_that("summary carries the diagnostics as values plus the curvature", {
  f <- mpl_fit()
  s <- summary(f)

  expect_s3_class(s, "summary.lppls_fit")
  expect_s3_class(s$curvature, "lppls_curvature")
  expect_equal(s$best_sse, min(f$fit[[2]]$value))
  expect_equal(s$mpl$mple, f$mpl_output$tc_hat_mpl)

  out <- capture.output(print(s))
  expect_true(any(grepl("Curvature:", out)))
  expect_true(any(grepl("Likelihood intervals:", out)))
  expect_true(any(grepl("LI\\(0.05\\):", out)))
})


test_that("summary degrades when the series is unavailable", {
  f <- mpl_fit()
  f$log_price <- NULL
  s <- summary(f)
  expect_null(s$curvature)
  out <- capture.output(print(s))
  expect_false(any(grepl("Curvature:", out)))   # the rest still prints
  expect_true(any(grepl("Estimate:", out)))
})


test_that("plot returns the stored components and errors on absent ones", {
  f <- mpl_fit(pp = TRUE)
  expect_identical(plot(f, "param"), f$param_plot)
  expect_error(plot(f, "fit"), "no fit plot")
  expect_error(plot(f, "nonsense"))
})


test_that("the basin plot labels each fit and is always available", {
  f <- mpl_fit()                       # no plots requested at fit time
  p <- plot(f, "basin")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(f$fit[[2]]))
  expect_true("basin" %in% names(p$data))
  ## With no split every fit belongs to the single cluster.
  expect_setequal(unique(p$data$basin), "single")

  ## Force two clusters and the labelling follows the reported split.
  f$fit[[2]]$value <- c(rep(1, 5), rep(10, nrow(f$fit[[2]]) - 5))
  d <- .fit_diagnostics(f)
  p2 <- plot(f, "basin")
  expect_true(d$split)
  expect_setequal(unique(p2$data$basin), c("good", "inferior"))
  expect_equal(sum(p2$data$basin == "inferior"), d$n_inferior)
})


test_that("plot defaults to the basin view", {
  f <- mpl_fit()
  expect_equal(plot(f)$labels$title, plot(f, "basin")$labels$title)
})


test_that("param_basin facets every parameter and shares the basin split", {
  f <- mpl_fit()
  f$fit[[2]]$value <- c(rep(1, 5), rep(10, nrow(f$fit[[2]]) - 5))
  d <- .fit_diagnostics(f)

  p <- plot(f, "param_basin")
  expect_s3_class(p, "ggplot")
  ## One facet per estimated quantity: everything in fit[[2]] but ID and tc.
  expect_setequal(unique(p$data$param),
                  setdiff(names(f$fit[[2]]), c("ID", "tc", "basin")))
  expect_equal(nrow(p$data),
               nrow(f$fit[[2]]) * length(unique(p$data$param)))

  ## The split is the one print() reports, not a separate calculation.
  expect_setequal(unique(p$data$basin), c("good", "inferior"))
  n_par <- length(unique(p$data$param))
  expect_equal(sum(p$data$basin == "inferior"), d$n_inferior * n_par)

  ## It carries the same subtitle as the basin view, and unlike "param" it does
  ## not need a plot to have been requested at fit time.
  expect_equal(p$labels$subtitle, plot(f, "basin")$labels$subtitle)
  expect_null(f$param_plot)
})
