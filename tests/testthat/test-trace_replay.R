# Tests for the reconstructed BFGS trace path.
#
# fit_lppls() does not record the optimiser's trajectory; it rebuilds it by
# re-running L-BFGS-B with maxit = 1, 2, 3, ... and stopping when the iterate
# stops changing. These tests pin the properties that stopping rule is meant to
# guarantee: every plotted point is distinct, the path ends at the optimum, and
# generating it never perturbs the fit.

## Small, fast fixture shared by the tests below.
trace_fixture <- function(num_searches = 1L, tp = c(1, 0, 0), seed = 1L) {
  set.seed(seed)
  n       <- 200L
  tc_true <- 220
  log_p   <- eval_lppls(seq_len(n), 8, -0.025, 0.0015, 0.001, tc_true, 0.6, 8) +
    rnorm(n) * 0.01
  list(
    log_p = log_p, t = seq_len(n),
    lower = c(0.1, 6, -0.03, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
    m_init = 0.3, o_init = 10,
    fit = suppressWarnings(fit_lppls(
      log_price = log_p, fh = 30, hold_out = 0,
      lower = c(0.1, 6, -0.03, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
      tc_init = tc_true, m_init = 0.3, o_init = 10,
      num_searches = num_searches, mode = "F2", tp = tp))
  )
}

## The plotted path, in order, as a data frame of coordinates.
trace_path <- function(p) {
  idx <- which(vapply(p$layers, function(l) class(l$geom)[1] == "GeomPath",
                      logical(1)))[1]
  d <- p$layers[[idx]]$data
  d[order(d$step), setdiff(names(d), "step"), drop = FALSE]
}


test_that("the trace contains no duplicated points", {
  skip_on_cran()

  f <- trace_fixture()
  path <- trace_path(f$fit$trace_plot_mo)

  ## More than just a start dot: the optimiser did move.
  expect_gt(nrow(path), 1L)

  ## The property the stopping rule exists to guarantee. The previous
  ## implementation sized the trace from optim()'s evaluation count, which
  ## over-ran the number of iterations by a factor of ~2-12 and left the tail
  ## of the path stacked on the optimum.
  expect_equal(anyDuplicated(path), 0L)

  ## Equivalently: consecutive points always differ.
  consecutive_same <- vapply(
    seq_len(nrow(path) - 1L),
    function(i) isTRUE(all.equal(unlist(path[i, ]), unlist(path[i + 1L, ]))),
    logical(1)
  )
  expect_false(any(consecutive_same))
})


test_that("the trace ends at the converged optimum", {
  skip_on_cran()

  f <- trace_fixture()
  path <- trace_path(f$fit$trace_plot_mo)

  bc   <- create_beta_calculator()
  sse2 <- function(par, tc, log_p, t) {
    b <- bc(log_p, t, tc, par[1], par[2])
    sum((log_p - eval_lppls(t, b["A"], b["B"], b["C1"], b["C2"],
                            tc, par[1], par[2], mode = 0))^2, na.rm = TRUE)
  }
  tc_hat <- f$fit$fit[[1]]$tc

  ## The replay is truncated at successive iteration limits, so its last point
  ## must equal an uncapped run from the same start.
  full <- stats::optim(
    c(f$m_init, f$o_init), sse2, tc = tc_hat, log_p = f$log_p, t = f$t,
    lower = f$lower[1:2], upper = f$upper[1:2], method = "L-BFGS-B")$par
  expect_equal(unname(unlist(path[nrow(path), c("m", "omega")])),
               unname(full), tolerance = 1e-8)

  ## Extending the limit past the end of the trace must not move it further,
  ## i.e. the path really has settled rather than been cut off.
  beyond <- stats::optim(
    c(f$m_init, f$o_init), sse2, tc = tc_hat, log_p = f$log_p, t = f$t,
    lower = f$lower[1:2], upper = f$upper[1:2], method = "L-BFGS-B",
    control = list(maxit = nrow(path) + 25L))$par
  expect_equal(unname(beyond), unname(full), tolerance = 1e-12)
})


test_that("the trace length matches the number of real optimiser steps", {
  skip_on_cran()

  f <- trace_fixture()
  path <- trace_path(f$fit$trace_plot_mo)

  bc   <- create_beta_calculator()
  sse2 <- function(par, tc, log_p, t) {
    b <- bc(log_p, t, tc, par[1], par[2])
    sum((log_p - eval_lppls(t, b["A"], b["B"], b["C1"], b["C2"],
                            tc, par[1], par[2], mode = 0))^2, na.rm = TRUE)
  }
  tc_hat <- f$fit$fit[[1]]$tc
  step_at <- function(k) stats::optim(
    c(f$m_init, f$o_init), sse2, tc = tc_hat, log_p = f$log_p, t = f$t,
    lower = f$lower[1:2], upper = f$upper[1:2], method = "L-BFGS-B",
    control = list(maxit = k))$par

  ## Count the distinct iterates independently, then compare. The plotted path
  ## prepends the starting point, so it holds one more.
  n_steps <- 0L
  prev <- NULL
  for (k in seq_len(100L)) {
    p_k <- step_at(k)
    if (!is.null(prev) && identical(p_k, prev)) break
    n_steps <- n_steps + 1L
    prev <- p_k
  }
  expect_equal(nrow(path), n_steps + 1L)
  expect_lt(n_steps, lpplsF:::.lppls_max_trace_steps)
})


test_that("generating the trace does not change the fit", {
  skip_on_cran()

  with_trace    <- trace_fixture(tp = c(1, 1, 1))$fit
  without_trace <- trace_fixture(tp = c(0, 0, 0))$fit

  expect_equal(with_trace$fit[[1]], without_trace$fit[[1]])
  expect_equal(with_trace$fit[[2]], without_trace$fit[[2]])
  expect_equal(with_trace$fit[[4]], without_trace$fit[[4]])

  ## and the trace plots are only attached when asked for
  expect_null(without_trace$trace_plot_mo)
  expect_null(without_trace$trace_plot_bm)
  expect_null(without_trace$trace_plot_bo)
})


test_that("all three trace selectors share one path", {
  skip_on_cran()

  fit <- trace_fixture(tp = c(1, 1, 1))$fit
  expect_s3_class(fit$trace_plot_mo, "ggplot")
  expect_s3_class(fit$trace_plot_bm, "ggplot")
  expect_s3_class(fit$trace_plot_bo, "ggplot")

  p_mo <- trace_path(fit$trace_plot_mo)
  p_bm <- trace_path(fit$trace_plot_bm)
  p_bo <- trace_path(fit$trace_plot_bo)

  ## One replay feeds all three, so they must agree in length ...
  expect_equal(nrow(p_bm), nrow(p_mo))
  expect_equal(nrow(p_bo), nrow(p_mo))

  ## ... and in the coordinates they have in common.
  expect_equal(p_bm$m, p_mo$m)
  expect_equal(p_bo$omega, p_mo$omega)

  ## The B column is present only for selectors 2 and 3.
  expect_true(all(c("m", "B") %in% names(p_bm)))
  expect_true(all(c("omega", "B") %in% names(p_bo)))
  expect_equal(p_bm$B, p_bo$B)

  ## Selectors 2 and 3 must be free of duplicates too.
  expect_equal(anyDuplicated(p_bm), 0L)
  expect_equal(anyDuplicated(p_bo), 0L)
})


test_that("the trace replay works when the best fit is a random restart", {
  skip_on_cran()

  fit <- trace_fixture(num_searches = 3L, tp = c(1, 0, 0), seed = 7L)$fit
  path <- trace_path(fit$trace_plot_mo)

  expect_gt(nrow(path), 1L)
  expect_equal(anyDuplicated(path), 0L)

  ## The start dot is the initial point of whichever search won: (m_init,
  ## o_init) for search 1, otherwise that search's seeded random draw.
  best_id <- fit$fit[[2]]$ID[1]
  set.seed(best_id)
  expected_start <- if (best_id == 1L) {
    c(0.3, 10)
  } else {
    c(stats::runif(1, 0.1, 0.9), stats::runif(1, 6, 13))
  }
  expect_equal(unname(unlist(path[1, c("m", "omega")])), expected_start)
})


test_that("the trace step cap is a sane positive bound", {
  expect_type(lpplsF:::.lppls_max_trace_steps, "integer")
  expect_gt(lpplsF:::.lppls_max_trace_steps, 1L)
})
