# Tests for lpplsF main fitting function

# Helper function to generate synthetic bubble data
generate_bubble_data <- function(n = 200, tc = 250, noise_sd = 0.01, seed = 123) {
  set.seed(seed)
  t <- 1:n
  true_params <- list(
    A = 4, B = -0.015, C1 = 0.002, C2 = 0.001,
    tc = tc, m = 0.5, omega = 9
  )
  log_p <- eval_lppls(t, true_params$A, true_params$B, true_params$C1,
                 true_params$C2, true_params$tc, true_params$m,
                 true_params$omega) + rnorm(n, 0, noise_sd)
  list(t = t, log_p = log_p, true_params = true_params)
}

test_that("fit_lppls validates mode parameter", {
  data <- generate_bubble_data(n = 50)
  expect_error(
    fit_lppls(log_price = data$log_p, mode = "invalid"),
    "mode must be one of: 'F1', 'F2', 'MPL'"
  )
})

test_that("fit_lppls validates bound vectors", {
  data <- generate_bubble_data(n = 50)
  expect_error(
    fit_lppls(log_price = data$log_p, lower = c(0.1, 6)),
    "lower and upper must be vectors of length 4"
  )
})

test_that("fit_lppls returns expected structure in F1 mode", {
  skip_on_cran()  # Skip on CRAN due to computation time

  data <- generate_bubble_data(n = 100, tc = 150)
  result <- fit_lppls(
    log_price = data$log_p,
    fh = 30,
    hold_out = 10,
    tc_init = 120,
    mode = "F1",
    num_searches = 5,
    fb = FALSE
  )

  # Check structure
  expect_type(result, "list")
  expect_true("fit" %in% names(result))
  expect_type(result$fit, "list")

  # Check fit[[1]] has all required parameters
  expected_params <- c("tc", "m", "omega", "A", "B", "C1", "C2", "D", "value")
  expect_true(all(expected_params %in% names(result$fit[[1]])))

  # Check fit[[2]] is a tibble with results
  expect_s3_class(result$fit[[2]], "tbl_df")
  expect_equal(nrow(result$fit[[2]]), 5)  # num_searches = 5
})

test_that("fit_lppls returns expected structure in F2 mode", {
  skip_on_cran()

  data <- generate_bubble_data(n = 100, tc = 150)
  result <- fit_lppls(
    log_price = data$log_p,
    fh = 20,
    hold_out = 10,
    tc_init = 110,
    mode = "F2",
    num_searches = 3,
    fb = FALSE
  )

  # Check structure
  expect_type(result, "list")
  expect_true("fit" %in% names(result))

  # Check fit[[1]] has all required parameters
  expected_params <- c("tc", "m", "omega", "A", "B", "C1", "C2", "D", "value")
  expect_true(all(expected_params %in% names(result$fit[[1]])))

  # Check fit[[4]] contains all tc iterations
  expect_type(result$fit[[4]], "list")
  expect_equal(length(result$fit[[4]]), 20)  # fh = 20
})

test_that("fit_lppls parameters are within specified bounds", {
  skip_on_cran()

  data <- generate_bubble_data(n = 100, tc = 150)
  lower <- c(0.1, 6, -1e14, 0.8)
  upper <- c(0.9, 13, -1e-14, 1e6)

  result <- fit_lppls(
    log_price = data$log_p,
    fh = 20,
    hold_out = 10,
    lower = lower,
    upper = upper,
    tc_init = 110,
    mode = "F2",
    num_searches = 3,
    fb = FALSE
  )

  fit <- result$fit[[1]]

  # Check m is in bounds
  expect_gte(fit$m, lower[1])
  expect_lte(fit$m, upper[1])

  # Check omega is in bounds
  expect_gte(fit$omega, lower[2])
  expect_lte(fit$omega, upper[2])

  # Check tc is in forecast horizon
  n_model <- 100 - 10  # n - hold_out
  expect_gte(fit$tc, n_model + 1)
  expect_lte(fit$tc, n_model + 20)  # fh = 20
})

test_that("fit_lppls recovers true parameters from clean data", {
  skip_on_cran()

  # Generate data with very low noise
  data <- generate_bubble_data(n = 150, tc = 200, noise_sd = 0.001)

  result <- fit_lppls(
    log_price = data$log_p,
    fh = 80,
    hold_out = 0,
    tc_init = 200,
    m_init = 0.5,
    o_init = 9,
    mode = "F2",
    num_searches = 10,
    fb = FALSE
  )

  fit <- result$fit[[1]]
  true <- data$true_params

  # Check parameters are close to true values
  # (allow some tolerance due to optimization)
  expect_equal(fit$tc, true$tc, tolerance = 20)
  expect_equal(fit$m, true$m, tolerance = 0.2)
  expect_equal(fit$omega, true$omega, tolerance = 2)
})

test_that("fit_lppls generates fit plot when requested", {
  skip_on_cran()

  data <- generate_bubble_data(n = 100, tc = 150)
  result <- fit_lppls(
    log_price = data$log_p,
    fh = 20,
    hold_out = 10,
    tc_init = 110,
    mode = "F2",
    num_searches = 3,
    fp = TRUE,
    fb = FALSE
  )

  expect_true("fit_plot" %in% names(result))
  expect_s3_class(result$fit_plot, "ggplot")
})

test_that("fit_lppls out_of_range_tracker is populated correctly", {
  skip_on_cran()

  data <- generate_bubble_data(n = 100, tc = 150)
  result <- fit_lppls(
    log_price = data$log_p,
    fh = 20,
    hold_out = 10,
    tc_init = 110,
    mode = "F2",
    num_searches = 3,
    fb = FALSE
  )

  # Check structure of tracker
  expect_type(result$out_of_range_tracker, "list")
  expect_true("B" %in% names(result$out_of_range_tracker))
  expect_true("D" %in% names(result$out_of_range_tracker))
})

test_that("fit_lppls handles edge case of short time series", {
  skip_on_cran()

  set.seed(42)
  log_p <- log(cumsum(1 + rnorm(30, 0.01, 0.01)))

  # Should run without error, though results may not be meaningful
  result <- fit_lppls(
    log_price = log_p,
    fh = 10,
    hold_out = 5,
    tc_init = 30,
    mode = "F1",
    num_searches = 2,
    fb = FALSE
  )

  expect_type(result, "list")
  expect_true("fit" %in% names(result))
})

# ---------------------------------------------------------------------------
# num_searches = 1 trace fidelity (P-6): the replay must start from
# (m_init, o_init), so the first plotted segment is the true first BFGS step.
# ---------------------------------------------------------------------------

test_that("num_searches = 1 trace path replays the real BFGS step from (m_init, o_init)", {
  skip_on_cran()

  set.seed(1)
  n       <- 200L
  tc_true <- 220
  log_p   <- eval_lppls(seq_len(n), 8, -0.025, 0.0015, 0.001, tc_true, 0.6, 8) +
    rnorm(n) * 0.01
  t       <- seq_len(n)
  lower   <- c(0.1, 6, -0.03, 0.8)
  upper   <- c(0.9, 13, -1e-14, 1e6)
  m_init  <- 0.3
  o_init  <- 10

  fit <- suppressWarnings(fit_lppls(
    log_price = log_p, fh = 30, hold_out = 0, lower = lower, upper = upper,
    tc_init = tc_true, m_init = m_init, o_init = o_init,
    num_searches = 1, mode = "F2", tp = c(1, 0, 0)))

  ## Pull the (m, omega) path out of the trace plot's GeomPath layer.
  ## Match the geom class exactly -- GeomContour inherits from GeomPath.
  p    <- fit$trace_plot_mo
  idx  <- which(vapply(p$layers, function(l) class(l$geom)[1] == "GeomPath",
                       logical(1)))[1]
  path <- p$layers[[idx]]$data
  path <- path[order(path$step), ]

  ## The red start dot is the true initial values ...
  expect_equal(unname(unlist(path[1, c("m", "omega")])), c(m_init, o_init))

  ## ... and the first segment is one L-BFGS-B step FROM those values, not from a
  ## random restart. Reconstruct the F2 objective the fitter minimises and replay
  ## a single iteration; this must equal the plotted second point.
  bc   <- create_beta_calculator()
  sse2 <- function(par, tc, log_p, t) {
    b <- bc(log_p, t, tc, par[1], par[2])
    sum((log_p - eval_lppls(t, b["A"], b["B"], b["C1"], b["C2"],
                            tc, par[1], par[2], mode = 0))^2, na.rm = TRUE)
  }
  tc_hat <- fit$fit[[1]]$tc
  step1  <- stats::optim(c(m_init, o_init), sse2, tc = tc_hat, log_p = log_p, t = t,
                         lower = lower[1:2], upper = upper[1:2],
                         method = "L-BFGS-B", control = list(maxit = 1))$par
  expect_equal(unname(unlist(path[2, c("m", "omega")])), unname(step1),
               tolerance = 1e-6)

  ## Sanity: the path converges to the reported estimate.
  expect_equal(unname(unlist(path[nrow(path), c("m", "omega")])),
               c(fit$fit[[1]]$m, fit$fit[[1]]$omega), tolerance = 1e-2)
})
