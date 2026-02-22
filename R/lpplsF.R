#' Fit LPPLS Model to Financial Time Series
#'
#' Fits the Log-Periodic Power Law Singularity (LPPLS) model to log-price
#' data using the calibration methods from Filimonov & Sornette (2013).
#' Supports multiple optimization modes and provides diagnostic outputs.
#'
#' @param time_ID Numeric vector of time indices (T1:T2 where T1 = 1).
#' @param log_price Numeric vector of log-prices (must be same length as time_ID).
#' @param fh Integer, forecast horizon length in time units (default 30).
#'   The search for tc extends from T2+1 to T2+fh.
#' @param hold_out Integer, number of observations to hold out from the end
#'   of the sample for validation (default 15).
#' @param lower Numeric vector of length 4, lower bounds for
#'   c(m, omega, B, damping). Default: c(0.1, 6, -1e14, 0.8).
#' @param upper Numeric vector of length 4, upper bounds for
#'   c(m, omega, B, damping). Default: c(0.9, 13, -1e-14, 1e6).
#' @param tc_init Numeric, initial value for critical time tc (default 1000).
#' @param m_init Numeric, initial value for power law exponent m (default 0.5).
#' @param o_init Numeric, initial value for frequency omega (default 13).
#' @param num_searches Integer, number of optimization runs with different
#'   starting values (default 20).
#' @param mode Character, optimization mode:
#'   \itemize{
#'     \item "F1": Optimize tc, m, omega simultaneously
#'     \item "F2": Optimize m, omega first (for each tc), then optimize tc
#'     \item "MPL": Modified Profile Likelihood with confidence intervals
#'   }
#' @param mpl_cutoff Numeric vector of length 3, cutoff levels for MPL
#'   confidence intervals (default c(0.05, 0.1, 0.5)).
#' @param mpl_plot Logical, whether to generate MPL plot (default FALSE).
#' @param fp Logical, whether to generate fit plot (default FALSE).
#' @param cp Logical, whether to generate contour plot (default FALSE).
#' @param sp Logical, whether to generate surface plot (default FALSE).
#' @param tp Numeric vector of length 3, trace plot selectors for
#'   c(m-omega, B-m, B-omega) (default c(0,0,0)).
#' @param tp_id Integer, which trace step to use for contour (default 1).
#' @param pp Logical, whether to generate parameter plot (default FALSE).
#' @param mp Logical, whether to generate matrix plot (default FALSE).
#' @param factr Numeric, convergence tolerance for L-BFGS-B (default 1e-08).
#' @param fb Logical, whether to print feedback during execution (default FALSE).
#'
#' @return A list containing:
#'   \describe{
#'     \item{fit}{List with fitted model results. In F1/F2 mode:
#'       fit[[1]] = best fit parameters;
#'       fit[[2]] = tibble of all fits sorted by SSE}
#'     \item{mpl_output}{(MPL mode only) List with likelihood intervals,
#'       R values, log-likelihoods}
#'     \item{mpl_plot}{(if requested) ggplot2 object for MPL plot}
#'     \item{fit_plot}{(if requested) ggplot2 object for fit plot}
#'     \item{contour_data}{(if requested) List with x, y, z for contour}
#'     \item{contour_plot}{(if requested) plotly object for contour plot}
#'     \item{surface_plot}{(if requested) plotly object for surface plot}
#'     \item{trace_plot_mo}{(if requested) ggplot2 object for m-omega trace}
#'     \item{trace_plot_bm}{(if requested) ggplot2 object for B-m trace}
#'     \item{trace_plot_bo}{(if requested) ggplot2 object for B-omega trace}
#'     \item{param_plot}{(if requested) ggplot2 object for parameter plot}
#'     \item{matrix_plot}{(if requested) ggplot2 object for matrix plot}
#'     \item{out_of_range_tracker}{List tracking which fits had out-of-range
#'       B or damping values}
#'   }
#'
#' @details
#' The LPPLS model is:
#' \deqn{y(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c - t)) + C_2 \sin(\omega \log(t_c - t))]}
#'
#' The model has 7 parameters: A, B, C1, C2 (linear) and tc, m, omega (nonlinear).
#' The Filimonov calibration exploits this structure by solving for the linear
#' parameters analytically given the nonlinear ones.
#'
#' **Mode F1**: Optimizes all nonlinear parameters (tc, m, omega) simultaneously
#' using L-BFGS-B with multiple random starting points.
#'
#' **Mode F2**: A two-stage procedure:
#' \enumerate{
#'   \item For each candidate tc value, optimize (m, omega)
#'   \item Select the tc that minimizes the overall SSE
#' }
#' This is more robust than F1 for ill-conditioned problems.
#'
#' **Mode MPL**: Modified Profile Likelihood inference from Filimonov et al. (2017).
#' Provides likelihood-based confidence intervals for tc.
#'
#' @references
#' Filimonov, V., & Sornette, D. (2013). A stable and robust calibration scheme
#' of the log-periodic power law model. Physica A, 392(17), 3698-3707.
#'
#' Filimonov, V., Demos, G., & Sornette, D. (2017). Modified profile likelihood
#' inference and interval forecast of the burst of financial bubbles.
#' Quantitative Finance, 17(8), 1167-1186.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic bubble data
#' set.seed(123)
#' t <- 1:200
#' true_params <- list(A = 4, B = -0.015, C1 = 0.002, C2 = 0.001,
#'                     tc = 250, m = 0.5, omega = 9)
#' log_p <- eval_lppls(t, true_params$A, true_params$B, true_params$C1,
#'                    true_params$C2, true_params$tc, true_params$m,
#'                    true_params$omega) + rnorm(200, 0, 0.01)
#'
#' # Fit model
#' result <- fit_lppls(time_ID = t, log_price = log_p,
#'                  fh = 100, hold_out = 0,
#'                  tc_init = 250, mode = "F2",
#'                  num_searches = 10, fp = TRUE)
#'
#' # View best fit parameters
#' print(result$fit[[1]])
#'
#' # View fit plot
#' print(result$fit_plot)
#' }
#'
#' @import ggplot2
#' @import plotly
#' @importFrom dplyr arrange filter
#' @importFrom tibble tibble
#' @importFrom tidyr gather
#' @importFrom stats optim runif
#' @export
fit_lppls <- function(
  time_ID,
  log_price,
  fh = 30,
  hold_out = 15,
  lower = c(0.1, 6, -1e14, 0.8),
  upper = c(0.9, 13, -1e-14, 1e6),
  tc_init = 1000,
  m_init = 0.5,
  o_init = 13,
  num_searches = 20,
  mode = "F2",
  mpl_plot = FALSE,
  mpl_cutoff = c(0.05, 0.1, 0.5),
  fp = FALSE,
  cp = FALSE,
  sp = FALSE,
  tp = c(0, 0, 0),
  tp_id = 1,
  pp = FALSE,
  mp = FALSE,
  factr = 1e-08,
  fb = FALSE
) {

  # Input validation
  if (length(time_ID) != length(log_price)) {
    stop("time_ID and log_price vectors must have the same length.")
  }

  if (!mode %in% c("F1", "F2", "MPL")) {
    stop("mode must be one of: 'F1', 'F2', 'MPL'")
  }

  if (length(lower) != 4 || length(upper) != 4) {
    stop("lower and upper must be vectors of length 4")
  }

  start_time <- Sys.time()

  # Separate modeling and holdout data
  n_total <- length(time_ID)
  n_model <- n_total - hold_out

  t <- time_ID[1:n_model]
  log_p <- log_price[1:n_model]
  n <- length(t)

  # Create beta calculator (uses symbolic math for efficiency)
  beta_calculator <- create_beta_calculator()

  # Internal SSE functions
  # SSE with all nonlinear parameters as vector
  SSE1 <- function(par, log_p, t) {
    tc <- par[1]
    m <- par[2]
    omega <- par[3]

    # Ensure tc > max(t)
    if (tc <= max(t)) return(Inf)

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(t, beta["A"], beta["B"], beta["C1"], beta["C2"],
                         tc, m, omega, mode = 0)
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  # SSE with fixed tc
  SSE2 <- function(par, tc, log_p, t) {
    m <- par[1]
    omega <- par[2]

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(t, beta["A"], beta["B"], beta["C1"], beta["C2"],
                         tc, m, omega, mode = 0)
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  # SSE with fixed m and omega
  SSE3 <- function(tc, m, omega, log_p, t) {
    if (tc <= max(t)) return(Inf)

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(t, beta["A"], beta["B"], beta["C1"], beta["C2"],
                         tc, m, omega, mode = 0)
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  # Initialize outputs
  fit <- list()
  mpl_output <- NULL
  mpl_plot_obj <- NULL
  fit_plot_obj <- NULL
  contour_data <- NULL
  contour_plot_obj <- NULL
  surface_plot_obj <- NULL
  trace_plot_mo <- NULL
  trace_plot_bm <- NULL
  trace_plot_bo <- NULL
  param_plot_obj <- NULL
  matrix_plot_obj <- NULL
  out_of_range_tracker <- list(B = list(), D = list())

  # =========================================================================
  # MODE F1: Simultaneous optimization of tc, m, omega
  # =========================================================================
  if (mode == "F1") {
    if (fb) message("Mode: F1 - Simultaneous optimization")

    # Initialize results tibble
    fit_tbl <- tibble::tibble(
      tc = numeric(num_searches),
      m = numeric(num_searches),
      omega = numeric(num_searches),
      A = numeric(num_searches),
      B = numeric(num_searches),
      C1 = numeric(num_searches),
      C2 = numeric(num_searches),
      D = numeric(num_searches),
      value = numeric(num_searches),
      ID = integer(num_searches)
    )

    for (i in seq_len(num_searches)) {
      if (fb) message(sprintf("Iteration %d/%d...", i, num_searches))

      # Set initial values
      if (i == 1) {
        init_par <- c(tc_init, m_init, o_init)
      } else {
        init_par <- c(
          stats::runif(1, n + 1, n + fh),
          stats::runif(1, lower[1], upper[1]),
          stats::runif(1, lower[2], upper[2])
        )
      }

      # Optimize
      opt_result <- stats::optim(
        par = init_par,
        fn = SSE1,
        log_p = log_p,
        t = t,
        lower = c(n + 1, lower[1], lower[2]),
        upper = c(n + fh, upper[1], upper[2]),
        method = "L-BFGS-B",
        control = list(factr = factr)
      )

      # Calculate linear coefficients
      beta_vals <- beta_calculator(log_p, t, opt_result$par[1],
                                   opt_result$par[2], opt_result$par[3])

      # Calculate damping
      damp <- calculate_damping(
        opt_result$par[2], beta_vals["B"],
        opt_result$par[3], beta_vals["C1"], beta_vals["C2"]
      )

      # Track out-of-range values
      if (beta_vals["B"] >= upper[3]) {
        out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] <- i
      }
      if (damp < lower[4]) {
        out_of_range_tracker$D[[length(out_of_range_tracker$D) + 1]] <- i
      }

      # Store results
      fit_tbl[i, ] <- list(
        tc = opt_result$par[1],
        m = opt_result$par[2],
        omega = opt_result$par[3],
        A = beta_vals["A"],
        B = beta_vals["B"],
        C1 = beta_vals["C1"],
        C2 = beta_vals["C2"],
        D = damp,
        value = opt_result$value,
        ID = i
      )
    }

    # Sort by SSE value
    fit_tbl <- dplyr::arrange(fit_tbl, value)
    fit[[2]] <- fit_tbl
    fit[[1]] <- as.list(fit_tbl[1, ])

    if (fb) message("F1 optimization complete")
  }

  # =========================================================================
  # MODE F2/MPL: Two-stage optimization
  # =========================================================================
  if (mode == "F2" || mode == "MPL") {
    if (fb) message("Mode: F2/MPL - Two-stage optimization")

    # Storage for results at each tc value
    fit2 <- list()
    fit2_filtered <- list()
    fit2_best_for_each_tc <- tibble::tibble(
      ID = integer(fh),
      value = numeric(fh),
      tc = numeric(fh),
      m = numeric(fh),
      omega = numeric(fh),
      A = numeric(fh),
      B = numeric(fh),
      C1 = numeric(fh),
      C2 = numeric(fh),
      D = numeric(fh)
    )
    fit2_best_for_each_tc_filtered <- tibble::tibble(
      ID = integer(0),
      value = numeric(0),
      tc = numeric(0),
      m = numeric(0),
      omega = numeric(0),
      A = numeric(0),
      B = numeric(0),
      C1 = numeric(0),
      C2 = numeric(0),
      D = numeric(0)
    )

    fit[[4]] <- list()

    # Loop over each candidate tc value
    for (k in seq_len(fh)) {
      tc_k <- n + k

      if (fb && k %% 5 == 1) {
        message(sprintf("Processing tc = %d (%d/%d)", tc_k, k, fh))
      }

      # Initialize results for this tc value
      fit2_tmp <- tibble::tibble(
        ID = integer(num_searches),
        value = numeric(num_searches),
        tc = numeric(num_searches),
        m = numeric(num_searches),
        omega = numeric(num_searches),
        A = numeric(num_searches),
        B = numeric(num_searches),
        C1 = numeric(num_searches),
        C2 = numeric(num_searches),
        D = numeric(num_searches)
      )

      for (i in seq_len(num_searches)) {
        # Set seed for reproducibility (matches trace plot)
        set.seed(i)

        # Initial values
        if (i == 1) {
          init_par <- c(m_init, o_init)
        } else {
          init_par <- c(
            stats::runif(1, lower[1], upper[1]),
            stats::runif(1, lower[2], upper[2])
          )
        }

        # Optimize m and omega for fixed tc
        opt_result <- stats::optim(
          par = init_par,
          fn = SSE2,
          tc = tc_k,
          log_p = log_p,
          t = t,
          lower = c(lower[1], lower[2]),
          upper = c(upper[1], upper[2]),
          method = "L-BFGS-B",
          control = list(factr = factr)
        )

        # Calculate linear coefficients
        beta_vals <- beta_calculator(log_p, t, tc_k,
                                     opt_result$par[1], opt_result$par[2])

        # Calculate damping
        damp <- calculate_damping(
          opt_result$par[1], beta_vals["B"],
          opt_result$par[2], beta_vals["C1"], beta_vals["C2"]
        )

        # Track out-of-range
        if (beta_vals["B"] >= upper[3]) {
          out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] <-
            list(tc_num = k, rand_iter_num = i)
        }
        if (damp < lower[4]) {
          out_of_range_tracker$D[[length(out_of_range_tracker$D) + 1]] <-
            list(tc_num = k, rand_iter_num = i)
        }

        # Store results
        fit2_tmp[i, ] <- list(
          ID = i,
          value = opt_result$value,
          tc = tc_k,
          m = opt_result$par[1],
          omega = opt_result$par[2],
          A = beta_vals["A"],
          B = beta_vals["B"],
          C1 = beta_vals["C1"],
          C2 = beta_vals["C2"],
          D = damp
        )
      }

      # Store all fits for this tc
      fit[[4]][[k]] <- fit2_tmp

      # Sort and filter
      fit2[[k]] <- dplyr::arrange(fit2_tmp, value)
      fit2_filtered[[k]] <- dplyr::filter(fit2_tmp, B < upper[3])
      fit2_filtered[[k]] <- dplyr::arrange(fit2_filtered[[k]], value)

      # Keep best fit for this tc
      fit2_best_for_each_tc[k, ] <- fit2[[k]][1, ]

      # Add to filtered list if any fits passed filter
      if (nrow(fit2_filtered[[k]]) > 0) {
        fit2_best_for_each_tc_filtered <- rbind(
          fit2_best_for_each_tc_filtered,
          fit2_filtered[[k]][1, ]
        )
      }
    }

    # Sort results across all tc values
    fit2_best_for_each_tc <- dplyr::arrange(fit2_best_for_each_tc, value)
    fit2_best_for_each_tc_filtered <- dplyr::arrange(
      fit2_best_for_each_tc_filtered, value
    )

    # Select best fit (prefer filtered if available)
    if (nrow(fit2_best_for_each_tc_filtered) == 0) {
      fit2_best <- fit2_best_for_each_tc[1, ]
      warning("No fits passed the filter. Using unfiltered best fit.")
    } else {
      fit2_best <- fit2_best_for_each_tc_filtered[1, ]
    }

    # Final optimization of tc given best m and omega
    if (fb) message("Optimizing tc given best m and omega...")

    opt_final <- stats::optim(
      par = tc_init,
      fn = SSE3,
      m = fit2_best$m,
      omega = fit2_best$omega,
      log_p = log_p,
      t = t,
      lower = n + 1,
      upper = n + fh,
      method = "Brent"
    )

    # Calculate final coefficients
    beta_vals <- beta_calculator(log_p, t, opt_final$par,
                                 fit2_best$m, fit2_best$omega)
    damp <- calculate_damping(
      fit2_best$m, beta_vals["B"],
      fit2_best$omega, beta_vals["C1"], beta_vals["C2"]
    )

    # Store final results
    fit[[1]] <- list(
      tc = opt_final$par,
      m = fit2_best$m,
      omega = fit2_best$omega,
      A = beta_vals["A"],
      B = beta_vals["B"],
      C1 = beta_vals["C1"],
      C2 = beta_vals["C2"],
      D = damp,
      value = opt_final$value
    )

    if (nrow(fit2_best_for_each_tc_filtered) > 0) {
      fit[[2]] <- fit2_best_for_each_tc_filtered
      fit[[3]] <- fit2_best_for_each_tc
    } else {
      fit[[2]] <- fit2_best_for_each_tc
    }

    if (fb) message("F2 optimization complete")
  }

  # =========================================================================
  # MPL: Modified Profile Likelihood (additional calculations)
  # =========================================================================
  if (mode == "MPL") {
    if (fb) message("Computing Modified Profile Likelihood...")

    mpl_output <- compute_mpl(
      fit = fit,
      log_p = log_p,
      t = t,
      fh = fh,
      cutoff = mpl_cutoff,
      fb = fb
    )

    if (mpl_plot) {
      mpl_plot_obj <- create_mpl_plot(
        mpl_output = mpl_output,
        fit = fit,
        n = n,
        fh = fh
      )
    }
  }

  # =========================================================================
  # Generate optional plots
  # =========================================================================

  # Fit plot
  if (fp) {
    fit_plot_obj <- create_fit_plot(
      fit = fit[[1]],
      time_ID = time_ID,
      log_price = log_price,
      n_model = n_model
    )
  }

  # Contour and surface plots
  if (cp || sp) {
    contour_result <- create_contour_data(
      fit = fit[[1]],
      log_p = log_p,
      t = t,
      lower = lower,
      upper = upper,
      beta_calculator = beta_calculator,
      fb = fb
    )
    contour_data <- contour_result

    if (cp) {
      contour_plot_obj <- create_contour_plot(contour_data, fit[[1]], lower)
    }
    if (sp) {
      surface_plot_obj <- create_surface_plot(contour_data, fit[[1]])
    }
  }

  # Parameter plot
  if (pp && length(fit) >= 2) {
    param_plot_obj <- create_param_plot(fit)
  }

  # Matrix plot
  if (mp && (mode == "F2" || mode == "MPL")) {
    matrix_plot_obj <- create_matrix_plot(
      n = n,
      fh = fh,
      lower = lower,
      upper = upper,
      log_p = log_p,
      t = t,
      SSE2_func = SSE2
    )
  }

  # Report timing
  end_time <- Sys.time()
  if (fb) message(sprintf("Total time: %.2f seconds", end_time - start_time))

  # Report out-of-range warnings
  num_estimates <- ifelse(mode == "F1", num_searches, fh * num_searches)
  if (length(out_of_range_tracker$B) > 0) {
    warning(sprintf("%d of %d estimates had B out of range.",
                    length(out_of_range_tracker$B), num_estimates))
  }
  if (length(out_of_range_tracker$D) > 0) {
    warning(sprintf("%d of %d estimates had D out of range.",
                    length(out_of_range_tracker$D), num_estimates))
  }

  # Return results
  list(
    fit = fit,
    mpl_output = mpl_output,
    mpl_plot = mpl_plot_obj,
    fit_plot = fit_plot_obj,
    contour_data = contour_data,
    contour_plot = contour_plot_obj,
    surface_plot = surface_plot_obj,
    trace_plot_mo = trace_plot_mo,
    trace_plot_bm = trace_plot_bm,
    trace_plot_bo = trace_plot_bo,
    param_plot = param_plot_obj,
    matrix_plot = matrix_plot_obj,
    out_of_range_tracker = out_of_range_tracker
  )
}
