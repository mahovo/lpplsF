#' Fit LPPLS Model to Financial Time Series
#'
#' Fits the Log-Periodic Power Law Singularity (LPPLS) model to log-price
#' data using the calibration methods from Filimonov & Sornette (2013).
#' Supports multiple optimization modes and provides diagnostic outputs.
#'
#' @param log_price Numeric vector of log-prices.
#' @param fh Integer, forecast horizon length in time units (default 30).
#'   The search for `tc` extends from `T2+1` to `T2+fh` where `fh>1`.
#' @param hold_out Integer, number of observations to hold out from the end
#'   of the sample for validation (default 0). That is, the number of data
#'   points after `T2`.
#' @param lower Numeric vector of length 4, lower bounds for
#'   `c(m, omega, B, D)`. Default: `c(0.1, 6, -1e14, 0.8)`.\cr
#'   Filimonov 2017:\cr
#'   \eqn{B < 0}, \eqn{D \geq 0.8}, \eqn{0.1 \leq m \leq 0.9}, \eqn{6 \leq \omega \leq 13}.\cr
#'   NOTE: For `trace_plot` mode 2 and 3 (`bm` and `bo`), lower for `B` should
#'   be set to the desired minimum value of the B-axis in the contour
#'   plotting.\cr
#'   Example:\cr
#'   \preformatted{lower = c(0.1, 6, -0.03, 0.8)}
#'   If the optimization algorithm goes below the lower limits at any point
#'   the lower limit for the B-axis will be extended automatically.
#' @param upper Numeric vector of length 4, upper bounds for
#'   `c(m, omega, B, D)`. Default: `c(0.9, 13, -1e-14, 1e6)`.
#' @param tc_init Numeric, initial value of critical time `tc` for optim()
#'   (default 1000). `tc_init` should be in \eqn{[T_2+1, T_3]}, where
#'   \eqn{T_3 = T_2 + s}, and \eqn{s} is the forecast horizon.
#' @param m_init Numeric, initial value of power law exponent `m` for optim()
#'   (default 0.5). `m_init` should be in [`lower[1]`, `upper[1]`].
#'   The algorithm will search for initial values of `m` and `omega` that give
#'   `B<0`. If such values are found, the given init values will be
#'   overwritten.
#' @param o_init Numeric, initial value of frequency `omega` for `optim()`
#'   (default 13). `o_init` should be in [`lower[2]`, `upper[2]`]. Default is
#'   13, which increases the probability of `B<0` (of course not for random
#'   initial values). The algorithm will search for initial values of `m` and
#'   `omega` that give `B<0`. If such values are found, the given init values
#'   will be overwritten.
#' @param num_searches Integer, number of times to repeat the optimization wrt
#'   `tc`, `m` and `omega`. First iteration uses given initial values. Any
#'   subsequent iterations use random initial values.
#' @param mode Character, optimization mode:
#'   \describe{
#'     \item{F1}{Optimize wrt `tc`, `m`, `omega` simultaneously}
#'     \item{F2}{Optimize wrt `m`, `omega` first (for each `tc`), then optimize
#'       wrt `tc`}
#'     \item{MPL}{Modified Profile Likelihood with likelihood intervals}
#'   }
#' @param mpl_cutoff Numeric vector of length 3, cutoff levels (\eqn{c}) for
#'   likelihood interval (default `c(0.05, 0.1, 0.5)`). See Filimonov2017
#'   equation (39). \eqn{c = 0.05} indicates that values of tc outside the
#'   likelihood interval have a probability of 0.05.
#' @param mpl_plot Logical, whether to generate MPL plot (default `FALSE`).
#' @param fp Logical, whether to generate fit plot (default `FALSE`).
#' @param cp Logical, whether to generate contour plot (default `FALSE`). `tc`
#'   value is fixed to the estimate with the lowest objective function value.
#' @param sp Logical, whether to generate surface plot (default `FALSE`). `tc`
#'   value is fixed to the estimate with the lowest objective function value.
#' @param tp Numeric vector of length 3, 0/1 selectors for trace plots.
#'   Plot 1: `m-omega`, plot 2: `B-m`, plot 3:`B-omega` (default `c(0, 0, 0)`).\cr
#'   Example: `c(1, 0, 1)` selects trace plot 1 and 3. Traces convergence of
#'   parameters during optimization of F2 objective function. Only works when
#'   `mode = "F2"`. The plotted trace path is the one for the value of `tc` that
#'   produced the smallest value of `SSE2`.\cr
#'   NOTE: Contour plots in trace plot 2 and 3 only apply to the first step in
#'   the trace. In a trace plot wrt `B` and `m`, `omega` is not fixed (and vice
#'   versa for `m` and `omega`).
#' @param tp_id Integer, which trace step to use for contour plot (default 1).
#'   If the indicated id is too high, 1 will be chosen.
#' @param pp Logical, whether to generate parameter plot (default `FALSE`).
#' @param mp Logical, whether to generate matrix plot (default `FALSE`).
#'   x-axis is `tc`, y-axis is `SSE`.
#'   3 columns of plots for smallest and biggest value of `m`, and value between
#'   the two.
#'   3 columns of plots for smallest and biggest value of `omega`, and value
#'   between the two.
#' @param factr Numeric, convergence tolerance for optim control list when
#'   using L-BFGS-B. See `optim()` documentation.
#' @param fb Logical, whether to print feedback during execution (default
#'   `FALSE`).
#'
#' @return A list containing:
#'   \describe{
#'     \item{fit}{List with fitted model results. In F1/F2 mode:\cr
#'       \describe{
#'         \item{mode = "F1"}{List of fits with random initial values.
#'           \describe{
#'             \item{`fit[[1]]`}{
#'               list of best fit coefficients optimized wrt `m` and `omega`.
#'             }
#'             \item{`fit[[2]]`}{
#'               Tibble of fits with random starting parameters.\cr
#'               Sorted by objective function value (best at top).\cr
#'               Output format:\cr
#'               \preformatted{
#'                 [
#'                   list(tc, m, omega, A, B, C1, C2, D, value_min, ID),
#'                   list(tc, m, omega, A, B, C1, C2, D, value, ID),
#'                   ...,
#'                   list(tc, m, omega, A, B, C1, C2, D, value_max, ID)
#'                 ]
#'               }
#'             }
#'             \item{`fit[[3]]`}{
#'               If no fits passed the filter: Only unfiltered fits returned.
#'             }
#'           }
#'           `ID` is an identifier indicating the order pre-sorting. These are
#'           the indexes reported in "out of range" warnings.\cr
#'         }
#'         \item{mode = "F2"}{List with two or three elements:
#'           \describe{
#'             \item{`fit[[1]]`}{
#'               List of best fit coefficients optimized wrt `m` and `omega`.
#'             }
#'             \item{`fit[[2]]`}{
#'               If at least one fit passed the filter: Tibble of best fits for
#'               objective function `F2` for each value of `tc`. Columns are:
#'               `c(ID, value, tc, m, omega, A, B, C1, C2, damp)`. Sorted by
#'               objective function value (best at top). For each value of `tc`
#'               the best fit is picked from `fit2_tmp`.\cr
#'               If no fits passed the filter: Only unfiltered fits returned.
#'             }
#'             \item{`fit[[3]]`}{
#'               If at least one fit passed the filter: Unfiltered fits returned.\cr
#'               If no fits passed the filter: Nothing is returned in `fit[[3]]`
#'             }
#'             \item{`fit[[4]]`}{
#'                  List of all tmp fits\cr
#'                  For each `tc` value, there is a tibble with all fits using
#'                  random starting points.
#'             }
#'           }
#'         }
#'       }
#'     }
#'     \item{mpl_output}{(MPL mode only) List with likelihood intervals,
#'       R values, log-likelihoods\cr
#'         `mpl = list(LI, R, LL, MLL)`\cr
#'         \describe{
#'              \item{LI}{vector of 3 elements, one for each value of `mpl_cutoff`.}
#'              \item{R}{number, relative likelihood}
#'              \item{LL}{vector of MPL log-likelihoods for each `tc` value}
#'              \item{MLL}{maximum of `LL`}
#'              \item{tc_hat_mpl}{`tc` for which `LL` takes it's maximum}
#'       }
#'     }
#'     \item{mpl_plot}{(if requested) ggplot2 object for Modified Profile
#'       Likelihood plot with likelihood intervals}
#'     \item{fit_plot}{(if requested) ggplot2 object for fit plot}
#'     \item{contour_data}{(if requested) List with `x`, `y` and `z` data for
#'       contour plot of SSE wrt `m` and `omega`.}
#'     \item{contour_plot}{(if requested) plotly object for contour plot}
#'     \item{surface_plot}{(if requested) plotly object for surface plot}
#'     \item{trace_plot_mo}{(if requested) ggplot2 object containing trace plot
#'       for `m` and `omega`}
#'     \item{trace_plot_bm}{(if requested) ggplot2 object containing trace plot
#'       for `B` and `m`}
#'     \item{trace_plot_bo}{(if requested) ggplot2 object for containing trace
#'       plot for `B` and `omega`}
#'     \item{param_plot}{(if requested) ggplot2 object for parameter plot}
#'     \item{matrix_plot}{(if requested) ggplot2 object for matrix plot}
#'     \item{out_of_range_tracker}{List tracking which fits had out-of-range
#'       `B` or damping values\cr
#'       \describe{
#'         \item{F1 mode}{List of random iteration IDs what resulted in
#'           out-of-range parameters}
#'         \item{F2 mode}{List of random tc value IDs and iteration IDs what
#'           resulted in out-of-range parameters.\cr
#'           `tc` value ID start with 1 for `T2+1`. That is, the first time step
#'            after the modelling period, or the first time step of the
#'            forecasting horizon.}
#'       }
#'     }
#'   }
#'
#' @details
#' The LPPLS model is:
#' \deqn{y(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c - t)) + C_2 \sin(\omega \log(t_c - t))]}
#'
#' The model has 7 parameters: \eqn{A}, \eqn{B}, \eqn{C_1}, \eqn{C_2} (linear) and \eqn{t_c}, \eqn{m}, \eqn{\omega} (nonlinear).
#' The Filimonov calibration exploits this structure by solving for the linear
#' parameters analytically given the nonlinear ones.
#'
#' **Mode F1**: Optimizes all nonlinear parameters `c(tc, m, omega)`
#' simultaneously using L-BFGS-B with multiple random starting points.
#'
#' **Mode F2**: A two-stage procedure:
#' \enumerate{
#'   \item For each candidate `tc` value, optimize `c(m, omega)`
#'   \item Select the `tc` that minimizes the overall SSE
#' }
#' This is more robust than F1 for ill-conditioned problems.
#'
#' **Mode MPL**: Modified Profile Likelihood inference from Filimonov et al. (2017).
#' Provides likelihood-based confidence intervals for `tc`.
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
#' result <- fit_lppls(log_price = log_p,
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
#' @importFrom dplyr arrange filter
#' @importFrom tibble tibble
#' @importFrom tidyr gather
#' @importFrom stats optim runif
#' @export
fit_lppls <- function(
  log_price,
  fh = 30,
  hold_out = 0,
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

  if (!mode %in% c("F1", "F2", "MPL")) {
    stop("mode must be one of: 'F1', 'F2', 'MPL'")
  }

  if (length(lower) != 4 || length(upper) != 4) {
    stop("lower and upper must be vectors of length 4")
  }

  start_time <- Sys.time()

  ## Separate modeling and holdout data
  n_total <- length(log_price)
  n <- n_total - hold_out

  time_id <- 1:n_total

  t <- time_id[1:n]
  log_p <- log_price[1:n]

  ## Create beta calculator (uses symbolic math for efficiency)
  beta_calculator <- create_beta_calculator()

  ## Internal Sum of Squared Errors functions
  ##
  ## SSE with all nonlinear parameters as vector
  ## par vector is {tc, m, omega}, used as initial parameter values in optim().
  SSE1 <- function(par, log_p, t) {
    tc <- par[1]
    m <- par[2]
    omega <- par[3]

    ## Ensure tc > max(t)
    if (tc <= max(t)) return(Inf)

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(
      t,
      beta["A"],
      beta["B"],
      beta["C1"],
      beta["C2"],
      tc,
      m,
      omega,
      mode = 0
    )
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  ## SSE with fixed tc
  ## par vector is {m, omega}, used as initial parameter values in optim().
  SSE2 <- function(par, tc, log_p, t) {
    m <- par[1]
    omega <- par[2]

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(
      t,
      beta["A"],
      beta["B"],
      beta["C1"],
      beta["C2"],
      tc,
      m,
      omega,
      mode = 0
    )
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  ## SSE with fixed m and omega
  ## par is tc, used as initial parameter value in optim().
  SSE3 <- function(tc, m, omega, log_p, t) {
    if (tc <= max(t)) return(Inf)

    beta <- beta_calculator(log_p, t, tc, m, omega)
    fitted <- eval_lppls(
      t,
      beta["A"],
      beta["B"],
      beta["C1"],
      beta["C2"],
      tc,
      m,
      omega,
      mode = 0
    )
    sum((log_p - fitted)^2, na.rm = TRUE)
  }

  ## Initialize outputs
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

  ## MODE F1: Simultaneous optimization of tc, m, omega ======================
  if (mode == "F1") {
    if (fb) message("Mode: F1 - Simultaneous optimization")

    ## Initialize results tibble
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

      ## Set initial values
      if (i == 1) {
        init_par <- c(tc_init, m_init, o_init)
      } else {
        init_par <- c(
          stats::runif(1, n + 1, n + fh),
          stats::runif(1, lower[1], upper[1]),
          stats::runif(1, lower[2], upper[2])
        )
      }

      ## Optimize
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

      ## Calculate linear coefficients
      beta_vals <- beta_calculator(
        log_p,
        t,
        opt_result$par[1],
        opt_result$par[2],
        opt_result$par[3]
      )

      ## Calculate damping
      damp <- calculate_damping(
        opt_result$par[2], beta_vals["B"],
        opt_result$par[3], beta_vals["C1"], beta_vals["C2"]
      )

      ## Track out-of-range values
      if (beta_vals["B"] >= upper[3]) {
        out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] <- i
      }
      if (damp < lower[4]) {
        out_of_range_tracker$D[[length(out_of_range_tracker$D) + 1]] <- i
      }

      ## Store results
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

    ## Sort by SSE value
    fit_tbl <- dplyr::arrange(fit_tbl, value)
    fit[[2]] <- fit_tbl
    fit[[1]] <- as.list(fit_tbl[1, ]) ## fit[[2]] is a tibble, fit is a list

    if (fb) message("F1 optimization complete")
  }

  ## MODE F2/MPL: Two-stage optimization =====================================
  if (mode == "F2" || mode == "MPL") { ## "F2"-mode output needed for "MPL"-mode
    if (fb) message("Mode: F2/MPL - Two-stage optimization")

    ## List of fits wrt m and omega for a single fixed value of tc.
    ## This is the combination of m and omega that minimizes SSE2.
    ## For each value of tc, fit2 contains a list of all fits with random
    ## initial values.
    ## For each value of tc_k, fit2 stores the contents of fit2_tmp in
    ## iteration k.
    fit2 <- list()

    fit2_filtered <- list()

    ## Initialize List with best fit optimized wrt m and omega for each value of
    ## tc. I.e. pick the best from fit2 for each tc.
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

    ## Initialize fit2_best_for_each_tc filtered for B < 0
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

    ## Loop over each candidate tc value
    for (k in seq_len(fh)) {
      tc_k <- n + k

      if (fb && k %% 5 == 1) {
        message(sprintf("Processing tc = %d (%d/%d)", tc_k, k, fh))
      }

      ## Initialize results for this tc value
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

      ## Use for trace plot
      opt2_counts <- list()

      for (i in seq_len(num_searches)) {
        ## Set seed for reproducibility (matches trace plot)
        set.seed(i)

        ## Initial values
        if (i == 1) {
          init_par <- c(m_init, o_init)
        } else {
          init_par <- c(
            stats::runif(1, lower[1], upper[1]),
            stats::runif(1, lower[2], upper[2])
          )
        }

        ## Optimize m and omega for fixed tc
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

        if (any(tp != 0)) opt2_counts[[i]] <- opt_result$counts[["function"]]

        ## Calculate linear coefficients
        beta_vals <- beta_calculator(
          log_p,
          t,
          tc_k,
          opt_result$par[1],
          opt_result$par[2]
        )

        ## Calculate damping
        damp <- calculate_damping(
          opt_result$par[1],
          beta_vals["B"],
          opt_result$par[2],
          beta_vals["C1"],
          beta_vals["C2"]
        )

        ## Track out-of-range
        if (beta_vals["B"] >= upper[3]) {
          out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] <-
            list(tc_num = k, rand_iter_num = i)
        }
        if (damp < lower[4]) {
          out_of_range_tracker$D[[length(out_of_range_tracker$D) + 1]] <-
            list(tc_num = k, rand_iter_num = i)
        }

        ## Add list of fitted coefficients to list of fits
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

      ## Store all fits for this tc
      fit[[4]][[k]] <- fit2_tmp

      ## Sort list by value of objective function F2.
      ## Keep all fits for tc_k in a list.
      ## Filter out results where B >= 0.
      ## Then sort by SSE value (smallest at top).
      fit2[[k]] <- dplyr::arrange(fit2_tmp, value)
      fit2_filtered[[k]] <- dplyr::filter(fit2_tmp, B < upper[3])
      fit2_filtered[[k]] <- dplyr::arrange(fit2_filtered[[k]], value)

      ## Keep best fit for this tc
      fit2_best_for_each_tc[k, ] <- fit2[[k]][1, ]

      ## Add to filtered list if any fits passed filter
      if (nrow(fit2_filtered[[k]]) > 0) {
        fit2_best_for_each_tc_filtered <- rbind(
          fit2_best_for_each_tc_filtered,
          fit2_filtered[[k]][1, ]
        )
      }
    }

    ## Sort results across all tc values
    fit2_best_for_each_tc <- dplyr::arrange(fit2_best_for_each_tc, value)
    fit2_best_for_each_tc_filtered <- dplyr::arrange(
      fit2_best_for_each_tc_filtered, value
    )

    ## Select best fit (prefer filtered if available)
    if (nrow(fit2_best_for_each_tc_filtered) == 0) {
      fit2_best <- fit2_best_for_each_tc[1, ]
      warning("No fits passed the filter. Using unfiltered best fit.")
    } else {
      fit2_best <- fit2_best_for_each_tc_filtered[1, ]
    }

    ## Final optimization of tc given best m and omega
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

    ## Calculate final coefficients
    beta_vals <- beta_calculator(
      log_p,
      t,
      opt_final$par,
      fit2_best$m,
      fit2_best$omega
    )

    damp <- calculate_damping(
      fit2_best$m,
      beta_vals["B"],
      fit2_best$omega,
      beta_vals["C1"],
      beta_vals["C2"]
    )

    ## Final results
    ## First element is the best fit wrt tc. This is the tc that minimizes SSE3.
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
      if (fb) message("No fits passed the filter. Only unfiltered fits returned.\n")
    }

    if (fb) message("F2 optimization complete")
  }

  ## MPL: Modified Profile Likelihood (additional calculations) ==============
  if (mode == "MPL") {
    if (fb) message("Computing Modified Profile Likelihood...")

    mpl_output <- compute_mpl(
      fit = fit,
      log_p = log_p,
      t = t,
      fh = fh,
      cutoff = mpl_cutoff,
      lower = lower,
      upper = upper,
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

  ## Generate optional plots =================================================

  ## Fit plot
  if (fp) {
    fit_plot_obj <- fit_plot(
      fit       = fit[[1]],
      time_ID   = time_id,
      log_price = log_price,
      mode      = 0,
      T2        = n            # n was the old n_model (the T2 marker position)
    )
  }

  ## Contour and surface plots
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

  ## Parameter plot
  if (pp && length(fit) >= 2) {
    param_plot_obj <- create_param_plot(fit)
  }

  ## Matrix plot
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

  ## Trace plots (F2/MPL mode only) ==========================================
  if (any(tp != 0)) {
    if (mode == "F2" || mode == "MPL") {
      tc_val <- fit[[1]]$tc
      need_B <- (tp[2] == 1 || tp[3] == 1)

      ## Replay the best search's optimization, tracing the parameter path.
      ## Start from the same initial point that search used (cf. the search loop):
      ## (m_init, o_init) for the first search, otherwise its random start. This
      ## makes the first plotted segment the real first BFGS step.
      s <- seq_len(opt2_counts[[fit2_best$ID]])
      set.seed(fit2_best$ID)
      init_replay <- if (fit2_best$ID == 1L) {
        c(m_init, o_init)
      } else {
        c(stats::runif(1, lower[1], upper[1]),
          stats::runif(1, lower[2], upper[2]))
      }

      opt_trace <- sapply(s, function(iter) {
        opt_i <- stats::optim(
          par     = init_replay,
          fn      = SSE2, tc = tc_val, log_p = log_p, t = t,
          lower   = c(lower[1], lower[2]), upper = c(upper[1], upper[2]),
          method  = "L-BFGS-B", control = list(maxit = iter)
        )$par
        if (need_B) {
          opt_i <- c(
            opt_i,
            unname(beta_calculator(log_p, t, tc_val, opt_i[1], opt_i[2])["B"])
          )
        }
        opt_i
      })

      ## Prepend the actual starting point (the red start dot)
      if (need_B) {
        init_B <- unname(
          beta_calculator(log_p, t, tc_val, init_replay[1], init_replay[2])["B"]
        )
        opt_trace <- matrix(
          c(init_replay[1], init_replay[2], init_B, opt_trace),
          byrow = FALSE,
          nrow = 3
        )
      } else {
        opt_trace <- matrix(c(init_replay[1], init_replay[2], opt_trace),
                            byrow = FALSE, nrow = 2)
      }

      ## Generate requested trace plots
      if (tp[1] == 1) {
        trace_plot_mo <- create_trace_plot(
          1, opt_trace, lower, upper, SSE2, beta_calculator,
          tc_val, log_p, t, tp_id, fb
        )
      }
      if (tp[2] == 1) {
        trace_plot_bm <- create_trace_plot(
          2, opt_trace, lower, upper, SSE2, beta_calculator,
          tc_val, log_p, t, tp_id, fb
        )
      }
      if (tp[3] == 1) {
        trace_plot_bo <- create_trace_plot(
          3, opt_trace, lower, upper, SSE2, beta_calculator,
          tc_val, log_p, t, tp_id, fb
        )
      }
    } else {
      warning("Trace plots are only available when mode = 'F2' or 'MPL'")
    }
  }

  ## Report timing
  end_time <- Sys.time()
  if (fb) message(sprintf("Total time: %.2f seconds", end_time - start_time))

  ## Report out-of-range warnings
  num_estimates <- ifelse(mode == "F1", num_searches, fh * num_searches)
  if (length(out_of_range_tracker$B) > 0) {
    warning(sprintf("%d of %d estimates had B out of range.",
                    length(out_of_range_tracker$B), num_estimates))
  }
  if (length(out_of_range_tracker$D) > 0) {
    warning(sprintf("%d of %d estimates had D out of range.",
                    length(out_of_range_tracker$D), num_estimates))
  }

  ## Return results
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
