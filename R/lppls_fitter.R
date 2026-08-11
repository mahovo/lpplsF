#' Cap on the reconstructed BFGS trace
#'
#' The trace plots rebuild the optimiser's path by re-running L-BFGS-B with an
#' increasing iteration limit until the iterate stops changing (see the trace
#' block of [fit_lppls()]). This bounds that search, so a path that never
#' settles cannot loop indefinitely. Well above anything observed: the SPX and
#' e008/e009 calibrations settle within 1-7 iterations.
#'
#' @keywords internal
#' @noRd
.lppls_max_trace_steps <- 200L


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
#' @param factr Numeric, convergence tolerance for the L-BFGS-B `control` list
#'   (default `1e7`, which is also [stats::optim()]'s own default). It is a
#'   **multiplier on machine epsilon**, not an absolute tolerance: the optimiser
#'   stops when the relative reduction in the objective falls below
#'   `factr * .Machine$double.eps`, so `1e7` asks for roughly `1e-9`.
#'
#'   Smaller is not automatically better. Any value below about `1e0` requests a
#'   tolerance finer than double precision can represent, which cannot be
#'   satisfied: L-BFGS-B then runs until its line search fails and returns
#'   convergence code 52 rather than converging. Earlier versions of this package
#'   defaulted to `1e-08`, which requests `2.2e-24` and had exactly that effect --
#'   about five times more objective evaluations than `1e7`, ending in
#'   line-search failure for two thirds of the candidate `tc` values, and
#'   reaching the same optimum. See `misc/bench_beta.Rmd` for the measurements.
#' @param fb Logical, whether to print feedback during execution (default
#'   `FALSE`).
#' @param beta_method Character, how to solve for the linear parameters
#'   \eqn{A, B, C_1, C_2} given \eqn{(t_c, m, \omega)}. All choices compute the
#'   same least-squares solution and differ only in speed and conditioning; the
#'   solver is evaluated once per objective evaluation, so it dominates the cost
#'   of a calibration.
#'   \describe{
#'     \item{`"crossprod"`}{(default) Normal equations via [crossprod()] and
#'       [solve()]. Fastest, and exact to well within anything interpretable
#'       over the admissible \eqn{0.1 \le m \le 0.9}{0.1 <= m <= 0.9}.}
#'     \item{`"qr"`}{Householder QR via [.lm.fit()]. Never forms \eqn{X'X}, so it
#'       does not square the condition number. About 30% slower and markedly
#'       more accurate; prefer it if `lower[1]` is widened toward
#'       \eqn{m \to 0}{m -> 0}.}
#'     \item{`"chol"`}{Cholesky of the normal equations.}
#'     \item{`"symengine"`}{The thesis implementation: a symbolic 4x4 solve
#'       compiled to one closure per coefficient. Retained so the thesis
#'       calibration can be reproduced and inspected directly -- it is roughly
#'       ten times slower than the default and the least accurate of the four.}
#'   }
#'   See `misc/bench_beta.Rmd` for the timings and precision measurements behind
#'   these characterisations.
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
#'               the best fit is picked from that `tc`'s table in `fit[[4]]`.\cr
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
#'     \item{log_price}{The input series, so diagnostics that need the data do
#'       not make the caller supply it again. The calibration slice is
#'       `log_price[seq_len(fit_args$n)]`.}
#'     \item{fit_args}{List recording the settings this calibration ran with --
#'       `mode`, `n` (calibration length), `fh`, `hold_out`, `lower`, `upper`,
#'       `num_searches` and `mpl_cutoff`. Diagnostics need these: whether
#'       `m = 0.9` is a boundary solution or an interior optimum cannot be told
#'       from the estimate alone.}
#'   }
#'
#'   The list has class `"lppls_fit"`. A calibration can look convincing as a
#'   plot and still not be worth interpreting, so the object carries its own
#'   diagnostics:
#'   \describe{
#'     \item{[print.lppls_fit()]}{Printing the fit reports the estimate together
#'       with the checks that decide whether to trust it -- bounds, optimisation
#'       basins, search filters and the extent of the modified profile
#'       likelihood. Warnings are listed only when something is wrong.}
#'     \item{[summary.lppls_fit()]}{The same as values rather than text, plus the
#'       curvature check and whether each likelihood interval is closed by the
#'       likelihood or by the edge of the `tc` grid.}
#'     \item{[plot.lppls_fit()]}{Any plot attached at fit time, or the basin
#'       views, drawn on demand.}
#'     \item{[lppls_curvature()]}{Whether the estimate sits at a genuine interior
#'       minimum, which is what the modified profile likelihood assumes.}
#'   }
#'
#' @details
#' The LPPLS model is:
#' \deqn{y(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c - t)) + C_2 \sin(\omega \log(t_c - t))]}{y(t) = A + (tc - t)^m [B + C1 cos(omega log(tc - t)) + C2 sin(omega log(tc - t))]}
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
#' @section Profiling `mode = "MPL"` on macOS:
#' Running a sampling profiler over an MPL calibration can **crash the R
#' session** ("Illegal instruction: 4") on macOS builds that use Apple's
#' Accelerate framework as their BLAS -- which is the default for the CRAN macOS
#' build of R. This affects [utils::Rprof()] and everything layered on it
#' (`profvis`, the RStudio profiler, `proftools`).
#'
#' The fault is an interaction between R's profiler and Accelerate's `DGEMM`,
#' not a fault in this package: it reproduces in plain R with no packages
#' loaded, via `Rprof(tempfile()); for (i in 1:2e5) crossprod(A, B)` for
#' matrices `A`, `B` of a few hundred rows or more. MPL mode is the only part of
#' lpplsF that reaches it, because the modified profile likelihood is the only
#' place that multiplies two matrices; every other product here is a
#' one-argument [crossprod()], a matrix-vector product or a decomposition, none
#' of which are affected. Fitting itself is unaffected -- the crash requires a
#' profiler to be running.
#'
#' To profile MPL mode, take R's own matrix products instead of the BLAS for the
#' duration of the profiling session:
#'
#' ```r
#' old <- options(matprod = "internal")
#' Rprof("prof.out")
#' fit <- fit_lppls(log_price, mode = "MPL", ...)
#' Rprof(NULL)
#' options(old)
#' summaryRprof("prof.out")
#' ```
#'
#' `matprod = "internal"` costs only a few percent on a whole MPL fit, since the
#' fit does relatively few matrix-matrix products. It does sum in a different
#' order, so estimates move in the last couple of digits; that is immaterial for
#' a profiling run but is a reason not to set it globally.
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
#' @seealso [print.lppls_fit()], [summary.lppls_fit()], [plot.lppls_fit()] and
#'   [lppls_curvature()] for diagnosing the returned fit; [lppls_rolling()] to
#'   recalibrate over a rolling window.
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
  factr = 1e7,
  fb = FALSE,
  beta_method = c("crossprod", "qr", "chol", "symengine")
) {

  beta_method <- match.arg(beta_method)

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

  ## Create the beta calculator once, with the requested solver.
  beta_calculator <- create_beta_calculator(beta_method)

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

    ## Accumulators for the tc loop.
    ##
    ## These are plain preallocated vectors, assembled into the documented
    ## tibbles once the loop has finished. The loop runs fh x num_searches
    ## times and everything in it other than optim() used to be tibble row
    ## assignment, two dplyr::arrange() calls and a dplyr::filter() per tc,
    ## plus an rbind() that regrew the filtered table on every iteration --
    ## quadratic in fh, and ~1.2 ms per tc in total.
    ##
    ## Column order is the documented one: ID, value, tc, m, omega, A, B, C1,
    ## C2, D. ID is an integer identifier giving the pre-sorting order.

    ## Best fit for each tc, one row per candidate (index k).
    b_ID <- integer(fh)
    b_value <- numeric(fh)
    b_tc <- numeric(fh)
    b_m <- numeric(fh)
    b_omega <- numeric(fh)
    b_A <- numeric(fh)
    b_B <- numeric(fh)
    b_C1 <- numeric(fh)
    b_C2 <- numeric(fh)
    b_D <- numeric(fh)

    ## The same, restricted to fits with B < upper[3]. At most one row per tc,
    ## so fh entries is always enough; n_keep records how many were used.
    f_ID <- integer(fh)
    f_value <- numeric(fh)
    f_tc <- numeric(fh)
    f_m <- numeric(fh)
    f_omega <- numeric(fh)
    f_A <- numeric(fh)
    f_B <- numeric(fh)
    f_C1 <- numeric(fh)
    f_C2 <- numeric(fh)
    f_D <- numeric(fh)
    n_keep <- 0L

    ## Scratch for the random searches at a single tc, overwritten each k.
    s_ID <- seq_len(num_searches)
    s_value <- numeric(num_searches)
    s_m <- numeric(num_searches)
    s_omega <- numeric(num_searches)
    s_A <- numeric(num_searches)
    s_B <- numeric(num_searches)
    s_C1 <- numeric(num_searches)
    s_C2 <- numeric(num_searches)
    s_D <- numeric(num_searches)

    fit[[4]] <- list()

    ## Loop over each candidate tc value
    for (k in seq_len(fh)) {
      tc_k <- n + k

      if (fb && k %% 5 == 1) {
        message(sprintf("Processing tc = %d (%d/%d)", tc_k, k, fh))
      }

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

        ## Record this search's fitted coefficients
        s_value[i] <- opt_result$value
        s_m[i] <- opt_result$par[1]
        s_omega[i] <- opt_result$par[2]
        s_A[i] <- beta_vals["A"]
        s_B[i] <- beta_vals["B"]
        s_C1[i] <- beta_vals["C1"]
        s_C2[i] <- beta_vals["C2"]
        s_D[i] <- damp
      }

      ## Store all fits for this tc, in search order (ID is the pre-sorting
      ## index, so this table is deliberately left unsorted).
      fit[[4]][[k]] <- tibble::tibble(
        ID = s_ID,
        value = s_value,
        tc = rep(as.numeric(tc_k), num_searches),
        m = s_m,
        omega = s_omega,
        A = s_A,
        B = s_B,
        C1 = s_C1,
        C2 = s_C2,
        D = s_D
      )

      ## Keep the best fit for this tc, i.e. the smallest objective value.
      ## which.min() returns the FIRST minimum, which is what sorting by
      ## `value` and taking the top row did (arrange() is a stable sort); it
      ## also skips NaN, as arrange() sorted them last. If every value is NaN
      ## neither leaves a candidate, and sorting would have left row 1 in
      ## place, so fall back to it.
      j <- which.min(s_value)
      if (length(j) == 0L) j <- 1L

      b_ID[k] <- s_ID[j]
      b_value[k] <- s_value[j]
      b_tc[k] <- tc_k
      b_m[k] <- s_m[j]
      b_omega[k] <- s_omega[j]
      b_A[k] <- s_A[j]
      b_B[k] <- s_B[j]
      b_C1[k] <- s_C1[j]
      b_C2[k] <- s_C2[j]
      b_D[k] <- s_D[j]

      ## The same, among fits that pass the B < upper[3] filter. which() drops
      ## NA comparisons exactly as dplyr::filter() dropped NA rows.
      keep <- which(s_B < upper[3])
      if (length(keep) > 0L) {
        jf <- which.min(s_value[keep])
        jf <- if (length(jf) == 0L) keep[1L] else keep[jf]

        n_keep <- n_keep + 1L
        f_ID[n_keep] <- s_ID[jf]
        f_value[n_keep] <- s_value[jf]
        f_tc[n_keep] <- tc_k
        f_m[n_keep] <- s_m[jf]
        f_omega[n_keep] <- s_omega[jf]
        f_A[n_keep] <- s_A[jf]
        f_B[n_keep] <- s_B[jf]
        f_C1[n_keep] <- s_C1[jf]
        f_C2[n_keep] <- s_C2[jf]
        f_D[n_keep] <- s_D[jf]
      }
    }

    ## Assemble the documented tables from the accumulators
    fit2_best_for_each_tc <- tibble::tibble(
      ID = b_ID, value = b_value, tc = b_tc, m = b_m, omega = b_omega,
      A = b_A, B = b_B, C1 = b_C1, C2 = b_C2, D = b_D
    )

    kept <- seq_len(n_keep)
    fit2_best_for_each_tc_filtered <- tibble::tibble(
      ID = f_ID[kept], value = f_value[kept], tc = f_tc[kept],
      m = f_m[kept], omega = f_omega[kept], A = f_A[kept], B = f_B[kept],
      C1 = f_C1[kept], C2 = f_C2[kept], D = f_D[kept]
    )

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
      fb = fb,
      beta_calculator = beta_calculator
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
      set.seed(fit2_best$ID)
      init_replay <- if (fit2_best$ID == 1L) {
        c(m_init, o_init)
      } else {
        c(stats::runif(1, lower[1], upper[1]),
          stats::runif(1, lower[2], upper[2]))
      }

      ## The path is reconstructed rather than recorded: L-BFGS-B is re-run from
      ## init_replay with maxit = 1, 2, 3, ..., so each run contributes the point
      ## the optimiser had reached after that many iterations.
      ##
      ## How many iterations that is cannot be read off the search: optim()
      ## reports counts of calls to fn and gr, not iterations, and those exceed
      ## the number of steps by a factor of roughly 2-12 here. So the length is
      ## discovered instead -- extend maxit until the iterate stops changing.
      ## The run that shows no change ends the loop and is not plotted, so the
      ## trace holds every distinct point and no duplicates.
      trace_step <- function(iter) {
        stats::optim(
          par     = init_replay,
          fn      = SSE2, tc = tc_val, log_p = log_p, t = t,
          lower   = c(lower[1], lower[2]), upper = c(upper[1], upper[2]),
          method  = "L-BFGS-B", control = list(maxit = iter)
        )$par
      }

      path <- vector("list", .lppls_max_trace_steps)
      n_step <- 0L
      prev <- NULL
      for (iter in seq_len(.lppls_max_trace_steps)) {
        p_i <- trace_step(iter)
        if (!is.null(prev) && identical(p_i, prev)) break
        n_step <- n_step + 1L
        path[[n_step]] <- p_i
        prev <- p_i
      }
      if (n_step == .lppls_max_trace_steps) {
        warning("BFGS trace still moving after ", .lppls_max_trace_steps,
                " iterations; the plotted path is truncated.")
      }

      opt_trace <- vapply(path[seq_len(n_step)], function(p_i) {
        if (need_B) {
          c(p_i, unname(beta_calculator(log_p, t, tc_val, p_i[1], p_i[2])["B"]))
        } else {
          p_i
        }
      }, numeric(if (need_B) 3L else 2L))

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
  structure(
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
      out_of_range_tracker = out_of_range_tracker,
      ## The input series, so diagnostics that need the data (curvature at the
      ## estimate, residuals) do not make the caller supply it again. The
      ## calibration slice is log_price[seq_len(fit_args$n)].
      log_price = log_price,
      fit_args = list(
        mode         = mode,
        n            = n,
        fh           = fh,
        hold_out     = hold_out,
        lower        = lower,
        upper        = upper,
        num_searches = num_searches,
        mpl_cutoff   = mpl_cutoff,
        beta_method  = beta_method
      )
    ),
    class = "lppls_fit"
  )
}


#' Diagnostic quantities behind `print.lppls_fit()`
#'
#' Kept separate from the printing so the checks can be tested on values instead
#' of on formatted output. Components are `NULL` when the fit does not carry what
#' they need -- `fit_args` is absent from objects made before it was recorded, and
#' without `lower`/`upper` no boundary check is possible.
#'
#' @param x An `"lppls_fit"` object.
#' @param gap_tol Share of the SSE spread the largest gap must reach before the
#'   fits count as splitting into a good and an inferior cluster. Default 0.3.
#' @param tol Relative tolerance for "at a bound", as in `.boundary_pars()`.
#'   Default 1e-3.
#'
#' @keywords internal
#' @noRd
.fit_diagnostics <- function(x, gap_tol = 0.3, tol = 1e-3) {
  a   <- x$fit_args
  est <- x$fit[[1]]
  tbl <- if (length(x$fit) >= 2L && is.data.frame(x$fit[[2]])) x$fit[[2]] else NULL

  pick <- function(v, empty) if (is.null(v)) empty else v
  d <- list(
    mode         = pick(a$mode, NA_character_),
    n            = pick(a$n, NA_integer_),
    fh           = pick(a$fh, NA_integer_),
    hold_out     = pick(a$hold_out, NA_integer_),
    num_searches = pick(a$num_searches, NA_integer_),
    estimate     = est,
    n_fits       = if (is.null(tbl)) NA_integer_ else nrow(tbl)
  )

  ## Boundary solutions: the point estimate, then the whole table of fits.
  bounds_known <- !is.null(a$lower) && !is.null(a$upper)
  d$pinned <- if (bounds_known) .boundary_pars(est, a$lower, a$upper, tol) else NULL
  d$boundary_counts <- if (bounds_known && !is.null(tbl)) {
    .boundary_counts(tbl, a$lower, a$upper, tol)
  } else {
    NULL
  }

  ## Optimisation basins. Distance from the best SSE will not do on its own: In
  ## "F2"/"MPL" the conditional SSE legitimately rises away from the optimal tc,
  ## so a perfectly good profile also has fits well above the minimum. What marks
  ## a failed search is a *gap* -- the sorted SSEs fall into a good cluster and a
  ## distinctly worse one. So measure the largest gap between consecutive sorted
  ## values as a share of the total spread: scale-free, bounded in [0, 1], and
  ## near 0 for a smooth profile whatever the noise level.
  if (!is.null(tbl) && sum(is.finite(tbl$value)) >= 3L) {
    v    <- sort(tbl$value[is.finite(tbl$value)])
    span <- diff(range(v))
    d$best_sse  <- v[1]
    d$worst_sse <- v[length(v)]
    if (span > 0) {
      gaps <- diff(v)
      i    <- which.max(gaps)
      d$basin_gap     <- gaps[i] / span
      d$n_inferior    <- length(v) - i
      d$frac_inferior <- d$n_inferior / length(v)
      d$split_at      <- c(v[i], v[i + 1L])
      d$gap_tol       <- gap_tol
      d$split         <- d$basin_gap >= gap_tol
    }
  }

  ## Search filters: B >= upper[3] and damping < lower[4].
  d$oor_B <- length(x$out_of_range_tracker$B)
  d$oor_D <- length(x$out_of_range_tracker$D)

  ## Modified profile likelihood, when it was computed.
  m <- x$mpl_output
  if (!is.null(m)) {
    ll   <- m$LL
    ok   <- which(!is.na(ll))
    mple <- m$tc_hat_mpl
    d$mpl <- list(
      defined = length(ok),
      total   = length(ll),
      mple    = mple,
      ## An argmax on the first or last defined grid point means the maximum is
      ## where the curve stops, not where the likelihood turns over.
      at_edge = if (length(ok) == 0L) NA_character_
                else if (which.max(ll) == max(ok)) "upper"
                else if (which.max(ll) == min(ok)) "lower"
                else NA_character_,
      cutoff  = a$mpl_cutoff,
      li      = m$LI,
      ## A bound sitting exactly on the estimate is a truncated interval: It was
      ## closed by the end of the defined region, not by a drop in likelihood.
      li_touches_mple = vapply(
        m$LI,
        function(li) !any(is.na(li)) && !is.na(mple) && any(li == mple),
        logical(1)
      )
    )

    ## Is each interval closed by the likelihood, or only by the grid?
    ## compute_mpl() takes range() of the tc whose R exceeds log(cutoff). If that
    ## set has holes the range spans them silently, so the reported interval
    ## contains grid points that are not themselves in it. That happens when the
    ## curve leaves the good basin rather than descending through the cutoff, and
    ## it is the difference between "narrow because the likelihood is sharp" and
    ## "narrow because the search only worked here".
    if (!is.null(a$mpl_cutoff) && length(a$mpl_cutoff) == length(m$LI)) {
      tc_grid <- seq_along(ll) + d$n
      d$mpl$li_check <- lapply(seq_along(m$LI), function(j) {
        ids <- which(log(a$mpl_cutoff[j]) < m$R)
        if (!length(ids)) {
          return(list(
            cutoff         = a$mpl_cutoff[j],
            n_in           = 0L,
            span           = NA_integer_,
            contiguous     = NA,
            bridged        = NA_integer_,
            all_good_basin = NA
          ))
        }
        span <- max(ids) - min(ids) + 1L
        ## Do the members coincide with the good basin? If every tc inside is a
        ## good-basin fit, the interval is tracing the search, not the curvature.
        good <- NA
        if (isTRUE(d$split) && !is.null(tbl) && !is.na(d$n)) {
          v_in <- tbl$value[match(tc_grid[ids], tbl$tc)]
          good <- all(!is.na(v_in) & v_in <= d$split_at[1])
        }
        list(
          cutoff         = a$mpl_cutoff[j],
          n_in           = length(ids),
          span           = span,
          contiguous     = length(ids) == span,
          bridged        = span - length(ids),
          all_good_basin = good
        )
      })
    }
  }

  d
}


#' Print an LPPLS calibration with its diagnostics
#'
#' Prints the point estimate together with the checks that decide whether it can
#' be trusted: whether the solution sits at an optimisation bound, how many of the
#' individual fits reached the best basin found, how many were filtered on `B` or
#' damping, and -- in `"MPL"` mode -- whether the modified profile likelihood is
#' defined across the whole grid and whether its likelihood intervals are closed
#' by the likelihood or merely by the edge of that grid.
#'
#' Findings that undermine the fit are repeated as a `Diagnostics` block, so a
#' degenerate calibration announces itself rather than having to be inferred from
#' the numbers.
#'
#' @param x An object of class `"lppls_fit"` from [fit_lppls()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [fit_lppls()], [summary.lppls_fit()], [lppls_curvature()]
#'
#' @export
print.lppls_fit <- function(x, ...) {
  .print_fit_diagnostics(.fit_diagnostics(x))
  invisible(x)
}


#' Format a diagnostic summary
#'
#' Shared by [print.lppls_fit()] and [print.summary.lppls_fit()]; the latter
#' passes a curvature report to append.
#'
#' @keywords internal
#' @noRd
.print_fit_diagnostics <- function(d, curv = NULL) {
  e    <- d$estimate
  unit <- if (identical(d$mode, "F1")) "restarts" else "conditional fits"

  cat(sprintf("LPPLS calibration (mode = %s)\n", d$mode))
  cat(sprintf("  Window:    n = %s   fh = %s   hold_out = %s   num_searches = %s\n",
              d$n, d$fh, d$hold_out, d$num_searches))
  cat(sprintf("  Estimate:  tc = %.2f   m = %.4f   omega = %.4f   D = %.3f   SSE = %.6g\n",
              e$tc, e$m, e$omega, e$D, e$value))

  if (is.null(d$boundary_counts)) {
    cat("  Bounds:    not recorded for this fit\n")
  } else {
    cat(sprintf("  Bounds:    estimate %s;  fits at a bound: m %d/%d, omega %d/%d\n",
                if (length(d$pinned)) {
                  paste0(paste(d$pinned, collapse = " and "), " at a bound")
                } else {
                  "interior"
                },
                d$boundary_counts[["m"]], d$n_fits,
                d$boundary_counts[["omega"]], d$n_fits))
  }

  if (!is.null(d$basin_gap)) {
    if (isTRUE(d$split)) {
      cat(sprintf(
        "  Basins:    SSE %.6g..%.6g;  largest gap %.2f of the spread splits off %d/%d %s (%.0f%%)\n",
        d$best_sse, d$worst_sse, d$basin_gap, d$n_inferior, d$n_fits, unit,
        100 * d$frac_inferior))
    } else {
      cat(sprintf(
        "  Basins:    SSE %.6g..%.6g;  no distinct inferior cluster (largest gap %.2f of the spread)\n",
        d$best_sse, d$worst_sse, d$basin_gap))
    }
  }

  cat(sprintf("  Filters:   B out of range %d   D out of range %d\n", d$oor_B, d$oor_D))

  if (!is.null(d$mpl)) {
    cat(sprintf("  MPL:       defined at %d/%d tc   MPLE = %s\n",
                d$mpl$defined, d$mpl$total, d$mpl$mple))
    if (length(d$mpl$li)) {
      lab <- if (length(d$mpl$cutoff) == length(d$mpl$li)) {
        sprintf("LI(%.2f)", d$mpl$cutoff)
      } else {
        sprintf("LI[%d]", seq_along(d$mpl$li))
      }
      cat("             ",
          paste(mapply(function(l, li) sprintf("%s [%s, %s]", l, li[1], li[2]),
                       lab, d$mpl$li),
                collapse = "   "),
          "\n", sep = "")
    }
  }

  flags <- character(0)
  if (!is.null(d$pinned) && length(d$pinned) > 0L) {
    flags <- c(flags, sprintf(
      "%s pinned at an optimisation bound - a boundary solution, so the modified profile likelihood and its intervals are unreliable",
      paste(d$pinned, collapse = " and ")))
  }
  if (isTRUE(d$split)) {
    flags <- c(flags, sprintf(
      "%d of %d %s (%.0f%%) sit in a distinct inferior basin (SSE >= %.6g against a best of %.6g) - raise num_searches",
      d$n_inferior, d$n_fits, unit, 100 * d$frac_inferior,
      d$split_at[2], d$best_sse))
  }
  if (!is.null(d$mpl)) {
    if (d$mpl$defined < d$mpl$total) {
      flags <- c(flags, sprintf(
        "MPL undefined at %d of %d tc (det(X'X - H) <= 0 there)",
        d$mpl$total - d$mpl$defined, d$mpl$total))
    }
    if (!is.na(d$mpl$at_edge)) {
      flags <- c(flags, sprintf(
        "MPLE is the %s edge of the defined region, not an interior maximum",
        d$mpl$at_edge))
    }
    if (any(d$mpl$li_touches_mple)) {
      flags <- c(flags,
        "a likelihood interval is bounded by the MPLE itself - closed by the grid, not by the likelihood")
    }
  }

  ## Intervals: closed by the likelihood, or bridged across a hole?
  if (!is.null(d$mpl$li_check)) {
    for (k in d$mpl$li_check) {
      if (isTRUE(k$contiguous) || is.na(k$contiguous)) next
      flags <- c(flags, sprintf(
        "LI(%.2f) spans %d grid points but only %d are above its cutoff - range() bridges %d that are not in the interval",
        k$cutoff, k$span, k$n_in, k$bridged))
    }
    if (any(vapply(d$mpl$li_check, function(k) isTRUE(k$all_good_basin), logical(1)))) {
      flags <- c(flags,
        "every tc inside an interval is a good-basin fit - the interval traces where the search succeeded, not the curvature of the likelihood")
    }
  }

  if (!is.null(curv)) {
    cat(sprintf("  Curvature: %s at the estimate  (det(X'X - H) = %+.3e;  profile Hessian eigenvalues %s)\n",
                curv$classification, curv$det_full,
                paste(sprintf("%+.3e", curv$eigen_profile), collapse = ", ")))
    if (!identical(curv$classification, "positive definite")) {
      flags <- c(flags, sprintf(
        "the profiled (m, omega) Hessian is %s, so det(X'X - H) <= 0 and the MPL is undefined at the estimate%s",
        curv$classification,
        if (length(d$pinned)) " - unsurprising, since a bound is active" else ""))
    }
  }

  if (length(flags)) {
    cat("\n  Diagnostics:\n")
    cat(paste0("    ! ", flags, collapse = "\n"), "\n", sep = "")
  }

  invisible(d)
}


#' Plot an LPPLS calibration
#'
#' Returns one of the plots attached to the fit, or draws the basin view. The
#' attached ones exist only if they were asked for at fit time (`fp`, `pp`, `mp`,
#' `tp`, `mpl_plot`, `cp`, `sp`); `"basin"` is always available, since it is
#' built from the per-fit table.
#'
#' `"basin"` plots each fit's SSE against its `tc`, split at the gap
#' [print.lppls_fit()] reports. When a search has failed on part of the grid the
#' two clusters separate visibly, which is invisible in the point estimate.
#'
#' `"param_basin"` is the per-`tc` parameter plot with the same split applied to
#' every parameter, not just the objective. It answers what `"basin"` cannot:
#' *which* parameters the inferior fits differ in -- typically `omega` jumping
#' between local optima while `m` stays pinned at its bound.
#'
#' @param x An object of class `"lppls_fit"` from [fit_lppls()].
#' @param which Which plot to return. `"basin"` is drawn on demand; the rest are
#'   returned from the fit.
#' @param ... Ignored.
#'
#' @return A `ggplot2` object (a `plotly` object for `"contour"`/`"surface"`).
#'
#' @seealso [fit_lppls()], [print.lppls_fit()]
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal
#' @export
plot.lppls_fit <- function(x, which = c("basin", "param_basin", "fit", "mpl",
                                        "param", "matrix", "trace_mo",
                                        "trace_bm", "trace_bo",
                                        "contour", "surface"), ...) {
  which <- match.arg(which)

  if (!which %in% c("basin", "param_basin")) {
    slot <- c(fit = "fit_plot", mpl = "mpl_plot", param = "param_plot",
              matrix = "matrix_plot", trace_mo = "trace_plot_mo",
              trace_bm = "trace_plot_bm", trace_bo = "trace_plot_bo",
              contour = "contour_plot", surface = "surface_plot")[[which]]
    p <- x[[slot]]
    if (is.null(p)) {
      stop("This fit has no ", which, " plot; re-run fit_lppls() asking for it.")
    }
    return(p)
  }

  tbl <- if (length(x$fit) >= 2L && is.data.frame(x$fit[[2]])) x$fit[[2]] else NULL
  if (is.null(tbl)) stop("This fit has no per-fit table to draw a basin plot from.")
  d <- .fit_diagnostics(x)

  tbl$basin <- if (isTRUE(d$split)) {
    ifelse(tbl$value <= d$split_at[1], "good", "inferior")
  } else {
    "single"
  }
  sub <- if (isTRUE(d$split)) {
    sprintf("largest gap %.2f of the spread; %d of %d fits in the inferior basin",
            d$basin_gap, d$n_inferior, nrow(tbl))
  } else {
    "no distinct inferior cluster"
  }

  pal <- c(good = "royalblue", inferior = "firebrick", single = "black")

  if (which == "basin") {
    return(
      ggplot2::ggplot(tbl, ggplot2::aes(x = tc, y = value, color = basin)) +
        ggplot2::geom_point(size = 0.6) +
        ggplot2::scale_color_manual(name = NULL, values = pal) +
        ggplot2::labs(x = "Critical Time (tc)", y = "SSE of the conditional fit",
                      title = "Optimisation basins", subtitle = sub) +
        ggplot2::theme_minimal() +
        overlay_legend_theme(inside = c(0.99, 0.99), justification = c(1, 1))
    )
  }

  ## param_basin: the per-tc parameter plot, split the same way, so the basins
  ## can be read off every parameter rather than only the objective.
  plot_df <- tidyr::gather(tbl, key = "param", value = "estimate",
                           -c(ID, tc, basin))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = tc, y = estimate, color = basin)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = x$fit[[1]]$tc, color = "red") +
    ggplot2::scale_color_manual(name = NULL, values = pal) +
    ggplot2::facet_wrap(~param, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Critical Time (tc)", y = NULL,
                  title = "Per-tc estimates by optimisation basin", subtitle = sub)
}


#' Curvature of the LPPLS objective at the estimate
#'
#' Tests whether the calibration sits at a genuine interior minimum. The
#' modified profile likelihood assumes it does: its correction term uses
#' \eqn{\det(X'X - H)}{det(X'X - H)}, the observed information, which is only guaranteed
#' positive at an interior stationary point.
#'
#' By the Schur complement, \eqn{\det(X'X - H)}{det(X'X - H)} factorises into the determinant
#' of the linear-parameter block (always positive) times the determinant of the
#' 2x2 Hessian of the \eqn{(m, \omega)} profile. So the sign of the full
#' determinant is decided entirely by the curvature in the two nonlinear
#' parameters, and that 2x2 block is what this function reports. A saddle there
#' is the signature of an estimate resting on an optimisation bound.
#'
#' @param x An object of class `"lppls_fit"` from [fit_lppls()].
#' @param log_price Optional numeric vector, the series the fit was calibrated
#'   on. Only needed for fits made before `fit_lppls()` stored it.
#'
#' @return An object of class `"lppls_curvature"`: a list with the full and
#'   profile determinants, the profile Hessian and its eigenvalues, the active
#'   bounds, and a `classification` of `"minimum"`, `"saddle"` or `"maximum"`.
#'
#' @seealso [fit_lppls()], [print.lppls_fit()]
#'
#' @export
lppls_curvature <- function(x, log_price = NULL) {
  if (!inherits(x, "lppls_fit")) {
    stop("'x' must be an object of class \"lppls_fit\" from fit_lppls().")
  }
  lp <- if (!is.null(log_price)) log_price else x$log_price
  if (is.null(lp)) {
    stop("The fit does not carry its input series (made before fit_lppls() ",
         "stored it); supply it with the 'log_price' argument.")
  }
  a <- x$fit_args
  n <- if (is.null(a$n)) min(x$fit[[2]]$tc) - 1 else a$n
  if (length(lp) < n) stop("'log_price' is shorter than the calibration window.")

  e     <- x$fit[[1]]
  t     <- seq_len(n)
  log_p <- lp[seq_len(n)]
  Psi   <- c(e$m, e$omega, e$A, e$B, e$C1, e$C2)

  X <- compute_X_matrix(Psi, e$tc, t)
  H <- compute_H_matrix(Psi, e$tc, log_p, t)
  M <- crossprod(X) - H

  ## Schur complement onto (m, omega): the Hessian of the profiled objective.
  lin <- M[3:6, 3:6, drop = FALSE]
  S   <- M[1:2, 1:2, drop = FALSE] -
         M[1:2, 3:6, drop = FALSE] %*% solve(lin) %*% M[3:6, 1:2, drop = FALSE]
  ev  <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values

  ## Curvature only. Positive definite does not imply the point is stationary:
  ## with a bound active the gradient need not vanish, and the two questions are
  ## reported separately.
  cls <- if (all(ev > 0)) {
    "positive definite"
  } else if (all(ev < 0)) {
    "negative definite"
  } else {
    "saddle"
  }
  pinned <- if (!is.null(a$lower) && !is.null(a$upper)) {
    .boundary_pars(e, a$lower, a$upper)
  } else {
    NULL
  }

  structure(
    list(
      det_full        = det(M),
      det_linear      = det(lin),
      det_profile     = det(S),
      hessian_profile = S,
      eigen_profile   = ev,
      classification  = cls,
      pinned          = pinned,
      tc              = e$tc,
      n               = n
    ),
    class = "lppls_curvature"
  )
}


#' @export
print.lppls_curvature <- function(x, ...) {
  cat("LPPLS curvature at the estimate\n")
  cat(sprintf("  tc = %.2f   n = %d\n", x$tc, x$n))
  cat(sprintf("  det(X'X - H)          = %+.4e   (MPL is defined only where this is > 0)\n",
              x$det_full))
  cat(sprintf("  det(linear block)     = %+.4e   (positive by construction)\n", x$det_linear))
  cat(sprintf("  det(profile Hessian)  = %+.4e   (sets the sign above)\n", x$det_profile))
  cat(sprintf("  eigenvalues (m, omega): %s\n",
              paste(sprintf("%+.4e", x$eigen_profile), collapse = "   ")))
  cat(sprintf("  => profile Hessian is %s%s\n", x$classification,
              switch(x$classification,
                     `positive definite` = ": the curvature condition behind the MPL holds",
                     `saddle`            = ": not a minimum in (m, omega), so det(X'X - H) < 0",
                     `negative definite` = ": a maximum in (m, omega)")))
  if (length(x$pinned)) {
    cat(sprintf("  note: %s sits at an optimisation bound, so the estimate is a constrained one -\n",
                paste(x$pinned, collapse = " and ")))
    cat("        the gradient need not vanish there, whatever the curvature says\n")
  }
  invisible(x)
}


#' Full diagnostic summary of an LPPLS calibration
#'
#' Everything [print.lppls_fit()] reports, plus the curvature check at the
#' estimate when the fit carries (or is given) its input series. The returned
#' object holds the diagnostics as values, so they can be used programmatically
#' rather than parsed out of printed text.
#'
#' @param object An object of class `"lppls_fit"` from [fit_lppls()].
#' @param log_price Optional numeric vector passed to [lppls_curvature()], for
#'   fits made before `fit_lppls()` stored the series.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.lppls_fit"`.
#'
#' @seealso [fit_lppls()], [lppls_curvature()]
#'
#' @export
summary.lppls_fit <- function(object, log_price = NULL, ...) {
  d <- .fit_diagnostics(object)
  d$curvature <- tryCatch(lppls_curvature(object, log_price = log_price),
                          error = function(e) NULL)
  structure(d, class = "summary.lppls_fit")
}


#' @export
print.summary.lppls_fit <- function(x, ...) {
  d <- x
  class(d) <- NULL
  .print_fit_diagnostics(d, curv = d$curvature)
  if (!is.null(d$mpl$li_check)) {
    cat("\n  Likelihood intervals:\n")
    for (k in d$mpl$li_check) {
      cat(sprintf("    LI(%.2f): %d grid points above the cutoff, interval spans %s%s\n",
                  k$cutoff, k$n_in,
                  if (is.na(k$span)) "-" else as.character(k$span),
                  if (isTRUE(k$contiguous)) "  (contiguous)"
                  else if (is.na(k$contiguous)) ""
                  else sprintf("  (bridges %d points below the cutoff)", k$bridged)))
    }
  }
  invisible(x)
}
