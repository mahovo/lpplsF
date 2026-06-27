#' Lagrange Regularization of the LPPLS Objective Function
#'
#' Determines the optimal calibration start time \eqn{T_1^*} by Lagrange
#' regularization of the normalized residual sum of squares, following
#' Demos & Sornette (2019). A sequence of LPPLS models is calibrated with a
#' fixed end time \eqn{T_2} and a rolling start time \eqn{T_1}. The normalized
#' residual sum of squares \eqn{\chi^2_{np}(T_1)} is corrected for the linear
#' bias induced by the shrinking number of observations, and \eqn{T_1^*} is the
#' minimizer of the regularized objective \eqn{\chi^2_{\lambda^*}(T_1)}.
#'
#' @details
#' For each start time \eqn{T_1} in `t1_seq`, the model is calibrated on the
#' window `log_price[T1:length(log_price)]`. The last `hold_out` observations of
#' every window are reserved for the forecast/validation region, so the
#' calibration region is \eqn{[T_1, T_2]} with the common end
#' \eqn{T_2 = \texttt{length(log\_price)} - \texttt{hold\_out}}, held fixed
#' across windows.
#'
#' The normalized residual sum of squares is
#' \deqn{\chi^2_{np}(T_1) = \frac{1}{(T_2 - T_1) - p} \sum_{i = T_1}^{T_2} r_i^2,}
#' where \eqn{r_i} is the residual of calibration observation \eqn{i} and \eqn{p}
#' is the number of model parameters. Setting `p = 0` reproduces the
#' \dQuote{p-correction omitted} form used in the thesis; `p = 7` (the default)
#' subtracts the full LPPLS parameter count
#' \eqn{\Psi = (t_c, m, \omega, A, B, C_1, C_2)}.
#'
#' The bias is modelled as a linear trend
#' \eqn{\chi^2_{np}(T_1) = \beta_0 + \beta_1 T_1} fitted by ordinary least
#' squares. The regularized objective is
#' \deqn{\chi^2_{\lambda^*}(T_1) = \chi^2_{np}(T_1) - \lambda^* (T_2 - T_1),
#'   \qquad \lambda^* = -\beta_1,}
#' which has zero trend in \eqn{T_1}. The optimal start time is
#' \eqn{T_1^* = \arg\min_{T_1} \chi^2_{\lambda^*}(T_1)}.
#'
#' Adding or removing a constant from the regularized objective does not change
#' its minimizer, so the choice of \eqn{T_2} in the \eqn{-\lambda^*(T_2 - T_1)}
#' term shifts the curve vertically but leaves \eqn{T_1^*} unchanged.
#'
#' @param log_price Numeric vector of log-prices for the full series. Windows
#'   are taken as `log_price[T1:length(log_price)]` for each `T1` in `t1_seq`.
#' @param t1_seq Integer vector of rolling start indices (absolute positions
#'   into `log_price`). All entries must be `< T2 = length(log_price) -
#'   hold_out`.
#' @param hold_out Integer, number of observations reserved from the end of the
#'   full series for the forecast/validation region (default 0). Fixes the
#'   common calibration end \eqn{T_2}.
#' @param fh Integer, forecast horizon passed to [fit_lppls()] (default 30).
#' @param p Integer, number of parameters subtracted in the \eqn{\chi^2_{np}}
#'   denominator (default 7). Use `p = 0` to omit the correction (thesis form).
#' @param tc_init Numeric or `NULL`. Initial critical time for [fit_lppls()].
#'   If `NULL` (default) a per-window value `n_cal + floor(fh / 2)` is used,
#'   matching the thesis, where `n_cal` is the window's calibration length.
#' @param fit_args Named list of additional arguments forwarded to
#'   [fit_lppls()] (e.g. `lower`, `upper`, `m_init`, `o_init`, `num_searches`,
#'   `mode`, `mpl_cutoff`, `factr`). The arguments `log_price`, `fh`,
#'   `hold_out`, `tc_init` and `fb` are controlled by `lppls_lagrange()` and are
#'   ignored if supplied here. If `mode` is not supplied it defaults to `"F2"`,
#'   which yields the same point estimate as `"MPL"` for the Lagrange
#'   computation.
#' @param fb Logical, whether to print per-window feedback (default `FALSE`).
#'
#' @return An object of class `"lppls_lagrange"`: a list with
#'   \describe{
#'     \item{table}{Data frame with columns `T1`, `n` (calibration length),
#'       `chi_sq_np` and `chi_sq_lambda`.}
#'     \item{fits}{List of best-fit parameter lists (`fit[[1]]`) per window.}
#'     \item{trend}{The fitted `lm` object for `chi_sq_np ~ T1`.}
#'     \item{intercept, slope}{OLS trend coefficients \eqn{\beta_0, \beta_1}.}
#'     \item{lambda}{The regularization multiplier \eqn{\lambda^* = -\beta_1}.}
#'     \item{T2}{The common calibration end \eqn{T_2}.}
#'     \item{p}{The parameter count used in the denominator.}
#'     \item{t1_star}{The optimal start time \eqn{T_1^*}.}
#'   }
#'
#' @references
#' Demos, G., & Sornette, D. (2019). Birth or burst of financial bubbles:
#' which one is easier to diagnose? Quantitative Finance, 17(5), 657-675.
#' (arXiv:1707.07162)
#'
#' @seealso [chi_sq_np_plot()], [chi_sq_lambda_plot()], [fit_lppls()]
#'
#' @examples
#' \dontrun{
#' set.seed(3630)
#' N <- 2000
#' tc <- N + 1
#' t <- seq_len(N)
#' log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
#'                     tc = tc, m = 0.6, omega = 8) + rnorm(N) * 0.05
#'
#' lag <- lppls_lagrange(
#'   log_price = log_p,
#'   t1_seq    = seq(from = 1, by = 10, length.out = 101),
#'   hold_out  = 100,
#'   fh        = 150,
#'   p         = 0,
#'   fit_args  = list(
#'     lower = c(0.1, 6, -0.03, 0.8),
#'     upper = c(0.9, 13, -1e-14, 1e6),
#'     m_init = 0.2, o_init = 12, num_searches = 1, mode = "MPL"
#'   )
#' )
#' lag$t1_star
#' print(chi_sq_np_plot(lag))
#' print(chi_sq_lambda_plot(lag))
#' }
#'
#' @importFrom stats lm coef
#' @export
lppls_lagrange <- function(
  log_price,
  t1_seq,
  hold_out = 0,
  fh = 30,
  p = 7,
  tc_init = NULL,
  fit_args = list(),
  fb = FALSE
) {

  ## Input validation ========================================================
  if (!is.numeric(log_price)) {
    stop("'log_price' must be a numeric vector")
  }
  if (!is.numeric(t1_seq) || length(t1_seq) < 2) {
    stop("'t1_seq' must be a numeric vector of length >= 2")
  }
  if (any(t1_seq < 1) || any(t1_seq != round(t1_seq))) {
    stop("'t1_seq' must contain positive integers")
  }

  n_full <- length(log_price)
  T2 <- n_full - hold_out

  if (any(t1_seq >= T2)) {
    stop("All entries of 't1_seq' must be < T2 = length(log_price) - hold_out")
  }

  ## Strip arguments that this function controls from fit_args ===============
  reserved <- c("log_price", "fh", "hold_out", "tc_init", "fb")
  if (length(intersect(names(fit_args), reserved)) > 0) {
    warning(
      "fit_args entries ", paste(intersect(names(fit_args), reserved),
                                 collapse = ", "),
      " are controlled by lppls_lagrange() and will be ignored."
    )
    fit_args <- fit_args[setdiff(names(fit_args), reserved)]
  }
  if (is.null(fit_args$mode)) fit_args$mode <- "F2"

  ## Roll over T1 ============================================================
  K <- length(t1_seq)
  fits <- vector("list", K)
  chi_np <- numeric(K)
  n_cal_vec <- integer(K)

  for (k in seq_len(K)) {
    T1 <- t1_seq[k]
    window_lp <- log_price[T1:n_full]
    n_total <- length(window_lp)
    n_cal <- n_total - hold_out
    n_cal_vec[k] <- n_cal

    denom <- n_cal - p
    if (denom <= 0) {
      stop(sprintf(
        "Non-positive denominator (n_cal - p = %d) at T1 = %d; reduce p or t1_seq range.",
        denom, T1
      ))
    }

    tc_init_k <- if (is.null(tc_init)) n_cal + floor(fh / 2) else tc_init

    if (fb) {
      message(sprintf("Window %d/%d: T1 = %d, n_cal = %d", k, K, T1, n_cal))
    }

    res_fit <- do.call(
      fit_lppls,
      c(
        list(
          log_price = window_lp,
          fh        = fh,
          hold_out  = hold_out,
          tc_init   = tc_init_k,
          fb        = fb
        ),
        fit_args
      )
    )

    best <- res_fit$fit[[1]]
    fits[[k]] <- best

    ## chi^2_np: residual sum of squares over the calibration region [T1, T2]
    ## only, normalized by (n_cal - p). Uses the Filimonov formulation (mode 0).
    t_cal <- seq_len(n_cal)
    sse_k <- SSE(best, window_lp[t_cal], t_cal, mode = 0)
    chi_np[k] <- sse_k / denom
  }

  ## Linear bias trend and regularization ====================================
  trend <- stats::lm(chi_np ~ t1_seq)
  beta0 <- unname(stats::coef(trend)[[1]])
  beta1 <- unname(stats::coef(trend)[[2]])
  lambda <- -beta1                                  ## lambda* = -beta1

  chi_lambda <- chi_np - lambda * (T2 - t1_seq)     ## chi^2_lambda*(T1)
  t1_star <- t1_seq[which.min(chi_lambda)]

  structure(
    list(
      table = data.frame(
        T1 = t1_seq,
        n = n_cal_vec,
        chi_sq_np = chi_np,
        chi_sq_lambda = chi_lambda
      ),
      fits = fits,
      trend = trend,
      intercept = beta0,
      slope = beta1,
      lambda = lambda,
      T2 = T2,
      p = p,
      t1_star = t1_star
    ),
    class = "lppls_lagrange"
  )
}


#' @export
print.lppls_lagrange <- function(x, ...) {
  cat("LPPLS Lagrange regularization\n")
  cat(sprintf("  Windows (T1):  %d  (from %d to %d)\n",
              nrow(x$table), min(x$table$T1), max(x$table$T1)))
  cat(sprintf("  T2 (cal. end): %d\n", x$T2))
  cat(sprintf("  p (denominator correction): %d\n", x$p))
  cat(sprintf("  Trend: chi^2_np = %.6g + (%.6g) * T1\n", x$intercept, x$slope))
  cat(sprintf("  lambda* = -beta1 = %.6g\n", x$lambda))
  cat(sprintf("  Optimal T1* = %d\n", x$t1_star))
  invisible(x)
}


#' Plot the Normalized Residual Sum of Squares with Bias Trend
#'
#' Plots \eqn{\chi^2_{np}(T_1)} against the rolling start time \eqn{T_1} with the
#' fitted linear bias trend overlaid (red line). Reproduces the
#' `*_chi_sq_np_plot` figure from the thesis.
#'
#' @param x An object of class `"lppls_lagrange"` from [lppls_lagrange()].
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_lagrange()], [chi_sq_lambda_plot()]
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline xlab ylab theme_minimal
#' @export
chi_sq_np_plot <- function(x) {
  if (!inherits(x, "lppls_lagrange")) {
    stop("'x' must be an object of class 'lppls_lagrange'")
  }

  ggplot2::ggplot(x$table, ggplot2::aes(x = T1, y = chi_sq_np)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      intercept = x$intercept,
      slope = x$slope,
      linetype = "solid",
      color = "red",
      linewidth = 0.5
    ) +
    ggplot2::xlab(quote(T[1])) +
    ggplot2::ylab(quote(paste(chi[np]^2, (T[1])))) +
    ggplot2::theme_minimal()
}


#' Plot the Lagrange-Regularized Objective with the Optimal Start Time
#'
#' Plots \eqn{\chi^2_{\lambda^*}(T_1)} against the rolling start time \eqn{T_1}
#' with a vertical line (green) at the optimal start time \eqn{T_1^*}.
#' Reproduces the `*_chi_sq_lambda_plot` figure from the thesis.
#'
#' @param x An object of class `"lppls_lagrange"` from [lppls_lagrange()].
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_lagrange()], [chi_sq_np_plot()]
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline xlab ylab theme_minimal
#' @export
chi_sq_lambda_plot <- function(x) {
  if (!inherits(x, "lppls_lagrange")) {
    stop("'x' must be an object of class 'lppls_lagrange'")
  }

  ggplot2::ggplot(x$table, ggplot2::aes(x = T1, y = chi_sq_lambda)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      xintercept = x$t1_star,
      linetype = "solid",
      color = "green3",
      linewidth = 0.5
    ) +
    ggplot2::xlab(quote(T[1])) +
    ggplot2::ylab(quote(paste(chi[lambda]^2, (T[1])))) +
    ggplot2::theme_minimal()
}
