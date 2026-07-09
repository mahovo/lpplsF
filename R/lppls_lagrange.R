#' Lagrange Regularization of a Rolling-Window LPPLS Calibration
#'
#' Determines the optimal calibration start time \eqn{T_1^*} by Lagrange
#' regularization of the normalized residual sum of squares, following
#' Demos & Sornette (2019). It is a post-processing step on a rolling-window
#' calibration ([lppls_rolling()]): the per-window residual sums of squares are
#' corrected for the linear bias induced by the shrinking number of
#' observations, and \eqn{T_1^*} is the minimizer of the regularized objective
#' \eqn{\chi^2_{\lambda^*}(T_1)}.
#'
#' @details
#' For the rolling calibration `x`, the normalized residual sum of squares is
#' \deqn{\chi^2_{np}(T_1) = \frac{1}{n_{cal}(T_1) - p} \sum_{i = T_1}^{T_2} r_i^2,}
#' where the sum is the calibration-region residual sum of squares
#' (`x$table$sse_cal`), \eqn{n_{cal}(T_1)} is the window's calibration length and
#' \eqn{p} is the number of parameters subtracted. `p = 0` reproduces the
#' \dQuote{p-correction omitted} form used by Demos & Sornette (2019) and in the
#' thesis; `p = 7` (the default) subtracts the full LPPLS parameter count
#' \eqn{\Psi = (t_c, m, \omega, A, B, C_1, C_2)}.
#'
#' The bias is modelled as a linear trend
#' \eqn{\chi^2_{np}(T_1) = \beta_0 + \beta_1 T_1} fitted by ordinary least
#' squares. The regularized objective is
#' \deqn{\chi^2_{\lambda^*}(T_1) = \chi^2_{np}(T_1) - \lambda^* (T_2 - T_1),
#'   \qquad \lambda^* = -\beta_1,}
#' which has zero trend in \eqn{T_1}. The optimal start time is
#' \eqn{T_1^* = \arg\min_{T_1} \chi^2_{\lambda^*}(T_1)}. Adding a constant to the
#' regularized objective does not change its minimizer, so the choice of
#' \eqn{T_2} in the \eqn{-\lambda^*(T_2 - T_1)} term shifts the curve vertically
#' but leaves \eqn{T_1^*} unchanged.
#'
#' @param x An object of class `"lppls_rolling"` from [lppls_rolling()].
#' @param p Integer, number of parameters subtracted in the \eqn{\chi^2_{np}}
#'   denominator (default 7). Use `p = 0` to omit the correction (thesis /
#'   Demos & Sornette (2019) form).
#'
#' @return An object of class `"lppls_lagrange"`: a list with
#'   \describe{
#'     \item{table}{Data frame with columns `T1`, `n` (calibration length),
#'       `chi_sq_np` and `chi_sq_lambda`.}
#'     \item{rolling}{The input [lppls_rolling()] object.}
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
#' @seealso [lppls_rolling()], [chi_sq_np_plot()], [chi_sq_lambda_plot()]
#'
#' @examples
#' \dontrun{
#' set.seed(3630)
#' N <- 2000
#' log_p <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
#'                     tc = N + 1, m = 0.6, omega = 8) + rnorm(N) * 0.05
#' roll <- lppls_rolling(
#'   log_price = log_p,
#'   t1_seq    = seq(from = 1, by = 10, length.out = 101),
#'   hold_out  = 100, fh = 150,
#'   fit_args  = list(lower = c(0.1, 6, -0.03, 0.8),
#'                    upper = c(0.9, 13, -1e-14, 1e6),
#'                    m_init = 0.2, o_init = 12, num_searches = 1)
#' )
#' lag <- lppls_lagrange(roll, p = 0)
#' lag$t1_star
#' print(chi_sq_np_plot(lag))
#' print(chi_sq_lambda_plot(lag))
#' }
#'
#' @importFrom stats lm coef
#' @export
lppls_lagrange <- function(x, p = 7) {

  if (!inherits(x, "lppls_rolling")) {
    stop("'x' must be an object of class 'lppls_rolling' (see lppls_rolling())")
  }

  tab <- x$table
  denom <- tab$n_cal - p
  if (any(denom <= 0)) {
    stop("Non-positive denominator (n_cal - p) for some window; reduce p.")
  }

  chi_np <- tab$sse_cal / denom
  T1 <- tab$T1
  T2 <- x$T2

  ## Linear bias trend and regularization
  trend <- stats::lm(chi_np ~ T1)
  beta0 <- unname(stats::coef(trend)[[1]])
  beta1 <- unname(stats::coef(trend)[[2]])
  lambda <- -beta1                                  ## lambda* = -beta1

  chi_lambda <- chi_np - lambda * (T2 - T1)         ## chi^2_lambda*(T1)
  t1_star <- T1[which.min(chi_lambda)]

  structure(
    list(
      table = data.frame(
        T1 = T1,
        n = tab$n_cal,
        chi_sq_np = chi_np,
        chi_sq_lambda = chi_lambda
      ),
      rolling = x,
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
#' @seealso [lppls_lagrange()], [lppls_rolling()], [chi_sq_lambda_plot()]
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
#' @seealso [lppls_lagrange()], [lppls_rolling()], [chi_sq_np_plot()]
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
