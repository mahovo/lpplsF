#' Log-Periodic Power Law Singularity (LPPLS) Model
#'
#' Computes the LPPLS model value for given parameters and time points.
#' The model describes the log-price dynamics of financial assets near
#' critical points (bubbles/crashes).
#'
#' @param t Numeric vector of time indices.
#' @param A Numeric scalar, the log-price at critical time.
#' @param B Numeric scalar, the amplitude of power law growth (typically negative for bubbles).
#' @param C1 Numeric scalar, amplitude of cosine oscillation.
#' @param C2 Numeric scalar, amplitude of sine oscillation.
#' @param tc Numeric scalar, the critical time (singularity point).
#' @param m Numeric scalar, the power law exponent (typically 0.1 to 0.9).
#' @param omega Numeric scalar, the log-periodic frequency (typically 6 to 13).
#' @param mode Integer, the model variant:
#'   \itemize{
#'     \item 0: Filimonov formulation (default)
#'     \item 1: First order expansion (Johansen 1999, eq. 19)
#'     \item 2: Second order expansion (Johansen 1999, eq. 20)
#'     \item 3: Third order expansion (Johansen 1999, eq. 22)
#'   }
#' @param T1 Numeric scalar, scale parameter for mode 2 and 3 (default 500).
#' @param T2 Numeric scalar, scale parameter for mode 3 (default 1990).
#' @param omega2 Numeric scalar, secondary frequency for mode 2 and 3 (default 0).
#' @param omega3 Numeric scalar, tertiary frequency for mode 3 (default 0).
#'
#' @return Numeric vector of LPPLS model values at each time point.
#'
#' @details
#' The standard LPPLS model (mode = 0) is:
#' \deqn{LPPLS(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c - t)) + C_2 \sin(\omega \log(t_c - t))]}
#'
#' This model captures both the power law growth towards a critical time and the

#' log-periodic oscillations that are characteristic of speculative bubbles.
#'
#' @references
#' Filimonov, V., & Sornette, D. (2013). A stable and robust calibration scheme
#' of the log-periodic power law model. Physica A: Statistical Mechanics and
#' its Applications.
#'
#' Johansen, A. (1999). Predicting Financial Crashes Using Discrete Scale
#' Invariance. Journal of Risk.
#'
#' @examples
#' # Generate time series
#' t <- 1:100
#'
#' # Generate LPPLS curve with known parameters
#' params <- list(A = 4, B = -0.015, C1 = 0.0015, C2 = 0.001,
#'                tc = 150, m = 0.5, omega = 9)
#' y <- eval_lppls(t, params$A, params$B, params$C1, params$C2,
#'                params$tc, params$m, params$omega)
#'
#' # Plot the result
#' plot(t, y, type = "l", main = "LPPLS Model")
#'
#' @export
eval_lppls <- function(
  t,
  A,
  B,
  C1,
  C2,
  tc,
  m,
  omega,
  mode = 0,
  T1 = 500,
  T2 = 1990,
  omega2 = 0,
  omega3 = 0
) {

  # Input validation
  if (!is.numeric(t)) {
    stop("'t' must be a numeric vector")
  }
  if (any(t >= tc)) {
    warning("Some time values are >= tc; these will produce NaN values")
  }

  if (mode == 0) {
    # Filimonov formulation
    d <- omega * log(tc - t)
    A + (tc - t)^m * (B + C1 * cos(d) + C2 * sin(d))

  } else if (mode == 1) {
    # First order expansion
    # Johansen 1999, equation (19)
    # C2 repurposed as phi
    A + (tc - t)^m * (B + C1 * cos(omega * log(tc - t) + C2))

  } else if (mode == 2) {
    # Second order expansion
    # Johansen 1999, equation (20)
    tau <- tc - t
    A + tau^m / sqrt(1 + (tau / T1)^(2 * m)) *
      (B + C1 * cos(omega * log(tau) +
                      (omega2 / (2 * m)) * log(1 + (tau / T1)^(2 * m)) + C2))

  } else if (mode == 3) {
    # Third order expansion
    # Johansen 1999, equation (22)
    tau <- tc - t
    A + (
      tau^m / sqrt(1 + (tau / T1)^(2 * m) + (tau / T2)^(4 * m)) *
        (B + C1 * cos(
          omega + log(tau) + (omega2 / (2 * m)) *
            log(1 + (tau / T1)^(2 * m)) +
            (omega3 / (4 * m)) *
            log(1 + (tau / T2)^(4 * m)) +
            C2)))

  } else {
    stop("LPPLS mode must be 0, 1, 2 or 3.")
  }
}


#' Sum of Squared Errors for LPPLS Model
#'
#' Computes the sum of squared residuals between observed log-prices and
#' the LPPLS model predictions.
#'
#' @param par A named list with elements: tc, m, omega, A, B, C1, C2
#' @param log_p Numeric vector of observed log-prices.
#' @param t Numeric vector of time indices.
#' @param mode Integer, the LPPLS model variant (default 0).
#'
#' @return Numeric scalar, the sum of squared errors.
#'
#' @examples
#' # Create parameter list
#' par <- list(A = 4, B = -0.015, C1 = 0.0015, C2 = 0.001,
#'             tc = 150, m = 0.5, omega = 9)
#'
#' # Generate synthetic data with noise
#' t <- 1:100
#' log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
#'                    par$tc, par$m, par$omega) + rnorm(100, 0, 0.01)
#'
#' # Calculate SSE
#' sse_value <- SSE(par, log_p, t)
#'
#' @export
SSE <- function(
  par,
  log_p,
  t,
  mode = 0
) {
  if (length(log_p) != length(t)) {
    stop("'log_p' and 't' must have the same length")
  }

  res <- log_p - eval_lppls(t, par$A, par$B, par$C1, par$C2,
                        par$tc, par$m, par$omega, mode)
  drop(sum(res^2, na.rm = TRUE))
}
