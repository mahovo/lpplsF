#' Log-Periodic Power Law Singularity (LPPLS) Model
#'
#' Computes the LPPLS model value for given parameters and time points.
#' The model describes the log-price dynamics of financial assets near
#' critical points (bubbles/crashes).
#'
#' @param t Numeric vector of time indices.
#' @param A Numeric scalar, the log-price at critical time.
#' @param B Numeric scalar, the amplitude of power law growth (typically
#'   negative for bubbles).
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
#' @param omega2 Numeric scalar, secondary frequency for mode 2 and 3 (default
#'   0).
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

  ## Input validation
  if (!is.numeric(t)) {
    stop("'t' must be a numeric vector")
  }
  if (any(t >= tc)) {
    warning("Some time values are >= tc; these will produce NaN values")
  }

  if (mode == 0) {
    ## Filimonov formulation
    d <- omega * log(tc - t)
    A + (tc - t)^m * (B + C1 * cos(d) + C2 * sin(d))

  } else if (mode == 1) {
    ## First order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale
    ## Invariance), equation (19).
    ## See also Johansen et al. (2000), equation (18).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## C2 instead of phi
    ## m instead of beta
    A + (tc - t)^m * (B + C1 * cos(omega * log(tc - t) + C2))

  } else if (mode == 2) {
    ## Second order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale
    ## Invariance), equation (20).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## A instead of A_2
    ## B instead of B_2
    ## C1 instead of C_2
    ## C2 instead of phi_2
    ## m instead of beta
    ## omega instead of omega_1
    tau <- tc - t
    A + tau^m / sqrt(1 + (tau / T1)^(2 * m)) *
      (B + C1 * cos(omega * log(tau) +
                      (omega2 / (2 * m)) * log(1 + (tau / T1)^(2 * m)) + C2))

  } else if (mode == 3) {
    ## Third order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale
    ## Invariance), equation (22).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## A instead of A_3
    ## B instead of B_3
    ## C1 instead of C_3
    ## C2 instead of phi_3
    ## m instead of beta
    ## omega instead of omega_1
    tau <- tc - t
    A + (
      tau^m / sqrt(1 + (tau / T1)^(2 * m) + (tau / T2)^(4 * m)) *
        (B + C1 * cos(
          omega * log(tau) + (omega2 / (2 * m)) *
            log(1 + (tau / T1)^(2 * m)) +
            (omega3 / (4 * m)) *
            log(1 + (tau / T2)^(4 * m)) +
            C2)))

  } else {
    stop("LPPLS mode must be 0, 1, 2 or 3.")
  }
}


#' LPPLS parameter names, grouped by whether eval_lppls() defaults them
#'
#' `core` are the seven parameters eval_lppls() has no default for (so a caller
#' must supply them); `higher` are the four higher-order parameters eval_lppls()
#' defaults (T1 = 500, T2 = 1990, omega2 = 0, omega3 = 0) and that only modes 2-3
#' use. Consumed by SSE() and contour_plot_sse() to validate `par` and to forward
#' the higher-order parameters when present.
#'
#' @keywords internal
#' @noRd
.lppls_pars <- list(
  core   = c("A", "B", "C1", "C2", "tc", "m", "omega"),
  higher = c("T1", "T2", "omega2", "omega3")
)


#' Higher-order LPPLS parameters used by a given mode
#'
#' Returns the higher-order parameter names eval_lppls() actually uses for
#' `mode`: none for modes 0-1, `c("T1", "omega2")` for mode 2, and
#' `c("T1", "T2", "omega2", "omega3")` for mode 3. Drives the "silent default"
#' warning in SSE().
#'
#' @param mode Integer, the LPPLS model mode (0-3).
#' @return Character vector of parameter names (empty for modes 0-1).
#'
#' @keywords internal
#' @noRd
.lppls_mode_higher <- function(mode) {
  switch(as.character(mode),
    "2" = c("T1", "omega2"),
    "3" = c("T1", "T2", "omega2", "omega3"),
    character(0))
}


#' Sum of squared LPPLS errors (compute core, unvalidated)
#'
#' The unchecked engine behind SSE(): evaluates the LPPLS model at the
#' parameters carried in `par` and returns the residual sum of squares.
#' Every eval_lppls() parameter present in `par` -- including T1, T2,
#' omega2, omega3 -- is forwarded by name; any absent one falls back
#' to eval_lppls()'s default. Called once per grid cell by contour_plot_sse()
#' (which validates up front) and by the validating wrapper SSE(), so it
#' performs no input checking itself.
#'
#' @param par Named list of LPPLS parameters (see SSE()).
#' @param log_p Numeric vector of observed log-prices.
#' @param t Numeric vector of time indices, same length as `log_p`.
#' @param mode Integer, the LPPLS model mode (0-3).
#' @return Numeric scalar, the residual sum of squares.
#'
#' @keywords internal
#' @noRd
.sse <- function(par, log_p, t, mode) {
  keep <- intersect(names(par), c(.lppls_pars$core, .lppls_pars$higher))
  drop(sum(
    (log_p - do.call(
      eval_lppls,
      c(list(t = t, mode = mode), par[keep])
    )
    )^2,
     na.rm = TRUE)
  )
}


#' Sum of Squared Errors for the LPPLS Model
#'
#' Computes the residual sum of squares between observed log-prices and the
#' LPPLS model evaluated at a given parameter set. This is the objective the
#' calibration minimizes and the quantity [contour_plot_sse()] scans over a
#' parameter grid.
#'
#' @details
#' `par` is a named list of LPPLS parameters. The seven core parameters `A`,
#' `B`,
#' `C1`, `C2`, `tc`, `m`, `omega` are required. For the higher-order variants
#' (`mode = 2` or `3`), `par` may also carry the scale and secondary-frequency
#' parameters `T1`, `T2`, `omega2`, `omega3`; each is forwarded to
#' [eval_lppls()] when present and otherwise falls back to that function's
#' default (`T1 = 500`, `T2 = 1990`, `omega2 = 0`, `omega3 = 0`). Requesting
#' a higher-order mode without supplying the parameters it uses issues a
#' warning and falls back to those defaults, since silently defaulting them
#' can give a misleading SSE.
#'
#' Elements of `par` not recognised by [eval_lppls()] are ignored. Time points
#' at or beyond `tc` produce `NaN` from [eval_lppls()] and are dropped from the
#' sum (`na.rm = TRUE`).
#'
#' @param par Named list of LPPLS parameters. Must contain the core parameters
#'   `A`, `B`, `C1`, `C2`, `tc`, `m`, `omega`; may additionally contain `T1`,
#'   `T2`, `omega2`, `omega3` for `mode = 2` or `3`.
#' @param log_p Numeric vector of observed log-prices.
#' @param t Numeric vector of time indices, the same length as `log_p`.
#' @param mode Integer, the LPPLS model variant (`0`-`3`; default `0`).
#'
#' @return Numeric scalar, the sum of squared errors.
#'
#' @seealso [eval_lppls()], [contour_plot_sse()]
#'
#' @examples
#' t <- 1:100
#' par <- list(A = 4, B = -0.015, C1 = 0.0015, C2 = 0.001,
#'             tc = 150, m = 0.5, omega = 9)
#' log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
#'                     par$tc, par$m, par$omega) + rnorm(100, 0, 0.01)
#' SSE(par, log_p, t)
#'
#' # Higher-order variant (mode 2): supply T1 and omega2 alongside the core set
#' par2 <- c(par, list(T1 = 300, omega2 = 5))
#' SSE(par2, log_p, t, mode = 2)
#'
#' @export
SSE <- function(par, log_p, t, mode = 0) {
  if (!is.list(par) || is.null(names(par))) {
    stop("'par' must be a named list")
  }
  if (!is.numeric(log_p) || !is.numeric(t)) {
    stop("'log_p' and 't' must be numeric")
  }
  if (length(log_p) != length(t)) {
    stop("'log_p' and 't' must have the same length")
  }
  if (!isTRUE(mode %in% 0:3)) {
    stop("'mode' must be 0, 1, 2 or 3")
  }
  miss <- setdiff(.lppls_pars$core, names(par))
  if (length(miss)) {
    stop(
      "'par' is missing required parameter(s): ",
      paste(miss, collapse = ", ")
    )
  }
  absent <- setdiff(.lppls_mode_higher(mode), names(par))
  if (length(absent))
    warning(
      "mode ", mode, " uses ", paste(.lppls_mode_higher(mode), collapse = ", "),
      "; '", paste(absent, collapse = "', '"),
      "' not supplied in 'par' \u2014 using eval_lppls() default(s)"
    )

  .sse(par, log_p, t, mode)
}
