#' @keywords internal
#' @noRd
#'
#' Helper functions for LPPLS model fitting.
#' These are internal functions not exported to users.

# Helper function: f(t, tc, m) = (tc - t)^m
lppls_f <- function(t, tc, m) {
  (tc - t)^m
}

# Helper function: g(t, tc, m, omega) = (tc - t)^m * cos(omega * log(tc - t))
lppls_g <- function(t, tc, m, omega) {
  (tc - t)^m * cos(omega * log(tc - t))
}

# Helper function: h(t, tc, m, omega) = (tc - t)^m * sin(omega * log(tc - t))
lppls_h <- function(t, tc, m, omega) {
  (tc - t)^m * sin(omega * log(tc - t))
}

# Sum functions for matrix construction
sum_f <- function(t, tc, m) {
  sum(lppls_f(t, tc, m))
}

sum_g <- function(t, tc, m, omega) {
  sum(lppls_g(t, tc, m, omega))
}

sum_h <- function(t, tc, m, omega) {
  sum(lppls_h(t, tc, m, omega))
}

sum_ff <- function(t, tc, m) {
  f_vals <- lppls_f(t, tc, m)
  sum(f_vals * f_vals)
}

sum_gg <- function(t, tc, m, omega) {
  g_vals <- lppls_g(t, tc, m, omega)
  sum(g_vals * g_vals)
}

sum_hh <- function(t, tc, m, omega) {
  h_vals <- lppls_h(t, tc, m, omega)
  sum(h_vals * h_vals)
}

sum_fg <- function(t, tc, m, omega) {
  sum(lppls_f(t, tc, m) * lppls_g(t, tc, m, omega))
}

sum_fh <- function(t, tc, m, omega) {
  sum(lppls_f(t, tc, m) * lppls_h(t, tc, m, omega))
}

sum_gh <- function(t, tc, m, omega) {
  sum(lppls_g(t, tc, m, omega) * lppls_h(t, tc, m, omega))
}

sum_yf <- function(log_p, t, tc, m) {
  sum(log_p * lppls_f(t, tc, m))
}

sum_yg <- function(log_p, t, tc, m, omega) {
  sum(log_p * lppls_g(t, tc, m, omega))
}

sum_yh <- function(log_p, t, tc, m, omega) {
  sum(log_p * lppls_h(t, tc, m, omega))
}


#' Calculate Damping Coefficient
#'
#' Computes the damping coefficient D which measures the strength of
#' log-periodic oscillations relative to the power law growth.
#'
#' @param m Numeric scalar, power law exponent.
#' @param B Numeric scalar, amplitude of power law growth.
#' @param omega Numeric scalar, log-periodic frequency.
#' @param C1 Numeric scalar, amplitude of cosine oscillation.
#' @param C2 Numeric scalar, amplitude of sine oscillation.
#'
#' @return Numeric scalar, the damping coefficient.
#'
#' @details
#' The damping coefficient is defined as:
#' \deqn{D = \frac{m |B|}{\omega \sqrt{C_1^2 + C_2^2}}}
#'
#' Higher values indicate stronger damping (less pronounced oscillations).
#' Filimonov (2017) suggests D >= 0.8 for valid bubble detection.
#'
#' @keywords internal
#' @noRd
calculate_damping <- function(m, B, omega, C1, C2) {
  m * abs(B) / (omega * sqrt(C1^2 + C2^2))
}


#' Create Beta Calculator Using Symbolic Math
#'
#' Creates a function that calculates the linear parameters (A, B, C1, C2)
#' given fixed nonlinear parameters (tc, m, omega) using symbolic matrix
#' inversion for efficiency.
#'
#' @return A function that takes (log_p, t, tc, m, omega) and returns
#'   a numeric vector c(A, B, C1, C2).
#'
#' @importFrom symengine Symbol Matrix Vector solve DoubleVisitor
#' @keywords internal
#' @noRd
create_beta_calculator <- function() {
  ## Symbolic variables for X'X matrix elements
  xx11 <- symengine::Symbol("xx11")
  xx12 <- symengine::Symbol("xx12")
  xx13 <- symengine::Symbol("xx13")
  xx14 <- symengine::Symbol("xx14")
  xx21 <- symengine::Symbol("xx21")
  xx22 <- symengine::Symbol("xx22")
  xx23 <- symengine::Symbol("xx23")
  xx24 <- symengine::Symbol("xx24")
  xx31 <- symengine::Symbol("xx31")
  xx32 <- symengine::Symbol("xx32")
  xx33 <- symengine::Symbol("xx33")
  xx34 <- symengine::Symbol("xx34")
  xx41 <- symengine::Symbol("xx41")
  xx42 <- symengine::Symbol("xx42")
  xx43 <- symengine::Symbol("xx43")
  xx44 <- symengine::Symbol("xx44")
  xy1 <- symengine::Symbol("xy1")
  xy2 <- symengine::Symbol("xy2")
  xy3 <- symengine::Symbol("xy3")
  xy4 <- symengine::Symbol("xy4")

  ## Build X'X matrix
  XX <- symengine::Matrix(
    c(xx11, xx12, xx13, xx14,
      xx21, xx22, xx23, xx24,
      xx31, xx32, xx33, xx34,
      xx41, xx42, xx43, xx44),
    byrow = TRUE,
    nrow = 4
  )

  ## Build X'y vector
  Xy <- symengine::Vector(xy1, xy2, xy3, xy4)

  ## Solve for beta = (X'X)^{-1} X'y
  beta <- symengine::solve(XX, Xy)

  ## Create efficient visitor functions for each parameter
  args <- c(xx11, xx12, xx13, xx14,
            xx21, xx22, xx23, xx24,
            xx31, xx32, xx33, xx34,
            xx41, xx42, xx43, xx44,
            xy1, xy2, xy3, xy4)

  A_func <- symengine::DoubleVisitor(beta[1], args = args)
  B_func <- symengine::DoubleVisitor(beta[2], args = args)
  C1_func <- symengine::DoubleVisitor(beta[3], args = args)
  C2_func <- symengine::DoubleVisitor(beta[4], args = args)

  ## Return the calculator function
  ## Calculates values of liniear parameters, given nonlinear estimates as inputs.
  function(log_p, t, tc, m, omega) {
    force(tc)
    force(m)
    force(omega)

    n <- length(t)

    ## Compute all matrix elements
    vals <- list(
      xx11 = n,
      xx12 = sum_f(t, tc, m),
      xx13 = sum_g(t, tc, m, omega),
      xx14 = sum_h(t, tc, m, omega),
      xx21 = sum_f(t, tc, m),
      xx22 = sum_ff(t, tc, m),
      xx23 = sum_fg(t, tc, m, omega),
      xx24 = sum_fh(t, tc, m, omega),
      xx31 = sum_g(t, tc, m, omega),
      xx32 = sum_fg(t, tc, m, omega),
      xx33 = sum_gg(t, tc, m, omega),
      xx34 = sum_gh(t, tc, m, omega),
      xx41 = sum_h(t, tc, m, omega),
      xx42 = sum_fh(t, tc, m, omega),
      xx43 = sum_gh(t, tc, m, omega),
      xx44 = sum_hh(t, tc, m, omega),
      xy1 = sum(log_p),
      xy2 = sum_yf(log_p, t, tc, m),
      xy3 = sum_yg(log_p, t, tc, m, omega),
      xy4 = sum_yh(log_p, t, tc, m, omega)
    )

    ## Calculate coefficients
    a <- drop(do.call(A_func, vals))
    b <- drop(do.call(B_func, vals))
    c1 <- drop(do.call(C1_func, vals))
    c2 <- drop(do.call(C2_func, vals))

    c(A = a, B = b, C1 = c1, C2 = c2)
  }
}


#' Validate LPPLS Parameters
#'
#' Checks if estimated parameters fall within acceptable ranges.
#'
#' @param params Named list with elements: tc, m, omega, A, B, C1, C2, D
#' @param lower Numeric vector of lower bounds for c(m, omega, B, D)
#' @param upper Numeric vector of upper bounds for c(m, omega, B, D)
#'
#' @return List with elements:
#'   \itemize{
#'     \item valid: logical, TRUE if all parameters in range
#'     \item B_valid: logical, TRUE if B < upper[3]
#'     \item D_valid: logical, TRUE if D >= lower[4]
#'     \item m_valid: logical, TRUE if lower[1] <= m <= upper[1]
#'     \item omega_valid: logical, TRUE if lower[2] <= omega <= upper[2]
#'   }
#'
#' @keywords internal
#' @noRd
validate_params <- function(params, lower, upper) {
  B_valid <- params$B < upper[3]
  D_valid <- params$D >= lower[4]
  m_valid <- params$m >= lower[1] && params$m <= upper[1]
  omega_valid <- params$omega >= lower[2] && params$omega <= upper[2]

  list(
    valid = B_valid && D_valid && m_valid && omega_valid,
    B_valid = B_valid,
    D_valid = D_valid,
    m_valid = m_valid,
    omega_valid = omega_valid
  )
}
