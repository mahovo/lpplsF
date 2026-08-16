#' lpplsF: Log-Periodic Power Law Singularity Model for Financial Bubbles
#'
#' Implementation of the LPPLS model for detecting and analysing financial
#' bubbles, based on the calibration methods of Filimonov & Sornette (2013)
#' and the Modified Profile Likelihood inference from Filimonov, Demos &
#' Sornette (2017).
#'
#' @section Fitting:
#' \describe{
#'   \item{\code{\link{fit_lppls}}}{Calibrate the model to a log-price series}
#'   \item{\code{\link{eval_lppls}}}{Evaluate the model at given parameters}
#'   \item{\code{\link{SSE}}}{The objective the calibration minimises}
#' }
#'
#' @section Diagnosing a calibration:
#' A calibration can look convincing as a plot and still not be worth
#' interpreting: the estimate may rest on an optimisation bound, the per-`tc`
#' searches may have settled in different local optima, and the modified profile
#' likelihood may be undefined over part of the `tc` grid. Printing the object
#' returned by \code{\link{fit_lppls}} reports on all of this, and the rest go
#' further.
#' \describe{
#'   \item{\code{\link[=print.lppls_fit]{print}}}{The estimate together with the
#'     checks that decide whether to trust it; warnings are listed only when
#'     something is wrong, so a clean calibration reports none}
#'   \item{\code{\link[=summary.lppls_fit]{summary}}}{The same diagnostics as
#'     values, plus the curvature check and whether each likelihood interval is
#'     closed by the likelihood or by the edge of the grid}
#'   \item{\code{\link[=plot.lppls_fit]{plot}}}{Any plot attached at fit time, or
#'     the basin views, which show whether the per-`tc` searches split into a
#'     good and an inferior cluster}
#'   \item{\code{\link{lppls_curvature}}}{Whether the estimate sits at a genuine
#'     interior minimum, which is what the modified profile likelihood assumes}
#'   \item{\code{\link{residual_qq_plot}}}{Calibration residuals against a normal
#'     or fitted Student-t reference}
#' }
#'
#' @section Rolling windows and Lagrange regularization:
#' Each of these returns an object that summarises itself when printed.
#' \describe{
#'   \item{\code{\link{lppls_rolling}}}{Recalibrate over a rolling start time
#'     \eqn{T_1}{T1}, holding the end of the window fixed}
#'   \item{\code{\link{lppls_lagrange}}}{Lagrange regularization of a rolling
#'     calibration, locating the onset \eqn{T_1^*}{T1*}}
#'   \item{\code{\link{lr_trend_test}}}{Whether the fitted bias trend is
#'     distinguishable from noise, by parametric bootstrap}
#'   \item{\code{\link{rolling_tc_plot}}}{\eqn{\hat{t}_c}{tc-hat} and its likelihood
#'     intervals against \eqn{T_1}{T1}}
#'   \item{\code{\link{rolling_param_plot}}}{Each estimated parameter against
#'     \eqn{T_1}{T1}}
#'   \item{\code{\link{chi_sq_np_plot}}, \code{\link{chi_sq_lambda_plot}}}{The
#'     Lagrange diagnostics \eqn{\chi^2_{np}(T_1)}{chi^2_np(T1)} and
#'     \eqn{\chi^2_{\lambda}(T_1)}{chi^2_lambda(T1)}}
#' }
#'
#' @section Visualising the model and its landscape:
#' \describe{
#'   \item{\code{\link{fit_plot}}}{A calibrated series with its fitted LPPLS curve}
#'   \item{\code{\link{contour_plot_sse}}}{SSE over a grid of two parameters}
#'   \item{\code{\link{sse_surface_plot}}}{The same landscape as an interactive
#'     or a static surface}
#' }
#'
#' @section Model Description:
#' The LPPLS model describes the log-price dynamics of financial assets
#' approaching a critical point (bubble peak or crash):
#'
#' \deqn{y(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c - t)) + C_2 \sin(\omega \log(t_c - t))]}{y(t) = A + (tc - t)^m [B + C1 cos(omega log(tc - t)) + C2 sin(omega log(tc - t))]}
#'
#' where:
#' \itemize{
#'   \item \eqn{t_c} is the critical time
#'   \item \eqn{m} is the power law exponent (typically 0.1-0.9)
#'   \item \eqn{\omega} is the log-periodic frequency (typically 6-13)
#'   \item \eqn{A, B, C_1, C_2} are linear parameters
#' }
#'
#' @section Optimization Modes:
#' \describe{
#'   \item{F1}{Simultaneous optimization of all nonlinear parameters}
#'   \item{F2}{Two-stage: optimize (m, omega) for each tc, then select best tc}
#'   \item{MPL}{Modified profile likelihood, with likelihood intervals for tc}
#' }
#'
#' @references
#' Filimonov, V., & Sornette, D. (2013). A stable and robust calibration scheme
#' of the log-periodic power law model. Physica A, 392(17), 3698-3707.
#'
#' Filimonov, V., Demos, G., & Sornette, D. (2017). Modified profile likelihood
#' inference and interval forecast of the burst of financial bubbles.
#' Quantitative Finance, 17(8), 1167-1186.
#'
#' @docType package
#' @name lpplsF-package
#' @aliases lpplsF-package
"_PACKAGE"


## Column names referred to bare inside ggplot2::aes() and dplyr verbs. R's code
## analysis cannot see that these are resolved against a data frame at evaluation
## time, so it reports them as undefined globals; declaring them here is the
## conventional remedy. Keep in sync with the aes() calls -- a stale entry is
## harmless, a missing one reappears as a check NOTE.
utils::globalVariables(c(
  ## fit tables (fit[[2]], fit[[4]]) and the parameter / basin plots
  "ID", "tc", "value", "m", "omega", "B", "estimate", "basin",
  ## rolling tables and their likelihood-interval columns
  "T1", "tc_hat_mpl", "LI5l", "LI5u", "LI10l", "LI10u", "LI50l", "LI50u",
  ## Lagrange diagnostics
  "chi_sq_np", "chi_sq_lambda",
  ## MPL, SSE-surface and trace plots
  "log_mpl", "xintercept", "series", "sse", "vals", "log_p", "y"
))
