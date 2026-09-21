#' S&P 500 Composite Index Daily Closing Level
#'
#' The daily closing level of the S&P 500 Composite Index, spanning the period
#' 2009-05-18 to 2020-12-31. This dataset is included as it is the series
#' calibrated in the worked S&P 500 example, whose modified profile likelihood
#' estimate places the critical time at 2020-03-03, at the onset of the
#' COVID-19 crash. The window is the calibrated one and no more: 2,928 trading
#' days, fitted on the first 2,428 with `hold_out = 500` and `fh = 500`.
#'
#' The series is supplied as **levels**, not returns, because the LPPLS model is
#' a model of the log-price path: [fit_lppls()] takes `log(close)` directly.
#'
#' The original analysis used the CRSP S&P 500 Composite series (variable
#' `spindx`, CRSP Stock Market Indexes, obtained through Wharton Research Data
#' Services), which is licensed and may not be redistributed. The series shipped
#' here is the corresponding Yahoo Finance series over the same trading days.
#' The two agree closely: the trading calendars are identical over the window,
#' the mean absolute relative difference in level is 1.6e-07, and the largest
#' single-day discrepancy is 0.39 index points, on 2019-08-12. Recalibrating on
#' this series leaves every reported quantity unchanged — the modified profile
#' likelihood estimate \eqn{\hat{t}_c = 2717}{tc-hat = 2717} (2020-03-03),
#' \eqn{m = 0.9}, \eqn{\omega = 6.1596}{omega = 6.1596}, and all three
#' likelihood intervals are identical; only the Brent-refined \eqn{t_c}{tc} at
#' the final step of the MPL computation differs, in its seventh decimal place.
#' That step refines thesis equation 38 at the
#' \eqn{(\hat{m}, \hat{\omega})}{(m-hat, omega-hat)} which minimise equation 39
#' on the integer grid.
#'
#' @format `spx`:
#' A data.frame containing 2928 observations in 2 columns:
#' \describe{
#'   \item{date}{The trading date, of class \code{Date}.}
#'   \item{close}{The closing level of the index.}
#' }
#'
#' @source
#' Yahoo Finance, ticker \code{^GSPC}, retrieved 2026-08-24 using
#' \code{quantmod::getSymbols()}. The closing level is used; for an index the
#' closing and adjusted-closing series are identical.
#'
#' @examples
#' # the log-price path the LPPLS model is fitted to
#' plot(spx$date, log(spx$close), type = "l",
#'      xlab = "", ylab = "log price", main = "S&P 500 Composite")
#'
#' # the calibration window used in the worked example
#' n_cal <- nrow(spx) - 500
#' spx$date[n_cal]
#'
#' @keywords datasets
"spx"
