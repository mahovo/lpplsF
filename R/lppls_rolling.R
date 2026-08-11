#' Rolling-Window LPPLS Calibration (rolling \eqn{T_1}, fixed \eqn{T_2})
#'
#' Calibrates a sequence of LPPLS models that share a common calibration end
#' \eqn{T_2} while the start time \eqn{T_1} rolls over `t1_seq`. This is the
#' multi-window primitive underlying both the window-diagnostic plots
#' ([rolling_tc_plot()], [rolling_param_plot()]) and Lagrange regularization
#' ([lppls_lagrange()]): the (expensive) calibrations are computed once here and
#' reused by all downstream analyses.
#'
#' @details
#' For each \eqn{T_1} in `t1_seq` the model is calibrated on
#' `log_price[T1:length(log_price)]` via [fit_lppls()]. The last `hold_out`
#' observations of every window are reserved for the forecast/validation region,
#' so the calibration region is the absolute interval \eqn{[T_1, T_2]} with the
#' common end \eqn{T_2 = \texttt{length(log\_price)} - \texttt{hold\_out}}{T2 = length(log_price) - hold_out} held
#' fixed across windows.
#'
#' Because `fit_lppls()` re-indexes each window to start at 1, all per-window
#' results (`tc`, the likelihood intervals and `tc_hat_mpl`) are shifted back to
#' absolute time by adding \eqn{T_1 - 1} before being returned.
#'
#' The likelihood-interval columns (`LI5l` ... `tc_hat_mpl`) are populated only
#' when the fit is run in `mode = "MPL"` (the default); under `mode = "F2"` they
#' are `NA` but the point estimates and `sse_cal` (hence Lagrange) are identical.
#'
#' @param log_price Numeric vector of log-prices for the full series.
#' @param t1_seq Integer vector of rolling start indices (absolute positions into
#'   `log_price`). All entries must be `< T2 = length(log_price) - hold_out`.
#' @param hold_out Integer, number of observations reserved from the end of the
#'   full series for the forecast/validation region (default 0). Fixes \eqn{T_2}.
#' @param fh Integer, forecast horizon passed to [fit_lppls()] (default 30).
#' @param tc_init Numeric or `NULL`. Initial critical time for [fit_lppls()]. If
#'   `NULL` (default) a per-window value `n_cal + floor(fh / 2)` is used, where
#'   `n_cal` is the window's calibration length.
#' @param fit_args Named list of additional arguments forwarded to [fit_lppls()]
#'   (e.g. `lower`, `upper`, `m_init`, `o_init`, `num_searches`, `mode`,
#'   `mpl_cutoff`, `factr`). The arguments `log_price`, `fh`, `hold_out`,
#'   `tc_init` and `fb` are controlled by `lppls_rolling()` and are ignored if
#'   supplied here. If `mode` is not supplied it defaults to `"MPL"`.
#' @param fb Logical, whether to print per-window feedback (default `FALSE`).
#'
#' @return An object of class `"lppls_rolling"`: a list with
#'   \describe{
#'     \item{table}{Data frame with one row per window: `T1`, `n_cal`, the
#'       absolute point estimates `tc`, `m`, `omega`, `A`, `B`, `C1`, `C2`, `D`,
#'       the calibration-region residual sum of squares `sse_cal`, the absolute
#'       likelihood-interval bounds `LI5l`/`LI5u`/`LI10l`/`LI10u`/`LI50l`/`LI50u`
#'       and the absolute modified-profile-likelihood estimate `tc_hat_mpl`.}
#'     \item{fits}{List of best-fit parameter lists (`fit[[1]]`) per window, in
#'       window-local time.}
#'     \item{t1_seq, T2, hold_out, fh, n_full}{The calibration setup.}
#'     \item{fit_args}{The effective list of arguments forwarded to [fit_lppls()]
#'       (after stripping the controlled ones and defaulting `mode`); recorded so
#'       [lr_trend_test()] can re-roll bootstrap replicates with the same settings.}
#'   }
#'
#' @section Profiling on macOS:
#' `mode` defaults to `"MPL"`, so profiling a rolling calibration on macOS hits
#' the R-profiler/Accelerate crash described under "Profiling `mode = \"MPL\"` on
#' macOS" in [fit_lppls()]. The same `options(matprod = "internal")` workaround
#' applies. Unprofiled runs are unaffected.
#'
#' @seealso [fit_lppls()], [lppls_lagrange()], [rolling_tc_plot()],
#'   [rolling_param_plot()]
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
#' print(rolling_tc_plot(roll))
#' }
#'
#' @export
lppls_rolling <- function(
  log_price,
  t1_seq,
  hold_out = 0,
  fh = 30,
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
      " are controlled by lppls_rolling() and will be ignored."
    )
    fit_args <- fit_args[setdiff(names(fit_args), reserved)]
  }
  if (is.null(fit_args$mode)) fit_args$mode <- "MPL"

  ## Roll over T1 ============================================================
  K <- length(t1_seq)
  fits <- vector("list", K)
  rows <- vector("list", K)

  ## Helper: offset-aware likelihood-interval bound (NA when not MPL)
  li_bound <- function(mpl, j, i, off) {
    if (is.null(mpl) || is.null(mpl$LI) || length(mpl$LI) < j) return(NA_real_)
    v <- mpl$LI[[j]]
    if (length(v) < i || is.na(v[i])) return(NA_real_)
    v[i] + off
  }

  for (k in seq_len(K)) {
    T1 <- t1_seq[k]
    window_lp <- log_price[T1:n_full]
    n_total <- length(window_lp)
    n_cal <- n_total - hold_out
    off <- T1 - 1                                  ## window-local -> absolute

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
    mpl <- res_fit$mpl_output

    ## Residual sum of squares over the calibration region [T1, T2] (absolute),
    ## i.e. window-local [1, n_cal]. Filimonov formulation (mode 0).
    t_cal <- seq_len(n_cal)
    sse_cal <- SSE(best, window_lp[t_cal], t_cal, mode = 0)

    rows[[k]] <- data.frame(
      T1         = T1,
      n_cal      = n_cal,
      tc         = unname(best$tc) + off,
      m          = unname(best$m),
      omega      = unname(best$omega),
      A          = unname(best$A),
      B          = unname(best$B),
      C1         = unname(best$C1),
      C2         = unname(best$C2),
      D          = unname(best$D),
      sse_cal    = sse_cal,
      LI5l       = li_bound(mpl, 1, 1, off),
      LI5u       = li_bound(mpl, 1, 2, off),
      LI10l      = li_bound(mpl, 2, 1, off),
      LI10u      = li_bound(mpl, 2, 2, off),
      LI50l      = li_bound(mpl, 3, 1, off),
      LI50u      = li_bound(mpl, 3, 2, off),
      tc_hat_mpl = if (is.null(mpl)) NA_real_ else unname(mpl$tc_hat_mpl) + off
    )
  }

  structure(
    list(
      table    = do.call(rbind, rows),
      fits     = fits,
      t1_seq   = t1_seq,
      T2       = T2,
      hold_out = hold_out,
      fh       = fh,
      n_full   = n_full,
      fit_args = fit_args
    ),
    class = "lppls_rolling"
  )
}


#' @export
print.lppls_rolling <- function(x, ...) {
  cat("LPPLS rolling-window calibration\n")
  cat(sprintf("  Windows (T1):  %d  (from %d to %d, step %d)\n",
              nrow(x$table), min(x$table$T1), max(x$table$T1),
              if (length(x$t1_seq) > 1) x$t1_seq[2] - x$t1_seq[1] else NA))
  cat(sprintf("  T2 (cal. end): %d\n", x$T2))
  cat(sprintf("  hold_out: %d   fh: %d\n", x$hold_out, x$fh))
  has_li <- any(!is.na(x$table$tc_hat_mpl))
  cat(sprintf("  Likelihood intervals: %s\n",
              if (has_li) "present (MPL)" else "absent (run with mode = 'MPL')"))
  invisible(x)
}
