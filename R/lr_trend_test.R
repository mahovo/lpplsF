## Internal: fit AR(1)-GARCH(1,1) to a residual vector (AR mean via arima, GARCH via
## tsgarch). Returns the coefficients needed to simulate fresh paths.
.fit_ar_garch <- function(r) {
  a     <- stats::arima(r, order = c(1, 0, 0), include.mean = FALSE)
  innov <- as.numeric(stats::residuals(a))
  xr    <- xts::xts(innov, order.by = as.Date("2000-01-01") + seq_along(innov) - 1L)
  spec  <- tsgarch::garch_modelspec(xr, model = "garch", constant = FALSE, order = c(1, 1))
  ## only the coefficients are used; suppress tsgarch's std-error ("NaNs produced") noise
  cf    <- suppressWarnings(stats::coef(tsgarch::estimate.tsgarch.spec(spec)))
  list(ar1 = unname(a$coef["ar1"]), om0 = unname(cf["omega"]),
       alpha1 = unname(cf["alpha1"]), beta1 = unname(cf["beta1"]))
}

## Internal: simulate one AR(1)-GARCH(1,1) path of length n from fitted params.
.sim_ar_garch <- function(n, p, seed = NULL, burn = 500) {
  if (!is.null(seed)) set.seed(seed)
  a1 <- p$alpha1; b1 <- p$beta1; om <- p$om0; ar <- p$ar1
  z <- stats::rnorm(n + burn); s2 <- numeric(n + burn); ee <- numeric(n + burn)
  s2[1] <- om / max(1e-6, 1 - a1 - b1); ee[1] <- sqrt(s2[1]) * z[1]
  for (t in 2:(n + burn)) {
    s2[t] <- om + a1 * ee[t - 1]^2 + b1 * s2[t - 1]
    ee[t] <- sqrt(max(s2[t], 0)) * z[t]
  }
  y <- numeric(n + burn)
  for (t in 2:(n + burn)) y[t] <- ar * y[t - 1] + ee[t]
  y[(burn + 1):(burn + n)]
}

#' Autocorrelation-robust significance test for the Lagrange \eqn{\chi^2_{np}}{chi^2_np} trend
#'
#' Assesses whether the decreasing "bias" trend that Lagrange regularization acts on
#' (the slope \eqn{\beta_1}{beta_1} of \eqn{\chi^2_{np}(T_1)}{chi^2_np(T1)} against the window start
#' \eqn{T_1}) is statistically distinguishable from what the estimation noise produces
#' on its own. Three tests are reported: a naive OLS \eqn{p}-value (invalid here,
#' because adjacent rolling windows overlap and \eqn{\chi^2_{np}(T_1)}{chi^2_np(T1)} is heavily
#' autocorrelated); a heteroskedasticity- and autocorrelation-consistent (Newey--West)
#' \eqn{p}-value; and a **parametric AR(1)--GARCH(1,1) bootstrap**, which is the
#' appropriate test.
#'
#' @details
#' The bootstrap fits an AR(1)--GARCH(1,1) model (via [stats::arima()] and the
#' \pkg{tsgarch} package) to the residuals of the full-window (\eqn{T_1 = 1}) LPPLS fit,
#' simulates fresh noise paths, adds them to the fitted LPPLS curve, and repeats the
#' rolling calibration to build the null distribution of \eqn{\beta_1}. The observed
#' slope is compared against this null (two-sided). A validity guard is applied: if the
#' fitted GARCH persistence \eqn{\alpha + \beta \ge 1} (IGARCH) the single-LPPLS null is
#' misspecified, so the bootstrap is skipped, a warning is issued and `p_boot` is `NA`.
#'
#' Requires the suggested packages \pkg{tsgarch}, \pkg{xts} and (for `hac = TRUE`)
#' \pkg{sandwich}.
#'
#' @param rolling An `"lppls_rolling"` object from [lppls_rolling()]; it must carry the
#'   `fit_args` slot (rebuild with a current version of the package if absent).
#' @param log_price Numeric vector of log-prices the `rolling` was built on; its length
#'   must equal `rolling$n_full`.
#' @param B Integer, number of bootstrap replicates (default `199`).
#' @param cores Integer, number of cores for [parallel::mclapply()] (default `1`;
#'   ignored on Windows).
#' @param hac Logical, whether to also compute the Newey--West HAC \eqn{p}-value
#'   (default `TRUE`).
#' @param seed Integer or `NULL`; a base seed that makes the bootstrap paths
#'   reproducible.
#'
#' @return An object of class `"lr_trend_test"`: a list with the observed slope `beta1`
#'   and `lambda` (\eqn{= -\beta_1}); the naive `p_ols` / `se_ols`; the HAC `p_hac` /
#'   `se_hac`; the bootstrap `p_boot`, standardized distance `z` and null vector `null`;
#'   the fitted noise model (`ar1`, `garch` = `c(omega, alpha1, beta1)`, and
#'   `persist` \eqn{= \alpha + \beta}) with its `valid` flag; the effective replicate
#'   count `B`, the window count `n_windows`, and the `call`.
#'
#' @seealso [lppls_lagrange()], [lppls_rolling()]
#'
#' @examples
#' \dontrun{
#' roll <- lppls_rolling(log_p, t1_seq = seq(1, 1001, 25), hold_out = 100, fh = 150,
#'                       fit_args = list(mode = "F2", num_searches = 1))
#' lr_trend_test(roll, log_price = log_p, B = 199, cores = 4, seed = 1)
#' }
#'
#' @export
lr_trend_test <- function(rolling, log_price, B = 199, cores = 1L,
                          hac = TRUE, seed = NULL) {

  if (!inherits(rolling, "lppls_rolling")) {
    stop("'rolling' must be an object of class 'lppls_rolling' (from lppls_rolling())")
  }
  if (!is.numeric(log_price) || length(log_price) != rolling$n_full) {
    stop(sprintf("'log_price' must be a numeric vector of length n_full = %d",
                 rolling$n_full))
  }
  fit_args <- rolling$fit_args
  if (is.null(fit_args)) {
    stop("'rolling' has no 'fit_args' slot; rebuild it with lppls_rolling().")
  }
  need <- c("tsgarch", "xts", if (isTRUE(hac)) "sandwich")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("lr_trend_test() needs package(s): ", paste(miss, collapse = ", "),
         if ("sandwich" %in% miss) " (or set hac = FALSE)." else ".")
  }

  ## observed bias trend (chi^2_np ~ T1) ------------------------------------
  lag_obs <- lppls_lagrange(rolling, p = 0)
  b1  <- lag_obs$slope
  fit <- lag_obs$trend
  sm  <- summary(fit)
  p_ols  <- unname(sm$coefficients["T1", "Pr(>|t|)"])
  se_ols <- unname(sm$coefficients["T1", "Std. Error"])
  p_hac <- se_hac <- NA_real_
  if (isTRUE(hac)) {
    se_hac <- sqrt(sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE)["T1", "T1"])
    p_hac  <- 2 * stats::pt(-abs(b1 / se_hac), df = fit$df.residual)
  }

  ## fit AR-GARCH to the residuals of the T1 = 1 fit ------------------------
  T2    <- rolling$T2
  th    <- rolling$fits[[1]]
  ff    <- eval_lppls(seq_len(rolling$n_full), th$A, th$B, th$C1, th$C2,
                      th$tc, th$m, th$omega)
  resid <- log_price[seq_len(T2)] - ff[seq_len(T2)]
  pars  <- .fit_ar_garch(resid)
  persist <- pars$alpha1 + pars$beta1
  valid   <- persist < 1

  ## Parametric bootstrap: fresh AR-GARCH paths -> re-roll -> null of beta1.
  ## Skipped when the fitted noise model is misspecified (IGARCH, alpha + beta >= 1):
  ## the null would be unreliable, so we warn and return NA for the bootstrap.
  if (!valid) {
    warning("fitted GARCH persistence alpha + beta = ", signif(persist, 3),
            " >= 1 (IGARCH): the single-LPPLS null is misspecified; the bootstrap ",
            "is skipped and its p-value is NA.")
    b1b <- numeric(0); p_boot <- NA_real_; z <- NA_real_
  } else {
    t1_seq <- rolling$t1_seq; hold_out <- rolling$hold_out; fh <- rolling$fh
    one_boot <- function(b) {
      bseed <- if (is.null(seed)) NULL else seed + b
      y <- ff + .sim_ar_garch(rolling$n_full, pars, seed = bseed)
      tryCatch(
        lppls_lagrange(
          lppls_rolling(y, t1_seq = t1_seq, hold_out = hold_out, fh = fh,
                        fit_args = fit_args),
          p = 0)$slope,
        error = function(e) NA_real_)
    }
    idx <- seq_len(B)
    b1b <- if (cores > 1L && .Platform$OS.type != "windows") {
      unlist(parallel::mclapply(idx, one_boot, mc.cores = cores))
    } else {
      vapply(idx, one_boot, numeric(1))
    }
    b1b    <- b1b[is.finite(b1b)]
    p_boot <- mean(abs(b1b) >= abs(b1))
    z      <- abs(b1) / stats::sd(b1b)
  }

  structure(
    list(
      beta1  = b1, lambda = -b1,
      p_ols  = p_ols,  se_ols = se_ols,
      p_hac  = p_hac,  se_hac = se_hac,
      p_boot = p_boot, z = z, null = b1b,
      ar1    = pars$ar1,
      garch  = c(omega = pars$om0, alpha1 = pars$alpha1, beta1 = pars$beta1),
      persist = persist, valid = valid,
      B = length(b1b), n_windows = nrow(rolling$table),
      call = match.call()
    ),
    class = "lr_trend_test"
  )
}

#' @export
print.lr_trend_test <- function(x, ...) {
  cat("Lagrange chi^2_np bias-trend test\n")
  cat(sprintf("  slope beta1 = %.3g   (lambda = %.3g)   over %d windows\n",
              x$beta1, x$lambda, x$n_windows))
  cat(sprintf("  OLS p       = %.3g%s\n", x$p_ols,
              if (!is.na(x$p_hac)) sprintf("      HAC (Newey-West) p = %.3g", x$p_hac) else ""))
  cat(sprintf("  AR(1)-GARCH(1,1): ar1 = %.2f, alpha + beta = %.3f  [%s]\n",
              x$ar1, x$persist, if (isTRUE(x$valid)) "valid" else "IGARCH — unreliable"))
  if (isTRUE(x$valid)) {
    cat(sprintf("  parametric bootstrap: p = %.3f, z = %.2f  (B = %d)\n",
                x$p_boot, x$z, x$B))
  } else {
    cat("  parametric bootstrap: skipped (misspecified IGARCH null)\n")
  }
  invisible(x)
}
