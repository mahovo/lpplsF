#' Check LPPLS Estimates Against the Admissibility Conditions
#'
#' A calibration is admissible when its estimate meets the conditions
#' [fit_lppls()] works with: `m` and `omega` inside their bounds, `B` below
#' `upper[3]`, and the damping `D` at least `lower[4]` (Filimonov, Demos &
#' Sornette, 2017). `fit_lppls()` bounds its search in `m` and `omega`, but
#' checks `B` and `D` only on the result of each search, and
#' [print.lppls_fit()] reports how many results were out of range.
#' `lppls_admissible()` applies the same conditions to any set of estimates: the
#' estimate of a fit, its per-`tc` fits, every window of a rolling calibration,
#' or a data frame of estimates.
#'
#' @details
#' The conditions are those of `fit_lppls()`: `lower[1] <= m <= upper[1]`,
#' `lower[2] <= omega <= upper[2]`, `B < upper[3]` and `D >= lower[4]`.
#' `lower[3]` and `upper[4]` are not used.
#'
#' The share of admissible windows in a rolling calibration,
#' `mean(lppls_admissible(roll)$admissible)`, measures how consistently a bubble
#' signal holds up as the window start moves.
#'
#' @param x An `"lppls_fit"` from [fit_lppls()], an `"lppls_rolling"` from
#'   [lppls_rolling()], or estimates with elements `m`, `omega`, `B` and `D`:
#'   a data frame with one row per estimate, or a list holding one estimate,
#'   such as `fit$fit[[1]]`.
#' @param which For an `"lppls_fit"`: `"estimate"` (the default) checks the
#'   reported estimate, `"fits"` every per-`tc` fit (every restart in mode
#'   `"F1"`), including those the `B` filter set aside.
#' @param lower,upper Numeric vectors of length 4, bounds for
#'   `c(m, omega, B, D)` as in [fit_lppls()]. By default the bounds the
#'   calibration was run with; for a data frame or a list, the defaults of
#'   `fit_lppls()`.
#' @param ... Not used.
#'
#' @return The estimates checked, as a data frame with one row per estimate and
#'   five logical columns added: `m_ok`, `omega_ok`, `B_ok`, `D_ok`, and
#'   `admissible`, which is `TRUE` when all four are. A missing value counts as
#'   out of range.
#'
#' @references
#' Filimonov, V., Demos, G., & Sornette, D. (2017). Modified profile likelihood
#' inference and interval forecast of the burst of financial bubbles.
#' *Quantitative Finance*, 17(8), 1167-1186. \doi{10.1080/14697688.2016.1276298}
#'
#' @seealso [fit_lppls()], [print.lppls_fit()], [lppls_rolling()]
#'
#' @examples
#' set.seed(4)
#' log_p <- eval_lppls(1:200, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
#'                     tc = 221, m = 0.6, omega = 8) + rnorm(200, sd = 0.01)
#' fit <- fit_lppls(log_p, fh = 30, mode = "MPL",
#'                  m_init = 0.5, o_init = 8, num_searches = 2)
#'
#' ## The reported estimate
#' lppls_admissible(fit)
#'
#' ## Every per-tc fit, against a stricter damping floor of 1
#' adm <- lppls_admissible(fit, which = "fits", lower = c(0.1, 6, -1e14, 1))
#' mean(adm$admissible)
#'
#' @export
lppls_admissible <- function(x, ...) {
  UseMethod("lppls_admissible")
}


#' @rdname lppls_admissible
#' @export
lppls_admissible.lppls_fit <- function(x, which = c("estimate", "fits"),
                                       lower = NULL, upper = NULL, ...) {
  which <- match.arg(which)
  a <- x$fit_args
  if ((is.null(lower) && is.null(a$lower)) || (is.null(upper) && is.null(a$upper))) {
    stop("This fit does not record its bounds; supply 'lower' and 'upper'.")
  }
  if (is.null(lower)) lower <- a$lower
  if (is.null(upper)) upper <- a$upper

  est <- if (which == "estimate") {
    x$fit[[1]]
  } else if (length(x$fit) >= 3L && is.data.frame(x$fit[[3]])) {
    ## Modes "F2" and "MPL" keep the fits that passed the B filter in fit[[2]]
    ## and all of them in fit[[3]] -- when any passed; otherwise fit[[2]] holds
    ## them all, as it does in mode "F1".
    x$fit[[3]]
  } else {
    x$fit[[2]]
  }
  .admissible(est, lower, upper)
}


#' @rdname lppls_admissible
#' @export
lppls_admissible.lppls_rolling <- function(x, lower = NULL, upper = NULL, ...) {
  ## Bounds not passed to lppls_rolling() were fit_lppls()'s defaults, which
  ## .admissible() falls back to.
  if (is.null(lower)) lower <- x$fit_args$lower
  if (is.null(upper)) upper <- x$fit_args$upper
  .admissible(x$table, lower, upper)
}


#' @rdname lppls_admissible
#' @export
lppls_admissible.default <- function(x, lower = NULL, upper = NULL, ...) {
  .admissible(x, lower, upper)
}


#' The checks behind lppls_admissible() and the out-of-range flag of
#' print.lppls_fit()
#'
#' @param est A data frame of estimates, or a list holding one estimate.
#' @param lower,upper Bounds for `c(m, omega, B, D)`; `NULL` falls back to
#'   fit_lppls()'s defaults.
#' @return `est` as a data frame, with the five logical columns added.
#'
#' @keywords internal
#' @noRd
.admissible <- function(est, lower = NULL, upper = NULL) {
  if (is.null(lower)) lower <- eval(formals(fit_lppls)$lower)
  if (is.null(upper)) upper <- eval(formals(fit_lppls)$upper)
  if (!is.numeric(lower) || length(lower) != 4L ||
      !is.numeric(upper) || length(upper) != 4L) {
    stop("'lower' and 'upper' must be numeric vectors of length 4")
  }
  if (!is.data.frame(est)) {
    if (!is.list(est)) {
      stop("'x' must be an lppls_fit, an lppls_rolling, a data frame or a ",
           "list holding one estimate")
    }
    est <- as.data.frame(lapply(est, unname))
  }
  miss <- setdiff(c("m", "omega", "B", "D"), names(est))
  if (length(miss)) {
    stop("'x' has no ", paste(miss, collapse = ", "), " to check")
  }

  ok <- function(z) !is.na(z) & z
  est$m_ok       <- ok(est$m >= lower[1] & est$m <= upper[1])
  est$omega_ok   <- ok(est$omega >= lower[2] & est$omega <= upper[2])
  est$B_ok       <- ok(est$B < upper[3])
  est$D_ok       <- ok(est$D >= lower[4])
  est$admissible <- est$m_ok & est$omega_ok & est$B_ok & est$D_ok
  est
}
