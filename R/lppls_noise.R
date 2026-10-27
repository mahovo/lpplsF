## Noise sweeps ---------------------------------------------------------------
##
## How much noise can a calibration take before its estimates stop meaning
## anything? The answer is not a property of the data but of the model, so it is
## studied on a known signal: take a deterministic LPPLS path, add noise at a
## range of levels, recalibrate at each, and watch the estimates move.
##
## `lppls_noise_sweep()` runs that experiment; `noise_param_plot()` and
## `noise_tc_plot()` draw it.


#' Calibrate an LPPLS signal across a range of noise levels
#'
#' Adds noise to a known LPPLS signal at each of several levels, recalibrates at
#' every level, and collects the estimates into one table. Because the signal is
#' known, the table shows directly how far each estimate drifts from truth as the
#' noise grows, and how the \eqn{t_c} likelihood intervals widen with it.
#'
#' The noise is supplied as a function, so the experiment is not limited to
#' Gaussian errors. Any generator of standardised innovations will do, which is
#' how a heavy-tailed or serially dependent alternative is compared against the
#' Gaussian case on the same signal.
#'
#' @param signal Numeric vector, the deterministic log-price path to perturb.
#'   Typically the output of [eval_lppls()].
#' @param sd_seq Numeric vector of noise levels. Each entry scales one draw from
#'   `noise`, and contributes one row to the returned table.
#' @param noise Function of one argument `n`, returning `n` standardised
#'   innovations. Defaults to [stats::rnorm]. The draw is multiplied by the
#'   current entry of `sd_seq`, so `noise` should have unit scale for `sd_seq` to
#'   read as a standard deviation.
#' @param fit_args Named list of arguments passed to [fit_lppls()]. Use
#'   `mode = "MPL"` for the likelihood-interval columns; without it those columns
#'   are `NA` and [noise_tc_plot()] has nothing to draw.
#' @param seed Optional integer. When supplied, the generator is re-seeded with
#'   it before *each* level, so every level sees the same innovations scaled
#'   differently and the comparison across levels is not confounded by the draw.
#'   The caller's random number stream is left as it was. When `NULL` (default)
#'   the levels draw independently, from the caller's stream.
#' @param fb Logical, report progress per level (default `FALSE`).
#'
#' @return A data.frame of class `"lppls_noise_sweep"`, one row per entry of
#'   `sd_seq`, with columns:
#'   \describe{
#'     \item{sd}{The noise level.}
#'     \item{tc, m, omega, A, B, C1, C2, D}{The parameter estimates.}
#'     \item{value}{The objective at the optimum.}
#'     \item{tc_hat_mpl}{The modified profile likelihood estimate of \eqn{t_c},
#'       `NA` unless `fit_args` asked for `mode = "MPL"`.}
#'     \item{LI5l, LI5u, ...}{Lower and upper bounds of each likelihood interval,
#'       named for its cutoff as a percentage, so a cutoff of 0.05 gives `LI5l`
#'       and `LI5u`.}
#'   }
#'
#' @seealso [noise_param_plot()] and [noise_tc_plot()] to draw the result;
#'   [fit_lppls()] for the calibration itself.
#'
#' @examples
#' # a short, coarse sweep; a real one uses more levels and a longer signal
#' N      <- 400
#' signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
#'                      C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
#' sweep <- lppls_noise_sweep(
#'   signal, sd_seq = c(0.01, 0.05),
#'   fit_args = list(fh = 20, hold_out = 20, num_searches = 1, mode = "F2"),
#'   seed = 3630
#' )
#' sweep[, c("sd", "tc", "m", "omega")]
#'
#' @export
lppls_noise_sweep <- function(signal, sd_seq, noise = stats::rnorm,
                              fit_args = list(), seed = NULL, fb = FALSE) {

  if (!is.numeric(signal) || length(signal) < 2) {
    stop("'signal' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (!is.numeric(sd_seq) || !length(sd_seq)) {
    stop("'sd_seq' must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(sd_seq < 0)) {
    stop("'sd_seq' must not contain negative values", call. = FALSE)
  }
  if (!is.function(noise)) {
    stop("'noise' must be a function of one argument, n", call. = FALSE)
  }
  if (!is.list(fit_args)) {
    stop("'fit_args' must be a list", call. = FALSE)
  }

  n_obs <- length(signal)
  rows  <- vector("list", length(sd_seq))

  ## Re-seeding resets the caller's random number stream; put it back on exit.
  if (!is.null(seed)) .preserve_seed()

  for (i in seq_along(sd_seq)) {
    if (fb) {
      message(sprintf("noise level %d/%d (sd = %g)", i, length(sd_seq), sd_seq[i]))
    }
    ## Re-seeding per level makes the innovations common across levels, so a
    ## change in the estimates is attributable to the noise scale alone.
    if (!is.null(seed)) set.seed(seed)
    y <- signal + noise(n_obs) * sd_seq[i]

    fit <- do.call(fit_lppls, c(list(log_price = y), fit_args))
    rows[[i]] <- .noise_sweep_row(fit, sd_seq[i])
  }

  ## Levels may differ in which likelihood intervals were defined; align columns.
  nms <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    missing <- setdiff(nms, names(r))
    for (m in missing) r[[m]] <- NA_real_
    r[nms]
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("lppls_noise_sweep", "data.frame"))
}


## One row of the sweep table from one fit. Likelihood-interval columns are named
## for their cutoff as a percentage (0.05 -> LI5l/LI5u), which is the convention
## lppls_rolling() already uses.
.noise_sweep_row <- function(fit, sd_value) {
  e <- fit$fit[[1]]
  g <- function(v) if (is.null(v)) NA_real_ else unname(as.numeric(v))

  row <- data.frame(
    sd    = sd_value,
    tc    = g(e$tc),    m  = g(e$m),  omega = g(e$omega),
    A     = g(e$A),     B  = g(e$B),  C1    = g(e$C1), C2 = g(e$C2),
    D     = g(e$D),     value = g(e$value),
    tc_hat_mpl = NA_real_
  )

  mpl <- fit$mpl_output
  if (!is.null(mpl)) {
    row$tc_hat_mpl <- g(mpl$tc_hat_mpl)
    cut <- fit$fit_args$mpl_cutoff
    li  <- mpl$LI
    if (!is.null(li) && length(li)) {
      lab <- if (!is.null(cut) && length(cut) == length(li)) {
        sprintf("%g", cut * 100)
      } else {
        as.character(seq_along(li))
      }
      for (j in seq_along(li)) {
        row[[paste0("LI", lab[j], "l")]] <- g(li[[j]][1])
        row[[paste0("LI", lab[j], "u")]] <- g(li[[j]][2])
      }
    }
  }
  row
}


## Accept either a sweep object or a bare data.frame carrying the same columns,
## so tables saved before the class existed still plot.
.check_noise_sweep <- function(x, need = "sd") {
  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame from lppls_noise_sweep()", call. = FALSE)
  }
  missing <- setdiff(need, names(x))
  if (length(missing)) {
    stop("'x' is missing the column(s): ", paste(missing, collapse = ", "),
         ". Was it produced by lppls_noise_sweep()?", call. = FALSE)
  }
  x
}


#' Parameter estimates against noise level
#'
#' Draws each estimated parameter from a noise sweep as a function of the noise
#' level, one panel per parameter. The log-periodic parameters are typically the
#' first to degrade, which is the point of looking.
#'
#' @param x A `"lppls_noise_sweep"` table from [lppls_noise_sweep()], or any
#'   data.frame with an `sd` column and the parameter columns.
#' @param params Character vector naming the parameters to draw. Defaults to all
#'   six LPPLS parameters.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_noise_sweep()], [noise_tc_plot()]
#'
#' @examples
#' N      <- 400
#' signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
#'                      C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
#' sweep <- lppls_noise_sweep(
#'   signal, sd_seq = c(0.01, 0.03, 0.05),
#'   fit_args = list(fh = 20, hold_out = 20, num_searches = 1, mode = "F2"),
#'   seed = 3630
#' )
#' noise_param_plot(sweep)
#'
#' @export
noise_param_plot <- function(x, params = c("m", "omega", "A", "B", "C1", "C2")) {
  x <- .check_noise_sweep(x, c("sd", params))

  long <- do.call(rbind, lapply(params, function(p) {
    data.frame(sd = x$sd, param = p, estimate = x[[p]], stringsAsFactors = FALSE)
  }))
  long$param <- factor(long$param, levels = params)

  ggplot2::ggplot(long, ggplot2::aes(sd, estimate)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~param, scales = "free_y", ncol = 1) +
    ggplot2::xlab("noise sd") +
    ggplot2::ylab("estimate") +
    ggplot2::theme_minimal()
}


#' Critical-time estimate and likelihood intervals against noise level
#'
#' Draws the \eqn{t_c} point estimate from a noise sweep, together with the
#' modified profile likelihood estimate and the nested likelihood-interval bands,
#' as a function of the noise level. The bands widen as the noise grows, which is
#' the honest reading of how much the critical time can be pinned down.
#'
#' @param x A `"lppls_noise_sweep"` table from [lppls_noise_sweep()], or any
#'   data.frame with an `sd` column, a `tc` column and paired `LI*l`/`LI*u`
#'   columns.
#' @param true_tc Optional numeric. When supplied, a dashed horizontal reference
#'   line is drawn there, which is useful on a simulated signal whose critical
#'   time is known.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_noise_sweep()], [noise_param_plot()]
#'
#' @examples
#' N      <- 400
#' signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
#'                      C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
#' sweep <- lppls_noise_sweep(
#'   signal, sd_seq = c(0.01, 0.03, 0.05),
#'   fit_args = list(fh = 20, hold_out = 20, num_searches = 1,
#'                   mode = "MPL", mpl_cutoff = c(0.05, 0.5)),
#'   seed = 3630
#' )
#' noise_tc_plot(sweep, true_tc = N + 1)
#'
#' @export
noise_tc_plot <- function(x, true_tc = NULL) {
  x <- .check_noise_sweep(x, c("sd", "tc"))

  ## Find the LI column pairs and order them widest cutoff last, so the narrow
  ## bands stay visible on top of the wide ones.
  lo_cols <- grep("^LI[0-9.]+l$", names(x), value = TRUE)
  pairs   <- lapply(lo_cols, function(lo) {
    hi <- sub("l$", "u", lo)
    if (!hi %in% names(x)) return(NULL)
    list(lo = lo, hi = hi,
         pct = as.numeric(sub("^LI([0-9.]+)l$", "\\1", lo)))
  })
  pairs <- Filter(Negate(is.null), pairs)
  if (!length(pairs)) {
    stop("No likelihood-interval columns found. Run lppls_noise_sweep() with ",
         "fit_args = list(mode = 'MPL', ...).", call. = FALSE)
  }
  pairs <- pairs[order(vapply(pairs, function(p) p$pct, numeric(1)))]

  bands <- do.call(rbind, lapply(pairs, function(p) {
    data.frame(sd = x$sd, lo = x[[p$lo]], hi = x[[p$hi]],
               band = sprintf("LI_%02g", p$pct), stringsAsFactors = FALSE)
  }))
  band_lv <- unique(bands$band)
  bands$band <- factor(bands$band, levels = band_lv)

  ## Match the palette lppls_rolling()'s tc plot uses, extended if needed.
  band_cols <- c("green", "orange", "red", "purple", "brown")
  cols <- stats::setNames(band_cols[seq_along(band_lv)], band_lv)
  cols <- c(cols, "MPLE" = "cyan", "PLE" = "black")

  p <- ggplot2::ggplot(x, ggplot2::aes(sd, tc)) +
    ggplot2::geom_ribbon(
      data = bands,
      mapping = ggplot2::aes(x = sd, ymin = lo, ymax = hi, color = band),
      inherit.aes = FALSE, fill = "grey50", alpha = 0.2, key_glyph = "path"
    ) +
    ggplot2::geom_line(ggplot2::aes(color = "PLE"), linewidth = 1)

  if ("tc_hat_mpl" %in% names(x) && !all(is.na(x$tc_hat_mpl))) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(y = tc_hat_mpl, color = "MPLE"), linewidth = 0.25)
  }
  if (!is.null(true_tc)) {
    p <- p + ggplot2::geom_hline(yintercept = true_tc, linetype = "dashed",
                                 color = "grey30")
  }

  p +
    ggplot2::scale_color_manual(name = NULL, values = cols,
                                breaks = names(cols)) +
    ggplot2::xlab("noise sd") +
    ggplot2::ylab(quote(hat(t)[c])) +
    ggplot2::theme_minimal()
}
