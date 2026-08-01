## Display label for an LPPLS parameter under a given mode.
## In the higher-order Landau modes (1-3) the coefficient slot `C2` is repurposed
## as the log-periodic phase, so it is labelled "phi"; every other parameter keeps
## its own name. Keeps axis labels correct for every caller of contour_plot_sse().
#' @noRd
.lppls_display_label <- function(param, mode = 0) {
  if (identical(param, "C2") && mode %in% 1:3) "phi" else param
}

## Default sweep range for a parameter, centred on `center`, used by
## sse_surface_plot() when the caller does not supply an explicit `span`.
#' @noRd
.lppls_default_range <- function(param, center, mode = 0) {
  sym <- function(v) { h <- if (abs(v) > 0) abs(v) else 1; v + c(-h, h) }
  if (param == "tc")                        return(center + c(-90, 90))
  if (param == "m")                         return(c(0.1, 0.9))
  if (param == "omega")                     return(c(1, 20))
  if (param == "C2" && mode %in% 1:3)       return(c(1, 20))   ## phase
  sym(center)                                                  ## linear coefficients, etc.
}

## Static base-R persp() SSE surface from a contour_plot_sse() grid.
#' @noRd
.lppls_persp_sse <- function(grid, xlab, ylab, theta = -35, phi = 22, ...) {
  op <- graphics::par(mar = c(0.4, 0.4, 2.2, 0.4))   # persp reserves unused label margins
  on.exit(graphics::par(op), add = TRUE)
  Z   <- t(grid$z)                                   # grid$z is [ny, nx] -> [nx, ny] for persp
  nrz <- nrow(Z); ncz <- ncol(Z)
  zf  <- (Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz]) / 4   # facet heights
  pal <- grDevices::hcl.colors(64, "viridis")
  graphics::persp(grid$x, grid$y, Z, theta = theta, phi = phi,
                  col = pal[cut(zf, 64, labels = FALSE)], border = NA, shade = 0.35,
                  xlab = xlab, ylab = ylab, zlab = "",
                  ticktype = "detailed", nticks = 4,
                  main = sprintf("SSE w.r.t. %s and %s", xlab, ylab), ...)
}

#' Plot the SSE landscape around a set of LPPLS parameter values
#'
#' Convenience wrapper around [contour_plot_sse()] that visualises the sum of
#' squared errors (SSE) of the LPPLS model over a grid of two chosen parameters,
#' holding the rest fixed at supplied values. Axis labels follow each parameter's
#' meaning under the given `mode` (in the higher-order Landau modes the `C2` slot
#' is the log-periodic phase and is labelled `phi`). By default an interactive
#' 3-D `plotly` surface is returned (no colour bar: the colours convey depth, and
#' the z-axis already is SSE); set `static = TRUE` for a base-R `persp()` render
#' suitable for print/PDF.
#'
#' @param log_price Numeric vector of observed log-prices.
#' @param t Numeric vector of time indices matching `log_price`
#'   (default `seq_along(log_price)`).
#' @param params Named list of LPPLS parameter values, e.g.
#'   `list(A =, B =, C1 =, C2 =, tc =, m =, omega =)`, plus any mode-specific
#'   parameters. Supplies the fixed values for the parameters that are not swept,
#'   and the centre value for those that are.
#' @param vary Character vector of length 2 naming the two parameters to sweep
#'   (for example `c("m", "omega")` or `c("tc", "C2")`).
#' @param span Optional named list giving the sweep range per varied parameter.
#'   Each entry is either a length-2 numeric `c(lower, upper)` (an explicit range)
#'   or a single numeric half-width `h` (range `centre` +/- `h`). Parameters not
#'   listed fall back to a per-parameter default: `tc`, centre +/- 90; `m`,
#'   `0.1`-`0.9`; `omega`, `1`-`20`; `C2` as a phase, `1`-`20`; other (linear)
#'   parameters, centre +/- |centre|.
#' @param mode Integer LPPLS mode (`0`-`3`); default `0`.
#' @param static Logical; `FALSE` (default) returns an interactive `plotly`
#'   surface, `TRUE` draws a static base-R `persp()` surface to the current device.
#' @param ... When `static = TRUE`, extra arguments passed to [graphics::persp()]
#'   (e.g. `theta`, `phi`).
#'
#' @return For `static = FALSE`, a `plotly` surface object. For `static = TRUE`,
#'   the SSE grid (`list(x, y, z)`) returned invisibly after drawing the surface.
#'
#' @seealso [contour_plot_sse()], [SSE()], [eval_lppls()]
#'
#' @examples
#' \donttest{
#' set.seed(3630)
#' n  <- 1000
#' lp <- vapply(seq_len(n), function(i)
#'   eval_lppls(i, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
#'              tc = n + 1, m = 0.6, omega = 8) + stats::rnorm(1) * 0.05,
#'   numeric(1))
#'
#' ## SSE over (m, omega), everything else fixed -- interactive surface.
#' sse_surface_plot(
#'   log_price = lp[1:910], t = 1:910,
#'   params = list(tc = 1001, m = 0.6, omega = 8, A = 8, B = -0.025,
#'                 C1 = 0.0015, C2 = 0.001),
#'   vary = c("m", "omega"))
#'
#' ## SSE over (tc, phi) in the 1st-order Landau mode, explicit ranges, static.
#' sse_surface_plot(
#'   log_price = lp[1:910], t = 1:910,
#'   params = list(tc = 1001, m = 0.6, omega = 8, A = 8, B = -0.025, C1 = 0.0015),
#'   vary = c("tc", "C2"), span = list(tc = c(911, 1090), C2 = c(1, 20)),
#'   mode = 1, static = TRUE)
#' }
#' @export
sse_surface_plot <- function(log_price,
                             t = seq_along(log_price),
                             params,
                             vary,
                             span = NULL,
                             mode = 0,
                             static = FALSE,
                             ...) {
  if (!is.list(params) || is.null(names(params)) || any(names(params) == ""))
    stop("'params' must be a fully named list of parameter values")
  if (!is.character(vary) || length(vary) != 2L || anyNA(vary))
    stop("'vary' must be a character vector of length 2")
  if (!is.null(span) && (!is.list(span) || is.null(names(span))))
    stop("'span' must be a named list (one entry per varied parameter) or NULL")

  ## Fixed parameters = everything supplied except the two being swept.
  par <- params[setdiff(names(params), vary)]

  ## Resolve a swept axis to a (lower, upper) range. A centre value in `params`
  ## is needed only when the range is not given explicitly as c(lower, upper).
  range_for <- function(p) {
    s <- if (is.null(span)) NULL else span[[p]]
    if (!is.null(s) && length(s) == 2L) return(sort(as.numeric(s)))   # explicit range
    if (!is.null(s) && length(s) != 1L)
      stop("span[['", p, "']] must have length 1 (half-width) or 2 (range)")
    if (!(p %in% names(params)))
      stop("swept parameter '", p, "' needs either a centre value in 'params' ",
           "or an explicit span = list(", p, " = c(lower, upper))")
    ctr <- params[[p]]
    if (is.null(s)) .lppls_default_range(p, ctr, mode)   # per-parameter default range
    else            ctr + c(-1, 1) * abs(s)              # half-width around the centre
  }
  r1 <- range_for(vary[1]); r2 <- range_for(vary[2])
  lower <- c(r1[1], r2[1]); upper <- c(r1[2], r2[2])

  vars <- list(x = vary[1], y = vary[2])
  xlab <- .lppls_display_label(vary[1], mode)
  ylab <- .lppls_display_label(vary[2], mode)

  if (!static) {
    return(contour_plot_sse(log_p = log_price, t = t, par = par, vars = vars,
                            lower = lower, upper = upper,
                            cp = FALSE, sp = TRUE, mode = mode)$surface_plot)
  }

  ## Static path: get just the SSE grid (no plotly built), then draw persp().
  grid <- contour_plot_sse(log_p = log_price, t = t, par = par, vars = vars,
                           lower = lower, upper = upper,
                           cp = FALSE, sp = FALSE, mode = mode)$contour_data
  .lppls_persp_sse(grid, xlab = xlab, ylab = ylab, ...)
  invisible(grid)
}
