# Tests for the trace plot's out-of-range regions.
#
# B and the damping D are not search coordinates: at each point of the
# (m, omega) surface they follow from the least-squares solution for the linear
# parameters. The trace plot shades where that solution would be out of range,
# B >= upper[3] or D < lower[4]. These tests drive create_trace_plot() directly,
# with bounds chosen to make each region present or absent.

## A small series, its objective and a two-point path.
trace_setup <- function() {
  set.seed(1)
  n     <- 200L
  log_p <- eval_lppls(seq_len(n), 8, -0.025, 0.0015, 0.001, 220, 0.6, 8) +
    rnorm(n) * 0.01
  bc <- create_beta_calculator()
  sc <- create_sse_calculator("crossprod", bc)
  list(log_p = log_p, t = seq_len(n), tc = 220, beta = bc,
       sse2 = function(par, tc, log_p, t) sc(log_p, t, tc, par[1], par[2]),
       path = matrix(c(0.3, 10, 0.6, 8), nrow = 2))
}

draw <- function(s, lower, upper) {
  create_trace_plot(s$path, lower, upper, s$sse2, s$beta, s$tc, s$log_p, s$t)
}

## Layers whose geom is exactly `cls` (GeomContour inherits from GeomPath).
layers_of <- function(p, cls) {
  Filter(function(l) class(l$geom)[1] == cls, p$layers)
}

## Bounds that put nothing out of range.
none_lower <- c(0.1, 6, -1, -Inf)
none_upper <- c(0.9, 13, Inf, Inf)


test_that("the shaded regions are exactly where B or D would be out of range", {
  skip_on_cran()

  s <- trace_setup()
  L <- draw(s, none_lower, none_upper)$data

  ## B and D at every lattice point, from least squares done here
  b <- t(mapply(function(m, w) s$beta(s$log_p, s$t, s$tc, m, w), L$m, L$omega))
  B <- b[, "B"]
  D <- L$m * abs(B) / (L$omega * sqrt(b[, "C1"]^2 + b[, "C2"]^2))

  ## Thresholds at the medians, so each region covers about half the lattice
  ub <- stats::median(B)
  ld <- stats::median(D)
  p  <- draw(s, c(0.1, 6, -1, ld), c(0.9, 13, ub, Inf))

  tiles <- layers_of(p, "GeomTile")
  expect_length(tiles, 2L)
  fills <- vapply(tiles, function(l) l$aes_params$fill, "")
  key   <- function(d) paste(d$m, d$omega)
  expect_setequal(key(tiles[[which(fills == "black")]]$data), key(L[B >= ub, ]))
  expect_setequal(key(tiles[[which(fills == "white")]]$data), key(L[D < ld, ]))

  ## Each region is edged at its threshold, and the edges carry the legend.
  expect_length(layers_of(p, "GeomContour"), 3L)   # the SSE contours + 2 edges
  built <- ggplot2::ggplot_build(p)
  expect_setequal(built$plot$scales$get_scales("colour")$get_limits(),
                  c("B out of range", "D out of range"))
})


test_that("nothing out of range draws no shading, no edges and no legend", {
  skip_on_cran()

  p <- draw(trace_setup(), none_lower, none_upper)
  expect_length(layers_of(p, "GeomTile"), 0L)
  expect_length(layers_of(p, "GeomContour"), 1L)   # only the SSE contours
  expect_null(p$scales$get_scales("colour"))
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("the shading leaves the surface and the path as they were", {
  skip_on_cran()

  s  <- trace_setup()
  p0 <- draw(s, none_lower, none_upper)
  L  <- p0$data
  p  <- draw(s, c(0.1, 6, -1, stats::median(L$D)),
             c(0.9, 13, stats::median(L$B), Inf))

  built_layer <- function(p, cls) {
    i <- which(vapply(p$layers, function(l) class(l$geom)[1] == cls, logical(1)))[1]
    ggplot2::ggplot_build(p)$data[[i]]
  }
  for (cls in c("GeomRaster", "GeomPath")) {
    expect_identical(built_layer(p, cls), built_layer(p0, cls))
  }
})


test_that("tp is a single TRUE or FALSE", {
  expect_error(fit_lppls(log_price = trace_setup()$log_p, tp = c(1, 0, 0)),
               "'tp' must be TRUE or FALSE")
  expect_error(fit_lppls(log_price = trace_setup()$log_p, tp = NA),
               "'tp' must be TRUE or FALSE")
})


test_that("plot() returns the trace plot as type \"trace\"", {
  skip_on_cran()

  s   <- trace_setup()
  fit <- suppressWarnings(fit_lppls(
    log_price = s$log_p, fh = 30, lower = c(0.1, 6, -0.03, 0.8),
    upper = c(0.9, 13, -1e-14, 1e6), m_init = 0.3, o_init = 10,
    num_searches = 1, mode = "F2", tp = TRUE))
  expect_identical(plot(fit, "trace"), fit$trace_plot)
  expect_error(plot(fit, "trace_bm"))
})
