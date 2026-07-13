# Tests for the exported fit_plot() (P-5): the unified fit-plot builder.
# It draws the observed series + fitted LPPLS curve (clipped to tc, so no NaN
# evaluation), a tc_hat marker (always), an optional T2 marker, and a merged
# colour/linetype legend. Same builder fit_lppls() attaches as $fit_plot.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_fit <- function(tc = 320) {
  list(A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001, tc = tc, m = 0.6, omega = 8)
}

## A clean, finite observed series: critical time beyond the window, so no NaN.
make_obs <- function(n = 300) {
  t <- seq_len(n)
  list(t = t, log_p = eval_lppls(t, 8, -0.025, 0.0015, 0.001,
                                 tc = n + 50, m = 0.6, omega = 8))
}

find_layer <- function(p, class_name) {
  idx <- which(vapply(p$layers, function(l) inherits(l$geom, class_name), logical(1)))
  if (length(idx) == 0L) stop(sprintf("No layer of geom class '%s' found", class_name))
  p$layers[[idx[1]]]
}

## The fitted-curve layer: the GeomLine carrying its own data frame (the observed
## line inherits the plot data, so its layer data is a waiver).
fitted_line <- function(p) {
  lines <- Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
  hit   <- Filter(function(l) is.data.frame(l$data) && nrow(l$data) > 0, lines)
  if (length(hit) == 0L) stop("no fitted-curve GeomLine layer found")
  hit[[1]]
}

## Rendered legend labels (robust across the internal scale machinery).
legend_labels <- function(p) {
  ggplot2::get_guide_data(ggplot2::ggplot_build(p), "colour")$.label
}

# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

test_that("fit_plot returns a ggplot with an observed line and a fitted line", {
  d <- make_obs()
  p <- fit_plot(make_fit(), d$t, d$log_p)
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(geoms == "GeomLine"), 2L)   # observed + fitted
  expect_false("GeomFunction" %in% geoms)      # clipped line, not geom_function
})

# ---------------------------------------------------------------------------
# Fitted curve is clipped to tc
# ---------------------------------------------------------------------------

test_that("the fitted curve is clipped to tc (stops at floor(tc) - 1, no NaN)", {
  d <- make_obs(300); fit <- make_fit(tc = 200)   # tc inside the observed window
  fl <- fitted_line(fit_plot(fit, d$t, d$log_p))
  expect_equal(min(fl$data$ID), min(d$t))
  expect_equal(max(fl$data$ID), floor(fit$tc) - 1L)
  expect_false(any(is.na(fl$data$log_p)))
})

# ---------------------------------------------------------------------------
# tc_hat marker (always) and optional T2 marker
# ---------------------------------------------------------------------------

test_that("without T2, the only vertical marker is tc_hat at fit$tc", {
  d <- make_obs(); fit <- make_fit(tc = 320)
  vl <- find_layer(fit_plot(fit, d$t, d$log_p), "GeomVline")
  expect_equal(nrow(vl$data), 1L)
  expect_equal(vl$data$series, "tc_hat")
  expect_equal(vl$data$xintercept, fit$tc)
})

test_that("with T2 supplied, both T2 and tc_hat markers are drawn", {
  d <- make_obs(); fit <- make_fit(tc = 320)
  vl <- find_layer(fit_plot(fit, d$t, d$log_p, T2 = 250), "GeomVline")
  expect_equal(nrow(vl$data), 2L)
  expect_setequal(vl$data$series, c("T2", "tc_hat"))
  expect_equal(vl$data$xintercept[vl$data$series == "T2"], 250)
  expect_equal(vl$data$xintercept[vl$data$series == "tc_hat"], fit$tc)
})

# ---------------------------------------------------------------------------
# Legend adapts to what is drawn
# ---------------------------------------------------------------------------

test_that("legend shows sim data / LPPLS fit / tc_hat, and T2 only when supplied", {
  d <- make_obs(); fit <- make_fit()
  expect_setequal(legend_labels(fit_plot(fit, d$t, d$log_p)),
                  c("sim data", "LPPLS fit", "tc_hat"))
  expect_setequal(legend_labels(fit_plot(fit, d$t, d$log_p, T2 = 250)),
                  c("sim data", "LPPLS fit", "T2", "tc_hat"))
})

test_that("colour and linetype are merged into a single legend", {
  d <- make_obs()
  p <- fit_plot(make_fit(), d$t, d$log_p)
  cs <- p$scales$get_scales("colour")
  ls <- p$scales$get_scales("linetype")
  expect_false(is.null(cs)); expect_false(is.null(ls))
  expect_identical(cs$name, ls$name)           # same (NULL) name => one guide
})

# ---------------------------------------------------------------------------
# mode threading + backward-compatible signature
# ---------------------------------------------------------------------------

test_that("mode is forwarded to eval_lppls (different modes give different curves)", {
  d <- make_obs(); fit <- make_fit(tc = 200)
  y0 <- fitted_line(fit_plot(fit, d$t, d$log_p, mode = 0))$data$log_p
  y1 <- fitted_line(fit_plot(fit, d$t, d$log_p, mode = 1))$data$log_p
  expect_false(isTRUE(all.equal(y0, y1)))
})

test_that("positional call fit_plot(fit, t, log_p) still works (backward compatible)", {
  d <- make_obs()
  expect_s3_class(fit_plot(make_fit(), d$t, d$log_p), "ggplot")
})

# ---------------------------------------------------------------------------
# Clean build (the CRAN win: no NaN warning)
# ---------------------------------------------------------------------------

test_that("fit_plot builds silently", {
  d <- make_obs(); fit <- make_fit()
  expect_silent(ggplot2::ggplot_build(fit_plot(fit, d$t, d$log_p)))
  expect_silent(ggplot2::ggplot_build(fit_plot(fit, d$t, d$log_p, T2 = 250)))
})

# ---------------------------------------------------------------------------
# Boundary-fit detection for the MPL (P-9): .boundary_pars() flags m / omega
# pinned at an optimisation bound, so compute_mpl() can warn that the modified
# profile likelihood and its likelihood intervals are unreliable (e.g. the SPX
# m = 0.9 boundary fit).
# ---------------------------------------------------------------------------

test_that(".boundary_pars flags m / omega pinned at an optimisation bound", {
  lower <- c(0.1, 6); upper <- c(0.9, 13)
  expect_identical(.boundary_pars(list(m = 0.9, omega = 8),  lower, upper), "m")           # m at upper
  expect_identical(.boundary_pars(list(m = 0.1, omega = 8),  lower, upper), "m")           # m at lower
  expect_identical(.boundary_pars(list(m = 0.6, omega = 13), lower, upper), "omega")       # omega at upper
  expect_identical(.boundary_pars(list(m = 0.9, omega = 6),  lower, upper), c("m", "omega"))
  expect_identical(.boundary_pars(list(m = 0.6, omega = 8),  lower, upper), character(0))   # interior
})
