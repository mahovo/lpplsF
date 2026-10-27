# Changelog

## lpplsF 0.1.0

- Initial CRAN submission.

### Features

- [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
  calibrates the LPPLS model in the Filimonov–Sornette parameterisation:
  by joint optimisation (`"F1"`), by profiling over a grid of critical
  times (`"F2"`), or by modified profile likelihood, with likelihood
  intervals for the critical time (`"MPL"`).
- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`lppls_curvature()`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md)
  report the diagnostics that decide whether a calibration can be
  trusted: boundary solutions, split optimisation basins, the `B` and
  damping filters, and where the modified profile likelihood is defined.
- [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md),
  [`lppls_lagrange()`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md)
  and
  [`lr_trend_test()`](https://mahovo.github.io/lpplsF/reference/lr_trend_test.md)
  date the start of a bubble by Lagrange regularisation, with an
  autocorrelation-robust test of the trend the regularisation removes.
- [`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md)
  recalibrates a known signal across a range of noise levels, with any
  noise generator.
- [`fit_plot()`](https://mahovo.github.io/lpplsF/reference/fit_plot.md),
  [`sse_surface_plot()`](https://mahovo.github.io/lpplsF/reference/sse_surface_plot.md),
  [`sse_contour_plot()`](https://mahovo.github.io/lpplsF/reference/sse_contour_plot.md),
  [`residual_qq_plot()`](https://mahovo.github.io/lpplsF/reference/residual_qq_plot.md)
  and the rolling, Lagrange and noise plots draw the fits, the SSE
  landscapes and the diagnostics.
- `spx` holds the S&P 500 daily close from 2009-05-18 to 2020-12-31.
