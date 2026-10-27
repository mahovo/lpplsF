# lpplsF 0.1.0

* Initial CRAN submission.

## Features

* `fit_lppls()` calibrates the LPPLS model in the Filimonov–Sornette
  parameterisation: by joint optimisation (`"F1"`), by profiling over a grid of
  critical times (`"F2"`), or by modified profile likelihood, with likelihood
  intervals for the critical time (`"MPL"`).
* `print()`, `summary()`, `plot()` and `lppls_curvature()` report the
  diagnostics that decide whether a calibration can be trusted: boundary
  solutions, split optimisation basins, the `B` and damping filters, and where
  the modified profile likelihood is defined.
* `lppls_rolling()`, `lppls_lagrange()` and `lr_trend_test()` date the start of
  a bubble by Lagrange regularisation, with an autocorrelation-robust test of
  the trend the regularisation removes.
* `lppls_noise_sweep()` recalibrates a known signal across a range of noise
  levels, with any noise generator.
* `fit_plot()`, `sse_surface_plot()`, `sse_contour_plot()`, `residual_qq_plot()`
  and the rolling, Lagrange and noise plots draw the fits, the SSE landscapes
  and the diagnostics.
* `spx` holds the S&P 500 daily close from 2009-05-18 to 2020-12-31.
