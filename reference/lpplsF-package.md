# lpplsF: Log-Periodic Power Law Singularity Model for Financial Bubbles

Implementation of the LPPLS model for detecting and analysing financial
bubbles, based on the calibration methods of Filimonov & Sornette (2013)
and the Modified Profile Likelihood inference from Filimonov, Demos &
Sornette (2017).

## Fitting

- [`fit_lppls`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md):

  Calibrate the model to a log-price series

- [`eval_lppls`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md):

  Evaluate the model at given parameters

- [`SSE`](https://mahovo.github.io/lpplsF/reference/SSE.md):

  The objective the calibration minimises

## Diagnosing a calibration

A calibration can look convincing as a plot and still not be worth
interpreting: the estimate may rest on an optimisation bound, the
per-`tc` searches may have settled in different local optima, and the
modified profile likelihood may be undefined over part of the `tc` grid.
Printing the object returned by
[`fit_lppls`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
reports on all of this, and the rest go further.

- [`print`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md):

  The estimate together with the checks that decide whether to trust it;
  warnings are listed only when something is wrong, so a clean
  calibration reports none

- [`summary`](https://mahovo.github.io/lpplsF/reference/summary.lppls_fit.md):

  The same diagnostics as values, plus the curvature check and whether
  each likelihood interval is closed by the likelihood or by the edge of
  the grid

- [`plot`](https://mahovo.github.io/lpplsF/reference/plot.lppls_fit.md):

  Any plot attached at fit time, or the basin views, which show whether
  the per-`tc` searches split into a good and an inferior cluster

- [`lppls_curvature`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md):

  Whether the estimate sits at a genuine interior minimum, which is what
  the modified profile likelihood assumes

- [`residual_qq_plot`](https://mahovo.github.io/lpplsF/reference/residual_qq_plot.md):

  Calibration residuals against a normal or fitted Student-t reference

## Rolling windows and Lagrange regularization

Each of these returns an object that summarises itself when printed.

- [`lppls_rolling`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md):

  Recalibrate over a rolling start time \\T_1\\, holding the end of the
  window fixed

- [`lppls_lagrange`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md):

  Lagrange regularization of a rolling calibration, locating the onset
  \\T_1^\*\\

- [`lr_trend_test`](https://mahovo.github.io/lpplsF/reference/lr_trend_test.md):

  Whether the fitted bias trend is distinguishable from noise, by
  parametric bootstrap

- [`rolling_tc_plot`](https://mahovo.github.io/lpplsF/reference/rolling_tc_plot.md):

  \\\hat{t}\_c\\ and its likelihood intervals against \\T_1\\

- [`rolling_param_plot`](https://mahovo.github.io/lpplsF/reference/rolling_param_plot.md):

  Each estimated parameter against \\T_1\\

- [`chi_sq_np_plot`](https://mahovo.github.io/lpplsF/reference/chi_sq_np_plot.md),
  [`chi_sq_lambda_plot`](https://mahovo.github.io/lpplsF/reference/chi_sq_lambda_plot.md):

  The Lagrange diagnostics \\\chi^2\_{np}(T_1)\\ and
  \\\chi^2\_{\lambda}(T_1)\\

## Visualising the model and its landscape

- [`fit_plot`](https://mahovo.github.io/lpplsF/reference/fit_plot.md):

  A calibrated series with its fitted LPPLS curve

- [`sse_contour_plot`](https://mahovo.github.io/lpplsF/reference/sse_contour_plot.md):

  SSE over a grid of two parameters

- [`sse_surface_plot`](https://mahovo.github.io/lpplsF/reference/sse_surface_plot.md):

  The same landscape as an interactive or a static surface

## Model Description

The LPPLS model describes the log-price dynamics of financial assets
approaching a critical point (bubble peak or crash):

\$\$y(t) = A + (t_c - t)^m \[B + C_1 \cos(\omega \log(t_c - t)) + C_2
\sin(\omega \log(t_c - t))\]\$\$

where:

- \\t_c\\ is the critical time

- \\m\\ is the power law exponent (typically 0.1-0.9)

- \\\omega\\ is the log-periodic frequency (typically 6-13)

- \\A, B, C_1, C_2\\ are linear parameters

## Optimization Modes

- F1:

  Simultaneous optimization of all nonlinear parameters

- F2:

  Two-stage: optimize (m, omega) for each tc, then select best tc

- MPL:

  Modified profile likelihood, with likelihood intervals for tc

## References

Filimonov, V., & Sornette, D. (2013). A stable and robust calibration
scheme of the log-periodic power law model. Physica A, 392(17),
3698-3707.

Filimonov, V., Demos, G., & Sornette, D. (2017). Modified profile
likelihood inference and interval forecast of the burst of financial
bubbles. Quantitative Finance, 17(8), 1167-1186.

## See also

Useful links:

- <https://mahovo.github.io/lpplsF/>

- <https://github.com/mahovo/lpplsF>

- Report bugs at <https://github.com/mahovo/lpplsF/issues>

## Author

**Maintainer**: Martin Hoshi Vognsen <martinvognsen@gmail.com>
([ORCID](https://orcid.org/0000-0003-4460-7763))
