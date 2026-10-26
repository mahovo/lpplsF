# Autocorrelation-robust significance test for the Lagrange \\\chi^2\_{np}\\ trend

Assesses whether the decreasing "bias" trend that Lagrange
regularization acts on (the slope \\\beta_1\\ of \\\chi^2\_{np}(T_1)\\
against the window start \\T_1\\) is statistically distinguishable from
what the estimation noise produces on its own. Three tests are reported:
a naive OLS \\p\\-value (invalid here, because adjacent rolling windows
overlap and \\\chi^2\_{np}(T_1)\\ is heavily autocorrelated); a
heteroskedasticity- and autocorrelation-consistent (Newey–West)
\\p\\-value; and a **parametric AR(1)–GARCH(1,1) bootstrap**, which is
the appropriate test.

## Usage

``` r
lr_trend_test(rolling, log_price, B = 199, cores = 1L, hac = TRUE, seed = NULL)
```

## Arguments

- rolling:

  An `"lppls_rolling"` object from
  [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md);
  it must carry the `fit_args` slot (rebuild with a current version of
  the package if absent).

- log_price:

  Numeric vector of log-prices the `rolling` was built on; its length
  must equal `rolling$n_full`.

- B:

  Integer, number of bootstrap replicates (default `199`).

- cores:

  Integer, number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`; ignored on Windows).

- hac:

  Logical, whether to also compute the Newey–West HAC \\p\\-value
  (default `TRUE`).

- seed:

  Integer or `NULL`; a base seed that makes the bootstrap paths
  reproducible.

## Value

An object of class `"lr_trend_test"`: a list with the observed slope
`beta1` and `lambda` (\\= -\beta_1\\); the naive `p_ols` / `se_ols`; the
HAC `p_hac` / `se_hac`; the bootstrap `p_boot`, standardized distance
`z` and null vector `null`; the fitted noise model (`ar1`, `garch` =
`c(omega, alpha1, beta1)`, and `persist` \\= \alpha + \beta\\) with its
`valid` flag; the effective replicate count `B`, the window count
`n_windows`, and the `call`.

## Details

The bootstrap fits an AR(1)–GARCH(1,1) model (via
[`stats::arima()`](https://rdrr.io/r/stats/arima.html) and the tsgarch
package) to the residuals of the full-window (\\T_1 = 1\\) LPPLS fit,
simulates fresh noise paths, adds them to the fitted LPPLS curve, and
repeats the rolling calibration to build the null distribution of
\\\beta_1\\. The observed slope is compared against this null
(two-sided). A validity guard is applied: if the fitted GARCH
persistence \\\alpha + \beta \ge 1\\ (IGARCH) the single-LPPLS null is
misspecified, so the bootstrap is skipped, a warning is issued and
`p_boot` is `NA`.

Requires the suggested packages tsgarch, xts and (for `hac = TRUE`)
sandwich.

## See also

[`lppls_lagrange()`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md),
[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md)

## Examples

``` r
if (FALSE) { # \dontrun{
roll <- lppls_rolling(log_p, t1_seq = seq(1, 1001, 25), hold_out = 100, fh = 150,
                      fit_args = list(mode = "F2", num_searches = 1))
lr_trend_test(roll, log_price = log_p, B = 199, cores = 4, seed = 1)
} # }
```
