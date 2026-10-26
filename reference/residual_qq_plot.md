# Residual QQ Plot Diagnostic

Quantile-quantile plot of calibration residuals against a normal or a
fitted Student-t reference distribution, with reference lines at chosen
quantiles. A general goodness-of-fit diagnostic for the residuals of an
LPPLS calibration.

## Usage

``` r
residual_qq_plot(residuals, distribution = c("normal", "t"), ref = c(0.1, 0.9))
```

## Arguments

- residuals:

  Numeric vector of residuals (length \>= 2).

- distribution:

  Reference distribution for the theoretical quantiles: `"normal"` (the
  default) or `"t"`, a Student-t distribution whose degrees of freedom
  are fitted to `residuals` with
  [`MASS::fitdistr()`](https://rdrr.io/pkg/MASS/man/fitdistr.html).

- ref:

  Numeric probabilities at which to draw the green reference lines
  (default `c(0.1, 0.9)`): horizontal lines at the corresponding sample
  quantiles of `residuals`, vertical lines at the theoretical quantiles.

## Value

A `ggplot2` object.

## See also

[`chi_sq_np_plot()`](https://mahovo.github.io/lpplsF/reference/chi_sq_np_plot.md),
[`rolling_tc_plot()`](https://mahovo.github.io/lpplsF/reference/rolling_tc_plot.md)

## Examples

``` r
set.seed(1)
residual_qq_plot(rnorm(200))

if (requireNamespace("MASS", quietly = TRUE))
  residual_qq_plot(rt(200, df = 4), distribution = "t")
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
```
