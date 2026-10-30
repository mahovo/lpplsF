# Plot tc Estimates and Likelihood Intervals Across Rolling T1

Plots the point estimate of the critical time \\t_c\\ (profile
likelihood, black; modified profile likelihood, cyan) with a
likelihood-interval band for each cutoff (by default 5\\ calibration
start time \\T_1\\. This is a general window-diagnostic – useful with or
without Lagrange regularization – and reproduces the `*_output_T1`
figure from the thesis. Requires the rolling calibration to have been
run in `mode = "MPL"`.

## Usage

``` r
rolling_tc_plot(x)
```

## Arguments

- x:

  An object of class `"lppls_rolling"` from
  [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md).

## Value

A `ggplot2` object.

## See also

[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md),
[`rolling_param_plot()`](https://mahovo.github.io/lpplsF/reference/rolling_param_plot.md)
