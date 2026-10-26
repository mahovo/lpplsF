# Plot LPPLS Parameter Estimates Across Rolling T1

Plots the estimated LPPLS parameters \\m, \omega, A, B, C_1, C_2\\ as a
function of the rolling calibration start time \\T_1\\, one facet per
parameter. A general window-diagnostic that reproduces the
`*_params_plot` figure from the thesis.

## Usage

``` r
rolling_param_plot(x)
```

## Arguments

- x:

  An object of class `"lppls_rolling"` from
  [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md).

## Value

A `ggplot2` object.

## See also

[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md),
[`rolling_tc_plot()`](https://mahovo.github.io/lpplsF/reference/rolling_tc_plot.md)
