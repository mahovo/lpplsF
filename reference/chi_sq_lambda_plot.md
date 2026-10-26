# Plot the Lagrange-Regularized Objective with the Optimal Start Time

Plots \\\chi^2\_{\lambda^\*}(T_1)\\ against the rolling start time
\\T_1\\ with a vertical line (green) at the optimal start time
\\T_1^\*\\. Reproduces the `*_chi_sq_lambda_plot` figure from the
thesis.

## Usage

``` r
chi_sq_lambda_plot(x)
```

## Arguments

- x:

  An object of class `"lppls_lagrange"` from
  [`lppls_lagrange()`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md).

## Value

A `ggplot2` object.

## See also

[`lppls_lagrange()`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md),
[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md),
[`chi_sq_np_plot()`](https://mahovo.github.io/lpplsF/reference/chi_sq_np_plot.md)
