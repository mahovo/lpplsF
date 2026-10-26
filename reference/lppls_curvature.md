# Curvature of the LPPLS objective at the estimate

Tests whether the calibration sits at a genuine interior minimum.

## Usage

``` r
lppls_curvature(x, log_price = NULL)
```

## Arguments

- x:

  An object of class `"lppls_fit"` from
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).

- log_price:

  Optional numeric vector, the series the fit was calibrated on. Only
  needed for fits made before
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
  stored it.

## Value

An object of class `"lppls_curvature"`: a list with the full and profile
determinants, the profile Hessian and its eigenvalues, the active
bounds, and a `classification` of `"minimum"`, `"saddle"` or
`"maximum"`.

## Details

The modified profile likelihood assumes the estimate is an interior
minimum: its correction uses \\\det(X'X - H)\\, the observed
information, which is only guaranteed positive at an interior stationary
point. \\X\\ holds the model's first derivatives with respect to the six
parameters, and \\H\\ its second derivatives, each weighted by the
residual at that observation and summed (thesis eq. 55). \\H\\ is thus a
residual-weighted curvature term rather than a Hessian in its own right.
`lppls_curvature()` tests the determinant directly, and factorises it to
show where the sign comes from.

By the Schur complement, \\\det(X'X - H)\\ factorises into the
determinant of the linear-parameter block times the determinant of the
\\2 \times 2\\ Hessian of the \\(m, \omega)\\ profile. The model is
linear in \\A, B, C_1, C_2\\, so their second derivatives vanish, \\H\\
contributes nothing to the linear block, and what remains there is a
crossproduct with a positive determinant. The sign is therefore decided
entirely by curvature in the two nonlinear parameters: the \\2 \times
2\\ block the function reports.

## See also

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md)
