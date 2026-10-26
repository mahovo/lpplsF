# Full diagnostic summary of an LPPLS calibration

Everything
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md)
reports, plus the curvature check at the estimate when the fit carries
(or is given) its input series. The returned object holds the
diagnostics as values, so they can be used programmatically rather than
parsed out of printed text.

## Usage

``` r
# S3 method for class 'lppls_fit'
summary(object, log_price = NULL, ...)
```

## Arguments

- object:

  An object of class `"lppls_fit"` from
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).

- log_price:

  Optional numeric vector passed to
  [`lppls_curvature()`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md),
  for fits made before
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
  stored the series.

- ...:

  Ignored.

## Value

An object of class `"summary.lppls_fit"`.

## See also

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
[`lppls_curvature()`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md)
