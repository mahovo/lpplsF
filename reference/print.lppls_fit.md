# Print an LPPLS calibration with its diagnostics

Prints the point estimate together with the checks that decide whether
it can be trusted: whether the solution sits at an optimisation bound,
how many of the individual fits reached the best basin found, how many
were filtered on `B` or damping, and – in `"MPL"` mode – whether the
modified profile likelihood is defined across the whole grid and whether
its likelihood intervals are closed by the likelihood or merely by the
edge of that grid.

## Usage

``` r
# S3 method for class 'lppls_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `"lppls_fit"` from
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).

- ...:

  Ignored.

## Value

`x`, invisibly.

## Details

Findings that undermine the fit are repeated as a `Diagnostics` block,
so a degenerate calibration announces itself rather than having to be
inferred from the numbers.

The report is wrapped at `getOption("width")`. A wrapped line continues
indented under the value it belongs to, so the labels stay in their
column.

## See also

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
[`summary.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/summary.lppls_fit.md),
[`lppls_curvature()`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md)
