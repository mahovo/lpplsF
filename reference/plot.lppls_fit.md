# Plot an LPPLS calibration

Returns one of the plots attached to the fit, or draws one on demand.
The attached ones exist only if they were asked for at fit time (`fp`,
`pp`, `mp`, `tp`, `mpl_plot`, `cp`, `sp`). `"basin"`, `"param_basin"`
and `"mpl_relative"` are drawn from what the fit stores anyway, so they
need nothing asked for in advance.

## Usage

``` r
# S3 method for class 'lppls_fit'
plot(
  x,
  which = c("basin", "param_basin", "fit", "mpl", "mpl_relative", "param", "matrix",
    "trace", "contour", "surface"),
  ...
)
```

## Arguments

- x:

  An object of class `"lppls_fit"` from
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).

- which:

  Which plot to return. `"basin"`, `"param_basin"` and `"mpl_relative"`
  are drawn on demand; the rest are returned from the fit.

- ...:

  Ignored.

## Value

A `ggplot2` object (a `plotly` object for `"contour"`/`"surface"`).

## Details

`"basin"` plots each fit's SSE against its `tc`, split at the gap
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md)
reports. When a search has failed on part of the grid the two clusters
separate visibly, which is invisible in the point estimate.

`"param_basin"` is the per-`tc` parameter plot with the same split
applied to every parameter, not just the objective. It answers what
`"basin"` cannot: *which* parameters the inferior fits differ in –
typically `omega` jumping between local optima while `m` stays pinned at
its bound.

`"mpl_relative"` draws the modified profile likelihood relative to its
maximum, \\R = \log \mathrm{MPL} - \max \log \mathrm{MPL}\\, zoomed to
the band the likelihood-interval cutoffs live in, with each cutoff drawn
as a horizontal line. The absolute level of the log MPL carries no
information, whereas \\R\\ is what the intervals test, so each interval
can be read off where the curve crosses its cutoff. It needs a fit made
with `mode = "MPL"`.

## See also

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md)
