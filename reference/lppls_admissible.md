# Check LPPLS Estimates Against the Admissibility Conditions

A calibration is admissible when its estimate meets the conditions
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
works with: `m` and `omega` inside their bounds, `B` below `upper[3]`,
and the damping `D` at least `lower[4]` (Filimonov, Demos & Sornette,
2017).
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
bounds its search in `m` and `omega`, but checks `B` and `D` only on the
result of each search, and
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md)
reports how many results were out of range. `lppls_admissible()` applies
the same conditions to any set of estimates: the estimate of a fit, its
per-`tc` fits, every window of a rolling calibration, or a data frame of
estimates.

## Usage

``` r
lppls_admissible(x, ...)

# S3 method for class 'lppls_fit'
lppls_admissible(
  x,
  which = c("estimate", "fits"),
  lower = NULL,
  upper = NULL,
  ...
)

# S3 method for class 'lppls_rolling'
lppls_admissible(x, lower = NULL, upper = NULL, ...)

# Default S3 method
lppls_admissible(x, lower = NULL, upper = NULL, ...)
```

## Arguments

- x:

  An `"lppls_fit"` from
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
  an `"lppls_rolling"` from
  [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md),
  or estimates with elements `m`, `omega`, `B` and `D`: a data frame
  with one row per estimate, or a list holding one estimate, such as
  `fit$fit[[1]]`.

- ...:

  Not used.

- which:

  For an `"lppls_fit"`: `"estimate"` (the default) checks the reported
  estimate, `"fits"` every per-`tc` fit (every restart in mode `"F1"`),
  including those the `B` filter set aside.

- lower, upper:

  Numeric vectors of length 4, bounds for `c(m, omega, B, D)` as in
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).
  By default the bounds the calibration was run with; for a data frame
  or a list, the defaults of
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).

## Value

The estimates checked, as a data frame with one row per estimate and
five logical columns added: `m_ok`, `omega_ok`, `B_ok`, `D_ok`, and
`admissible`, which is `TRUE` when all four are. A missing value counts
as out of range.

## Details

The conditions are those of
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md):
`lower[1] <= m <= upper[1]`, `lower[2] <= omega <= upper[2]`,
`B < upper[3]` and `D >= lower[4]`. `lower[3]` and `upper[4]` are not
used.

The share of admissible windows in a rolling calibration,
`mean(lppls_admissible(roll)$admissible)`, measures how consistently a
bubble signal holds up as the window start moves.

## References

Filimonov, V., Demos, G., & Sornette, D. (2017). Modified profile
likelihood inference and interval forecast of the burst of financial
bubbles. *Quantitative Finance*, 17(8), 1167-1186.
[doi:10.1080/14697688.2016.1276298](https://doi.org/10.1080/14697688.2016.1276298)

## See also

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
[`print.lppls_fit()`](https://mahovo.github.io/lpplsF/reference/print.lppls_fit.md),
[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md)

## Examples

``` r
set.seed(4)
log_p <- eval_lppls(1:200, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                    tc = 221, m = 0.6, omega = 8) + rnorm(200, sd = 0.01)
fit <- fit_lppls(log_p, fh = 30, mode = "MPL",
                 m_init = 0.5, o_init = 8, num_searches = 2)

## The reported estimate
lppls_admissible(fit)
#>         tc         m    omega        A          B            C1         C2
#> 1 222.9766 0.5937427 8.258585 8.010661 -0.0261428 -0.0006679849 0.00169013
#>          D     value m_ok omega_ok B_ok D_ok admissible
#> 1 1.034206 0.0178573 TRUE     TRUE TRUE TRUE       TRUE

## Every per-tc fit, against a stricter damping floor of 1
adm <- lppls_admissible(fit, which = "fits", lower = c(0.1, 6, -1e14, 1))
mean(adm$admissible)
#> [1] 1
```
