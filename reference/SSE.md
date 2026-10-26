# Sum of Squared Errors for the LPPLS Model

Computes the residual sum of squares between observed log-prices and the
LPPLS model evaluated at a given parameter set. This is the objective
the calibration minimizes and the quantity
[`sse_contour_plot()`](https://mahovo.github.io/lpplsF/reference/sse_contour_plot.md)
scans over a parameter grid.

## Usage

``` r
SSE(par, log_p, t, mode = 0)
```

## Arguments

- par:

  Named list of LPPLS parameters. Must contain the core parameters `A`,
  `B`, `C1`, `C2`, `tc`, `m`, `omega`; may additionally contain `T1`,
  `T2`, `omega2`, `omega3` for `mode = 2` or `3`.

- log_p:

  Numeric vector of observed log-prices.

- t:

  Numeric vector of time indices, the same length as `log_p`.

- mode:

  Integer, the LPPLS model variant (`0`-`3`; default `0`).

## Value

Numeric scalar, the sum of squared errors.

## Details

`par` is a named list of LPPLS parameters. The seven core parameters
`A`, `B`, `C1`, `C2`, `tc`, `m`, `omega` are required. For the
higher-order variants (`mode = 2` or `3`), `par` may also carry the
scale and secondary-frequency parameters `T1`, `T2`, `omega2`, `omega3`;
each is forwarded to
[`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md)
when present and otherwise falls back to that function's default
(`T1 = 500`, `T2 = 1990`, `omega2 = 0`, `omega3 = 0`). Requesting a
higher-order mode without supplying the parameters it uses issues a
warning and falls back to those defaults, since silently defaulting them
can give a misleading SSE.

Elements of `par` not recognised by
[`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md)
are ignored. Time points at or beyond `tc` produce `NaN` from
[`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md)
and are dropped from the sum (`na.rm = TRUE`).

## See also

[`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md),
[`sse_contour_plot()`](https://mahovo.github.io/lpplsF/reference/sse_contour_plot.md)

## Examples

``` r
t <- 1:100
par <- list(A = 4, B = -0.015, C1 = 0.0015, C2 = 0.001,
            tc = 150, m = 0.5, omega = 9)
log_p <- eval_lppls(t, par$A, par$B, par$C1, par$C2,
                    par$tc, par$m, par$omega) + rnorm(100, 0, 0.01)
SSE(par, log_p, t)
#> [1] 0.01098369

# Higher-order variant (mode 2): supply T1 and omega2 alongside the core set
par2 <- c(par, list(T1 = 300, omega2 = 5))
SSE(par2, log_p, t, mode = 2)
#> [1] 0.07333604
```
