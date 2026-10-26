# Critical-time estimate and likelihood intervals against noise level

Draws the \\t_c\\ point estimate from a noise sweep, together with the
modified profile likelihood estimate and the nested likelihood-interval
bands, as a function of the noise level. The bands widen as the noise
grows, which is the honest reading of how much the critical time can be
pinned down.

## Usage

``` r
noise_tc_plot(x, true_tc = NULL)
```

## Arguments

- x:

  A `"lppls_noise_sweep"` table from
  [`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md),
  or any data.frame with an `sd` column, a `tc` column and paired
  `LI*l`/`LI*u` columns.

- true_tc:

  Optional numeric. When supplied, a dashed horizontal reference line is
  drawn there, which is useful on a simulated signal whose critical time
  is known.

## Value

A `ggplot2` object.

## See also

[`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md),
[`noise_param_plot()`](https://mahovo.github.io/lpplsF/reference/noise_param_plot.md)

## Examples

``` r
N      <- 400
signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
                     C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
sweep <- lppls_noise_sweep(
  signal, sd_seq = c(0.01, 0.03, 0.05),
  fit_args = list(fh = 20, hold_out = 20, num_searches = 1,
                  mode = "MPL", mpl_cutoff = c(0.05, 0.5)),
  seed = 3630
)
noise_tc_plot(sweep, true_tc = N + 1)

```
