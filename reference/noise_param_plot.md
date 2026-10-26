# Parameter estimates against noise level

Draws each estimated parameter from a noise sweep as a function of the
noise level, one panel per parameter. The log-periodic parameters are
typically the first to degrade, which is the point of looking.

## Usage

``` r
noise_param_plot(x, params = c("m", "omega", "A", "B", "C1", "C2"))
```

## Arguments

- x:

  A `"lppls_noise_sweep"` table from
  [`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md),
  or any data.frame with an `sd` column and the parameter columns.

- params:

  Character vector naming the parameters to draw. Defaults to all six
  LPPLS parameters.

## Value

A `ggplot2` object.

## See also

[`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md),
[`noise_tc_plot()`](https://mahovo.github.io/lpplsF/reference/noise_tc_plot.md)

## Examples

``` r
N      <- 400
signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
                     C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
sweep <- lppls_noise_sweep(
  signal, sd_seq = c(0.01, 0.03, 0.05),
  fit_args = list(fh = 20, hold_out = 20, num_searches = 1, mode = "F2"),
  seed = 3630
)
noise_param_plot(sweep)

```
