# Calibrate an LPPLS signal across a range of noise levels

Adds noise to a known LPPLS signal at each of several levels,
recalibrates at every level, and collects the estimates into one table.
Because the signal is known, the table shows directly how far each
estimate drifts from truth as the noise grows, and how the \\t_c\\
likelihood intervals widen with it.

## Usage

``` r
lppls_noise_sweep(
  signal,
  sd_seq,
  noise = stats::rnorm,
  fit_args = list(),
  seed = NULL,
  fb = FALSE
)
```

## Arguments

- signal:

  Numeric vector, the deterministic log-price path to perturb. Typically
  the output of
  [`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md).

- sd_seq:

  Numeric vector of noise levels. Each entry scales one draw from
  `noise`, and contributes one row to the returned table.

- noise:

  Function of one argument `n`, returning `n` standardised innovations.
  Defaults to [stats::rnorm](https://rdrr.io/r/stats/Normal.html). The
  draw is multiplied by the current entry of `sd_seq`, so `noise` should
  have unit scale for `sd_seq` to read as a standard deviation.

- fit_args:

  Named list of arguments passed to
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md).
  Use `mode = "MPL"` for the likelihood-interval columns; without it
  those columns are `NA` and
  [`noise_tc_plot()`](https://mahovo.github.io/lpplsF/reference/noise_tc_plot.md)
  has nothing to draw.

- seed:

  Optional integer. When supplied, the generator is re-seeded with it
  before *each* level, so every level sees the same innovations scaled
  differently and the comparison across levels is not confounded by the
  draw. When `NULL` (default) the levels draw independently.

- fb:

  Logical, report progress per level (default `FALSE`).

## Value

A data.frame of class `"lppls_noise_sweep"`, one row per entry of
`sd_seq`, with columns:

- sd:

  The noise level.

- tc, m, omega, A, B, C1, C2, D:

  The parameter estimates.

- value:

  The objective at the optimum.

- tc_hat_mpl:

  The modified profile likelihood estimate of \\t_c\\, `NA` unless
  `fit_args` asked for `mode = "MPL"`.

- LI5l, LI5u, ...:

  Lower and upper bounds of each likelihood interval, named for its
  cutoff as a percentage, so a cutoff of 0.05 gives `LI5l` and `LI5u`.

## Details

The noise is supplied as a function, so the experiment is not limited to
Gaussian errors. Any generator of standardised innovations will do,
which is how a heavy-tailed or serially dependent alternative is
compared against the Gaussian case on the same signal.

## See also

[`noise_param_plot()`](https://mahovo.github.io/lpplsF/reference/noise_param_plot.md)
and
[`noise_tc_plot()`](https://mahovo.github.io/lpplsF/reference/noise_tc_plot.md)
to draw the result;
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
for the calibration itself.

## Examples

``` r
# a short, coarse sweep; a real one uses more levels and a longer signal
N      <- 400
signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015,
                     C2 = 0.001, tc = N + 1, m = 0.6, omega = 8)
sweep <- lppls_noise_sweep(
  signal, sd_seq = c(0.01, 0.05),
  fit_args = list(fh = 20, hold_out = 20, num_searches = 1, mode = "F2"),
  seed = 3630
)
sweep[, c("sd", "tc", "m", "omega")]
#>     sd  tc         m    omega
#> 1 0.01 400 0.6055168 8.003865
#> 2 0.05 400 0.6151187 8.237078
```
