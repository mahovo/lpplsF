# The \`lpplsF\` package: A guided tour

`lpplsF` fits the **log-periodic power-law singularity (LPPLS)** model
of Sornette and co-authors — a model of the faster-than-exponential,
log-periodically decorated price growth that can precede a
financial-bubble crash. This document is a **guided tour of the
package**: every calibration, diagnostic and plot it offers, worked
through on a series of examples that build on one another.

**How to read it.** Each section builds on the previous one, from a
single calibration on simulated data, through the optimizer’s behaviour
and its sensitivity to noise, to rolling-window bubble-onset detection
and finally a real market series. You can read it front to back, or jump
to the section matching the function you need.

**Theory, and where to find it.** This guide is about *using* the
package; it assumes the model rather than deriving it. The accompanying
master’s thesis sets out the theory behind every method used here: the
model’s derivation, the calibration scheme, the modified profile
likelihood and the Lagrange regularization. Each section below names the
thesis figure it corresponds to, so you can read the two side by side.

**Reproducibility.** The single-model calibrations run live every time
the document is knitted. The multi-fit noise and rolling-window tables
(tens of minutes to hours of CPU) are shipped precomputed and loaded
instead. The specification for each is given in its own section (see
[Reproducibility and the cache](#reproducibility-and-the-cache)), and
the dual-regime rolling calibration can be re-run in place by knitting
with

``` r

rmarkdown::render("lpplsF_guide.Rmd", params = list(recompute = TRUE))
```

**Interactive surfaces.** The SSE landscapes are interactive `plotly`
surfaces (drag to rotate). `sse_surface_plot(..., static = TRUE)` draws
the same grids as static base-R
[`persp()`](https://rdrr.io/r/graphics/persp.html) surfaces, which is
how the package’s PDF overview,
[`vignette("lpplsF")`](https://mahovo.github.io/lpplsF/articles/lpplsF.md),
shows them.

## The model and its calibration landscape

The core generator is
[`eval_lppls()`](https://mahovo.github.io/lpplsF/reference/eval_lppls.md).
A single LPPLS trajectory is a deterministic function of a critical time
\\t_c\\, the power-law exponent \\m\\, the log-periodic angular
frequency \\\omega\\, and linear coefficients \\A, B, C_1, C_2\\:

``` r

## ten points spread across a series whose critical time is tc = 1001
round(eval_lppls(seq(1, 1000, length.out = 10),
                 A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                 tc = 1001, m = 0.6, omega = 8), 4)
```

    ##  [1] 6.3886 6.4297 6.5710 6.7978 6.9676 6.9962 7.1630 7.3707 7.6021 7.9765

Calibration means minimizing the sum of squared errors (SSE) between the
model and observed log-prices over the nonlinear parameters \\(t_c, m,
\omega)\\; the linear parameters are profiled out. Before fitting, it
helps to *see* that SSE landscape.
[`sse_surface_plot()`](https://mahovo.github.io/lpplsF/reference/sse_surface_plot.md)
evaluates the SSE over a grid of two parameters, holding the rest fixed,
and draws that surface. Give the values of the parameters you are not
scanning, name the two you are, and give each a range:

``` r

sse_surface_plot(
  log_p, t_obs,
  params = list(tc = 1001, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001),
  vary   = c("m", "omega"),
  span   = list(m = c(0.1, 0.9), omega = c(6, 13)),
  mode   = 0
)
```

Every surface in this section is produced that way.
[`sse_contour_plot()`](https://mahovo.github.io/lpplsF/reference/sse_contour_plot.md),
which it wraps, returns the grid itself and can draw a contour plot of
it as well.

### Simulated LPPLS(0) data

*Thesis Figure 10.* A single LPPLS(0) process (the Filimonov–Sornette
parameterization, `mode = 0`) of length \\N = 1000\\ plus Gaussian noise
(\\\sigma = 0.05\\); the SSE is evaluated on the first 910 observations
(the critical time \\t_c = 1001\\ lies beyond the fitting window).

| \\t_c\\ | \\m\\ | \\\omega\\ | \\A\\ |  \\B\\ | \\C_1\\ | \\C_2\\ |
|--------:|------:|-----------:|------:|-------:|--------:|--------:|
|    1001 |   0.6 |          8 |     8 | -0.025 |  0.0015 |   0.001 |

``` r

set.seed(3630)
sigma <- 0.05
N     <- 1000
tc    <- N + 1                      # 1001
sim   <- vapply(1:N, function(i) {
  eval_lppls(i, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
             tc = tc, m = 0.6, omega = 8) + rnorm(1) * sigma
}, numeric(1))

log_p <- sim[1:910]                 # fitting window
t_obs <- (1:N)[1:910]
```

#### SSE over \\m\\ and \\\omega\\ — thesis Figure 10(a)

Flat/low for small \\m\\, a weak depression near \\m = 0.6\\, then a
steep rise for larger \\m\\; only small oscillating variation along
\\\omega\\.

#### SSE over \\t_c\\ and \\m\\ — thesis Figure 10(b)

Flat for small \\m\\, rising strongly for larger \\m\\; the variation
along \\t_c\\ is small for small \\m\\.

#### SSE over \\t_c\\ and \\\omega\\ — thesis Figure 10(c)

Many strong oscillations along the \\\omega\\ direction — the source of
the multi-modal, hard-to-optimize likelihood that motivates the profile
methods below.

### Higher-order Landau variants

*Thesis Figure 11.* The package also implements the 1st- and 2nd-order
Landau expansions (`mode = 1`, `mode = 2`). In those modes \\C_2\\ plays
the role of the phase \\\phi\\, so a “\\\phi\\ axis” is produced by
scanning \\C_2\\. The axes below are labelled \\\phi\\ accordingly.

#### \\t_c\\ and \\\phi\\, 1st-order Landau — thesis Figure 11(a)

This surface uses the **\\N = 1000\\ mode-0 data from the previous
section**, examined under the mode-1 \\(t_c, \phi)\\ landscape — modes 0
and 1 are equivalent parameterizations of the same 1st-order Landau
LPPL. The phase periodicity (\\2\pi\\) produces three local minima
across \\\phi \in \[1, 20\]\\, plus a boundary solution at small
\\\phi\\.

#### A 2nd-order Landau simulation

The 2nd-order model adds a pre-critical regime of length \\\Delta t =
500\\ with secondary frequency \\\Delta\omega = 0\\. SSE is evaluated on
the first 1910 of \\N = 2000\\ observations.

| \\t_c\\ | \\m\\ | \\\omega\\ | \\A\\ |  \\B\\ | \\C_1\\ | \\C_2\\ | \\\Delta t\\ |
|--------:|------:|-----------:|------:|-------:|--------:|--------:|-------------:|
|    2001 |   0.6 |          8 |     8 | -0.025 |  0.0015 |   0.001 |          500 |

``` r

set.seed(3630)
N2  <- 2000
tc2 <- N2 + 1
sim2 <- vapply(1:N2, function(i) {
  eval_lppls(i, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
             tc = tc2, m = 0.6, omega = 8, mode = 2, T1 = 500) + rnorm(1) * sigma
}, numeric(1))
log_p2 <- sim2[1:1910]
t_obs2 <- (1:N2)[1:1910]
```

#### \\t_c\\ and \\\phi\\, 2nd-order Landau — thesis Figure 11(b)

Note the wider \\t_c\\ axis: the pre-critical regime lets parameters
wander outside the usual LPPLS filter values.

#### \\t_c\\ and \\\omega\\, 2nd-order Landau — thesis Figure 11(c)

Even more rugged than the Filimonov variant: narrow, deep ravines and
strong oscillations along \\\omega\\.

## Fitting a single model

[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
is the workhorse. With `mode = "MPL"` it calibrates by **modified
profile likelihood**: for each candidate \\t_c\\ in the search window it
profiles out \\(m, \omega)\\ (and the linear parameters), then builds a
likelihood interval for \\t_c\\. We calibrate two simulated series: one
with Gaussian noise (`M1`), one with AR(1)-GARCH(1,1) noise (`M2`). Each
holds out the last 200 points, fits on \\\[1, 1800\]\\, and searches
\\t_c\\ over \\\[1801, 2100\]\\; the true \\t_c = 2001\\.

``` r

## One-row tidy estimate table from a fit_lppls() result.
est_row <- function(fit) {
  e <- fit$fit[[1]]
  data.frame(tc = e$tc, m = e$m, omega = e$omega,
             A = unname(e$A), B = unname(e$B),
             C1 = unname(e$C1), C2 = unname(e$C2), D = unname(e$D))
}

## The calibration shared by M1 and M2. fp/pp/tp ask fit_lppls() to also attach
## the fit plot, the per-tc parameter plot and the BFGS trace plot to the
## returned object, so they are available later without refitting.
fit_call <- function(log_price) {
  fit_lppls(
    log_price = log_price, fh = 300, hold_out = 200,
    lower = c(0.1, 6, -0.03, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
    tc_init = 1950, m_init = 0.2, o_init = 12, num_searches = 1,
    mode = "MPL", mpl_cutoff = c(0.05, 0.1, 0.5), beta_method = "crossprod",
    fp = TRUE, pp = TRUE, tp = TRUE
  )
}
```

The Gaussian series is regenerated from `set.seed(3630)`. The GARCH
series ships with this guide as `garch_series.rds`: it was drawn with
`rugarch`, whose simulated paths depend on the package version, so it is
stored rather than re-simulated.

``` r

set.seed(3630)
Nf <- 2000
gauss_series <- vapply(1:Nf, function(i) {
  eval_lppls(i, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
             tc = Nf + 1, m = 0.6, omega = 8) + rnorm(1) * sigma
}, numeric(1))
garch_series <- readRDS(dpath("garch_series.rds"))$log_p
```

``` r

## 300 profile-likelihood optimizations per series
fit_gauss <- fit_call(gauss_series)
fit_garch <- fit_call(garch_series)
```

### The fit plot

In each plot the blue line is the simulated series, the red curve the
estimated LPPLS process (up to the estimated \\t_c\\), the green dashed
line marks \\T_2 = 1800\\, and the red dashed line the estimated
\\t_c\\. *Thesis Figure 13; the estimates are those of thesis Table 3.*

#### M1 — Gaussian noise — thesis Figure 13(a)

``` r

fit_gauss$fit_plot
```

![Simulated log-price series with Gaussian noise over 2,000 time steps
in blue, closely tracked by the fitted LPPLS curve in red. Dashed lines
mark T2 = 1800 in green and the estimated critical time in red, near the
end of the series.](lpplsF_guide_files/figure-html/figfit-gauss-1.png)

``` r

knitr::kable(est_row(fit_gauss), digits = 5,
             caption = "M1 (Gaussian noise) estimate; tc_hat = 1996.")
```

|      tc |       m |   omega |       A |       B |      C1 |    C2 |       D |
|--------:|--------:|--------:|--------:|--------:|--------:|------:|--------:|
| 1996.01 | 0.60502 | 7.95536 | 7.97948 | -0.0239 | 0.00169 | 4e-04 | 1.04878 |

M1 (Gaussian noise) estimate; tc_hat = 1996. {.table}

#### M2 — AR(1)-GARCH(1,1) noise — thesis Figure 13(b)

``` r

fit_garch$fit_plot
```

![The simulated series with AR(1)-GARCH(1,1) noise in blue and its
fitted LPPLS curve in red, with dashed lines at T2 = 1800 in green and
at the estimated critical time in red, a little earlier than in the
Gaussian case.](lpplsF_guide_files/figure-html/figfit-garch-1.png)

``` r

knitr::kable(est_row(fit_garch), digits = 5,
             caption = "M2 (AR(1)-GARCH(1,1) noise) estimate; tc_hat = 1980.")
```

|      tc |       m |   omega |       A |        B |      C1 |      C2 |       D |
|--------:|--------:|--------:|--------:|---------:|--------:|--------:|--------:|
| 1979.96 | 0.56747 | 7.93666 | 8.04315 | -0.03272 | 0.00239 | 0.00011 | 0.97822 |

M2 (AR(1)-GARCH(1,1) noise) estimate; tc_hat = 1980. {.table}

### Reading the fit: `print()` diagnostics

A calibration can look perfectly good as a plot and still be
untrustworthy. The [`print()`](https://rdrr.io/r/base/print.html) method
for
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
reports the estimate together with the checks that decide whether it
means anything — this is the fastest way to find out whether a fit is
worth interpreting.

``` r

fit_gauss
```

    ## LPPLS calibration (mode = MPL)
    ##   Window:    n = 1800   fh = 300   hold_out = 200   num_searches = 1
    ##   Estimate:  tc = 1996.01   m = 0.6050   omega = 7.9554   D = 1.049
    ##              SSE = 4.24025
    ##   Bounds:    estimate interior;  fits at a bound: m 0/300, omega 0/300
    ##   Basins:    SSE 4.24025..4.88259;
    ##              no distinct inferior cluster (largest gap 0.02 of the spread)
    ##   Filters:   B out of range 0   D out of range 0
    ##   MPL:       defined at 300/300 tc   MPLE = 1996
    ##              LI(0.05) [1970, 2023]   LI(0.10) [1973, 2019]
    ##              LI(0.50) [1984, 2009]

Line by line:

- **Bounds.** \\(m, \omega)\\ are optimised inside a box. If the optimum
  sits *on* a wall, the SSE is still falling in that direction and the
  estimate is a constrained one — the regularity conditions behind the
  modified profile likelihood do not hold, and its intervals mean
  nothing. The check is reported twice: for the point estimate, and as a
  count across all the per-\\t_c\\ conditional fits. The second matters
  on its own, because an interior point estimate can still rest on a
  profile that is pinned nearly everywhere.
- **Basins.** Each candidate \\t_c\\ gets its own \\(m, \omega)\\
  search, and with few restarts some land in a worse local optimum than
  others. Distance from the best SSE will not detect this — the
  conditional SSE *should* rise away from the optimal \\t_c\\, so a
  healthy profile also has fits well above the minimum. What betrays a
  failed search is a **gap**: the sorted SSEs falling into a good
  cluster and a distinctly worse one. The reported number is the largest
  gap between consecutive sorted values as a share of the total spread,
  which is scale-free and near zero for a smooth profile.
- **Filters.** How many candidate fits were rejected for \\B \ge 0\\ or
  for damping, \\D\\, below its floor. These are the two constraints
  [`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
  reports internally, otherwise visible only as warnings while the fit
  runs.
- **MPL.** How much of the \\t_c\\ grid the modified profile likelihood
  is defined on. It is `NA` wherever \\\det(X'X - H) \le 0\\, which at a
  constrained optimum carries no guarantee of sign.

A `Diagnostics:` block of `!` lines appears only when something is
wrong, so a clean fit prints none — as here. [The S&P 500
example](#what-print-says-about-this-one) shows one that prints six.

``` r

fit_garch
```

    ## LPPLS calibration (mode = MPL)
    ##   Window:    n = 1800   fh = 300   hold_out = 200   num_searches = 1
    ##   Estimate:  tc = 1979.96   m = 0.5675   omega = 7.9367   D = 0.978
    ##              SSE = 4.4327
    ##   Bounds:    estimate interior;  fits at a bound: m 0/300, omega 0/300
    ##   Basins:    SSE 4.4327..5.62348;
    ##              no distinct inferior cluster (largest gap 0.01 of the spread)
    ##   Filters:   B out of range 0   D out of range 0
    ##   MPL:       defined at 300/300 tc   MPLE = 1980
    ##              LI(0.05) [1957, 2004]   LI(0.10) [1959, 2001]
    ##              LI(0.50) [1969, 1991]

### Per-\\t_c\\ parameter estimates — thesis Figure 14

Passing `pp = TRUE` to
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
attaches a parameter plot: each estimated parameter (and the objective)
as a function of the candidate \\t_c\\, with the objective minimum
(\\\hat t_c\\) marked. \\C_1\\ and \\C_2\\ trace near-sinusoidal curves
that distort far from the true \\t_c\\ — more so under GARCH noise. No
extra calibration runs; these reuse the fits above.

#### M1 — Gaussian noise — thesis Figure 14(a)

``` r

fit_gauss$param_plot
```

![Eight stacked panels, one each for A, B, C1, C2, D, m, omega and the
objective value, plotted against the candidate critical time from 1800
to 2100, with a red vertical line at the estimate. A, B, D, m and omega
change smoothly, C1 and C2 oscillate, and the objective is lowest at the
red line.](lpplsF_guide_files/figure-html/param-gauss-1.png)

#### M2 — AR(1)-GARCH(1,1) noise — thesis Figure 14(b)

``` r

fit_garch$param_plot
```

![The same eight panels for the AR(1)-GARCH(1,1) series: each estimate
against the candidate critical time, with a red vertical line at the
estimate. C1 and C2 again oscillate, more distorted than under Gaussian
noise.](lpplsF_guide_files/figure-html/param-garch-1.png)

## Reading the optimizer: BFGS traces

*Thesis Figure 15.* With `tp = TRUE`,
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
attaches a trace plot. It shows the SSE landscape over \\(m, \omega)\\
at the estimated \\t_c\\, overlaid with the BFGS path — red dot at the
start \\(m\_{\text{init}}, \omega\_{\text{init}}) = (0.2, 12)\\, green
diamond at the estimate.

\\B\\ and the damping \\D\\ are not search coordinates. At every point
of the landscape they follow from the least-squares solution for the
linear parameters, so the plot shows them as regions instead: shaded
where that solution would be out of range, black for \\B \ge 0\\ and
white for \\D\\ below its floor of 0.8, each edged with a dashed line.
The conditions are checked only where a search ends, so a path may cross
a shaded region on its way.

### \\(m, \omega)\\ topology — thesis Figure 15

The governing topology for calibration: a single interior minimum BFGS
finds easily. The GARCH minimum sits slightly further from the true
\\(0.6, 8)\\. In both, the white lobe at small \\m\\ is where the
damping would fall below 0.8, and the Gaussian path touches its edge on
the way to the estimate. \\B\\ is negative across the whole box, so
nothing is shaded black.

``` r

fit_gauss$trace_plot   # Gauss
```

![Heat map of the SSE over m from 0.1 to 0.9 and omega from 6 to 13,
with white contour lines around a single interior basin. An orange
optimiser path runs from the red start point at m = 0.2, omega = 12 to a
green diamond at the estimate near m = 0.6, omega = 8. A white shaded
lobe at small m, between omega 7.5 and 10, marks where the damping would
fall below 0.8; the path touches its
edge.](lpplsF_guide_files/figure-html/trace-mo-gauss-1.png)

``` r

fit_garch$trace_plot   # GARCH
```

![The same SSE heat map over m and omega for the AR(1)-GARCH(1,1)
series: a single basin, with the optimiser path from the red start point
to a green diamond at the estimate, slightly further from the true m =
0.6, omega = 8. A similar white lobe at small m marks where the damping
would fall below
0.8.](lpplsF_guide_files/figure-html/trace-mo-garch-1.png)

## Robustness to noise

*Thesis Figures 16, 17 and 18.* How much noise can a calibration take
before its estimates stop meaning anything? That is a property of the
model rather than of any particular series, so it is studied on a known
signal: take a deterministic LPPLS path, add noise at a range of levels,
recalibrate at each, and watch the estimates move.

[`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md)
runs that experiment. It takes the clean signal, a vector of noise
levels, and the arguments to pass to
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md),
and returns one row per level:

``` r

N      <- 2000
signal <- eval_lppls(seq_len(N), A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
                     tc = N + 1, m = 0.6, omega = 8)

sweep <- lppls_noise_sweep(
  signal,
  sd_seq   = seq(0.005, 0.100, by = 0.005),      # 20 noise levels
  fit_args = list(fh = 150, hold_out = 100,
                  lower = c(0.1, 6, -0.03, 0.8),
                  upper = c(0.9, 13, -1e-14, 1e6),
                  m_init = 0.2, o_init = 12, num_searches = 1,
                  mode = "MPL", mpl_cutoff = c(0.05, 0.1, 0.5)),
  seed = 3630
)
```

The `noise` argument decides what kind of perturbation is studied. It
defaults to `rnorm`, and takes any function of `n` returning
standardised innovations, so a heavy-tailed or serially dependent
alternative can be compared against the Gaussian case on the very same
signal:

``` r

## Student-t innovations, rescaled to unit variance so `sd_seq` still reads as
## a standard deviation
heavy <- function(n) rt(n, df = 4) / sqrt(4 / 2)
lppls_noise_sweep(signal, sd_seq = seq(0.005, 0.100, by = 0.005),
                  noise = heavy, fit_args = fit_args, seed = 3630)
```

The four sweeps below are of that shape: M1 and M2 under Gaussian noise,
and the same pair under AR(1)-GARCH(1,1) noise, where M2 is M1 with
\\C_1, C_2\\ negated (opposite phase). That is **80 calibrations**, so
the tables are shipped precomputed (`noise_*.rds`) and loaded here; see
[Reproducibility and the cache](#reproducibility-and-the-cache).

### Parameter estimates vs noise level — thesis Figure 16

[`noise_param_plot()`](https://mahovo.github.io/lpplsF/reference/noise_param_plot.md)
draws each estimate against the noise level, one panel per parameter.
\\\omega, C_1, C_2\\ (the log-periodic parameters) are the most
noise-sensitive.

#### Gaussian noise, M1

``` r

noise_param_plot(noise_m1_gauss)
```

![Six stacked panels showing the estimates of m, omega, A, B, C1 and C2
against the noise standard deviation, from 0.005 to 0.1, for M1 under
Gaussian noise. Each drifts steadily as the noise grows; C1 roughly
halves and C2 changes sign, while m barely
moves.](lpplsF_guide_files/figure-html/noise-sd-gauss-1.png)

#### AR(1)-GARCH(1,1) noise, M1

``` r

noise_param_plot(noise_m1_garch)
```

![The same six panels for M1 under AR(1)-GARCH(1,1) noise: each estimate
against the noise standard deviation, with the log-periodic parameters
omega, C1 and C2 again the most
sensitive.](lpplsF_guide_files/figure-html/noise-sd-garch-1.png)

### \\t_c\\ likelihood intervals vs noise level

[`noise_tc_plot()`](https://mahovo.github.io/lpplsF/reference/noise_tc_plot.md)
draws the \\t_c\\ estimate and the MPLE with the nested
likelihood-interval bands. Pass `true_tc` on a simulated signal and it
marks the true value, here \\t_c = 2001\\.

#### Gaussian noise: M1, then M2 at opposite phase — thesis Figure 17

``` r

noise_tc_plot(noise_m1_gauss, true_tc = 2001)
```

![The estimated critical time against the noise standard deviation for
M1 under Gaussian noise. The PLE and MPLE lines coincide and fall from
about 2000 to about 1978 as the noise grows, inside nested
likelihood-interval bands that widen with it. A dashed line marks the
true value, 2001, which stays inside the widest
band.](lpplsF_guide_files/figure-html/noise-li-m1-gauss-1.png)

![The same plot for M2, at opposite phase. Under Gaussian noise the PLE
and MPLE rise from about 2002 to about 2018 as the noise grows, inside
likelihood-interval bands that widen with it, while the true value,
2001, stays just below the lower edge of the innermost
band.](lpplsF_guide_files/figure-html/noise-li-m2-gauss-1.png)

#### AR(1)-GARCH(1,1) noise: M1, then M2 at opposite phase — thesis Figure 18

![The estimated critical time and its likelihood-interval bands against
the noise standard deviation for M1 under AR(1)-GARCH(1,1) noise, with a
dashed line at the true value,
2001.](lpplsF_guide_files/figure-html/noise-li-m1-garch-1.png)

![The estimated critical time and its likelihood-interval bands against
the noise standard deviation for M2 under AR(1)-GARCH(1,1) noise, with a
dashed line at the true value,
2001.](lpplsF_guide_files/figure-html/noise-li-m2-garch-1.png)

## Rolling windows and Lagrange regularization

*Thesis Figures 21 to 23.* Lagrange regularization (Demos & Sornette,
2019) tries to find the onset of a bubble by rolling the calibration
start \\T_1\\ forward (fixed end \\T_2\\) and looking for the window
where the normalized residual \\\chi^2\_{np}(T_1)\\ stops falling. The
package splits this into a reusable primitive and cheap post-processing:

- [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md)
  runs one calibration per \\T_1\\ (the expensive part);
- [`lppls_lagrange()`](https://mahovo.github.io/lpplsF/reference/lppls_lagrange.md),
  [`rolling_tc_plot()`](https://mahovo.github.io/lpplsF/reference/rolling_tc_plot.md),
  [`rolling_param_plot()`](https://mahovo.github.io/lpplsF/reference/rolling_param_plot.md),
  [`chi_sq_np_plot()`](https://mahovo.github.io/lpplsF/reference/chi_sq_np_plot.md),
  [`chi_sq_lambda_plot()`](https://mahovo.github.io/lpplsF/reference/chi_sq_lambda_plot.md)
  are all post-processing on it.

In a pure LPPL process every window is already LPPL, so \\\chi^2\_{np}\\
is flat and regularization has nothing to correct. So the example is a
**dual-regime** series consisting of a flat “ranging” segment for \\t \<
600\\ that switches to an LPPL process at \\t = 600\\
(`tc = 2001, m = 0.6, omega = 8, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001`;
Gaussian noise \\\sigma = 0.05\\, seed 3630). The rolling approaches the
onset from below: as \\T_1\\ climbs toward 600 each window sheds more of
the pre-LPPL data and the fit tightens, so \\\chi^2\_{np}\\ decreases;
once \\T_1\\ passes 600 the window is pure LPPL again and the trend
flattens. The rolling (\\T_1 = 1 \to 1001\\ step 10, \\T_2 = 1900\\; 101
windows) is shipped precomputed as `lagrange_rolling_regime_gauss.rds`.

``` r

## The dual-regime series: flat at the LPPL value at t = 600 for t < 600, then
## the LPPL process itself, with N(0, 0.05^2) noise throughout. In this example
## the two segments draw their noise from separate seeds: `seed` for the flat
## prefix, `seed + 1` for the bubble. So both set.seed() calls are needed to
## reproduce it. (This is not really important, it just happens to be the way
## the figures for the thesis were generated.)
N_rg <- 2000
build_rg <- function(seed, onset = 600) {
  bub <- eval_lppls(onset:N_rg, 8, -0.025, 0.0015, 0.001, N_rg + 1, 0.6, 8)
  set.seed(seed + 1)
  c(rep(bub[1], onset - 1), bub) +
    c({ set.seed(seed); rnorm(onset - 1) * 0.05 }, rnorm(N_rg - onset + 1) * 0.05)
}

## Regenerated on both paths: it is instant, and the residual guard below needs
## the series even when the calibration comes from the cache.
lp_rg <- build_rg(3630)
```

``` r

## 101 windows, each a 150-candidate tc scan (MPL: the likelihood
## intervals rolling_tc_plot() needs).
roll_rg <- lppls_rolling(
  lp_rg,
  t1_seq = seq(1, 1001, by = 10),
  hold_out = 100,
  fh = 150,
  fit_args = list(
    lower = c(0.1, 6, -0.03, 0.8),
    upper = c(0.9, 13, -1e-14, 1e6),
    m_init = 0.2,
    o_init = 12,
    num_searches = 1,
    mode = "MPL",
    mpl_cutoff = c(0.05, 0.1, 0.5),
    factr = 1e7)
  )
```

``` r

lag_rg <- lppls_lagrange(roll_rg, p = 0)
lag_rg
```

    ## LPPLS Lagrange regularization
    ##   Windows (T1):  101  (from 1 to 1001)
    ##   T2 (cal. end): 1900
    ##   p (denominator correction): 0
    ##   Trend: chi^2_np = 0.0054087 + (-3.9517e-06) * T1
    ##   lambda* = -beta1 = 3.9517e-06
    ##   Optimal T1* = 501

Optimal \\T_1^\* =\\ **501** (true onset 600); \\\lambda =\\ 3.95^{-6}
\\\> 0\\ — a genuine bias trend (a pure LPPL series would give \\\lambda
\approx 0\\).

[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md)
keeps the fit from every window. To examine the one that Lagrange
regularization selected, locate it in the sequence and pull out the
quantities that describe it.

``` r

k_rg  <- which(roll_rg$t1_seq == lag_rg$t1_star)  # position of the optimal window
T1_rg <- roll_rg$t1_seq[k_rg]                     # its calibration start, T1*
nc_rg <- roll_rg$table$n_cal[k_rg]                # how many points it calibrated on
T2_rg <- T1_rg + nc_rg - 1                        # its calibration end, T2
f_rg  <- roll_rg$fits[[k_rg]]                     # the fit itself
```

### The generated series and its fit — thesis Figure 21

Blue the generated series, red the fitted LPPLS process up to \\\hat
t_c\\, dashed markers at the calibration bounds and \\\hat t_c\\. The
series spans two regimes and the fit is calibrated on only the later
part, so the red curve starts at \\T_1^\*\\ — the flat regime is shown
but never fitted.
[`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md)
stores the fit in window-local time (\\t = 1\\ at \\T_1\\), so \\\hat
t_c\\ is shifted by the window offset to draw it on the global clock.

``` r

f_glob    <- f_rg
f_glob$tc <- f_rg$tc + (T1_rg - 1)   # window-local -> global clock

fit_plot(f_glob, seq_along(lp_rg), lp_rg, T1 = T1_rg, T2 = T2_rg)
```

![The dual-regime series in blue, flat with noise up to t = 600 and then
rising as an LPPLS bubble to t = 2000. The fitted LPPLS curve in red
starts at the selected T1, marked by an orange dashed line, and ends at
the estimated critical time, marked in red; a green dashed line marks T2
= 1900.](lpplsF_guide_files/figure-html/regime-gauss-fit-1.png)

Give
[`fit_plot()`](https://mahovo.github.io/lpplsF/reference/fit_plot.md)
both window bounds and it draws the curve only where the fit applies.

### Rolling-\\T_1\\ estimates — thesis Figure 22

How the \\t_c\\ estimate and the six parameters move as \\T_1\\ rolls
across the onset. For \\T_1\\ well below 600 the window still contains
ranging data and the estimates wander; as \\T_1\\ approaches and passes
600 the window becomes predominantly, then purely, LPPL and \\\hat t_c\\
settles near the true 2001.

``` r

rolling_tc_plot(roll_rg)
```

![The estimated critical time against the window start T1, from 1 to
1001, with nested likelihood-interval bands. It jumps between about 1930
and 2050 while the windows still hold much of the flat regime, then
climbs and settles near 1990 once T1 passes about 600, where the widest
bands contain the true value,
2001.](lpplsF_guide_files/figure-html/regime-gauss-tc-1.png)

``` r

rolling_param_plot(roll_rg)
```

![Six stacked panels showing A, B, C1, C2, m and omega against the
window start T1. All six jump erratically for T1 below about 270, then
settle: A near 8, B just below 0, C1 and C2 close to 0, m rising to
about 0.6, and omega rising to about 8 by T1 =
600.](lpplsF_guide_files/figure-html/regime-gauss-params-1.png)

### Lagrange regularization — thesis Figure 23

\\\chi^2\_{np}(T_1)\\ with the OLS bias trend (red), then the
regularized objective \\\chi^2\_{\lambda^\*}(T_1)\\ with the optimal
\\T_1^\*\\ (green):

``` r

chi_sq_np_plot(lag_rg)
```

![Scatter of the normalised residual chi-squared against T1 with a red
fitted linear trend sloping downward. The points fall from about 0.006
at T1 = 0 to about 0.0024 by T1 = 570, scattered in two tiers below T1 =
270, and stay flat after
that.](lpplsF_guide_files/figure-html/regime-gauss-chi-np-1.png)

``` r

chi_sq_lambda_plot(lag_rg)
```

![Scatter of the regularised objective against T1, V-shaped with its
minimum near T1 = 500, where a green vertical line marks the optimal
T1\*. Below T1 = 270 the points scatter in two
tiers.](lpplsF_guide_files/figure-html/regime-gauss-chi-lambda-1.png)

> **Is that decrease real, or a noise artifact?** That question, and the
> autocorrelation-robust
> [`lr_trend_test()`](https://mahovo.github.io/lpplsF/reference/lr_trend_test.md)
> that answers it, are taken up in a separate investigation,
> [`lagrange_inference_report.Rmd`](https://github.com/mahovo/lpplsF/blob/main/misc/lagrange_inference_report.Rmd),
> in the package’s `misc/` directory on GitHub. It is not built with the
> package.

### Residuals at \\T_1^\*\\

The residual QQ plot for the window at \\T_1^\*\\. The residuals are
near-normal.

``` r

residual_qq_plot(res_rg)
```

![Normal QQ plot of the residuals of the window at T1\*, with a red
reference line and green lines at the 10% and 90% quantiles. The points
lie on the line except for slight departures in the extreme
tails.](lpplsF_guide_files/figure-html/regime-gauss-res-1.png)

## Empirical example: the S&P 500

*Thesis Figures 25 to 30.* The package ships the S&P 500 series as the
dataset `spx` (see
[`?spx`](https://mahovo.github.io/lpplsF/reference/spx.md)), covering
exactly the calibrated window, 2009-05-18 to 2020-12-31. We fit on
\\\[1, 2428\]\\ (to 2019-01-08) with `hold_out = 500`, `fh = 500`.

``` r

## `spx` is lazy-loaded with the package: 2,928 daily closes, 2009-05-18..2020-12-31.
df_spx <- data.frame(Date = spx$date, logpris = log(spx$close))
df_spx$ID <- seq_len(nrow(df_spx))
n_cal_spx <- nrow(df_spx) - 500    # 2428
```

> **A note on the series.** The thesis calibrated the CRSP S&P 500
> Composite series (`spindx`, obtained through Wharton Research Data
> Services). That data is licensed and cannot be redistributed, so the
> package ships the corresponding Yahoo Finance series (`^GSPC`) for the
> same trading days. The two agree to 1.6e-07 on average, and
> recalibrating on the shipped series leaves every quantity reported
> below unchanged: \\\hat t_c = 2717\\ (2020-03-03), \\m = 0.9\\,
> \\\omega = 6.1596\\, and all three likelihood intervals identical.
> Only the Brent-refined \\t_c\\ at the final step of the MPL
> computation moves, in its seventh decimal. That step refines thesis
> equation 38 at the \\(\hat m, \hat\omega)\\ which minimise equation 39
> on the integer grid.

``` r

## fh = 500 profile-likelihood optimizations
spx_fit <- fit_lppls(
  log_price = df_spx$logpris, fh = 500, hold_out = 500,
  lower = c(0.1, 6, -0.5, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
  tc_init = n_cal_spx + 250, m_init = 0.2, o_init = 12, num_searches = 1,
  mode = "MPL", mpl_cutoff = c(0.05, 0.1, 0.5),
  dates = df_spx$Date,          # calendar labels for the 1..n index
  style = "empirical",          # grey data; the forecast past T2 in red
  fp = TRUE, mpl_plot = TRUE, tp = TRUE, pp = TRUE)
```

``` r

e <- spx_fit$fit[[1]]
```

The point estimate: PLE \\t_c =\\ **2020-03-03**, with \\m =\\ 0.9
(pinned at the upper boundary) and \\\omega =\\ 6.16.

#### What `print()` says about this one

[Reading the fit](#reading-the-fit-print-diagnostics) introduced the
diagnostics on a healthy calibration, which printed no warnings. This is
the contrasting case, and it is worth reading before any of the figures
below:

``` r

spx_fit
```

    ## LPPLS calibration (mode = MPL)
    ##   Window:    n = 2428   fh = 500   hold_out = 500   num_searches = 1
    ##   Estimate:  tc = 2717.01 (2020-03-03)   m = 0.9000   omega = 6.1596   D = 1.545
    ##              SSE = 5.28678
    ##   Bounds:    estimate m at a bound;  fits at a bound: m 500/500, omega 169/500
    ##   Basins:    SSE 5.28678..8.88078;  largest gap 0.56 of the spread splits off
    ##              394/500 conditional fits (79%)
    ##   Filters:   B out of range 0   D out of range 0
    ##   MPL:       defined at 455/500 tc   MPLE = 2717 (2020-03-03)
    ##              LI(0.05) [2696, 2757]   (2020-01-31 .. 2020-04-29)
    ##              LI(0.10) [2696, 2757]   (2020-01-31 .. 2020-04-29)
    ##              LI(0.50) [2696, 2744]   (2020-01-31 .. 2020-04-09)
    ## 
    ##   Diagnostics:
    ##     ! m pinned at an optimisation bound - a boundary solution, so the modified
    ##       profile likelihood and its intervals are unreliable
    ##     ! 394 of 500 conditional fits (79%) sit in a distinct inferior basin (SSE >=
    ##       7.99365 against a best of 5.28678) - raise num_searches
    ##     ! MPL undefined at 45 of 500 tc (det(X'X - H) <= 0 there)
    ##     ! LI(0.05) spans 62 grid points but only 58 are above its cutoff - range()
    ##       bridges 4 that are not in the interval
    ##     ! LI(0.10) spans 62 grid points but only 58 are above its cutoff - range()
    ##       bridges 4 that are not in the interval
    ##     ! every tc inside an interval is a good-basin fit - the interval traces
    ##       where the search succeeded, not the curvature of the likelihood

Six flags, and every caveat in the rest of this section follows from
them. The estimate is a **boundary solution**, and not marginally —
\\m\\ is at its upper bound in the point estimate *and* in every
conditional fit, with \\\partial\\\mathrm{SSE}/\partial m \< 0\\ there,
so the constraint is genuinely active and the regularity conditions
behind the modified profile likelihood do not hold. Most of the
conditional fits sit in a **distinctly inferior basin**, which is why
the MPL curve below is two-tiered rather than smooth, and why
`num_searches = 1` is the wrong setting for this series. And the MPL is
**undefined on part of the grid**, so the curve has gaps the likelihood
intervals cannot see across.

The last three are about the likelihood intervals themselves. LI(0.05)
and LI(0.10) each **bridge grid points** that are not above their
cutoff, and every \\t_c\\ inside an interval is a **good-basin fit**, so
the intervals trace where the search succeeded rather than the curvature
of the likelihood. [The modified profile
likelihood](#modified-profile-likelihood-thesis-figure-26) below takes
this up.

None of this is visible in the point estimate, which looks entirely
ordinary.

It would be easy to read the warnings above as a defect in the data, in
the settings, or in the package, and to go looking for a configuration
that silences them. That would be the wrong lesson. Boundary solutions,
split basins and an undefined modified profile likelihood are
*properties of the LPPLS model itself*: the objective is non-convex, the
log-periodic parameters are only weakly identified, and the regularity
conditions the modified profile likelihood assumes are not guaranteed to
hold at a constrained optimum.

The master’s thesis for which this package was written takes up exactly
this question: When an LPPLS calibration is worth believing, and what to
do when it is not. The S&P 500 case above is its Chapter 9.

### Fit plot — thesis Figure 25

Because the calibration was given `dates`, the attached `$fit_plot` is
already drawn on a calendar axis. Pass `dates` to
\[[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)\]
and every plot and printout it produces speaks in dates; omit it and
everything falls back to the 1..n time index.

It was also given `style = "empirical"`. An observed price series is
given and cannot be otherwise, so it is drawn in grey, and the fitted
curve splits at \\T_2\\: blue where data supports it, red beyond, where
the curve is a forecast. The default `style = "simulated"` draws a
generated series in colour with a single fitted line, which is what the
earlier sections use.

``` r

plot(spx_fit, "fit")
```

![The S&P 500 log price from 2009 to 2020 in grey, with the fitted LPPLS
curve in blue up to T2 in January 2019 and its forecast continuation in
red up to the estimated critical time in March 2020. Dashed vertical
lines mark T2 in blue and the critical time in red; the index crashes
around the critical time and then
recovers.](lpplsF_guide_files/figure-html/spx-fitplot-1.png)

The same plot can be drawn from a bare estimate with
[`fit_plot()`](https://mahovo.github.io/lpplsF/reference/fit_plot.md),
which takes `dates` directly:

``` r

fit_plot(spx_fit$fit[[1]], seq_len(nrow(df_spx)), df_spx$logpris,
         T2 = n_cal_spx, dates = df_spx$Date, style = "empirical")
```

### Modified profile likelihood — thesis Figure 26

The MPLE coincides with the PLE, at 2020-03-03. The `MPL:` line of the
printout above gives that date next to the index, and each likelihood
interval in both forms. The modified profile likelihood is `NA` wherever
\\\det(X'X - H) \le 0\\. At a constrained optimum that determinant
carries no guarantee of sign, which is why the same line counted gaps in
the grid.

> **Diagnostic caveat.** The point estimate agrees with the PLE, but the
> bands are not inference — the `Bounds:` and `Basins:` lines above
> already say why. Two consequences the printed summary does not spell
> out. First, the basin split is what makes the curve **two-tiered**
> rather than smooth: the two tiers are \\\approx 500\\ log-likelihood
> units apart, so the plot below is not one curve but two.
>
> Second, and the reason to distrust the bands specifically: the
> likelihood intervals here trace **where the optimiser succeeded**, not
> the curvature of the likelihood. Every grid point inside LI_05 is a
> good-basin fit and every point outside it is not; its lower edge is
> exactly where the good-basin run containing the PLE begins, and its
> upper edge is the end of the *next* such run, which
> [`range()`](https://rdrr.io/r/base/range.html) bridges across the gap
> between them. LI_05 and LI_10 coincide for the same reason. The
> relative-scale plot below makes this visible: only LI_50 has a
> crossing. **Read the MPLE, not the bands.**

``` r

spx_fit$mpl_plot
```

![The log modified profile likelihood of the S&P 500 fit against
candidate critical times from about 2460 to 2930, in two tiers roughly
500 units apart, with gaps: an upper tier peaking in a short plateau
near 7383 around the estimate, and a long lower tier between about 6770
and 6890. Vertical lines mark the PLE, the MPLE and the bounds of the
5%, 10% and 50% likelihood
intervals.](lpplsF_guide_files/figure-html/spx-mplplot-1.png)

#### The same curve on a relative scale

The absolute level of the log MPL carries no information. Its magnitude
comes almost entirely from the variance term of the modified profile
likelihood (thesis eq. 53), \\-\tfrac{n - p - 2}{2}\log
\hat{s}\_{t_c}\\. A well-fitting model has a small residual variance,
here \\\hat{s} = \sigma^2 \approx 0.0022\\, whose logarithm is large and
negative. Multiplying by \\-(n - p - 2)/2 \approx -1210\\ flips the sign
and scales it by the sample size, which gives \\\approx 7417\\. The
determinant ratio contributes only \\\approx -34\\, so the curve peaks
near \\7383\\. The formula also drops additive constants that do not
depend on \\t_c\\. What the likelihood intervals actually test is the
*difference* \\R = \log\mathrm{MPL} - \max\log\mathrm{MPL}\\ against
\\\log(\text{cutoff})\\, and the fit already stores it as
`mpl_output$R`.

Plotting \\R\\ and restricting the view to the band the cutoffs live in
provides a better view of the optimum. `plot(fit, "mpl_relative")` draws
that view from the stored `mpl_output`, so it needs nothing requested at
fit time. The horizontal lines are the cutoffs themselves, so each
interval can be read off where the curve crosses one.

``` r

plot(spx_fit, "mpl_relative")
```

![The same curve on a relative scale, log MPL minus its maximum, with
horizontal lines at the cutoffs log 0.5, log 0.1 and log 0.05. The curve
exists only between about 2696 and 2757: it peaks at 0 at the estimate
and falls to about -1.15, so it crosses only the log 0.5
line.](lpplsF_guide_files/figure-html/spx-mplplot-rel-1.png)

Only `LI_50` has a crossing: the curve descends from the peak to about
\\-1.15\\ and then stops, so it never reaches \\\log 0.10\\ or \\\log
0.05\\. Those two intervals end where the good basin ends.

### The rest of the diagnostics

[`print()`](https://rdrr.io/r/base/print.html) is the quick triage.
Three further tools answer the questions it raises: *where* the basins
are, whether the estimate sits at a genuine minimum, and whether the
intervals are closed by the likelihood or by the grid.

#### The basin split, drawn

`plot(fit, "basin")` shows every conditional fit’s SSE against its
\\t_c\\, split at the gap [`print()`](https://rdrr.io/r/base/print.html)
reported. It needs nothing that was not already computed, so it works
even on a fit that requested no plots.

``` r

plot(spx_fit, "basin")
```

![SSE of each conditional fit against its candidate critical time,
coloured by basin: short runs of good-basin fits in blue with SSE
between about 5.3 and 6, and inferior-basin fits in red with SSE between
about 8 and 9 covering most of the
grid.](lpplsF_guide_files/figure-html/spx-basin-1.png)

The same view of the well-behaved calibration from [Fitting a single
model](#fitting-a-single-model), for contrast — one cloud, rising
smoothly away from the optimal \\t_c\\, which is what a profile should
look like:

``` r

plot(fit_gauss, "basin")
```

![The same view for the Gaussian-noise calibration: a single smooth
curve of conditional fits with no distinct inferior cluster, falling
from about 4.9 at tc = 1800 to a minimum near 2000 and rising again
toward 2100.](lpplsF_guide_files/figure-html/gauss-basin-1.png)

#### What the basins are made of

The view above shows *that* the fits split, but not what they split on.
`"param_basin"` applies the same division to every estimated parameter
rather than only the objective, which is where the mechanism becomes
visible.

``` r

plot(spx_fit, "param_basin")
```

![Eight stacked panels showing A, B, C1, C2, D, m, omega and the
objective value of every conditional fit against its candidate critical
time, coloured by basin, with a red vertical line at the estimate. The
basins separate most clearly in omega, near 6 for the good fits and
between 10 and 13 for the inferior ones, while m is 0.9 in both; the
inferior fits' C1 and C2 oscillate rapidly and the good ones trace
smooth arcs.](lpplsF_guide_files/figure-html/spx-param-basin-1.png)

``` r

## The split is the one summary() reports; here it is per parameter.
cut_at <- summary(spx_fit)$split_at[1]
pt     <- spx_fit$fit[[2]]
pt$basin <- ifelse(pt$value <= cut_at, "good", "inferior")
rng <- function(z) sprintf("%.3f - %.3f", min(z), max(z))
knitr::kable(
  aggregate(cbind(m, omega, D) ~ basin, data = pt, FUN = rng),
  caption = "Range of each nonlinear estimate within each basin.")
```

| basin    | m             | omega           | D             |
|:---------|:--------------|:----------------|:--------------|
| good     | 0.900 - 0.900 | 6.000 - 6.392   | 1.412 - 1.548 |
| inferior | 0.900 - 0.900 | 10.251 - 13.000 | 2.363 - 3.365 |

Range of each nonlinear estimate within each basin. {.table}

The two basins are separated almost entirely by \\\omega\\: The good one
sits near the bottom of the \\\[6, 13\]\\ search box, the inferior one
near the top. \\m\\ is at \\0.9\\ in **both**, so it plays no part in
the split at all.

That separates two failures which the SSE view alone runs together. The
pinned \\m\\ is a **boundary** problem, and it affects every fit on the
grid. The two tiers are an **\\\omega\\ multimodality** problem, and it
affects only which optimum a single L-BFGS-B start happens to reach — a
`num_searches` setting, not a property of the data. The `C1`/`C2` facets
tell the same story: The inferior fits oscillate through a
high-frequency regime while the good ones trace a smooth low-frequency
arc.

#### Curvature at the estimate

The modified profile likelihood assumes the estimate is an interior
minimum: Its correction uses \\\det(X'X - H)\\, the observed
information, which is only guaranteed positive at an interior stationary
point. \\X\\ holds the model’s first derivatives with respect to the six
parameters, and \\H\\ its second derivatives, each weighted by the
residual at that observation and summed (thesis eq. 55). \\H\\ is thus a
residual-weighted curvature term rather than a Hessian in its own right.
[`lppls_curvature()`](https://mahovo.github.io/lpplsF/reference/lppls_curvature.md)
tests the determinant directly, and factorises it to show where the sign
comes from.

By the Schur complement, \\\det(X'X - H)\\ factorises into the
determinant of the linear-parameter block times the determinant of the
\\2 \times 2\\ Hessian of the \\(m, \omega)\\ profile. The model is
linear in \\A, B, C_1, C_2\\, so their second derivatives vanish, \\H\\
contributes nothing to the linear block, and what remains there is a
crossproduct with a positive determinant. The sign is therefore decided
entirely by curvature in the two nonlinear parameters: the \\2 \times
2\\ block the function reports.

``` r

lppls_curvature(spx_fit)
```

    ## LPPLS curvature at the estimate
    ##   tc = 2717.01   n = 2428
    ##   det(X'X - H)          = +6.5784e+29   (MPL is defined only where this is > 0)
    ##   det(linear block)     = +1.5668e+29   (positive by construction)
    ##   det(profile Hessian)  = +4.1986e+00   (sets the sign above)
    ##   eigenvalues (m, omega): +7.5551e+00   +5.5573e-01
    ##   => profile Hessian is positive definite: the curvature condition behind the
    ##      MPL holds
    ##   note: m sits at an optimisation bound, so the estimate is a constrained one -
    ##         the gradient need not vanish there, whatever the curvature says

Read the classification and the bound note as two separate facts. The
curvature *is* positive definite here, so the \\t_c\\ the MPL reports is
not an artefact of a degenerate Hessian. But \\m\\ rests on its bound,
so the gradient need not vanish: the estimate is a constrained one
whatever the curvature says. It is the combination of sound curvature
and an active constraint that makes the point estimate usable and the
intervals not.

#### Everything at once: `summary()`

[`summary()`](https://rdrr.io/r/base/summary.html) returns the
diagnostics as values rather than text, so they can be tested against
instead of parsed. Printing it gives the
[`print()`](https://rdrr.io/r/base/print.html) report plus the curvature
line and an interval breakdown.

``` r

s <- summary(spx_fit)
s$basin_gap
```

    ## [1] 0.5611065

``` r

s$mpl$li_check[[1]]
```

    ## $cutoff
    ## [1] 0.05
    ## 
    ## $n_in
    ## [1] 58
    ## 
    ## $span
    ## [1] 62
    ## 
    ## $contiguous
    ## [1] FALSE
    ## 
    ## $bridged
    ## [1] 4
    ## 
    ## $all_good_basin
    ## [1] TRUE

That last list is the interval integrity check.
[`fit_lppls()`](https://mahovo.github.io/lpplsF/reference/fit_lppls.md)
forms each interval as [`range()`](https://rdrr.io/r/base/range.html) of
the \\t_c\\ whose \\R\\ exceeds \\\log(\text{cutoff})\\ — and
[`range()`](https://rdrr.io/r/base/range.html) does not care whether
that set has holes. When it does, the reported interval contains grid
points that are not themselves in it:

``` r

s
```

    ## LPPLS calibration (mode = MPL)
    ##   Window:    n = 2428   fh = 500   hold_out = 500   num_searches = 1
    ##   Estimate:  tc = 2717.01 (2020-03-03)   m = 0.9000   omega = 6.1596   D = 1.545
    ##              SSE = 5.28678
    ##   Bounds:    estimate m at a bound;  fits at a bound: m 500/500, omega 169/500
    ##   Basins:    SSE 5.28678..8.88078;  largest gap 0.56 of the spread splits off
    ##              394/500 conditional fits (79%)
    ##   Filters:   B out of range 0   D out of range 0
    ##   MPL:       defined at 455/500 tc   MPLE = 2717 (2020-03-03)
    ##              LI(0.05) [2696, 2757]   (2020-01-31 .. 2020-04-29)
    ##              LI(0.10) [2696, 2757]   (2020-01-31 .. 2020-04-29)
    ##              LI(0.50) [2696, 2744]   (2020-01-31 .. 2020-04-09)
    ##   Curvature: positive definite at the estimate  (det(X'X - H) = +6.578e+29;
    ##              profile Hessian eigenvalues +7.555e+00, +5.557e-01)
    ## 
    ##   Diagnostics:
    ##     ! m pinned at an optimisation bound - a boundary solution, so the modified
    ##       profile likelihood and its intervals are unreliable
    ##     ! 394 of 500 conditional fits (79%) sit in a distinct inferior basin (SSE >=
    ##       7.99365 against a best of 5.28678) - raise num_searches
    ##     ! MPL undefined at 45 of 500 tc (det(X'X - H) <= 0 there)
    ##     ! LI(0.05) spans 62 grid points but only 58 are above its cutoff - range()
    ##       bridges 4 that are not in the interval
    ##     ! LI(0.10) spans 62 grid points but only 58 are above its cutoff - range()
    ##       bridges 4 that are not in the interval
    ##     ! every tc inside an interval is a good-basin fit - the interval traces
    ##       where the search succeeded, not the curvature of the likelihood
    ## 
    ##   Likelihood intervals:
    ##     LI(0.05): 58 grid points above the cutoff, interval spans 62
    ##               (bridges 4 points below the cutoff)
    ##     LI(0.10): 58 grid points above the cutoff, interval spans 62
    ##               (bridges 4 points below the cutoff)
    ##     LI(0.50): 49 grid points above the cutoff, interval spans 49  (contiguous)

`LI(0.05)` and `LI(0.10)` each span more grid points than are above
their cutoff: the interval bridges the gap between two runs of
good-basin fits. `LI(0.50)` is contiguous. This is the same conclusion
the relative-scale plot above shows visually, arrived at from the
numbers.

### Fit diagnostics — thesis Figure 30

The BFGS path slides to the \\(m, \omega)\\ boundary corner; the
per-\\t_c\\ parameter panel shows \\m\\ pinned at 0.9 with jumps where
it hits the boundary.

``` r

spx_fit$trace_plot
```

![Heat map of the SSE over m and omega for the S&P 500 fit, falling
steadily toward large m and small omega. The orange optimiser path runs
straight from the red start point at m = 0.2, omega = 12 to a green
diamond in the corner at m = 0.9, omega near 6, on the edge of the
search box. A small white lobe in the bottom-left corner, at small m and
omega, marks where the damping would fall below
0.8.](lpplsF_guide_files/figure-html/spx-trace-1.png)

``` r

spx_fit$param_plot
```

![Eight stacked panels showing each estimate of the S&P 500 calibration
against the candidate critical time, with a red vertical line at the
estimate. m sits at its upper bound of 0.9 for every candidate, and the
other parameters jump between two branches, one for each optimisation
basin.](lpplsF_guide_files/figure-html/spx-param-1.png)

### Residual QQ diagnostics — thesis Figure 27

Residuals of the SPX calibration (\\t \< \hat t_c\\) against a normal
and a fitted Student-t reference, both via
[`residual_qq_plot()`](https://mahovo.github.io/lpplsF/reference/residual_qq_plot.md).

``` r

residual_qq_plot(res_spx)
```

![Normal QQ plot of the S&P 500 calibration residuals with a red
reference line and green lines at the 10% and 90% quantiles. The points
follow the line in the middle but bend away from it at both ends, below
it in the lower tail and above it in the upper: heavier tails than the
normal.](lpplsF_guide_files/figure-html/spx-qq-norm-1.png)

``` r

residual_qq_plot(res_spx, "t")
```

![The same residuals against a fitted Student-t reference. The points
follow the line through the middle but flatten at both ends, above it in
the lower tail and below it in the upper: lighter tails than the fitted
t.](lpplsF_guide_files/figure-html/spx-qq-t-1.png)

### Lagrange regularization on the S&P 500 — thesis Figure 29

Rolling \\T_1 = 1 \to 501\\ step 25 (21 windows), \\T_2 = 2428\\ fixed;
shipped precomputed as `lagrange_rolling_spx.rds` (21-window MPL fit).
It is the `spx_fit` call above rolled over \\T_1\\, so window 1 is
`spx_fit` itself:

``` r

roll_spx <- lppls_rolling(
  df_spx$logpris, t1_seq = seq(1, 501, 25), hold_out = 500, fh = 500,
  fit_args = list(lower = c(0.1, 6, -0.5, 0.8), upper = c(0.9, 13, -1e-14, 1e6),
                  m_init = 0.2, o_init = 12, num_searches = 1, mode = "MPL",
                  mpl_cutoff = c(0.05, 0.1, 0.5)))
```

``` r

roll_spx <- readRDS(dpath("lagrange_rolling_spx.rds"))
lag_spx  <- lppls_lagrange(roll_spx, p = 0)
lag_spx
```

    ## LPPLS Lagrange regularization
    ##   Windows (T1):  21  (from 1 to 501)
    ##   T2 (cal. end): 2428
    ##   p (denominator correction): 0
    ##   Trend: chi^2_np = 0.00191588 + (-1.77776e-06) * T1
    ##   lambda* = -beta1 = 1.77776e-06
    ##   Optimal T1* = 76

Optimal \\T_1^\* =\\ **76** (= 2009-09-02). Here \\\beta_1 \< 0\\,
\\\lambda \> 0\\, so the Demos bias trend is genuinely present.

``` r

chi_sq_np_plot(lag_spx)
```

![Scatter of the normalised residual chi-squared against T1 for the 21
S&P 500 rolling windows, T1 from 1 to 501. It falls from about 0.0022 to
about 0.0011, with a plateau between T1 = 100 and 275, and a red fitted
linear trend slopes downward through
it.](lpplsF_guide_files/figure-html/spxroll-chi-np-1.png)

``` r

chi_sq_lambda_plot(lag_spx)
```

![Scatter of the regularised objective against T1 for the S&P 500
rolling windows: irregular, lowest at T1 = 76, where a green vertical
line marks the optimal T1\*, with a second, slightly higher dip near T1
= 350.](lpplsF_guide_files/figure-html/spxroll-chi-lambda-1.png)

``` r

rolling_tc_plot(roll_spx)   # PLE saturates at the tc bound
```

![The estimated critical time against T1 for the S&P 500 rolling
windows, with nested likelihood-interval bands. The PLE and MPLE climb
from about 2717 at T1 = 1 to the top of the candidate grid, about 2928,
by T1 = 125 and stay there, while the bands close in on the same
limit.](lpplsF_guide_files/figure-html/spxroll-output-T1-1.png)

## Reproducibility and the cache

This document separates three tiers of computation:

1.  **Live (always run).** The [SSE
    surfaces](#the-model-and-its-calibration-landscape), all data
    simulation and the three single-model calibrations `fit_gauss`,
    `fit_garch` and `spx_fit` are recomputed on every knit; the
    calibration code is shown in full above.
2.  **Cached, re-runnable in place.** The dual-regime rolling
    calibration is cached as `lagrange_rolling_regime_gauss.rds`. Knit
    with `params = list(recompute = TRUE)` to re-run it and refresh the
    cache.
3.  **Precomputed multi-fit tables (reload-only).** The 20-level noise
    experiments (`noise_*.rds`) and the S&P 500 rolling calibration
    (`lagrange_rolling_spx.rds`) are shipped precomputed. The
    specification for each is given in its own section, so they can be
    regenerated with
    [`lppls_noise_sweep()`](https://mahovo.github.io/lpplsF/reference/lppls_noise_sweep.md)
    and
    [`lppls_rolling()`](https://mahovo.github.io/lpplsF/reference/lppls_rolling.md).

The only input that cannot be regenerated deterministically is the
AR(1)-GARCH(1,1) series. `rugarch`’s simulated paths depend on the
package version, so it is always loaded from `garch_series.rds`.

    ## R version 4.6.2 (2026-10-29)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_4.0.3 lpplsF_0.1.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] plotly_4.12.1       sass_0.4.10         generics_0.1.4     
    ##  [4] tidyr_1.3.2         digest_0.6.39       magrittr_2.0.5     
    ##  [7] evaluate_1.0.5      grid_4.6.2          RColorBrewer_1.1-3 
    ## [10] fastmap_1.2.0       jsonlite_2.0.0      httr_1.4.9         
    ## [13] purrr_1.2.2         crosstalk_1.2.2     viridisLite_0.4.3  
    ## [16] scales_1.4.0        isoband_0.3.0       textshaping_1.0.5  
    ## [19] jquerylib_0.1.4     cli_3.6.6           rlang_1.3.0        
    ## [22] symengine_0.2.14    withr_3.0.3         cachem_1.1.0       
    ## [25] yaml_2.3.12         otel_0.2.0          tools_4.6.2        
    ## [28] dplyr_1.2.1         vctrs_0.7.3         R6_2.6.1           
    ## [31] lifecycle_1.0.5     fs_2.1.0            htmlwidgets_1.6.4  
    ## [34] MASS_7.3-66         ragg_1.5.2          pkgconfig_2.0.3    
    ## [37] desc_1.4.3          pkgdown_2.2.1       pillar_1.11.1      
    ## [40] bslib_0.12.0        gtable_0.3.6        glue_1.8.1         
    ## [43] data.table_1.18.6.1 Rcpp_1.1.2          systemfonts_1.3.2  
    ## [46] xfun_0.62           tibble_3.3.1        tidyselect_1.2.1   
    ## [49] knitr_1.52          farver_2.1.2        htmltools_0.5.9    
    ## [52] rmarkdown_2.32      labeling_0.4.3      compiler_4.6.2     
    ## [55] S7_0.2.2
