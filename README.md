# lpplsF: Log-Periodic Power Law Singularity Model for Financial Bubbles

<!-- badges: start -->
[![R-CMD-check](https://github.com/mahovo/lpplsF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mahovo/lpplsF/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`lpplsF` implements the Log-Periodic Power Law Singularity (LPPLS) model for detecting and analysing financial bubbles. The package is based on:

1. **Filimonov & Sornette (2013)**: A stable and robust calibration scheme of the log-periodic power law model. *Physica A: Statistical Mechanics and its Applications*.

2. **Filimonov, Demos & Sornette (2017)**: Modified profile likelihood inference and interval forecast of the burst of financial bubbles. *Quantitative Finance*.

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("mahovo/lpplsF")

# Once on CRAN
# install.packages("lpplsF")
```

## Quick Start

```r
library(lpplsF)

# Generate synthetic bubble data
set.seed(123)
t <- 1:200
true_params <- list(A = 4, B = -0.015, C1 = 0.002, C2 = 0.001,
                    tc = 250, m = 0.5, omega = 9)
log_p <- eval_lppls(t, true_params$A, true_params$B, true_params$C1,
                    true_params$C2, true_params$tc, true_params$m,
                    true_params$omega) + rnorm(200, 0, 0.01)

# Fit the LPPLS model
result <- fit_lppls(
  time_ID = t,
  log_price = log_p,
  fh = 100,           # Forecast horizon
  hold_out = 0,       # No holdout data
  tc_init = 250,      # Initial guess for critical time
  mode = "F2",        # Two-stage optimization
  num_searches = 20,  # Number of random starts
  fp = TRUE           # Generate fit plot
)

# View results
print(result$fit[[1]])

# Display fit plot
print(result$fit_plot)
```

## Model Description

The LPPLS model describes log-price dynamics approaching a critical point:

$$y(t) = A + (t_c - t)^m [B + C_1 \cos(\omega \log(t_c-t)) + C_2 \sin(\omega \log(t_c-t))]$$


where:  

- $t_c$ is the critical time (bubble peak)  
- $m$ is the power law exponent (typically 0.1-0.9)  
- $\omega$ (omega) is the log-periodic frequency (typically 6-13)  
- $A$, $B$, $C_1$, $C_2$ are linear parameters  

The damping coefficient `D = m|B| / (ω√(C1² + C2²))` measures the strength of oscillations relative to power law growth. Filimonov (2017) suggests D ≥ 0.8 for valid bubble detection.

## Optimization Modes

| Mode | Description |
|------|-------------|
| `F1` | Simultaneous optimization of tc, m, ω |
| `F2` | Two-stage: optimize (m, ω) for each tc, then select best tc |
| `MPL` | Modified Profile Likelihood with confidence intervals for tc |

## Main Functions

| Function | Description |
|----------|-------------|
| `fit_lppls()` | Main fitting function |
| `eval_lppls()` | Compute LPPLS model values |
| `SSE()` | Compute sum of squared errors |
| `fit_plot()` | Generate fit visualization |
| `contour_plot_sse()` | Generate SSE contour plots |

## Requirements

- R ≥ 4.0.0
- ggplot2
- dplyr
- tibble
- tidyr
- symengine
- plotly

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md) for development workflow using Doom Emacs.

## License

MIT License

## Citation

If you use this package in your research, please cite:

```bibtex
@Manual{lpplsF,
  title = {lpplsF: Log-Periodic Power Law Singularity Model for Financial Bubbles},
  author = {Martin Hoshi Vognsen},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/mahovo/lpplsF}
}
```
