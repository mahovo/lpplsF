# Create Custom Contour Plot for SSE

Generate a contour plot of the sum of squared errors as a function of
two parameters, with the remaining parameters fixed.

## Usage

``` r
sse_contour_plot(
  log_p,
  t,
  par,
  vars,
  lower = c(0, -0.001),
  upper = c(1, 0.001),
  cp = TRUE,
  sp = FALSE,
  fb = FALSE,
  point = c(NA, NA),
  mode = 0
)
```

## Arguments

- log_p:

  Numeric vector of log-prices.

- t:

  Numeric vector of time indices.

- par:

  Named list of fixed parameter values.

- vars:

  Named list specifying x and y variables as character strings.

- lower:

  Numeric vector of length 2, lower bounds for x and y axes.

- upper:

  Numeric vector of length 2, upper bounds for x and y axes.

- cp:

  Logical, whether to generate contour plot (default TRUE).

- sp:

  Logical, whether to generate surface plot (default FALSE).

- fb:

  Logical, whether to print feedback (default FALSE).

- point:

  Numeric vector of length 2, optional point to highlight.

- mode:

  Integer, LPPLS model mode (default 0).

## Value

A list with three elements:

- contour_plot:

  A `ggplot2` contour plot, or `NULL` when `cp = FALSE`.

- surface_plot:

  A `plotly` 3-D surface, or `NULL` when `sp = FALSE`.

- contour_data:

  The SSE grid itself, always returned: a list with the grid values `x`
  and `y` of the two parameters and the matrix `z` of SSE values, with
  `z[j, i]` the SSE at `(x[i], y[j])`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create contour plot for m vs omega
par <- list(tc = 250, A = 4, B = -0.015, C1 = 0.002, C2 = 0.001)
vars <- list(x = "m", y = "omega")
result <- sse_contour_plot(log_p, t, par, vars,
                           lower = c(0.1, 6), upper = c(0.9, 13))
print(result$contour_plot)
} # }
```
