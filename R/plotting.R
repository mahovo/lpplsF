#' Which nonlinear parameters are pinned at an optimisation bound
#'
#' Returns the names (`"m"` / `"omega"`) of the nonlinear parameters sitting at
#' (within `tol` of) their `lower`/`upper` bound. A non-empty result flags a
#' boundary solution, for which the information matrix is (near-)singular and the
#' modified profile likelihood / likelihood intervals are unreliable.
#'
#' @keywords internal
#' @noRd
.boundary_pars <- function(fit, lower, upper, tol = 1e-3) {
  at <- function(x, lo, hi) (x - lo) <= tol * (hi - lo) || (hi - x) <= tol * (hi - lo)
  pinned <- character(0)
  if (at(fit$m,     lower[1], upper[1])) pinned <- c(pinned, "m")
  if (at(fit$omega, lower[2], upper[2])) pinned <- c(pinned, "omega")
  pinned
}


#' How many fits in a per-fit table sit at an optimisation bound
#'
#' Vectorised companion to `.boundary_pars()`, applied to a `fit[[2]]` table: in
#' `"F2"`/`"MPL"` mode its rows are the per-`tc` conditional fits, in `"F1"` mode
#' the random restarts. A point estimate away from the bounds can still rest on a
#' profile that is pinned almost everywhere, which is what this counts.
#'
#' @return Named integer vector, `c(m = , omega = )`.
#'
#' @keywords internal
#' @noRd
.boundary_counts <- function(tbl, lower, upper, tol = 1e-3) {
  at <- function(x, lo, hi) (x - lo) <= tol * (hi - lo) | (hi - x) <= tol * (hi - lo)
  c(m     = sum(at(tbl$m,     lower[1], upper[1]), na.rm = TRUE),
    omega = sum(at(tbl$omega, lower[2], upper[2]), na.rm = TRUE))
}





#' Modified Profile Likelihood Functions and Plotting
#'
#' @keywords internal
#' @noRd

#' Compute Modified Profile Likelihood
#'
#' Calculates the Modified Profile Likelihood and confidence intervals
#' for the critical time tc, following Filimonov et al. (2017).
#'
#' @param fit List of fitted model results from fit_lppls().
#' @param log_p Numeric vector of log-prices.
#' @param t Numeric vector of time indices.
#' @param fh Integer, forecast horizon.
#' @param cutoff Numeric vector of cutoff levels for confidence intervals.
#' @param lower,upper Numeric bounds passed to fit_lppls(); used to flag a
#'   boundary solution (m or omega pinned at a bound) for which the MPL is
#'   unreliable.
#' @param fb Logical, whether to print feedback.
#' @param beta_calculator Function from [create_beta_calculator()]. Passed in
#'   rather than built here so the solver choice is honoured and the calculator
#'   is constructed once per fit instead of once per candidate `tc`.
#'
#' @return List with:
#'   \itemize{
#'     \item LI: List of likelihood intervals for each cutoff
#'     \item R: Numeric vector of relative likelihood values
#'     \item LL: Numeric vector of log-likelihood values
#'     \item MLL: Maximum log-likelihood
#'     \item tc_hat_mpl: tc value at maximum likelihood
#'   }
#'
#' @keywords internal
#' @noRd
compute_mpl <- function(fit, log_p, t, fh, cutoff, lower, upper, fb = FALSE,
                        beta_calculator = create_beta_calculator()) {
  n <- length(t)
  log_cutoff <- log(cutoff)

  ## Initialize results
  R_tbl <- data.frame(R = numeric(fh), LL = numeric(fh))

  ## Get profile likelihood estimate
  par_hat <- c(
    fit[[1]]$tc,
    fit[[1]]$m,
    fit[[1]]$omega,
    fit[[1]]$A,
    fit[[1]]$B,
    fit[[1]]$C1,
    fit[[1]]$C2
  )

  pinned <- .boundary_pars(fit[[1]], lower, upper)
  if (length(pinned) > 0) {
    warning("LPPLS fit is a boundary solution (", paste(pinned, collapse = ", "),
            " at the optimisation bound); the modified profile likelihood and its ",
            "likelihood intervals are unreliable for this fit.")
  }

  ## Compute log-likelihood for each tc value
  for (i in seq_len(fh)) {
    if (fb && i %% 10 == 1) {
      message(sprintf("Computing LL for tc = %d to %d", n + i, min(n + i + 9, n + fh)))
    }

    ## Get fit for tc = n + i
    fit_i <- fit[[2]] %>% dplyr::filter(tc == n + i)

    if (nrow(fit_i) == 0) {
      ## Use unfiltered fits if filtered is empty
      if (length(fit) >= 3) {
        fit_i <- fit[[3]] %>% dplyr::filter(tc == n + i)
      }
    }

    if (nrow(fit_i) > 0) {
      par_hat_tc <- c(
        fit_i$m[1],
        fit_i$omega[1],
        fit_i$A[1],
        fit_i$B[1],
        fit_i$C1[1],
        fit_i$C2[1]
      )

      R_tbl$LL[i] <- compute_mpl_loglik(
        tc = n + i,
        tc_hat = par_hat[1],
        Psi_hat_tc = par_hat_tc,
        Psi_hat = par_hat[2:7],
        log_p = log_p,
        t = t,
        beta_calc = beta_calculator
      )
    } else {
      R_tbl$LL[i] <- NA
    }
  }

  ## If the modified profile likelihood is undefined for every tc (e.g. an
  ## indefinite information matrix at a boundary fit, so det1 <= 0 throughout),
  ## return NA estimates rather than failing.
  if (all(is.na(R_tbl$LL))) {
    warning("MPL log-likelihood is NA for all tc (degenerate/boundary fit); returning NA estimates.")
    return(list(
      LI = replicate(length(cutoff), c(NA_real_, NA_real_), simplify = FALSE),
      R = R_tbl$R, LL = R_tbl$LL, MLL = NA_real_, tc_hat_mpl = NA_real_
    ))
  }

  ## Calculate maximum and relative likelihood
  MLL <- max(R_tbl$LL, na.rm = TRUE)
  R_tbl$R <- R_tbl$LL - MLL

  tc_hat_mpl <- n + which.max(R_tbl$LL)

  ## Compute likelihood intervals
  tc_range <- (n + 1):(n + fh)
  LI <- list()

  for (j in seq_along(cutoff)) {
    LI_IDs <- which(log_cutoff[j] < R_tbl$R)
    LI_vect <- tc_range[LI_IDs]

    if (length(LI_vect) > 0) {
      LI[[j]] <- range(LI_vect)
    } else {
      warning(sprintf("No values of tc in likelihood interval for cutoff %.2f", cutoff[j]))
      LI[[j]] <- c(NA, NA)
    }
  }

  list(
    LI = LI,
    R = R_tbl$R,
    LL = R_tbl$LL,
    MLL = MLL,
    tc_hat_mpl = tc_hat_mpl
  )
}


#' Compute MPL Log-Likelihood for Single tc Value
#'
#' @keywords internal
#' @noRd
compute_mpl_loglik <- function(tc, tc_hat, Psi_hat_tc, Psi_hat, log_p, t,
                               beta_calc = create_beta_calculator()) {
  n <- length(t)
  p <- 6

  beta_vals <- beta_calc(log_p, t, tc, Psi_hat_tc[1], Psi_hat_tc[2])

  ## Compute fitted values and SSE
  fitted <- eval_lppls(
    t,
    beta_vals["A"],
    beta_vals["B"],
    beta_vals["C1"],
    beta_vals["C2"],
    tc, Psi_hat_tc[1],
    Psi_hat_tc[2],
    mode = 0
  )
  s_tc <- sum((log_p - fitted)^2, na.rm = TRUE) / n

  ## Compute matrices
  X_hat_tc <- compute_X_matrix(Psi_hat_tc, tc, t)
  X_hat <- compute_X_matrix(Psi_hat, tc_hat, t)
  H_hat_tc <- compute_H_matrix(Psi_hat_tc, tc, log_p, t)

  ## Compute log-likelihood
  ## One-argument crossprod() for X'X: bit-identical to crossprod(X, X) but it
  ## takes the symmetric-rank-update BLAS path, so it is cheaper and returns an
  ## exactly symmetric matrix -- which matters here, since det(X'X - H) is only
  ## the observed information when both terms are symmetric. The cross term
  ## below is a genuine two-matrix product and has no such shortcut.
  XtX_tc <- crossprod(X_hat_tc)
  XtX_cross <- crossprod(X_hat, X_hat_tc)

  det1 <- det(XtX_tc - H_hat_tc)
  det2 <- det(XtX_cross)

  if (is.na(det1) || is.na(det2) || det2 == 0 || det1 <= 0) {
    return(NA)
  }

  log(sqrt(abs(det1)) / abs(det2)) - ((n - p - 2) / 2) * log(s_tc)
}


#' Compute X Matrix for MPL
#'
#' Gradient matrix at parameter estimates.
#' Filimonov (2017), equation (36).
#'
#' @keywords internal
#' @noRd
compute_X_matrix <- function(Psi, tc, t) {
  n <- length(t)
  p <- 6 ## length(Psi) = 6
  X <- matrix(0, nrow = n, ncol = p)


  ## Psi = (m, omega, a, b, c1, c2)
  m <- Psi[1]
  omega <- Psi[2]
  # a <- Psi[3]  ## unused in gradient
  b <- Psi[4]
  c1 <- Psi[5]
  c2 <- Psi[6]

  ## The basis, for every t_i at once. Column-wise rather than row-wise: this
  ## is called 2 x fh times per MPL fit, and the loop it replaces cost ~7.5x
  ## more for bit-identical output.
  tau <- tc - t
  log_tau <- log(tau)
  tau_m <- tau^m
  cos_term <- cos(omega * log_tau)
  sin_term <- sin(omega * log_tau)

  ## Gradient, one column per parameter
  ## Filimonov2017, equation (B16)

  ## d LPPLS(t_i; tc, psi) / d m
  X[, 1] <- tau_m * log_tau * (b + c1 * cos_term + c2 * sin_term)
  ## d LPPLS(t_i; tc, psi) / d omega
  X[, 2] <- tau_m * log_tau * (-c1 * sin_term + c2 * cos_term)
  ## d LPPLS(t_i; tc, psi) / d a
  X[, 3] <- 1
  ## d LPPLS(t_i; tc, psi) / d b
  X[, 4] <- tau_m
  ## d LPPLS(t_i; tc, psi) / d c1
  X[, 5] <- tau_m * cos_term
  ## d LPPLS(t_i; tc, psi) / d c2
  X[, 6] <- tau_m * sin_term

  X
}


#' Compute H Matrix for MPL
#'
#' H matrix for modified profile likelihood.
#' Filimonov (2017), equation (37).
#'
#' H holds residual-weighted second derivatives, so it is symmetric; both
#' triangles must be filled or `det(X'X - H)` is not the observed information.
#'
#' @keywords internal
#' @noRd
compute_H_matrix <- function(Psi, tc, log_p, t) {
  p <- 6 ## length(Psi) = 6
  H <- matrix(0, nrow = p, ncol = p)

  m <- Psi[1]
  omega <- Psi[2]
  a <- Psi[3]
  b <- Psi[4]
  c1 <- Psi[5]
  c2 <- Psi[6]

  ## The basis, for every t_i at once; each entry of H is then one sum over t
  ## rather than an accumulation in an R loop. Called once per candidate tc,
  ## and ~10.4x cheaper than the loop it replaces, for bit-identical output.
  tau <- tc - t
  log_tau <- log(tau)
  tau_m <- tau^m
  cos_term <- cos(omega * log_tau)
  sin_term <- sin(omega * log_tau)

  ## LPPLS value
  lppls_val <- a + tau_m * (b + c1 * cos_term + c2 * sin_term)
  res <- log_p - lppls_val

  ## Second derivatives, summed over t
  ## (m, m)
  H[1, 1] <- sum(res * tau_m * log_tau^2 * (b + c1 * cos_term + c2 * sin_term))
  ## (m, omega)
  H[1, 2] <- sum(res * tau_m * log_tau^2 * (-c1 * sin_term + c2 * cos_term))
  ## (m, b)
  H[1, 4] <- sum(res * tau_m * log_tau)
  ## (m, c1)
  H[1, 5] <- sum(res * tau_m * log_tau * cos_term)
  ## (m, c2)
  H[1, 6] <- sum(res * tau_m * log_tau * sin_term)
  ## (omega, omega)
  H[2, 2] <- sum(res * (-1) * tau_m * log_tau^2 * (c1 * cos_term + c2 * sin_term))
  ## (omega, c1)
  H[2, 5] <- sum(res * (-1) * tau_m * log_tau * sin_term)
  ## (omega, c2)
  H[2, 6] <- sum(res * tau_m * log_tau * cos_term)

  ## Symmetry

  ## (omega, m)
  H[2, 1] <- H[1, 2]

  ## (b, m)
  H[4, 1] <- H[1, 4]

  ## (c1, m)
  H[5, 1] <- H[1, 5]

  ## (c2, m)
  H[6, 1] <- H[1, 6]

  ## (c1, omega)
  H[5, 2] <- H[2, 5]

  ## (c2, omega)
  H[6, 2] <- H[2, 6]

  ## The rest is zero

  H
}


#' Markers shared by the MPL plots
#'
#' Collects the vertical markers into one data frame so colour and linetype map
#' to a single legend: the PLE / MPLE point estimates and the LI bounds (each
#' interval contributes two lines but one legend entry). Returned with the
#' colour and linetype of every series, so the absolute and the relative MPL
#' plot draw them identically.
#'
#' Each interval is named for its cutoff as a percentage, `LI_05` for 0.05, as
#' in [noise_tc_plot()]. Colours and linetypes go to the intervals in order, so
#' the default cutoffs keep their green, orange and red.
#'
#' @param cutoff Numeric, the likelihood-interval cutoffs `mpl_output` was
#'   computed with.
#'
#' @keywords internal
#' @noRd
.mpl_markers <- function(mpl_output, fit, cutoff = c(0.05, 0.1, 0.5)) {
  vmarks <- data.frame(
    xintercept = c(fit[[1]]$tc, mpl_output$tc_hat_mpl),
    series = c("PLE", "MPLE")
  )
  ## Cutoffs that do not match the intervals (a fit made without recording
  ## them) can only be numbered, as print.lppls_fit() does.
  li_names <- if (length(cutoff) == length(mpl_output$LI)) {
    sprintf("LI_%02g", 100 * cutoff)
  } else {
    sprintf("LI[%d]", seq_along(mpl_output$LI))
  }
  k <- length(li_names)
  for (i in seq_along(mpl_output$LI)) {
    li <- mpl_output$LI[[i]]
    if (!any(is.na(li))) {
      vmarks <- rbind(
        vmarks,
        data.frame(xintercept = c(li[1], li[2]), series = li_names[i])
      )
    }
  }

  list(
    vmarks   = vmarks,
    li_names = li_names,
    levels   = c("PLE", "MPLE", li_names),
    cols     = c(PLE = "blue", MPLE = "blue",
                 stats::setNames(rep_len(c("green3", "orange", "red", "purple",
                                           "brown", "darkcyan"), k), li_names)),
    ltys     = c(PLE = "solid", MPLE = "dashed",
                 stats::setNames(rep_len(c("dashed", "dotted", "dotdash",
                                           "longdash", "twodash", "dashed"), k),
                                 li_names))
  )
}


#' Create MPL Plot
#'
#' @keywords internal
#' @noRd
create_mpl_plot <- function(mpl_output, fit, n, fh, cutoff = c(0.05, 0.1, 0.5)) {
  tc_seq <- (n + 1):(n + fh)

  df_plot <- data.frame(
    tc = tc_seq,
    log_mpl = mpl_output$LL
  )

  mk      <- .mpl_markers(mpl_output, fit, cutoff)
  vmarks  <- mk$vmarks
  present <- mk$levels[mk$levels %in% vmarks$series]
  vmarks$series <- factor(vmarks$series, levels = present)

  ggplot2::ggplot(df_plot, ggplot2::aes(x = tc, y = log_mpl)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      data = vmarks,
      mapping = ggplot2::aes(xintercept = xintercept,
                             color = series, linetype = series),
      key_glyph = "path"
    ) +
    ggplot2::scale_color_manual(name = NULL, values = mk$cols, breaks = present) +
    ggplot2::scale_linetype_manual(name = NULL, values = mk$ltys, breaks = present) +
    ggplot2::labs(
      x = "Critical Time (tc)",
      y = "Log Modified Profile Likelihood",
      title = "Modified Profile Likelihood"
    ) +
    ggplot2::theme_minimal() +
    overlay_legend_theme(
      inside = c(0.99, 0.99),
      justification = c(1, 1)
    )
}


#' Create the relative MPL plot
#'
#' The modified profile likelihood relative to its maximum,
#' `R = log MPL - max log MPL`, zoomed to the band the likelihood-interval
#' cutoffs live in. Its absolute level carries no information, whereas `R` is
#' what the intervals test: each interval is where `R` stays above
#' `log(cutoff)`. The cutoffs are drawn as horizontal lines, so an interval can
#' be read off where the curve crosses one, and a cutoff the curve never reaches
#' still shows that.
#'
#' @param cutoff Numeric, the likelihood-interval cutoffs of the fit.
#' @param pad Numeric, how far below the lowest cutoff the view extends.
#'
#' @keywords internal
#' @noRd
create_mpl_relative_plot <- function(mpl_output, fit, n, cutoff, pad = 1.2) {
  mk       <- .mpl_markers(mpl_output, fit, cutoff)
  vmarks   <- mk$vmarks
  h_series <- mk$li_names[seq_along(cutoff)]
  present  <- mk$levels[mk$levels %in% c(vmarks$series, h_series)]
  vmarks$series <- factor(vmarks$series, levels = present)
  hmarks <- data.frame(yintercept = log(cutoff),
                       series = factor(h_series, levels = present))

  df_plot <- data.frame(tc = (n + 1):(n + length(mpl_output$R)),
                        R  = mpl_output$R)

  ggplot2::ggplot(df_plot, ggplot2::aes(x = tc, y = R)) +
    ggplot2::geom_point(size = 0.5, na.rm = TRUE) +
    ggplot2::geom_hline(
      data = hmarks,
      mapping = ggplot2::aes(yintercept = yintercept,
                             color = series, linetype = series),
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      data = vmarks,
      mapping = ggplot2::aes(xintercept = xintercept,
                             color = series, linetype = series),
      key_glyph = "path"
    ) +
    ggplot2::scale_color_manual(name = NULL, values = mk$cols, breaks = present) +
    ggplot2::scale_linetype_manual(name = NULL, values = mk$ltys, breaks = present) +
    ggplot2::coord_cartesian(ylim = c(min(log(cutoff)) - pad, 0.4)) +
    ggplot2::labs(
      x = "Critical Time (tc)",
      y = "Relative log MPL",
      title = "Modified Profile Likelihood (relative)"
    ) +
    ggplot2::theme_minimal() +
    overlay_legend_theme(
      inside = c(0.99, 0.99),
      justification = c(1, 1)
    )
}


#' Compact semi-transparent legend overlaid in the panel's upper-left corner
#'
#' A theme fragment that places the legend inside the plotting region (top-left)
#' instead of in the right margin, so it does not steal horizontal space.
#'
#' @keywords internal
#' @noRd
overlay_legend_theme <- function(
    inside = c(0.01, 0.99),
    justification = c(0, 1)
  ) {
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = inside,
    legend.justification.inside = justification,
    legend.title = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.6), color = NA),
    legend.key = ggplot2::element_rect(fill = NA, color = NA),
    legend.margin = ggplot2::margin(2, 4, 2, 2)
  )
}


#' Create Contour Data
#'
#' Generate SSE data for contour plot of SSE wrt m and omega.
#'
#' @keywords internal
#' @noRd
create_contour_data <- function(fit, log_p, t, lower, upper, sse_calculator,
                                fb = FALSE) {
  if (fb) message("Generating contour plot data...")

  tc_val <- fit$tc

  x_contour <- seq(lower[1], upper[1], length.out = 101)  ## m
  y_contour <- seq(lower[2], upper[2], length.out = 101)  ## omega

  ## The same objective the search minimised, so the surface and the reported
  ## optimum are computed the same way.
  sse_func <- function(m_val, omega_val) {
    sse_calculator(log_p, t, tc_val, m_val, omega_val)
  }

  z_contour <- outer(x_contour, y_contour, Vectorize(sse_func))
  z_contour <- t(z_contour)  ## Transpose for correct orientation, as expected by
                             ## contour plot

  list(x = x_contour, y = y_contour, z = z_contour, tc = tc_val)
}


#' Create Contour Plot
#'
#' @keywords internal
#' @noRd
create_contour_plot <- function(contour_data, fit, lower) {
  z_coord <- contour_data$z[
    which.min(abs(contour_data$y - fit$omega)),
    which.min(abs(contour_data$x - fit$m))
  ]

  plotly::plot_ly(
    x = contour_data$x,
    y = contour_data$y,
    z = contour_data$z,
    type = "contour",
    contours = list(showlabels = FALSE, exponentformat = "e"),
    ncontours = 10,
    colorbar = list(exponentformat = "e", title = "SSE")
  ) %>%
    plotly::layout(
      title = sprintf("SSE Contour (tc = %.2f)", contour_data$tc),
      xaxis = list(title = "m"),
      yaxis = list(title = "omega")
    ) %>%
    plotly::add_trace(
      x = fit$m,
      y = fit$omega,
      type = "scatter",
      mode = "markers",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      x = c(fit$m, fit$m, lower[1]),
      y = c(lower[2], fit$omega, fit$omega),
      type = "scatter",
      mode = "lines",
      line = list(dash = "dash", width = 1, color = "orange"),
      showlegend = FALSE
    ) %>%
    plotly::add_annotations(
      x = fit$m,
      y = fit$omega,
      xanchor = "left",
      yanchor = "bottom",
      text = sprintf("(%.3f, %.2f, %.4e)", fit$m, fit$omega, z_coord),
      font = list(color = "white", size = 10),
      showarrow = FALSE
    )
}


#' Create Surface Plot
#'
#' @keywords internal
#' @noRd
create_surface_plot <- function(contour_data, fit) {
  plotly::plot_ly(
    x = contour_data$x,
    y = contour_data$y,
    z = contour_data$z
  ) %>%
    plotly::add_surface(
      contours = list(
        z = list(
          show = TRUE,
          usecolormap = TRUE,
          highlightcolor = "#ff0000",
          project = list(z = TRUE)
        )
      ),
      colorbar = list(exponentformat = "e", title = "SSE")
    ) %>%
    plotly::layout(
      title = sprintf("SSE Surface (tc = %.2f)", contour_data$tc),
      scene = list(
        camera = list(eye = list(x = 1.1, y = -1.3, z = 0.1)),
        xaxis = list(title = "m"),
        yaxis = list(title = "omega"),
        zaxis = list(title = "SSE")
      )
    )
}


#' Create Parameter Plot
#'
#' @keywords internal
#' @noRd
create_param_plot <- function(fit) {
  plot_df <- tidyr::gather(
    data = fit[[2]],
    key = "param",
    value = "estimate",
    -c(ID, tc)
  )

  tc_hat <- fit[[1]]$tc

  ggplot2::ggplot(plot_df, ggplot2::aes(x = tc, y = estimate)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = tc_hat, color = "red") +
    ggplot2::facet_wrap(~param, scales = "free_y", ncol = 1) # +
    # ggplot2::theme_minimal()
}


#' Create Matrix Plot
#'
#' @keywords internal
#' @noRd
create_matrix_plot <- function(n, fh, lower, upper, log_p, t, SSE2_func) {
  tc_seq <- (n + 1):(n + fh)
  m_seq <- c(lower[1], (lower[1] + upper[1]) / 2, upper[1])
  omega_seq <- c(lower[2], (lower[2] + upper[2]) / 2, upper[2])

  results <- expand.grid(
    tc = tc_seq,
    m = m_seq,
    omega = omega_seq
  )
  results$sse <- NA

  for (i in seq_len(nrow(results))) {
    results$sse[i] <- SSE2_func(
      c(results$m[i], results$omega[i]),
      results$tc[i],
      log_p,
      t
    )
  }

  ggplot2::ggplot(results, ggplot2::aes(x = tc, y = sse)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(omega ~ m, labeller = "label_both") +
    ggplot2::theme_minimal()
}


#' Create Trace Plot
#'
#' Visualizes the L-BFGS-B optimization path over the SSE surface in
#' `(m, omega)`, the objective the search minimises: at each point the linear
#' parameters are solved by least squares. `B` and the damping `D` are not
#' search coordinates; they follow from that least-squares solution, so they are
#' drawn as regions. Shaded are the parts of the surface where the solution would
#' be out of range, `B >= upper[3]` in black and `D < lower[4]` in white, each
#' edged with a dashed line.
#'
#' @param opt_trace Matrix with 2 rows (m, omega) and one column per step.
#' @param lower Numeric vector of length 4, parameter lower bounds.
#' @param upper Numeric vector of length 4, parameter upper bounds.
#' @param SSE2_func SSE2 function (fixed tc) from fit_lppls.
#' @param beta_calc Beta calculator function from fit_lppls.
#' @param tc_val Numeric, the tc value for the trace.
#' @param log_p Numeric vector of log-prices.
#' @param t Numeric vector of time indices.
#' @param fb Logical, whether to print feedback.
#'
#' @return A ggplot2 object.
#'
#' @keywords internal
#' @noRd
create_trace_plot <- function(
    opt_trace,
    lower,
    upper,
    SSE2_func,
    beta_calc,
    tc_val,
    log_p,
    t,
    fb = FALSE
  ) {
  if (fb) message("Generating trace plot...")

  lattice_dim <- 100

  ## m vs omega lattice
  lattice <- as.matrix(expand.grid(
    m = seq(
      min(lower[1], min(opt_trace[1, ])),
      max(upper[1], max(opt_trace[1, ])),
      length.out = lattice_dim
    ),
    omega = seq(
      min(lower[2], min(opt_trace[2, ])),
      max(upper[2], max(opt_trace[2, ])),
      length.out = lattice_dim
    )
  ))
  lattice_vals <- apply(lattice, 1, SSE2_func,
                        tc = tc_val, log_p = log_p, t = t)

  lattice <- as.data.frame(lattice)
  lattice$vals <- lattice_vals

  ## The least-squares B and damping D at each lattice point. fit_lppls() checks
  ## both on the result of a search rather than bounding the search, and reports
  ## a fit with B >= upper[3] or D < lower[4] as out of range.
  beta <- vapply(seq_len(nrow(lattice)), function(i) {
    beta_calc(log_p, t, tc_val, lattice$m[i], lattice$omega[i])
  }, numeric(4))
  lattice$B <- beta["B", ]
  lattice$D <- calculate_damping(lattice$m, lattice$B, lattice$omega,
                                 beta["C1", ], beta["C2", ])

  regions <- list(
    "B out of range" = list(var = "B", at = upper[3], fill = "black",
                            inside = which(lattice$B >= upper[3])),
    "D out of range" = list(var = "D", at = lower[4], fill = "white",
                            inside = which(lattice$D < lower[4]))
  )

  opt_df <- data.frame(m = opt_trace[1, ], omega = opt_trace[2, ])
  opt_df$step <- seq_len(nrow(opt_df))

  p <- ggplot2::ggplot(lattice, ggplot2::aes(x = m, y = omega)) +
    ggplot2::geom_raster(ggplot2::aes(fill = vals)) +
    ggplot2::geom_contour(ggplot2::aes(z = vals), col = "white")

  ## Shade each non-empty region and edge it, unless it fills the whole lattice
  ## and so has no edge to draw. The edge carries the legend.
  dm <- diff(range(lattice$m)) / (lattice_dim - 1)
  dw <- diff(range(lattice$omega)) / (lattice_dim - 1)
  edged <- character(0)
  for (r in names(regions)) {
    g <- regions[[r]]
    if (length(g$inside) == 0L) next
    p <- p + ggplot2::geom_tile(data = lattice[g$inside, ], width = dm,
                                height = dw, fill = g$fill, alpha = 0.4)
    if (length(g$inside) < nrow(lattice)) {
      p <- p + ggplot2::geom_contour(
        ggplot2::aes(z = !!as.name(g$var), colour = !!r),
        breaks = g$at, linetype = "dashed", linewidth = 0.7
      )
      edged <- c(edged, r)
    }
  }
  if (length(edged) > 0L) {
    p <- p +
      ggplot2::scale_colour_manual(
        name = NULL, values = vapply(regions, `[[`, "", "fill")
      ) +
      ggplot2::theme(legend.key = ggplot2::element_rect(fill = "grey50"))
  }

  p <- p +
    ggplot2::geom_path(data = opt_df, col = "orange") +
    ggplot2::geom_point(data = opt_df, col = "orange") +
    ggplot2::geom_point(
      data = utils::head(opt_df, 1),
      shape = 21,
      fill = "red",
      color = "orange",
      size = 3
    ) +
    ggplot2::geom_point(
      data = utils::tail(opt_df, 1),
      shape = 23,
      fill = "green",
      color = "orange",
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::coord_cartesian(expand = FALSE)

  if (fb) message("Trace plot generated")
  p
}


#' Plot a Fitted LPPLS Model
#'
#' Draws the observed log-price series together with the fitted LPPLS curve, a
#' vertical marker at the estimated critical time \eqn{\hat{t}_c}{tc-hat}, and---when the
#' calibration end is supplied via `T2`---a marker at \eqn{T_2}. Colour and line
#' type are merged into a single legend (by default `sim data` / `LPPLS fit` / `T2` /
#' `tc_hat`). This is the same builder [fit_lppls()] attaches as its `fit_plot`
#' component (called there with `T2` set to the calibration-window end).
#'
#' @param fit A named list of point estimates with elements `A`, `B`, `C1`,
#'   `C2`, `tc`, `m`, `omega` (e.g. `fit_lppls(...)$fit[[1]]`).
#' @param time_ID Numeric vector of time indices for `log_price`.
#' @param log_price Numeric vector of observed log-prices.
#' @param mode Integer, the LPPLS model variant passed to [eval_lppls()] for the
#'   fitted curve (default 0, the Filimonov parameterization that [fit_lppls()]
#'   calibrates).
#' @param T2 Numeric or `NULL`. The calibration-window end. When supplied, a
#'   dashed marker is drawn there and added to the legend; when `NULL` (default)
#'   only the \eqn{\hat{t}_c}{tc-hat} marker is shown.
#' @param T1 Numeric or `NULL`. The calibration-window *start*. When supplied, a
#'   dashed marker is drawn there and the fitted curve begins at `T1` rather than
#'   at the first plotted observation. This is what a rolling calibration needs:
#'   the data before `T1` is drawn but was never fitted, so a curve extending
#'   back over it would misrepresent the fit. When `NULL` (default) the curve
#'   starts at the first plotted point.
#' @param dates `NULL` (default) or a `Date` vector the same length as
#'   `log_price`, strictly increasing. When given, the x axis is drawn as a
#'   calendar axis instead of the 1..n time index: the series, the fitted curve
#'   and both markers are all mapped through it. A forecast
#'   \eqn{\hat{t}_c}{tc-hat} lying past the end of the series is extrapolated at
#'   the series' median spacing, so treat such a date as approximate.
#' @param style One of `"simulated"` (default) or `"empirical"`, selecting what
#'   the colours are asked to say.
#'
#'   `"simulated"` draws the series in the foreground colour and the fitted curve
#'   as a single line, with the window markers in the legend. It suits a series
#'   the caller generated, where the markers are part of the specification being
#'   checked.
#'
#'   `"empirical"` draws the series in grey, the way an inactive control is: an
#'   observed price series is given and cannot be otherwise. The fitted curve
#'   splits at `T2` into the part supported by data and the part beyond it, and
#'   only the latter is red, because past `T2` the curve is a forecast and the
#'   colour is the warning. The markers stay out of the legend, leaving three
#'   entries. The split needs points on both sides of `T2`; without it, or with
#'   \eqn{\hat{t}_c}{tc-hat} inside the calibration window, a single line is
#'   drawn instead.
#' @param data_label `NULL` (default) or a string naming the observed series in
#'   the legend. The default follows `style`: `"sim data"` for `"simulated"`,
#'   `"data"` for `"empirical"`.
#' @param linewidth Width of the series and fitted lines, passed to
#'   [ggplot2::geom_line()]. Defaults to `0.5`, ggplot2's own default; a long
#'   daily series usually reads better thinner.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [fit_lppls()], [eval_lppls()]
#'
#' @examples
#' set.seed(123)
#' t <- 1:200
#' log_p <- eval_lppls(t, 4, -0.015, 0.002, 0.001, 250, 0.5, 9) +
#'   rnorm(200, 0, 0.01)
#' est <- list(A = 4, B = -0.015, C1 = 0.002, C2 = 0.001,
#'             tc = 250, m = 0.5, omega = 9)
#'
#' # default: the 1..n time index on the x axis
#' fit_plot(est, t, log_p)
#'
#' # the same fit on a calendar axis
#' fit_plot(est, t, log_p, dates = as.Date("2020-01-01") + t - 1)
#'
#' # an observed series: grey data, and the forecast past T2 marked in red
#' fit_plot(est, t, log_p, T2 = 180, style = "empirical",
#'          dates = as.Date("2020-01-01") + t - 1)
#'
#' @export
fit_plot <- function(fit, time_ID, log_price, mode = 0, T2 = NULL, T1 = NULL,
                     dates = NULL, style = c("simulated", "empirical"),
                     data_label = NULL, linewidth = 0.5) {
  style <- match.arg(style)
  dates <- .check_dates(dates, length(log_price))

  ## With `dates`, every x position is mapped through the calendar: the series,
  ## the fitted curve and the vertical markers. `tc` may sit past the end of the
  ## series -- a forecast critical time -- in which case its date is extrapolated
  ## (see .ids_to_dates()).
  xmap <- if (is.null(dates)) identity else function(i) .ids_to_dates(i, dates)

  ## The two styles differ in what the colours are asked to say.
  ##
  ## "simulated": the series is one the caller generated, so it is drawn in the
  ## foreground colour and the fitted curve is a single line. The window markers
  ## join the legend, because on a synthetic example they are part of the
  ## specification the reader is checking.
  ##
  ## "empirical": the series is given and cannot be otherwise, so it is drawn in
  ## grey, the way an inactive control is. The fitted curve splits at T2 into the
  ## part supported by data and the part beyond it, and only the latter is red:
  ## past T2 the curve is a forecast, and the colour is the warning. The markers
  ## stay out of the legend, which leaves three entries carrying the distinction
  ## that matters.
  if (is.null(data_label)) {
    data_label <- if (style == "empirical") "data" else "sim data"
  }
  empirical <- identical(style, "empirical")

  plot_data <- data.frame(ID = xmap(time_ID), log_p = log_price)

  ## The fitted curve spans the calibration window. Without T1 it starts at the
  ## first plotted point, which is right for a fit calibrated from the beginning
  ## of the series; with T1 it starts there instead, which is what a rolling
  ## calibration needs -- the earlier data is shown but was never fitted.
  fit_from <- if (is.null(T1)) min(time_ID) else T1
  fit_to   <- floor(fit$tc) - 1L
  curve_at <- function(ids) data.frame(
    ID    = xmap(ids),
    log_p = eval_lppls(ids, fit$A, fit$B, fit$C1, fit$C2,
                       fit$tc, fit$m, fit$omega, mode = mode)
  )

  ## Split the curve at T2 only when there is something on both sides of it.
  ## Without T2, or with the whole curve inside the calibration window, there is
  ## no forecast to distinguish and one line is the honest drawing.
  split_at <- if (empirical && !is.null(T2) && T2 > fit_from && T2 < fit_to) T2 else NULL

  fit_lab  <- if (empirical) "fit" else "LPPLS fit"
  pred_lab <- "prediction"

  lvls <- c(data_label, fit_lab,
            if (!is.null(split_at)) pred_lab,
            if (!empirical) c(if (!is.null(T1)) "T1", if (!is.null(T2)) "T2", "tc_hat"))
  cols <- stats::setNames(
    c(if (empirical) "darkgray" else "royalblue1",
      if (empirical) "royalblue" else "red",
      "red", "orange", "green", "red"),
    c(data_label, fit_lab, pred_lab, "T1", "T2", "tc_hat"))
  ltys <- stats::setNames(
    c("solid", "solid", "solid", "dashed", "dashed", "dashed"),
    c(data_label, fit_lab, pred_lab, "T1", "T2", "tc_hat"))

  x_lab <- if (is.null(dates)) "Time Index" else "Date"

  ## In the empirical style the markers are not in the legend, so every entry is
  ## a solid line and mapping linetype would add a scale carrying no information.
  ln_aes <- function(lab) {
    if (empirical) ggplot2::aes(color = lab)
    else ggplot2::aes(color = lab, linetype = lab)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(ID, log_p)) +
    ggplot2::geom_line(ln_aes(data_label), linewidth = linewidth)

  if (is.null(split_at)) {
    p <- p + ggplot2::geom_line(
      data = curve_at(seq.int(fit_from, fit_to)),
      mapping = ln_aes(fit_lab), linewidth = linewidth)
  } else {
    ## The segments share the point at T2, so the two lines meet without a gap.
    p <- p +
      ggplot2::geom_line(
        data = curve_at(seq.int(fit_from, split_at)),
        mapping = ln_aes(fit_lab), linewidth = linewidth) +
      ggplot2::geom_line(
        data = curve_at(seq.int(split_at, fit_to)),
        mapping = ln_aes(pred_lab), linewidth = linewidth)
  }

  if (empirical) {
    ## Markers as plain rules, outside the legend.
    if (!is.null(T1)) {
      p <- p + ggplot2::geom_vline(xintercept = xmap(T1), color = "orange",
                                   linetype = "dashed")
    }
    if (!is.null(T2)) {
      p <- p + ggplot2::geom_vline(xintercept = xmap(T2), color = "royalblue",
                                   linetype = "dashed")
    }
    p <- p + ggplot2::geom_vline(xintercept = xmap(fit$tc), color = "red",
                                 linetype = "dashed")
  } else {
    marks <- data.frame(xintercept = xmap(fit$tc), series = "tc_hat")
    if (!is.null(T2)) {
      marks <- rbind(data.frame(xintercept = xmap(T2), series = "T2"), marks)
    }
    if (!is.null(T1)) {
      marks <- rbind(data.frame(xintercept = xmap(T1), series = "T1"), marks)
    }
    p <- p + ggplot2::geom_vline(
      data = marks,
      mapping = ggplot2::aes(xintercept = xintercept, color = series,
                             linetype = series),
      inherit.aes = FALSE, key_glyph = "path")
  }

  p <- p + ggplot2::scale_color_manual(name = NULL, values = cols, breaks = lvls)
  if (!empirical) {
    p <- p + ggplot2::scale_linetype_manual(name = NULL, values = ltys,
                                            breaks = lvls)
  }

  p +
    ggplot2::labs(x = x_lab, y = "Log Price") +
    ggplot2::theme_minimal() +
    overlay_legend_theme()
}


#' Create Custom Contour Plot for SSE
#'
#' Generate a contour plot of the sum of squared errors as a function of
#' two parameters, with the remaining parameters fixed.
#'
#' @param log_p Numeric vector of log-prices.
#' @param t Numeric vector of time indices.
#' @param par Named list of fixed parameter values.
#' @param vars Named list specifying x and y variables as character strings.
#' @param lower Numeric vector of length 2, lower bounds for x and y axes.
#' @param upper Numeric vector of length 2, upper bounds for x and y axes.
#' @param cp Logical, whether to generate contour plot (default TRUE).
#' @param sp Logical, whether to generate surface plot (default FALSE).
#' @param fb Logical, whether to print feedback (default FALSE).
#' @param point Numeric vector of length 2, optional point to highlight.
#' @param mode Integer, LPPLS model mode (default 0).
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{contour_plot}{A `ggplot2` contour plot, or `NULL` when `cp = FALSE`.}
#'     \item{surface_plot}{A `plotly` 3-D surface, or `NULL` when `sp = FALSE`.}
#'     \item{contour_data}{The SSE grid itself, always returned: a list with the
#'       grid values `x` and `y` of the two parameters and the matrix `z` of SSE
#'       values, with `z[j, i]` the SSE at `(x[i], y[j])`.}
#'   }
#'
#' @examples
#' ## A simulated bubble whose critical time is t = 221
#' set.seed(4)
#' t     <- 1:200
#' log_p <- eval_lppls(t, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001,
#'                     tc = 221, m = 0.6, omega = 8) + rnorm(200, sd = 0.01)
#'
#' ## The SSE over (m, omega), with the other parameters held at their true
#' ## values
#' fixed  <- list(tc = 221, A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001)
#' result <- sse_contour_plot(log_p, t, fixed, vars = list(x = "m", y = "omega"),
#'                            lower = c(0.1, 6), upper = c(0.9, 13))
#' result$contour_plot
#' dim(result$contour_data$z)   # the grid itself
#'
#' @export
sse_contour_plot <- function(
    log_p,
    t,
    par,
    vars,
    lower = c(0.0, -0.001),
    upper = c(1.0, 0.001),
    cp = TRUE,
    sp = FALSE,
    fb = FALSE,
    point = c(NA, NA),
    mode = 0
  ) {

  if (length(log_p) != length(t)) {
    stop("Length of log_p and t must be the same")
  }

  all_pars <- c(.lppls_pars$core, .lppls_pars$higher)
  if (!is.list(par) || (length(par) && is.null(names(par)))) {
    stop("'par' must be a named list")
  }
  if (!is.list(vars) || !all(c("x", "y") %in% names(vars)) ||
      !is.character(vars$x) || !is.character(vars$y) ||
      length(vars$x) != 1L || length(vars$y) != 1L) {
    stop("'vars' must be a list with single character elements 'x' and 'y'")
  }
  if (!all(c(vars$x, vars$y) %in% all_pars)) {
    stop("'vars$x'/'vars$y' must name eval_lppls parameters: ",
        paste(all_pars, collapse = ", "))
  }
  if (any(c(vars$x, vars$y) %in% names(par)))                        ## can't both fix and scan a param
    stop("scanned variable(s) also fixed in 'par': ",
         paste(intersect(c(vars$x, vars$y), names(par)), collapse = ", "))
  if (!is.numeric(lower) || !is.numeric(upper) ||
        length(lower) != 2L || length(upper) != 2L) {
    stop("'lower' and 'upper' must be numeric vectors of length 2")
  }
  if (any(lower >= upper)) {
    stop("each 'lower' bound must be strictly less than its 'upper'")
  }
  if (!isTRUE(mode %in% 0:3)) stop("'mode' must be 0, 1, 2 or 3")
  miss <- setdiff(.lppls_pars$core, c(names(par), vars$x, vars$y))   ## fixed + scanned cover the core
  if (length(miss)) {
    stop(
      "core parameter(s) missing from 'par' + 'vars': ",
      paste(miss, collapse = ", "))
  }
  absent <- setdiff(.lppls_mode_higher(mode), c(names(par), vars$x, vars$y))
  if (length(absent))
    warning("mode ", mode, " uses ", paste(.lppls_mode_higher(mode), collapse = ", "),
            "; '", paste(absent, collapse = "', '"),
            "' not in 'par'/'vars' \u2014 using eval_lppls() default(s)")

  if (fb) message("Generating contour plot data...")

  contour_plot_out <- NULL
  surface_plot_out <- NULL

  ## Axis labels follow each parameter's meaning under `mode`: in the higher-order
  ## Landau modes the C2 slot is the log-periodic phase, so it is shown as "phi".
  xlab <- .lppls_display_label(vars$x, mode)
  ylab <- .lppls_display_label(vars$y, mode)

  x_contour <- seq(lower[1], upper[1], length.out = 101)

  x_contour <- seq(lower[1], upper[1], length.out = 101)
  y_contour <- seq(lower[2], upper[2], length.out = 101)

  sse_contour <- function(x, y) {
    xy_list <- list(x, y)
    names(xy_list) <- c(vars$x, vars$y)
    .sse(c(par, xy_list), log_p, t, mode)     # was: SSE(par = c(par, xy_list), ...)
  }

  z_contour <- t(outer(x_contour, y_contour, Vectorize(sse_contour)))
  contour_data <- list(x = x_contour, y = y_contour, z = z_contour)

  if (fb) message("Contour data generated")

  if (cp) {
    if (fb) message("Generating contour plot...")

    contour_plot_out <- plotly::plot_ly(
      x = contour_data$x,
      y = contour_data$y,
      z = contour_data$z,
      type = "contour",
      contours = list(showlabels = FALSE, exponentformat = "e"),
      ncontours = 10,
      colorbar = list(exponentformat = "e", title = "SSE")
    ) %>%
    plotly::layout(
        title = sprintf("SSE w.r.t. %s and %s", xlab, ylab),
        xaxis = list(title = xlab),
        yaxis = list(title = ylab)
      )

    if (!is.na(sum(point))) {
      z_coord <- sse_contour(point[1], point[2])
      contour_plot_out <- contour_plot_out %>%
        plotly::add_trace(
          x = point[1],
          y = point[2],
          type = "scatter",
          mode = "markers",
          showlegend = FALSE
        ) %>%
        plotly::add_annotations(
          x = point[1],
          y = point[2],
          xanchor = "middle",
          yanchor = "bottom",
          text = sprintf("(%.4f, %.4f, %.4e)", point[1], point[2], z_coord),
          font = list(color = "white", size = 10),
          showarrow = FALSE
        )
    }

    if (fb) message("Contour plot generated")
  }

  if (sp) {
    if (fb) message("Generating surface plot...")

    surface_plot_out <- plotly::plot_ly(
      x = contour_data$x,
      y = contour_data$y,
      z = contour_data$z
    ) %>%
      plotly::add_surface(
        contours = list(
          z = list(
            show = TRUE,
            usecolormap = TRUE,
            highlightcolor = "#ff0000",
            project = list(z = TRUE)
          )
        ),
        showscale = FALSE          ## z-axis already is SSE; drop the redundant colorbar
      ) %>%
      plotly::layout(
        title = sprintf("SSE w.r.t. %s and %s", xlab, ylab),
        scene = list(
          camera = list(eye = list(x = 1.1, y = -1.3, z = 0.1)),
          xaxis = list(title = xlab),
          yaxis = list(title = ylab),
          zaxis = list(title = "SSE")
        )
      )

    if (fb) message("Surface plot generated")
  }

  list(
    contour_plot = contour_plot_out,
    surface_plot = surface_plot_out,
    contour_data = contour_data
  )
}


#' Plot tc Estimates and Likelihood Intervals Across Rolling T1
#'
#' Plots the point estimate of the critical time \eqn{t_c} (profile likelihood,
#' black; modified profile likelihood, cyan) with a likelihood-interval band for
#' each cutoff (by default 5\%, 10\% and 50\%), as a function of the rolling
#' calibration start time
#' \eqn{T_1}. This is a general window-diagnostic -- useful with or without
#' Lagrange regularization -- and reproduces the `*_output_T1` figure from the
#' thesis. Requires the rolling calibration to have been run in `mode = "MPL"`.
#'
#' @param x An object of class `"lppls_rolling"` from [lppls_rolling()].
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_rolling()], [rolling_param_plot()]
#'
#' @export
rolling_tc_plot <- function(x) {
  if (!inherits(x, "lppls_rolling")) {
    stop("'x' must be an object of class 'lppls_rolling' (see lppls_rolling())")
  }
  if (all(is.na(x$table$tc_hat_mpl))) {
    stop("No likelihood intervals available: run lppls_rolling() with mode = 'MPL'.")
  }

  tab <- x$table

  ## The interval columns are named for their cutoff, LI5l/LI5u for 0.05; find
  ## the pairs by name, as noise_tc_plot() does, and take them in cutoff order.
  lo_cols <- grep("^LI[0-9.]+l$", names(tab), value = TRUE)
  pairs   <- lapply(lo_cols, function(lo) {
    hi <- sub("l$", "u", lo)
    if (!hi %in% names(tab)) return(NULL)
    list(lo = lo, hi = hi, pct = as.numeric(sub("^LI([0-9.]+)l$", "\\1", lo)))
  })
  pairs <- Filter(Negate(is.null), pairs)
  pairs <- pairs[order(vapply(pairs, function(p) p$pct, numeric(1)))]
  bands <- vapply(pairs, function(p) sprintf("LI_%02g", p$pct), character(1))

  ## Colours go to the bands in order, so the default three keep theirs.
  cols <- c(stats::setNames(rep_len(c("green", "orange", "red", "purple", "brown",
                                      "darkcyan"), length(bands)), bands),
            "MPLE" = "cyan", "PLE" = "black")

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = T1, y = tc))
  for (i in seq_along(pairs)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = !!as.name(pairs[[i]]$lo), ymax = !!as.name(pairs[[i]]$hi),
                   color = !!bands[i]),
      fill = "grey50",
      alpha = 0.2,
      key_glyph = "path"
    )
  }

  p +
    ggplot2::geom_line(ggplot2::aes(color = "PLE"), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = tc_hat_mpl, color = "MPLE"), linewidth = 0.25) +
    ggplot2::scale_color_manual(name = NULL, values = cols, breaks = names(cols)) +
    ggplot2::xlab(quote(T[1])) +
    ggplot2::ylab(quote(hat(t)[c])) +
    ggplot2::theme_minimal() +
    overlay_legend_theme(
      inside = c(0.99, 0.99),
      justification = c(1, 1)
    )
}


#' Plot LPPLS Parameter Estimates Across Rolling T1
#'
#' Plots the estimated LPPLS parameters \eqn{m, \omega, A, B, C_1, C_2} as a
#' function of the rolling calibration start time \eqn{T_1}, one facet per
#' parameter. A general window-diagnostic that reproduces the `*_params_plot`
#' figure from the thesis.
#'
#' @param x An object of class `"lppls_rolling"` from [lppls_rolling()].
#'
#' @return A `ggplot2` object.
#'
#' @seealso [lppls_rolling()], [rolling_tc_plot()]
#'
#' @importFrom tidyr gather
#' @export
rolling_param_plot <- function(x) {
  if (!inherits(x, "lppls_rolling")) {
    stop("'x' must be an object of class 'lppls_rolling' (see lppls_rolling())")
  }

  d <- tidyr::gather(
    x$table[, c("m", "omega", "A", "B", "C1", "C2", "T1")],
    key = "param", value = "estimate", -T1
  )

  ggplot2::ggplot(d, ggplot2::aes(x = T1, y = estimate)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~param, scales = "free_y", ncol = 1) +
    ggplot2::xlab(quote(T[1])) # +
    # ggplot2::theme_minimal()
}


#' Residual QQ Plot Diagnostic
#'
#' Quantile-quantile plot of calibration residuals against a normal or a fitted
#' Student-t reference distribution, with reference lines at chosen quantiles. A
#' general goodness-of-fit diagnostic for the residuals of an LPPLS calibration.
#'
#' @param residuals Numeric vector of residuals (length >= 2).
#' @param distribution Reference distribution for the theoretical quantiles:
#'   `"normal"` (the default) or `"t"`, a Student-t distribution whose degrees of
#'   freedom are fitted to `residuals` with [MASS::fitdistr()].
#' @param ref Numeric probabilities at which to draw the green reference lines
#'   (default `c(0.1, 0.9)`): horizontal lines at the corresponding sample
#'   quantiles of `residuals`, vertical lines at the theoretical quantiles.
#'
#' @return A `ggplot2` object.
#'
#' @seealso [chi_sq_np_plot()], [rolling_tc_plot()]
#'
#' @importFrom ggplot2 ggplot aes stat_qq stat_qq_line geom_hline geom_vline labs theme_minimal
#' @importFrom stats quantile qnorm qt
#' @export
#'
#' @examples
#' set.seed(1)
#' residual_qq_plot(rnorm(200))
#' if (requireNamespace("MASS", quietly = TRUE))
#'   residual_qq_plot(rt(200, df = 4), distribution = "t")
residual_qq_plot <- function(residuals, distribution = c("normal", "t"),
                             ref = c(0.1, 0.9)) {
  distribution <- match.arg(distribution)
  if (!is.numeric(residuals) || length(residuals) < 2L) {
    stop("'residuals' must be a numeric vector of length >= 2.")
  }

  d     <- data.frame(y = residuals)
  y_ref <- stats::quantile(residuals, ref, names = FALSE)

  if (distribution == "normal") {
    x_ref     <- stats::qnorm(ref)
    qq_layers <- list(
      ggplot2::stat_qq(size = 0.4),
      ggplot2::stat_qq_line(color = "red")
    )
  } else {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is needed for distribution = \"t\"; please install it.")
    }
    ## fitdistr()'s optimiser tries negative scales on its way to the estimate,
    ## and dt() then warns "NaNs produced". Those warnings say nothing about the
    ## fit, so muffle them, and only them.
    t_fit <- withCallingHandlers(
      MASS::fitdistr(residuals, "t"),
      warning = function(w) {
        if (grepl("NaNs produced", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    nu        <- t_fit$estimate[["df"]]
    dp        <- list(df = nu)
    x_ref     <- stats::qt(ref, df = nu)
    qq_layers <- list(
      ggplot2::stat_qq(distribution = stats::qt, dparams = dp, size = 0.4),
      ggplot2::stat_qq_line(distribution = stats::qt, dparams = dp, color = "red")
    )
  }

  ggplot2::ggplot(d, ggplot2::aes(sample = y)) +
    qq_layers +
    ggplot2::geom_hline(yintercept = y_ref, color = "green3") +
    ggplot2::geom_vline(xintercept = x_ref, color = "green3") +
    ggplot2::labs(x = "theoretical quantiles", y = "sample quantiles") +
    ggplot2::theme_minimal()
}
