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


#' Create MPL Plot
#'
#' @keywords internal
#' @noRd
create_mpl_plot <- function(mpl_output, fit, n, fh) {
  tc_seq <- (n + 1):(n + fh)

  df_plot <- data.frame(
    tc = tc_seq,
    log_mpl = mpl_output$LL
  )

  ## Collect the vertical markers into one data frame so colour and linetype
  ## map to a single legend: the PLE / MPLE point estimates and the LI bounds
  ## (each interval contributes two lines but one legend entry).
  vmarks <- data.frame(
    xintercept = c(fit[[1]]$tc, mpl_output$tc_hat_mpl),
    series = c("PLE", "MPLE")
  )
  li_names <- c("LI_05", "LI_10", "LI_50")
  for (i in seq_along(mpl_output$LI)) {
    li <- mpl_output$LI[[i]]
    if (!any(is.na(li))) {
      vmarks <- rbind(
        vmarks,
        data.frame(xintercept = c(li[1], li[2]), series = li_names[i])
      )
    }
  }

  all_lvls <- c("PLE", "MPLE", "LI_05", "LI_10", "LI_50")
  present <- all_lvls[all_lvls %in% vmarks$series]
  vmarks$series <- factor(vmarks$series, levels = present)

  cols <- c(PLE = "blue", MPLE = "blue", LI_05 = "green3",
            LI_10 = "orange", LI_50 = "red")
  ltys <- c(PLE = "solid", MPLE = "dashed", LI_05 = "dashed",
            LI_10 = "dotted", LI_50 = "dotdash")

  ggplot2::ggplot(df_plot, ggplot2::aes(x = tc, y = log_mpl)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      data = vmarks,
      mapping = ggplot2::aes(xintercept = xintercept,
                             color = series, linetype = series),
      key_glyph = "path"
    ) +
    ggplot2::scale_color_manual(name = NULL, values = cols, breaks = present) +
    ggplot2::scale_linetype_manual(name = NULL, values = ltys, breaks = present) +
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
#' Visualizes the L-BFGS-B optimization path over an SSE contour surface.
#' Three selectors: 1 = (m, omega), 2 = (B, m), 3 = (B, omega).
#'
#' @param selector Integer (1, 2, or 3) selecting which parameter pair to plot.
#' @param opt_trace Matrix with 2 or 3 rows (m, omega, [B]) and one column per step.
#' @param lower Numeric vector of length 4, parameter lower bounds.
#' @param upper Numeric vector of length 4, parameter upper bounds.
#' @param SSE2_func SSE2 function (fixed tc) from fit_lppls.
#' @param beta_calc Beta calculator function from fit_lppls.
#' @param tc_val Numeric, the tc value for the trace.
#' @param log_p Numeric vector of log-prices.
#' @param t Numeric vector of time indices.
#' @param tp_id Integer, which trace step to use for contour in selectors 2/3.
#' @param fb Logical, whether to print feedback.
#'
#' @return A ggplot2 object.
#'
#' @keywords internal
#' @noRd
create_trace_plot <- function(
    selector,
    opt_trace,
    lower,
    upper,
    SSE2_func,
    beta_calc,
    tc_val,
    log_p,
    t,
    tp_id,
    fb = FALSE
  ) {
  if (fb) message(sprintf("Generating trace plot %d...", selector))

  lattice_dim <- 100

  if (selector == 1) {
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

  } else if (selector == 2) {
    ## B vs m lattice
    lattice <- as.matrix(expand.grid(
      m = seq(
        min(lower[1], min(opt_trace[1, ])),
        max(upper[1], max(opt_trace[1, ])),
        length.out = lattice_dim
      ),
      B = seq(
        min(lower[3], min(opt_trace[3, ])),
        max(upper[3], max(opt_trace[3, ])),
        length.out = lattice_dim
      )
    ))

    beta_tp <- beta_calc(log_p, t, tc_val, opt_trace[1, 1], opt_trace[2, 1])

    if (tp_id < 1 || tp_id > ncol(opt_trace)) {
      tp_id <- 1
      warning("tp_id for (B,m) trace plot is <1 or too high. Replaced with tp_id=1\n")
    }

    SSE_tp <- function(par, a, c1, c2, tc, omega, log_p, t) {
      SSE(
        par = list(
          A = a,
          B = par[2],
          C1 = c1,
          C2 = c2,
          tc = tc,
          m = par[1],
          omega = omega
        ),
        log_p = log_p,
        t = t
      )
    }

    lattice_vals <- apply(
      lattice,
      1,
      SSE_tp,
      a = unname(beta_tp["A"]),
      c1 = unname(beta_tp["C1"]),
      c2 = unname(beta_tp["C2"]),
      tc = tc_val,
      omega = opt_trace[2, tp_id],
      log_p = log_p,
      t = t
    )

  } else if (selector == 3) {
    ## B vs omega lattice
    lattice <- as.matrix(expand.grid(
      omega = seq(
        min(lower[2], min(opt_trace[2, ])),
        max(upper[2], max(opt_trace[2, ])),
        length.out = lattice_dim
      ),
      B = seq(
        min(lower[3], min(opt_trace[3, ])),
        max(upper[3], max(opt_trace[3, ])),
        length.out = lattice_dim
      )
    ))

    beta_tp <- beta_calc(log_p, t, tc_val, opt_trace[1, 1], opt_trace[2, 1])

    if (tp_id > ncol(opt_trace)) {
      tp_id <- 1
      warning("tp_id for (B,omega) trace plot too high. Replaced with tp_id=1\n")
    }

    SSE_tp <- function(par, a, c1, c2, tc, m, log_p, t) {
      SSE(
        par = list(
          A = a,
          B = par[2],
          C1 = c1,
          C2 = c2,
          tc = tc,
          m = m,
          omega = par[1]
        ),
        log_p = log_p,
        t = t
      )
    }

    lattice_vals <- apply(
      lattice,
      1,
      SSE_tp,
      a = unname(beta_tp["A"]),
      c1 = unname(beta_tp["C1"]),
      c2 = unname(beta_tp["C2"]),
      tc = tc_val,
      m = opt_trace[1, tp_id],
      log_p = log_p,
      t = t
    )
  }

  lattice <- as.data.frame(lattice)
  lattice$vals <- lattice_vals

  if (selector == 1) {
    opt_df <- data.frame(m = opt_trace[1, ], omega = opt_trace[2, ])
    opt_df$step <- seq_len(nrow(opt_df))

    p <- ggplot2::ggplot(lattice, ggplot2::aes(x = m, y = omega)) +
      ggplot2::geom_raster(ggplot2::aes(fill = vals)) +
      ggplot2::geom_contour(ggplot2::aes(z = vals), col = "white") +
      ggplot2::geom_path(data = opt_df, col = "orange") +
      ggplot2::geom_point(data = opt_df, col = "orange") +
      ggplot2::geom_point(
        data = head(opt_df, 1),
        shape = 21,
        fill = "red",
        color = "orange",
        size = 3
      ) +
      ggplot2::geom_point(
        data = tail(opt_df, 1),
        shape = 23,
        fill = "green",
        color = "orange",
        size = 3
      ) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_cartesian(expand = FALSE)

  } else if (selector == 2) {
    opt_df <- data.frame(m = opt_trace[1, ], B = opt_trace[3, ])
    opt_df$step <- seq_len(nrow(opt_df))

    p <- ggplot2::ggplot(lattice, ggplot2::aes(x = m, y = B)) +
      ggplot2::geom_raster(ggplot2::aes(fill = vals)) +
      ggplot2::geom_contour(ggplot2::aes(z = vals), col = "white") +
      ggplot2::geom_path(data = opt_df, col = "orange") +
      ggplot2::geom_point(data = opt_df, col = "orange") +
      ggplot2::geom_point(
        data = head(opt_df, 1),
        shape = 21,
        fill = "red",
        color = "orange",
        size = 3
      ) +
      ggplot2::geom_point(
        data = tail(opt_df, 1),
        shape = 23,
        fill = "green",
        color = "orange",
        size = 3
      ) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_cartesian(expand = FALSE)

  } else if (selector == 3) {
    opt_df <- data.frame(omega = opt_trace[2, ], B = opt_trace[3, ])
    opt_df$step <- seq_len(nrow(opt_df))

    p <- ggplot2::ggplot(lattice, ggplot2::aes(x = omega, y = B)) +
      ggplot2::geom_raster(ggplot2::aes(fill = vals)) +
      ggplot2::geom_contour(ggplot2::aes(z = vals), col = "white") +
      ggplot2::geom_path(data = opt_df, col = "orange") +
      ggplot2::geom_point(data = opt_df, col = "orange") +
      ggplot2::geom_point(
        data = head(opt_df, 1),
        shape = 21,
        fill = "red",
        color = "orange",
        size = 3
      ) +
      ggplot2::geom_point(
        data = tail(opt_df, 1),
        shape = 23,
        fill = "green",
        color = "orange",
        size = 3
      ) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_cartesian(expand = FALSE)
  }

  if (fb) message("Trace plot generated")
  p
}


#' Plot a Fitted LPPLS Model
#'
#' Draws the observed log-price series together with the fitted LPPLS curve, a
#' vertical marker at the estimated critical time \eqn{\hat{t}_c}{tc-hat}, and---when the
#' calibration end is supplied via `T2`---a marker at \eqn{T_2}. Colour and line
#' type are merged into a single legend (`sim data` / `LPPLS fit` / `T2` /
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
#'
#' @return A `ggplot2` object.
#'
#' @seealso [fit_lppls()], [eval_lppls()]
#'
#' @examples
#' \dontrun{
#' fit <- fit_lppls(log_price = log_p, mode = "F2")
#' # public convenience call: only the tc_hat marker
#' fit_plot(fit$fit[[1]], seq_along(log_p), log_p)
#' # with the calibration end marked, as fit_lppls() draws it internally
#' fit_plot(fit$fit[[1]], seq_along(log_p), log_p, T2 = length(log_p) - 200)
#' }
#'
#' @export
fit_plot <- function(fit, time_ID, log_price, mode = 0, T2 = NULL) {
  plot_data <- data.frame(ID = time_ID, log_p = log_price)

  # tc_hat always; T2 only when the caller supplies it
  marks <- data.frame(xintercept = fit$tc, series = "tc_hat")
  if (!is.null(T2)) marks <- rbind(data.frame(xintercept = T2, series = "T2"), marks)

  lvls <- c("sim data", "LPPLS fit", if (!is.null(T2)) "T2", "tc_hat")
  cols <- c("sim data"="royalblue1", "LPPLS fit"="red", "T2"="green", "tc_hat"="red")
  ltys <- c("sim data"="solid", "LPPLS fit"="solid", "T2"="dashed", "tc_hat"="dashed")

  fit_t    <- seq.int(min(time_ID), floor(fit$tc) - 1L)
  fit_line <- data.frame(
    ID    = fit_t,
    log_p = eval_lppls(fit_t, fit$A, fit$B, fit$C1, fit$C2,
                       fit$tc, fit$m, fit$omega, mode = mode)
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(ID, log_p)) +
    ggplot2::geom_line(ggplot2::aes(color = "sim data", linetype = "sim data")) +
    ggplot2::geom_line(
      data = fit_line,
      mapping = ggplot2::aes(color = "LPPLS fit", linetype = "LPPLS fit")
    ) +
    ggplot2::geom_vline(
      data = marks,
      mapping = ggplot2::aes(xintercept = xintercept, color = series, linetype = series),
      inherit.aes = FALSE, key_glyph = "path") +
    ggplot2::scale_color_manual(name = NULL, values = cols, breaks = lvls) +
    ggplot2::scale_linetype_manual(name = NULL, values = ltys, breaks = lvls) +
    ggplot2::labs(x = "Time Index", y = "Log Price") +
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
#' @return List with elements: contour_plot, surface_plot (both may be NULL).
#'
#' @examples
#' \dontrun{
#' # Create contour plot for m vs omega
#' par <- list(tc = 250, A = 4, B = -0.015, C1 = 0.002, C2 = 0.001)
#' vars <- list(x = "m", y = "omega")
#' result <- contour_plot_sse(log_p, t, par, vars,
#'                            lower = c(0.1, 6), upper = c(0.9, 13))
#' print(result$contour_plot)
#' }
#'
#' @export
contour_plot_sse <- function(
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
#' black; modified profile likelihood, cyan) with the 5\%, 10\% and 50\%
#' likelihood-interval bands, as a function of the rolling calibration start time
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
  cols <- c("LI_05" = "green", "LI_10" = "orange", "LI_50" = "red",
            "MPLE" = "cyan", "PLE" = "black")

  ggplot2::ggplot(tab, ggplot2::aes(x = T1, y = tc)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = LI5l, ymax = LI5u, color = "LI_05"),
      fill = "grey50",
      alpha = 0.2,
      key_glyph = "path"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = LI10l, ymax = LI10u, color = "LI_10"),
      fill = "grey50",
      alpha = 0.2,
      key_glyph = "path"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = LI50l, ymax = LI50u, color = "LI_50"),
      fill = "grey50",
      alpha = 0.2,
      key_glyph = "path"
    ) +
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
    nu        <- MASS::fitdistr(residuals, "t")$estimate[["df"]]
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
