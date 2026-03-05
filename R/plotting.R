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
#' @param fb Logical, whether to print feedback.
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
compute_mpl <- function(fit, log_p, t, fh, cutoff, fb = FALSE) {
  n <- length(t)
  log_cutoff <- log(cutoff)

  # Initialize results
  R_tbl <- data.frame(R = numeric(fh), LL = numeric(fh))

  # Get profile likelihood estimate
  par_hat <- c(
    fit[[1]]$tc,
    fit[[1]]$m,
    fit[[1]]$omega,
    fit[[1]]$A,
    fit[[1]]$B,
    fit[[1]]$C1,
    fit[[1]]$C2
  )

  # Compute log-likelihood for each tc value
  for (i in seq_len(fh)) {
    if (fb && i %% 10 == 1) {
      message(sprintf("Computing LL for tc = %d to %d", n + i, min(n + i + 9, n + fh)))
    }

    # Get fit for tc = n + i
    fit_i <- fit[[2]] %>% dplyr::filter(tc == n + i)

    if (nrow(fit_i) == 0) {
      # Use unfiltered fits if filtered is empty
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
        t = t
      )
    } else {
      R_tbl$LL[i] <- NA
    }
  }

  # Calculate maximum and relative likelihood
  MLL <- max(R_tbl$LL, na.rm = TRUE)
  R_tbl$R <- R_tbl$LL - MLL

  tc_hat_mpl <- n + which.max(R_tbl$LL)

  # Compute likelihood intervals
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
compute_mpl_loglik <- function(tc, tc_hat, Psi_hat_tc, Psi_hat, log_p, t) {
  n <- length(t)
  p <- 6

  # Get SSE for this tc
  m_hat <- Psi_hat[1]
  omega_hat <- Psi_hat[2]

  # Create beta calculator
  beta_calc <- create_beta_calculator()
  beta_vals <- beta_calc(log_p, t, tc, Psi_hat_tc[1], Psi_hat_tc[2])

  # Compute fitted values and SSE
  fitted <- eval_lppls(t, beta_vals["A"], beta_vals["B"], beta_vals["C1"], beta_vals["C2"],
                       tc, Psi_hat_tc[1], Psi_hat_tc[2], mode = 0)
  s_tc <- sum((log_p - fitted)^2, na.rm = TRUE) / n

  # Compute matrices
  X_hat_tc <- compute_X_matrix(Psi_hat_tc, tc, t)
  X_hat <- compute_X_matrix(Psi_hat, tc_hat, t)
  H_hat_tc <- compute_H_matrix(Psi_hat_tc, tc, log_p, t)

  # Compute log-likelihood
  XtX_tc <- crossprod(X_hat_tc, X_hat_tc)
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

  for (i in seq_len(n)) {
    tau <- tc - t[i]
    log_tau <- log(tau)
    tau_m <- tau^m
    cos_term <- cos(omega * log_tau)
    sin_term <- sin(omega * log_tau)

    ## Gradient for single t_i
    ## Filimonov2017, equation (B16)

    ## d LPPLS(t_i; tc, psi) / d m
    X[i, 1] <- tau_m * log_tau * (b + c1 * cos_term + c2 * sin_term)
    ## d LPPLS(t_i; tc, psi) / d omega
    X[i, 2] <- tau_m * log_tau * (-c1 * sin_term + c2 * cos_term)
    ## d LPPLS(t_i; tc, psi) / d a
    X[i, 3] <- 1
    ## d LPPLS(t_i; tc, psi) / d b
    X[i, 4] <- tau_m
    ## d LPPLS(t_i; tc, psi) / d c1
    X[i, 5] <- tau_m * cos_term
    ## d LPPLS(t_i; tc, psi) / d c2
    X[i, 6] <- tau_m * sin_term
  }

  X
}


#' Compute H Matrix for MPL
#'
#' H matrix for modified profile likelihood.
#' Filimonov (2017), equation (37).
#'
#' @keywords internal
#' @noRd
compute_H_matrix <- function(Psi, tc, log_p, t) {
  n <- length(t)
  p <- 6 ## length(Psi) = 6
  H <- matrix(0, nrow = p, ncol = p)

  m <- Psi[1]
  omega <- Psi[2]
  a <- Psi[3]
  b <- Psi[4]
  c1 <- Psi[5]
  c2 <- Psi[6]

  for (i in seq_len(n)) {
    tau <- tc - t[i]
    log_tau <- log(tau)
    tau_m <- tau^m
    cos_term <- cos(omega * log_tau)
    sin_term <- sin(omega * log_tau)

    # LPPLS value
    lppls_val <- a + tau_m * (b + c1 * cos_term + c2 * sin_term)
    res <- log_p[i] - lppls_val

    # Second derivatives
    # (m, m)
    H[1, 1] <- H[1, 1] + res * tau_m * log_tau^2 * (b + c1 * cos_term + c2 * sin_term)
    # (m, omega)
    H[1, 2] <- H[1, 2] + res * tau_m * log_tau^2 * (-c1 * sin_term + c2 * cos_term)
    # (m, b)
    H[1, 4] <- H[1, 4] + res * tau_m * log_tau
    # (m, c1)
    H[1, 5] <- H[1, 5] + res * tau_m * log_tau * cos_term
    # (m, c2)
    H[1, 6] <- H[1, 6] + res * tau_m * log_tau * sin_term
    # (omega, omega)
    H[2, 2] <- H[2, 2] + res * (-1) * tau_m * log_tau^2 * (c1 * cos_term + c2 * sin_term)
    # (omega, c1)
    H[2, 5] <- H[2, 5] + res * (-1) * tau_m * log_tau * sin_term
    # (omega, c2)
    H[2, 6] <- H[2, 6] + res * tau_m * log_tau * cos_term
  }

  ## Symmetry

  ## (omega, m)
  H[2, 1] <- H[1, 2]

  ## (b, m)
  H[4, 1] <- H[1, 4]

  ## (c1, m)
  H[5, 1] <- H[1, 5]

  ## (c2, m)
  H[6, 1] <- H[1, 6]

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

  tc_hat <- fit[[1]]$tc
  tc_mpl <- mpl_output$tc_hat_mpl

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = tc, y = log_mpl)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = tc_hat, color = "blue") +
    ggplot2::geom_vline(xintercept = tc_mpl, color = "blue", linetype = "dashed") +
    ggplot2::labs(
      x = "Critical Time (tc)",
      y = "Log Modified Profile Likelihood",
      title = "Modified Profile Likelihood"
    ) +
    ggplot2::theme_minimal()

  # Add likelihood interval lines
  colors <- c("green3", "orange", "red")
  linetypes <- c("dashed", "dotted", "dotdash")

  for (i in seq_along(mpl_output$LI)) {
    li <- mpl_output$LI[[i]]
    if (!any(is.na(li))) {
      p <- p +
        ggplot2::geom_vline(xintercept = li[1], color = colors[i], linetype = linetypes[i]) +
        ggplot2::geom_vline(xintercept = li[2], color = colors[i], linetype = linetypes[i])
    }
  }

  p
}


#' Create Fit Plot
#'
#' @keywords internal
#' @noRd
create_fit_plot <- function(fit, time_id, log_price, n_model) {

  ## Full data set in [T1, T2 + hold_out]
  plot_data <- data.frame(ID = time_id, log_p = log_price)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = ID, y = log_p)) +
    ggplot2::geom_line(color = "royalblue1") +
    ggplot2::geom_function(
      fun = eval_lppls,
      n = floor(fit$tc) - 1,
      args = list(
        A = fit$A,
        B = fit$B,
        C1 = fit$C1,
        C2 = fit$C2,
        tc = fit$tc,
        m = fit$m,
        omega = fit$omega
      ),
      color = "red"
    ) +
    ggplot2::geom_vline(xintercept = n_model, color = "green", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = fit$tc, color = "red", linetype = "dashed") +
    ggplot2::labs(x = "Time Index", y = "Log Price") +
    ggplot2::theme_minimal()
}


#' Create Contour Data
#'
#' Generate SSE data for contour plot of SSE wrt m and omega.
#'
#' @keywords internal
#' @noRd
create_contour_data <- function(fit, log_p, t, lower, upper, beta_calculator, fb = FALSE) {
  if (fb) message("Generating contour plot data...")

  tc_val <- fit$tc

  x_contour <- seq(lower[1], upper[1], length.out = 101)  ## m
  y_contour <- seq(lower[2], upper[2], length.out = 101)  ## omega

  sse_func <- function(m_val, omega_val) {
    beta <- beta_calculator(log_p, t, tc_val, m_val, omega_val)
    fitted <- eval_lppls(t, beta["A"], beta["B"], beta["C1"], beta["C2"],
                         tc_val, m_val, omega_val, mode = 0)
    sum((log_p - fitted)^2, na.rm = TRUE)
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
    ggplot2::facet_wrap(~param, scales = "free_y", ncol = 1) +
    ggplot2::theme_minimal()
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
create_trace_plot <- function(selector, opt_trace, lower, upper,
                               SSE2_func, beta_calc,
                               tc_val, log_p, t, tp_id, fb = FALSE) {
  if (fb) message(sprintf("Generating trace plot %d...", selector))

  lattice_dim <- 100

  if (selector == 1) {
    # m vs omega lattice
    lattice <- as.matrix(expand.grid(
      m = seq(min(lower[1], min(opt_trace[1, ])),
              max(upper[1], max(opt_trace[1, ])), length.out = lattice_dim),
      omega = seq(min(lower[2], min(opt_trace[2, ])),
                  max(upper[2], max(opt_trace[2, ])), length.out = lattice_dim)
    ))
    lattice_vals <- apply(lattice, 1, SSE2_func,
                          tc = tc_val, log_p = log_p, t = t)

  } else if (selector == 2) {
    # B vs m lattice
    lattice <- as.matrix(expand.grid(
      m = seq(min(lower[1], min(opt_trace[1, ])),
              max(upper[1], max(opt_trace[1, ])), length.out = lattice_dim),
      B = seq(min(lower[3], min(opt_trace[3, ])),
              max(upper[3], max(opt_trace[3, ])), length.out = lattice_dim)
    ))

    beta_tp <- beta_calc(log_p, t, tc_val, opt_trace[1, 1], opt_trace[2, 1])

    if (tp_id < 1 || tp_id > ncol(opt_trace)) {
      tp_id <- 1
      warning("tp_id for (B,m) trace plot is <1 or too high. Replaced with tp_id=1")
    }

    SSE_tp <- function(par, a, c1, c2, tc, omega, log_p, t) {
      SSE(par = list(A = a, B = par[2], C1 = c1, C2 = c2,
                     tc = tc, m = par[1], omega = omega),
          log_p = log_p, t = t)
    }

    lattice_vals <- apply(lattice, 1, SSE_tp,
                          a = unname(beta_tp["A"]),
                          c1 = unname(beta_tp["C1"]),
                          c2 = unname(beta_tp["C2"]),
                          tc = tc_val,
                          omega = opt_trace[2, tp_id],
                          log_p = log_p, t = t)

  } else if (selector == 3) {
    # B vs omega lattice
    lattice <- as.matrix(expand.grid(
      omega = seq(min(lower[2], min(opt_trace[2, ])),
                  max(upper[2], max(opt_trace[2, ])), length.out = lattice_dim),
      B = seq(min(lower[3], min(opt_trace[3, ])),
              max(upper[3], max(opt_trace[3, ])), length.out = lattice_dim)
    ))

    beta_tp <- beta_calc(log_p, t, tc_val, opt_trace[1, 1], opt_trace[2, 1])

    if (tp_id > ncol(opt_trace)) {
      tp_id <- 1
      warning("tp_id for (B,omega) trace plot too high. Replaced with tp_id=1")
    }

    SSE_tp <- function(par, a, c1, c2, tc, m, log_p, t) {
      SSE(par = list(A = a, B = par[2], C1 = c1, C2 = c2,
                     tc = tc, m = m, omega = par[1]),
          log_p = log_p, t = t)
    }

    lattice_vals <- apply(lattice, 1, SSE_tp,
                          a = unname(beta_tp["A"]),
                          c1 = unname(beta_tp["C1"]),
                          c2 = unname(beta_tp["C2"]),
                          tc = tc_val,
                          m = opt_trace[1, tp_id],
                          log_p = log_p, t = t)
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
      ggplot2::geom_point(data = head(opt_df, 1), shape = 21,
                          fill = "red", color = "orange", size = 3) +
      ggplot2::geom_point(data = tail(opt_df, 1), shape = 23,
                          fill = "green", color = "orange", size = 3) +
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
      ggplot2::geom_point(data = head(opt_df, 1), shape = 21,
                          fill = "red", color = "orange", size = 3) +
      ggplot2::geom_point(data = tail(opt_df, 1), shape = 23,
                          fill = "green", color = "orange", size = 3) +
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
      ggplot2::geom_point(data = head(opt_df, 1), shape = 21,
                          fill = "red", color = "orange", size = 3) +
      ggplot2::geom_point(data = tail(opt_df, 1), shape = 23,
                          fill = "green", color = "orange", size = 3) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_cartesian(expand = FALSE)
  }

  if (fb) message("Trace plot generated")
  p
}


#' Create Fit Plot (Exported Version)
#'
#' Generate a plot showing the fitted LPPLS model against observed data.
#'
#' @param fit List of parameter values (A, B, C1, C2, tc, m, omega).
#' @param time_ID Numeric vector of time indices.
#' @param log_price Numeric vector of log-prices.
#' @param mode Integer, LPPLS model mode (default 0).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # After fitting a model
#' result <- fit_lppls(log_price = log_p, mode = "F2")
#' p <- fit_plot(result$fit[[1]], t, log_p)
#' print(p)
#' }
#'
#' @export
fit_plot <- function(fit, time_ID, log_price, mode = 0) {
  plot_data <- data.frame(ID = time_ID, log_p = log_price)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = ID, y = log_p)) +
    ggplot2::geom_line(linewidth = 0.5, color = "royalblue1") +
    ggplot2::geom_function(
      fun = eval_lppls,
      args = list(
        A = fit$A,
        B = fit$B,
        C1 = fit$C1,
        C2 = fit$C2,
        tc = fit$tc,
        m = fit$m,
        omega = fit$omega,
        mode = mode
      ),
      color = "red"
    ) +
    ggplot2::labs(x = "Time Index", y = "Log Price") +
    ggplot2::theme_minimal()
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
contour_plot_sse <- function(log_p, t, par, vars,
                              lower = c(0.0, -0.001), upper = c(1.0, 0.001),
                              cp = TRUE, sp = FALSE, fb = FALSE,
                              point = c(NA, NA), mode = 0) {

  if (length(log_p) != length(t)) {
    stop("Length of log_p and t must be the same")
  }

  if (fb) message("Generating contour plot data...")

  contour_plot_out <- NULL
  surface_plot_out <- NULL

  x_contour <- seq(lower[1], upper[1], length.out = 101)
  y_contour <- seq(lower[2], upper[2], length.out = 101)

  sse_contour <- function(x, y) {
    xy_list <- list(x, y)
    names(xy_list) <- c(vars$x, vars$y)
    SSE(par = c(par, xy_list), log_p = log_p, t = t, mode = mode)
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
        title = sprintf("SSE w.r.t. %s and %s", vars$x, vars$y),
        xaxis = list(title = vars$x),
        yaxis = list(title = vars$y)
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
        colorbar = list(exponentformat = "e", title = "SSE")
      ) %>%
      plotly::layout(
        title = sprintf("SSE w.r.t. %s and %s", vars$x, vars$y),
        scene = list(
          camera = list(eye = list(x = 1.1, y = -1.3, z = 0.1)),
          xaxis = list(title = vars$x),
          yaxis = list(title = vars$y),
          zaxis = list(title = "SSE")
        )
      )

    if (fb) message("Surface plot generated")
  }

  list(contour_plot = contour_plot_out, surface_plot = surface_plot_out)
}
