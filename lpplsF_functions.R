##############################################
## Implementation of the LPPLS model for analysing financial bubbles.
## By Martin Hoshi Vognsen
##
## This implementation is based on the following articles:
##
## 1) A stable and robust calibration scheme of the log-periodic power law model
## Filimonov, V., Sornette, D.
## Physica A: Statistical Mechanics and its Applications, 2013
##
## 2) Modified profile likelihood inference and interval forecast of the burst of financial bubbles
## Filimonov, V., Demos, G., Sornette, D.
## Quantitative Finance, 2017
##############################################

## Setup ====

library(ggplot2)
library(tidyverse)
library(symengine) ## Solve symbolic equations
library(plotly) ## Contour plot
library(rlist) ## Sort list of fits
#library(Rcpp)
#library(optimParallel)

## Functions ====

## Calculate complete model based on input time ID and price vectors.

## Inputs:
## 1) time_ID, time index
##         vector, T1:T2 where T1 = 1
## 2) log_price, price
##         vector, must be same length as t
## 3) fh, length of forecast horizon. fh = s means forecast period is [T2+1, T2+s], s>1.
##      Unit is indexes.
##      Default is 90 time units.
## 4) hold_out, number of time units for holdout data. 
##         That is, the number of data points after T2.
## 5) lower, lower limits of parameter filter
##      {m, omega, B, D}
##      Filimonov 2017: B < 0, D >= 0.8, 0.1 <= m <= 0.9, 6 <= omega <= 13 
##      Damping, D = (m * abs(B) / (omega * sqrt(C1^2 + C2^2))
##      Defaults: lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6)
##      NOTE: For trace_plot mode 2 and 3 (bm and bo), lower for B should be set to
##            the desired minimum value of the B-axis in the contour plot.
##            Example:
##            lower = c(0.1, 6, -0.03, 0.8)
##            If the optimization algorithm goes below the lower limits at any point
##            the lower limit for the B-axis will be extended automatically.
## 6) upper, upper limits of parameter filter
## 7) tc_init, initial tc value for optim().
##         tc_init should be in [T2+1, T2+s]
## 8) m_init, initial m value for optim()
##         m_init should be in [lower[1], upper[1]]
##         The algorithm will search for initial values of m and omega that give B<0.
##         If such values are found, the given init values will be overwritten.
## 9) o_init, initial omega value for optim()
##        Omega_init should be in [lower[2], upper[2]]
##        Default is 13, which increases the probability of B<0
##        (of course not for random initial values).
##        The algorithm will search for initial values of m and omega that give B<0.
##        If such values are found, the given init values will be overwritten.
## 10) num_searches, number of times to repeat the optimization wrt tc, m and omega.
##          First iteration uses given initial values.
##          Any subsequent iterations use random initial values.
## 11) mode (see also under output below), 
##          mode = "F1", optimize wrt tc, m and omega simultaneously
##          mode = "F2", optimize wrt m and omega simultaneously, then wrt tc
##          mode = "MPL", Modified Profile Likelihood with likelihood intervals and plot.
## 12) mpl_cutoff, cutoff level (c) for likelihood interval. See Filimonov2017 equation (39).
##               vector with 3 elements, for 3 levels.
##               Default is c = c(0.05, 0.1, 0.5).
##               c = 0.05 indicates that values of tc outside the likelihood interval
##               have a probability of 0.05.
##               Note that the intervals for level c = 0.05 can be extremely wide.
## 12) mpl_plot, boolean for Modified Profile Likelihood plot with likelihood intervals.
## 13) cp, boolean for contour plot.
##          tc value is fixed to the estimate with the lowest objective function value.
## 14) sp, boolean for contour plot.
##          tc value is fixed to the estimate with the lowest objective function value.
## 15) tp, boolean for vector trace plot selector.
##          Any combination of three plots are selected with a vector:
##          c(selector 1 on/off, selector 2 on/off, selector 3 on/off).
##          selector 1: m, omega
##          selector 2: B, m
##          selector 3: B, omega
##          Example
##          c(1, 0, 1) selects trace plot 1 and 3.
##          Traces convergence of parameters during optimization of F2 objective function.
##          Only works when mode = "F2".
##          The plotted trace path is the one for the value of tc that produced the
##          smallest value of SSE2.
##          NOTE: Contour plots in trace_plot 2 and 3 only apply to the first step in the trace.
##                In a trace plot wrt B and m, omega is not fixed. (And vice versa for m and omega.)
## 16) tp_id, integer
##          Indicate which step of the trace to use for contour plot.
##          If the indicated id is too high, 1 will be chosen.
## 17) fp, boolean for fit plot.
## 18) mp, boolean for matrix plot.
##              x-axis is tc, y-axis is sse.
##              3 columns of plots for smallest and biggest value of m, and value between the two
##              3 columns of plots for smallest and biggest value of omega, and value between the two
## 19) factr, factr argument for optim control list when using L-BFGS-B. See optim documentation.
## 20) fb, boolean for printing feedback during execution.



## Outputs:
## 1)
## mode = F1:
##    fit, list of fits with random initial values.
##         fit[[1]]:
##              list of best fit coefficients optimized wrt m and omega.
##         fit[[2]]:
##              Tibble of fits with random starting parameters.
##              Sorted by objective function value (best at top).
##              Output format:
##              [
##                   list(tc, m, omega, A, B, C1, C2, D, value_min, ID),
##                   list(tc, m, omega, A, B, C1, C2, D, value, ID),
##                   ...,
##                   list(tc, m, omega, A, B, C1, C2, D, value_max, ID)
##              ]
##         fit[[3]]:
##              If no fits passed the filter: Only unfiltered fits returned.
##    ID is an identifier indicating the order pre-sorting.
##    These are the indexes reported in "out of range" warnings.
## mode = F2:
##    fit, list with two or three elements:
##         fit[[1]]:
##              list of best fit coefficients optimized wrt m and omega.
##         fit[[2]]:
##              If at least one fit passed the filter:
##                   Tibble of best fits for objective function F2 for each value of tc.
##                   Columns are: {ID, value, tc, m, omega, A, B, C1, C2, damp}
##                   Sorted by objective function value (best at top).
##                   For each value of tc the best fit is picked from fit2_tmp.
##              If no fits passed the filter: 
##                   Only unfiltered fits returned.
##         fit[[3]]:
##              If at least one fit passed the filter:
##                   Unfiltered fits returned.
##              If no fits passed the filter: 
##                   Nothing is returned in fit[[3]]
##         fit[[4]]:
##              List of all tmp fits
##              For each tc value, there is a tibble with all fits using random starting points.
## 2) mpl_output
##         mpl = list(LI, R, LL, MLL)
##              LI, vector of 3 elements, one for each value of mpl_cutoff.
##              R, number, relative likelihood.
##              LL, vector of MPL log-likelihoods for each tc value.
##              MLL, maximum of LL.
##              tc_hat_mpl, tc for which LL takes it's maximum.
## 3) mpl_plot, Modified Profile Likelihood plot with likelihood intervals.
## 4) contour_data, list of x, y and z data for contour plot of SSE wrt m and omega.
## 5) fit_plot, object containing fit plot
## 6) contour_plot, object containing contour plot
## 7) surface_plot, object containing surface plot
## 8) trace_plot_mo, object containing trace plot for m and omega
## 9) trace_plot_bm, object containing trace plot for B and m
## 10) trace_plot_bo, object containing trace plot for B and omega
## 11) matrix_plot, object containing matrix plot
## 12) out_of_range_tracker
##         F1 mode: list of random iteration IDs what resulted in out-of-range parameters.
##         F2 mode: list of random tc value IDs and iteration IDs what resulted in 
##                  out-of-range parameters.
##                  tc value ID start with 1 for T2+1. That is, the first time step
##                  after the modelling period, or the first time step of the
##                  forecasting horizon.

lpplsF <- function(time_ID, log_price, fh = 30, hold_out = 15, lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6), tc_init = 1000, m_init = 0.5, o_init = 13, num_searches = 20, mode = "F2", mpl_plot = 0, mpl_cutoff = c(0.05, 0.1, 0.5), fp = 0, cp = 0, sp = 0, tp = c(0,0,0), tp_id = 1, pp = 0, mp = 0, factr = 1e-08, fb = 0) {
  if(length(time_ID) != length(log_price)) {
    warning(paste0("t and log_p vectors must be same length.\n"))
    return(print("Execution aborted..."))
  }
  start_time <- Sys.time()
  
  ## Only using [T1, T2] for modelling
  t <- head(time_ID, length(time_ID) - hold_out)
  log_p <- head(log_price, length(log_price) - hold_out)
  
  ## Functions ----
  
  LPPLS <- function(t, A, B, C1, C2, tc, m, omega) {
    d <- omega * log(tc - t)
    A + (tc - t)^m * (B + C1 * cos(d) + C2 * sin(d))
  }
  
  n <- function(t) {
    length(t)
  }
  y <- function(log_p) {
    log_p
  }
  f <- function(t, tc, m) {
    (tc - t)^m
  }
  g <- function(t, tc, m, omega) {
    (tc - t)^m * cos(omega * log(tc - t))
  }
  h <- function(t, tc, m, omega) {
    (tc - t)^m * sin(omega * log(tc - t))
  }
  sum_f <- function(t, tc, m) {
    sum(f(t, tc, m))
  }
  sum_g <- function(t, tc, m, omega) {
    sum(g(t, tc, m, omega))
  }
  sum_h <- function(t, tc, m, omega) {
    sum(h(t, tc, m, omega))
  }
  sum_ff <- function(t, tc, m) {
    sum(f(t, tc, m) * f(t, tc, m))
  }
  sum_gg  <- function(t, tc, m, omega) {
    sum(g(t, tc, m, omega) * g(t, tc, m, omega))
  }
  sum_hh  <- function(t, tc, m, omega) {
    sum(h(t, tc, m, omega) * h(t, tc, m, omega))
  }
  sum_fg <- function(t, tc, m, omega) {
    sum(f(t, tc, m) * g(t, tc, m, omega))
  }
  sum_fh <- function(t, tc, m, omega) {
    sum(f(t, tc, m) * h(t, tc, m, omega))
  }
  sum_gh <- function(t, tc, m, omega) {
    sum(g(t, tc, m, omega) * h(t, tc, m, omega))
  }
  sum_y <- function(log_p) {
    sum(y(log_p))
  }
  sum_yf <- function(log_p, t, tc, m) {
    sum(y(log_p) * f(t, tc, m))
  }
  sum_yg <- function(log_p, t, tc, m, omega) {
    sum(y(log_p) * g(t, tc, m, omega))
  }
  sum_yh <- function(log_p, t, tc, m, omega) {
    sum(y(log_p) * h(t, tc, m, omega))
  }
  
  ## Symbols ====
  ## Elements in XX matrix
  xx11 <- Symbol("xx11")
  xx12 <- Symbol("xx12")
  xx13 <- Symbol("xx13")
  xx14 <- Symbol("xx14")
  xx21 <- Symbol("xx21")
  xx22 <- Symbol("xx22")
  xx23 <- Symbol("xx23")
  xx24 <- Symbol("xx24")
  xx31 <- Symbol("xx31")
  xx32 <- Symbol("xx32")
  xx33 <- Symbol("xx33")
  xx34 <- Symbol("xx34")
  xx41 <- Symbol("xx41")
  xx42 <- Symbol("xx42")
  xx43 <- Symbol("xx43")
  xx44 <- Symbol("xx44")
  xy1 <- Symbol("xy1")
  xy2 <- Symbol("xy2")
  xy3 <- Symbol("xy3")
  xy4 <- Symbol("xy4")
  
  XX <- symengine::Matrix(
    c(xx11, xx12, xx13, xx14,
      xx21, xx22, xx23, xx24,
      xx31, xx32, xx33, xx34,
      xx41, xx42, xx43, xx44
    ),
    byrow = TRUE,
    nrow = 4
  )
  
  Xy <- symengine::Vector(xy1, xy2, xy3, xy4)
  
  beta <- symengine::solve(XX, Xy)

  A <- DoubleVisitor(beta[1], args = c(
    xx11, xx12, xx13, xx14,
    xx21, xx22, xx23, xx24,
    xx31, xx32, xx33, xx34,
    xx41, xx42, xx43, xx44,
    xy1, xy2, xy3, xy4
  ))
  B <- DoubleVisitor(beta[2], args = c(
    xx11, xx12, xx13, xx14,
    xx21, xx22, xx23, xx24,
    xx31, xx32, xx33, xx34,
    xx41, xx42, xx43, xx44,
    xy1, xy2, xy3, xy4
  ))
  C1 <- DoubleVisitor(beta[3], args = c(
    xx11, xx12, xx13, xx14,
    xx21, xx22, xx23, xx24,
    xx31, xx32, xx33, xx34,
    xx41, xx42, xx43, xx44,
    xy1, xy2, xy3, xy4
  ))
  C2 <- DoubleVisitor(beta[4], args = c(
    xx11, xx12, xx13, xx14,
    xx21, xx22, xx23, xx24,
    xx31, xx32, xx33, xx34,
    xx41, xx42, xx43, xx44,
    xy1, xy2, xy3, xy4
  ))
  
  ## beta_calculator() calculates values of liniear parameters, given nonlinear estimates as inputs.
  beta_calculator <- function(log_p, t, tc, m, omega) {
    force(tc)
    force(m)
    force(omega)
    a <- drop(A(
      xx11 = n(t),
      xx12 = sum_f(t, tc, m),
      xx13 = sum_g(t, tc, m, omega),
      xx14 = sum_h(t, tc, m, omega),
      xx21 = sum_f(t, tc, m),
      xx22 = sum_ff(t, tc, m),
      xx23 = sum_fg(t, tc, m, omega),
      xx24 = sum_fh(t, tc, m, omega),
      xx31 = sum_g(t, tc, m, omega),
      xx32 = sum_fg(t, tc, m, omega),
      xx33 = sum_gg(t, tc, m, omega),
      xx34 = sum_gh(t, tc, m, omega),
      xx41 = sum_h(t, tc, m, omega),
      xx42 = sum_fh(t, tc, m, omega),
      xx43 = sum_gh(t, tc, m, omega),
      xx44 = sum_hh(t, tc, m, omega),
      xy1 = sum_y(log_p),
      xy2 = sum_yf(log_p, t, tc, m),
      xy3 = sum_yg(log_p, t, tc, m, omega),
      xy4 = sum_yh(log_p, t, tc, m, omega)
    ))
    b <- drop(B(
      xx11 = n(t),
      xx12 = sum_f(t, tc, m),
      xx13 = sum_g(t, tc, m, omega),
      xx14 = sum_h(t, tc, m, omega),
      xx21 = sum_f(t, tc, m),
      xx22 = sum_ff(t, tc, m),
      xx23 = sum_fg(t, tc, m, omega),
      xx24 = sum_fh(t, tc, m, omega),
      xx31 = sum_g(t, tc, m, omega),
      xx32 = sum_fg(t, tc, m, omega),
      xx33 = sum_gg(t, tc, m, omega),
      xx34 = sum_gh(t, tc, m, omega),
      xx41 = sum_h(t, tc, m, omega),
      xx42 = sum_fh(t, tc, m, omega),
      xx43 = sum_gh(t, tc, m, omega),
      xx44 = sum_hh(t, tc, m, omega),
      xy1 = sum_y(log_p),
      xy2 = sum_yf(log_p, t, tc, m),
      xy3 = sum_yg(log_p, t, tc, m, omega),
      xy4 = sum_yh(log_p, t, tc, m, omega)
    ))
    c1 <- drop(C1(
      xx11 = n(t),
      xx12 = sum_f(t, tc, m),
      xx13 = sum_g(t, tc, m, omega),
      xx14 = sum_h(t, tc, m, omega),
      xx21 = sum_f(t, tc, m),
      xx22 = sum_ff(t, tc, m),
      xx23 = sum_fg(t, tc, m, omega),
      xx24 = sum_fh(t, tc, m, omega),
      xx31 = sum_g(t, tc, m, omega),
      xx32 = sum_fg(t, tc, m, omega),
      xx33 = sum_gg(t, tc, m, omega),
      xx34 = sum_gh(t, tc, m, omega),
      xx41 = sum_h(t, tc, m, omega),
      xx42 = sum_fh(t, tc, m, omega),
      xx43 = sum_gh(t, tc, m, omega),
      xx44 = sum_hh(t, tc, m, omega),
      xy1 = sum_y(log_p),
      xy2 = sum_yf(log_p, t, tc, m),
      xy3 = sum_yg(log_p, t, tc, m, omega),
      xy4 = sum_yh(log_p, t, tc, m, omega)
    ))
    c2 <- drop(C2(
      xx11 = n(t),
      xx12 = sum_f(t, tc, m),
      xx13 = sum_g(t, tc, m, omega),
      xx14 = sum_h(t, tc, m, omega),
      xx21 = sum_f(t, tc, m),
      xx22 = sum_ff(t, tc, m),
      xx23 = sum_fg(t, tc, m, omega),
      xx24 = sum_fh(t, tc, m, omega),
      xx31 = sum_g(t, tc, m, omega),
      xx32 = sum_fg(t, tc, m, omega),
      xx33 = sum_gg(t, tc, m, omega),
      xx34 = sum_gh(t, tc, m, omega),
      xx41 = sum_h(t, tc, m, omega),
      xx42 = sum_fh(t, tc, m, omega),
      xx43 = sum_gh(t, tc, m, omega),
      xx44 = sum_hh(t, tc, m, omega),
      xy1 = sum_y(log_p),
      xy2 = sum_yf(log_p, t, tc, m),
      xy3 = sum_yg(log_p, t, tc, m, omega),
      xy4 = sum_yh(log_p, t, tc, m, omega)
    ))
    c(a, b, c1, c2)
  }
  
  ## res(tc, m, omega, log_p, t) calculates residuals.
  res <- function(tc, m, omega, log_p, t) {
    force(tc)
    force(m)
    force(omega)
    t <- t[t < tc]
    beta <- beta_calculator(log_p, t, tc, m, omega)
    a <- beta[1]
    b <- beta[2]
    c1 <- beta[3]
    c2 <- beta[4]
    log_p - (
      a + (tc - t)^m * (b + c1 * cos(omega * log(tc - t)) + c2 * sin(omega * log(tc - t)))
    )
  }
  
  ## Note: SSE-functions are used to optimize wrt nonlinear parameters, given linear parameters.
  ## Sum of squared errors.
  ## par vector is {tc, m, omega}, used as initial parameter values in optim().
  SSE1 <- function(par, log_p, t) {
    tc <- par[1]
    m <- par[2]
    omega <- par[3]
    drop(sum(res(tc, m, omega, log_p, t)^2, na.rm = TRUE))
  }
  
  ## SSE2, SSE2(par, tc, log_p, t) function.
  ## Sum of squared errors with fixed tc.
  ## par vector is {m, omega}, used as initial parameter values in optim().
  SSE2 <- function(par, tc, log_p, t) {
    m <- par[1]
    omega <- par[2]
    drop(sum(res(tc, m, omega, log_p, t)^2, na.rm = TRUE))
  }
  
  ## SSE3, SSE3(tc, m, omega, log_p, t) function.
  ## Sum of squared errors with estimated m and omega.
  ## par is tc, used as initial parameter value in optim().
  SSE3 <- function(tc, m, omega, log_p, t) {
    drop(sum(res(tc, m, omega, log_p, t)^2, na.rm = TRUE))
  }

  ## mode F1 ====
  if(mode == "F1") {
    if(fb){print("mode: F1")}
    ## -> Calculations ----
    
    fit <- list() ## Initialize output list
    ## Initialize tibble for output fit[[2]]
    fit_tbl <- tibble(tc = 0, m = 0, omega = 0, A = 0, B = 0, C1 = 0, C2 = 0, D = 0, value = 0, ID = 0)
    
    out_of_range_tracker <- list(B = list(), D = list())
    if(fb){print("Iteration 1...")}
    
    ## L-BFGS-B (box-constrained BFGS)
    ## Constrained to avoid NaN's from tc = t.
    opt_tmp <- optim(par = c(tc_init, m_init, o_init),
                     fn = SSE1,
                     log_p = log_p,
                     t = t,
                     lower = c(n(t) + 1, lower[1], lower[2]),
                     upper = c(n(t) + fh, upper[1], upper[2]),
                     method = "L-BFGS-B",
                     control = list(factr = factr)
    )
    
    
    # Nelder-Mead (no constraints)
    # opt_tmp <- optim(par = c(tc_init, m_init, o_init),
    #                  fn = SSE1,
    #                  log_p = log_p,
    #                  t = t,
    #                  method = "Nelder-Mead"
    # )
    
    
    ## Calculate linear koefficients
    beta_vals <- beta_calculator(log_p, t, opt_tmp$par[[1]], opt_tmp$par[[2]], opt_tmp$par[[3]])
    #if(beta_vals[[2]] >= upper[3]) {warning(paste0("B out of range (", 1, ").\n"))}
    if(beta_vals[[2]] >= upper[3]) {out_of_range_tracker$B[length(out_of_range_tracker$B) + 1] = 1}
    
    ## Calculate damping (D)
    ## m * abs(B) / (omega * sqrt(C1^2 + C2^2))
    damp <- opt_tmp$par[[2]] * abs(beta_vals[2]) / (opt_tmp$par[[3]] * sqrt(beta_vals[3]^2 + beta_vals[4]^2))
    #if(damp <= lower[4]) {warning(paste0("D out of range (", 1, ").\n"))}
    if(damp <= lower[4]) {out_of_range_tracker$D[length(out_of_range_tracker$D) + 1] = 1}
    
    
    ## Add list of fitted coefficients to list of fits
    
    fit_tbl[1, ] <- list(tc = opt_tmp$par[[1]], m = opt_tmp$par[[2]], omega = opt_tmp$par[[3]], A = beta_vals[1], B = beta_vals[2], C1 = beta_vals[3], C2 = beta_vals[4], D = damp, value = opt_tmp$value, ID = 1)
    
    if(fb){print("Iteration 1 done")}
    
    ## Make a number of additional fits with random starting points
    if(num_searches > 1) {
      for(i in 2:num_searches) {
        if(fb){print(paste0("Iteration ",i,"..."))}
        
        ## L-BFGS-B (box-constrained BFGS)
        opt_tmp <- optim(par = c(runif(1, n(t) + 1 , n(t) + fh), runif(1, lower[1], upper[1]), runif(1, lower[2], upper[2])),
                         fn = SSE1,
                         log_p = log_p,
                         t = t,
                         lower = c(n(t) + 1, lower[1], lower[2]),
                         upper = c(n(t) + fh, upper[1], upper[2]),
                         method = "L-BFGS-B",
                         control = list(factr = factr)
        )
        
        ## Nelder-Mead (no constraints)
        # opt_tmp <- optim(par = c(tc_init, m_init, o_init),
        #                  fn = SSE1,
        #                  log_p = log_p,
        #                  t = t,
        #                  method = "Nelder-Mead"
        # )
        
        
        # if(opt_tmp$value < fit$value){
        #   fit <- opt_tmp
        # }
        ## Calculate linear koefficients
        beta_vals <- beta_calculator(log_p, t, opt_tmp$par[[1]], opt_tmp$par[[2]], opt_tmp$par[[3]])
        #if(beta_vals[[2]] >= upper[3]) {warning(paste0("B out of range (", i, ").\n"))}
        if(beta_vals[[2]] >= upper[3]) {out_of_range_tracker$B[length(out_of_range_tracker$B) + 1] = i}
        
        damp <- opt_tmp$par[[2]] * abs(beta_vals[2]) / (opt_tmp$par[[3]] * sqrt(beta_vals[3]^2 + beta_vals[4]^2))
        #if(damp <= lower[4]) {warning(paste0("D out of range (", i, ").\n"))}
        if(damp <= lower[4]) {out_of_range_tracker$D[length(out_of_range_tracker$D) + 1] = i}
        
        ## Add list of fitted coefficients to list of fits
        fit_tbl[i, ] <- list(tc = opt_tmp$par[[1]], m = opt_tmp$par[[2]], omega = opt_tmp$par[[3]], A = beta_vals[1], B = beta_vals[2], C1 = beta_vals[3], C2 = beta_vals[4], D = damp, value = opt_tmp$value, ID = i)
        if(fb){print(paste0("Iteration ",i," done"))}
      }
    }
    
    ## Sorted list by value of objective function
    if(fb){print(paste0("Sorting list..."))}
    fit[[2]] <- arrange(fit_tbl, value)
    fit[[1]] <- as.list(fit[[2]][1, ]) ## fit[[2]] is a tibble, fit is a list
  }
  
  
  
  ## mode F2 ====
  if(mode == "F2" || mode == "MPL") { ## "F2"-mode output needed for "MPL"-mode
    if(fb){print("mode: F2")}
    
    ## Initialize out-of-range tracker
    out_of_range_tracker <- list(B = list(), D = list())
    
    ## List of fits wrt m and omega for a single fixed value of tc.
    ## This is the combination of m and omega that minimizes SSE2.
    
    fit2 <- list() ## For each value of tc, fit2 contains a list of all fits
    ## with random initial values.
    ## For each value of tc_k, fit2 stores the contents of fit2_tmp
    ## in iteration k.
    
    fit2_filtered <- list() ## fit2 filtered
    
    fit2_best_for_each_tc  <- tibble(ID = 0, value = 0, tc = 0, m = 0, omega = 0, A = 0, B = 0, C1 = 0, C2 = 0, D = 0)
    ## Initialize list with best fit optimized wrt m and omega for each value.
    ## of tc. I.e. pick the best from fit2 for each tc.
    ## (Values are dummies. Doesn't work with NA og NULL.)
    
    ## Initialize fit2_best_for_each_tc filtered for B < 0
    fit2_best_for_each_tc_filtered  <- tibble(ID = 0, value = 0, tc = 0, m = 0, omega = 0, A = 0, B = 0, C1 = 0, C2 = 0, D = 0)
    
    fit <- list()  ## First element is the best fit wrt tc. 
    ## This is the tc that minimizes SSE3.
    ## Second element is fit2_best_for_each_tc.
    
    fit2_best_row <- 1
    
    fit[[4]] <- list()
    
    for(k in 1:fh) {
      #3 Optimize m and omega with fixed tc for each value of tc in [T2+1, T2+fh]:
      tc_k <- n(t) + k
      
      ## -> Calculations ----
      
      if(fb){print(paste0("tc_", k, ", Iteration 1..."))}
      
      ## L-BFGS-B
      opt2_counts <- list()
      set.seed(1) ## Set seed to match random values for trace_plot optimization
      opt_tmp <- optim(par = c(m_init, o_init),
                       fn = SSE2,
                       tc = tc_k,
                       log_p = log_p,
                       t = t,
                       lower = c(lower[1], lower[2]),
                       upper = c(upper[1], upper[2]),
                       method = "L-BFGS-B",
                       control = list(factr = factr)
      )
      if(tp[1] == 1 || tp[2] == 1 || tp[3] == 1) {opt2_counts[1] <- opt_tmp$counts["function"]} ## Use for trace plot
      
      ## Nelder-Mead (no constraints)
      # opt_tmp <- optim(par = c(m_init, o_init),
      #                  fn = SSE2,
      #                  tc = tc_k,
      #                  log_p = log_p,
      #                  t = t,
      #                  method = "Nelder-Mead"
      # )
      
      
      
      ## Calculate linear koefficients
      beta_vals <- beta_calculator(log_p, t, tc_k, opt_tmp$par[[1]], opt_tmp$par[[2]])

      #if(beta_vals[[2]] >= upper[3]) {warning(paste0("F2, iteration 1: B out of range (", 1, ").\n"))}
      if(beta_vals[[2]] >= upper[3]) {out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] = list(tc_num = k, rand_iter_num = 1)}
     
      ## Calculate damping (D)
      ## m * abs(B) / (omega * sqrt(C1^2 + C2^2))
      damp <- opt_tmp$par[[1]] * abs(beta_vals[2]) / (opt_tmp$par[[2]] * sqrt(beta_vals[3]^2 + beta_vals[4]^2))
      #if(damp <= lower[4]) {warning(paste0("F2, iteration 1: D out of range (", 1, ").\n"))}
      if(damp <= lower[4]) {out_of_range_tracker$B[[length(out_of_range_tracker$D) + 1]] = list(tc_num = k, rand_iter_num = 1)}

      ## Add list of fitted coefficients to list of fits
      fit2_tmp <- tibble(ID = 1, value = opt_tmp$value, tc = tc_k, m = opt_tmp$par[[1]], omega = opt_tmp$par[[2]], A = beta_vals[1], B = beta_vals[2], C1 = beta_vals[3], C2 = beta_vals[4], D = damp)
      
      tmp_fits2 <- list()
      
      if(fb){print(paste0("tc_", k, ", Iteration 1 done"))}
      
      ## Make a number of additional fits with random starting points
      if(num_searches > 1) {
        for(i in 2:num_searches) {
          if(fb){print(paste0("tc_", k, ", Iteration ", i, "..."))}
          
          ## Nelder-Mead (no constraints)
          # opt_tmp <- optim(par = c(runif(1, lower[1], upper[1]), runif(1, lower[2], upper[2])),
          #                  fn = SSE2,
          #                  tc = tc_k,
          #                  log_p = log_p,
          #                  t = t,
          #                  method = "Nelder-Mead"
          # )
          
          
          ## L-BFGS-B (constrained)
          set.seed(i) ## Same seed as used in trace_plot
          opt_tmp <- optim(par = c(runif(1, lower[1], upper[1]), runif(1, lower[2], upper[2])), ## Random init values for tc, m and omega
                           fn = SSE2,
                           tc = tc_k,
                           log_p = log_p,
                           t = t,
                           lower = c(lower[1], lower[2]),
                           upper = c(upper[1], upper[2]),
                           method = "L-BFGS-B",
                           control = list(factr = factr)
          )
          if(tp[1] == 1 || tp[2] == 1 || tp[3] == 1) {opt2_counts[i] <- opt_tmp$counts["function"]} ## Use for trace plot
          

          
          # if(opt_tmp$value < fit$value){
          #   fit <- opt_tmp
          # }
          ## Calculate linear koefficients
          beta_vals <- beta_calculator(log_p, t, tc_k, opt_tmp$par[[1]], opt_tmp$par[[2]])
          #if(beta_vals[[2]] >= upper[3]) {warning(paste0("F2, iteration >1: B out of range (", i, ").\n"))}
          if(beta_vals[[2]] >= upper[3]) {out_of_range_tracker$B[[length(out_of_range_tracker$B) + 1]] = list(tc_num = k, rand_iter_num = i)}
          
          damp <- opt_tmp$par[[1]] * abs(beta_vals[2]) / (opt_tmp$par[[2]] * sqrt(beta_vals[3]^2 + beta_vals[4]^2))
          #if(damp <= lower[4]) {warning(paste0("F2, iteration >1: D out of range (", i, ").\n"))}
          if(damp <= lower[4]) {out_of_range_tracker$B[[length(out_of_range_tracker$D) + 1]] = list(tc_num = k, rand_iter_num = i)}
          
          ## Add list of fitted coefficients to list of fits for objective function F2
          fit2_tmp[i, ]<- list(ID = i, value = opt_tmp$value, tc = tc_k, m = opt_tmp$par[[1]], omega = opt_tmp$par[[2]], A = beta_vals[1], B = beta_vals[2], C1 = beta_vals[3], C2 = beta_vals[4], D = damp)
          
        }
        tmp_fits2[[k]] <- fit2_tmp
        if(fb){print(paste0("tc_", k, ", Iteration ", i, " done"))}
      }
      fit[[4]][[k]] <- tmp_fits2
      
      ## Sort list by value of objective function F2.
      ## Keep all fits for tc_k in a list.
      if(fb){print(paste0("Sorting list of all fits for tc_", k))}


      ## Filter out results where B > 0.
      ## Then sort by SSE value (smallest at top).
      fit2[[k]] <- fit2_tmp %>% arrange(value)
      fit2_filtered[[k]] <- fit2_tmp %>% filter(B < upper[3]) %>% arrange(value)
      
      ## Keep best fit for each value of tc
      if(fb){print(paste0("Saving best fit for tc_", k))}

      ## Skip tc values where no fits have B < 0
      if(nrow(fit2_filtered[[k]]) > 0) {
        fit2_best_for_each_tc_filtered[fit2_best_row, ] <- fit2_filtered[[k]][1, ]
        fit2_best_row = fit2_best_row + 1
      }
      
      fit2_best_for_each_tc[k, ] <- fit2[[k]][1, ]
    }

    ## Sort list of best fits for each value of tc
    if(fb){print(paste0("Sorting list of best fits for each value of tc"))}

    fit2_best_for_each_tc <- arrange(fit2_best_for_each_tc, value)
    fit2_best_for_each_tc_filtered <- arrange(fit2_best_for_each_tc_filtered, value)
    
    ## Final chosen fit is the one at the top of the list
    ## If no fits passed filter, best fit is unfiltered:
    if(nrow(fit2_best_for_each_tc_filtered) == 1 && sum(fit2_best_for_each_tc_filtered) == 0) {
      fit2_best <- fit2_best_for_each_tc[1, ]
    } else { ## Else best fit is filtered
      fit2_best <- fit2_best_for_each_tc_filtered[1, ]
    }
    
    ## Optimize wrt tc given estimated m and omega
    if(fb){print(paste0("Optimizing SSE3..."))}

    ## Brent
    opt_tmp <- optim(par = tc_init,
                     fn = SSE3,
                     m = fit2_best$m,
                     omega = fit2_best$omega,
                     log_p = log_p,
                     t = t,
                     lower = n(t) + 1,
                     upper = n(t) + fh,
                     method = "Brent" ## Uses optimize() for single var opt
    )
    
    ## Gauss-Newton
    # opt_tmp <- nls.lm(
    #   par = tc_init,
    #   fn = SSE3,
    #   m = fit2_best$m,
    #   omega = fit2_best$omega,
    #   log_p = log_p,
    #   t = t,
    #   control = nls.control(maxiter = 500, tol = 1e-10, minFactor = 1/4096)
    # )
    
    ## Nelder-Mead (no constraints)
    # opt_tmp <- optim(par = tc_init, 
    #                  fn = SSE3, 
    #                  m = fit2_best$m,
    #                  omega = fit2_best$omega,
    #                  log_p = log_p,
    #                  t = t,
    #                  method = "Nelder-Mead"
    # )
    
    if(fb){print(paste0("Finished optimizing SSE3"))}
    
    ## Calculate linear koefficients
    if(fb){print(paste0("Calculating linear koefficient..."))}
    beta_vals <- beta_calculator(log_p, t, opt_tmp$par[[1]], fit2_best$m, fit2_best$omega)

    if(fb){print(paste0("Calculating damp..."))}
    damp <- fit2_best$m * abs(beta_vals[2]) / (fit2_best$omega * sqrt(beta_vals[3]^2 + beta_vals[4]^2))

    ## Final fit
    if(fb){print(paste0("Saving list of final fit parameters..."))}

    fit[[1]] <- list(tc = opt_tmp$par[[1]], m = fit2_best$m, omega = fit2_best$omega, A = beta_vals[1], B = beta_vals[2], C1 = beta_vals[3], C2 = beta_vals[4], D = damp, value = opt_tmp$value)
    
    if(nrow(fit2_best_for_each_tc_filtered) == 1 && sum(fit2_best_for_each_tc_filtered) == 0) {
      ## Add list of fits for objective function F2 for each value of tc.
      fit[[2]] <- fit2_best_for_each_tc
    } else {
      ## Add list of fits for objective function F2 for each value of tc.
      fit[[2]] <- fit2_best_for_each_tc_filtered
      fit[[3]] <- fit2_best_for_each_tc
    }
    if(nrow(fit2_best_for_each_tc_filtered) == 1 && sum(fit2_best_for_each_tc_filtered) == 0) {
      warning("No fits passed the filter. Only unfiltered fits returned.\n")
    }
    if(fb){print(paste0("Model fit done"))}
  }
  
  ## mode MPL ====
  ## Modified Profile Likelihood
  if(mode == "MPL") {
    
    if(fb){print(paste0("Mode: MPL"))}
    if(fb){print(paste0("Calculating likelihood intervals..."))}
    ## Gradient for single t_i
    ## Filimonov2017, equation (B16)
    grad_LPPLS <- function(Psi, t, tc) {
      ## par = (m, omega, a, b, c1, c2)
      m = Psi[1]; omega = Psi[2]; b = Psi[4]; c1 = Psi[5]; c2 = Psi[6];
      grad = rep(0, 6) ## length(Psi) = 6
      tau = tc - t
      
      ## d LPPLS(t_i; tc, psi) / d m
      grad[1] = tau^m * log(tau) * (b + c1 * cos(omega * log(tau)) + c2 * sin(omega * log(tau)))
      
      ## d LPPLS(t_i; tc, psi) / d omega
      grad[2] = tau^m * log(tau) * ((-1) * c1 * sin(omega * log(tau)) + c2 * cos(omega * log(tau)))
      
      ## d LPPLS(t_i; tc, psi) / d a
      grad[3] = 1
      
      ## d LPPLS(t_i; tc, psi) / d b
      grad[4] = tau^m
      
      ## d LPPLS(t_i; tc, psi) / d c1
      grad[5] = tau^m * cos(omega * log(tau))
      
      ## d LPPLS(t_i; tc, psi) / d c2
      grad[6] = tau^m * sin(omega * log(tau))
      
      grad
    }
    
    ## X matrix
    ## Filimonov2017, equation (36)
    X_MPL <- function(Psi, tc, t) {
      ## Psi = (m, omega, a, b, c1, c2)
      n = length(t)
      p = 6 ## length(Psi) = 6
      grad_i <- rep(0, p)
      X <- matrix(0, nrow = n, ncol = p)
      
      for(i in 1:n) {
        grad_i = grad_LPPLS(Psi, t[i], tc)
        for(j in 1:p) {
          X[i, j] = grad_i[j]
        }
      }
      X
    }
    
    ## H matrix
    ## Filimonov2017, equation (37)
    H_MPL <- function(Psi, tc, log_p, t) {
      ## Psi = (m, omega, a, b, c1, c2)
      n = length(t)
      p = 6 ## length(Psi) = 6
      m = Psi[1]; omega = Psi[2]; a = Psi[3]; b = Psi[4]; c1 = Psi[5]; c2 = Psi[6];
      LPPLS <- rep(0, n)
      H <- matrix(0, nrow = p, ncol = p)
      
      for(i in 1:n) {
        tau = tc - t[i]
        LPPLS[i] = a + tau^m * (b + c1 * cos(omega * log(tau)) + c2 * sin(omega * log(tau)));
        res = (log_p[i] - LPPLS[i])
        
        ## Second derivative of LPPLS wrt:
        ## (m, m)
        H[1,1] = H[1,1] + res * tau^m * log(tau)^2 * (b + c1 * cos(omega * log(tau)) + c2 * sin(omega * log(tau)))
        ## (m, omega)
        H[1,2] = H[1,2] + res * tau^m * log(tau)^2 * ((-1) * c1 * sin(omega * log(tau)) + c2 * cos(omega * log(tau)))
        
        ## (m, b)
        H[1,4] = H[1,4] + res * tau^m * log(tau)
        ## (m, c1)
        H[1,5] = H[1,5] + res * tau^m  * log(tau) * cos(omega * log(tau))
        ## (m, c2)
        H[1,6] = H[1,6] + res * tau^m  * log(tau) * sin(omega * log(tau))
        
        ## (omega, omega)
        H[2,2] = H[2,2] + res * (-1) * tau^m * log(tau)^2 * (c1 * cos(omega * log(tau)) + c2 * sin(omega * log(tau)))
        
        ## (omega, c1)
        H[2,5] = H[2,5] + res * (-1) * tau^m * log(tau) * sin(omega * log(tau))
        ## (omega, )c2
        H[2,6] = H[2,6] + res * tau^m * log(tau) * cos(omega * log(tau))
      }
      
      ## (omega, m)
      H[2,1] = H[1,2]
      
      ## (b, m)
      H[4,1] = H[1,4]
      
      ## (c1, m)
      H[5,1] = H[1,5]
      
      ## (c2, m)
      H[6,1] = H[1,6]
      
      ## The rest is zero
      
      H
    }
    
    ## Negativ modified profile log likelihood
    ## Only works for Filimonov equation (mode 0)
    ## Filimonov2017, equation 38
    ## Psi_hat = (m, omega, a, b, c1, c2)
    ## Negative because optim() minimizes as default (however, not using optim..... See below)
    LL_MPL_neg <- function(tc, tc_hat, Psi_hat_tc, Psi_hat, log_p, t) {
      n <- length(t)
      p <- length(Psi_hat)
      m_hat <- Psi_hat[1]
      omega_hat <- Psi_hat[2]
      s_tc <- SSE3(tc, m_hat, omega_hat, log_p, t)/n ## s_tc = F2 = SSE3
      X_hat_tc <- X_MPL(Psi_hat_tc, tc, t)
      X_hat <- X_MPL(Psi_hat, tc_hat, t)
      H_hat_tc <- H_MPL(Psi_hat_tc, tc, log_p, t)
      
    
      mpll <- log(
        sqrt(
          abs(
            det(
              crossprod(X_hat_tc, X_hat_tc) - H_hat_tc
            )
          )
        ) / abs(
          det(
            crossprod(X_hat, X_hat_tc)
          )
        )
      #) - ((n + p + 2)/2) * log(s_tc) ## Fejl
      ) - ((n - p - 2)/2) * log(s_tc)
      
      -mpll
    }
    
    ## Likelihood Interval
    ## Filimonov2017, equation (39)
    ##
    ## Likelihood Interval for modified profile likelihood.
    
    ## Filimonov2017, equation (39).
    
    ## Inputs:
    ## 1) cutoff, Cutoff for the likelihood interval.
    ##         log(cutoff) is used, and R is defined as R = log L - log L_max.
    ##         The likelihood interval (LI) is defined by log(cutoff) < R < 0.
    ## 2) fit, takes full F2 fit object from Filimonov().
    ## 4) log_p
    ## 5) t
    ## 6) fh, length of forecast horizon.
    ## 7) fb, boolean for printing feedback.
    ##
    ## par = c(tc, Psi_hat_tc), Parameters estimated for a given value of tc.
    ##              Psi_hat = c(m_hat_tc, omega_hat_tc, A_hat_tc, B_hat_tc, C1_hat_tc, C2_hat_tc)
    ## par_hat = c(tc_hat, Psi_hat): Parameters estimated for ML estimate of tc.
    ##              Psi_hat = c(m_hat, omega_hat, A_hat, B_hat, C1_hat, C2_hat)
    ##
    ## Outputs:
    ## 1) LI, vector of two interval endpoints.
    ## 2) LI_vect, vector of all points in the interval.
    ## 3) R, vector of R valuesfor each value of tc in (n + 1, n + fh).
    ## LL, vector of LL values for each value of tc in (n + 1, n + fh).
    ## 4) MLL, vector of LL values for each value of tc in (n + 1, n + fh).
    ## 5) LI_IDs, ID's of the elements in the R vector that match R > log(cutoff)
    
    
    LI_MPL  <- function(cutoff, fit, log_p, t, fh, fb) {
      
      n <- length(t)
      log_cutoff <- log(cutoff)
      
      ## Initialize tibble
      ## One row for each value of tc
      R_MPL_tbl <- tibble(R = rep(0, fh), LL = rep(0, fh))
      
      par_hat <-  as.numeric(fit[[1]])[1:7] ## Profile likelihood estimate from F2
      
      for(i in 1:fh) {
        ## Get fit for tc=N+i
        if(fb && i%%10 == 1) {print(paste0("Calculating LL for tc = ", n + i, "...", n + i + 9))}
        par_hat_tc <- as.numeric(fit[[2]] %>% filter(tc == n + i))[3:9]
        R_MPL_tbl$LL[i] <- -LL_MPL_neg(tc = n + i, tc_hat = par_hat[1], Psi_hat_tc = par_hat_tc[2:7], Psi_hat = par_hat[2:7], log_p, t)
      }
      
      MLL <- max(R_MPL_tbl$LL, na.rm = TRUE) ## Calculate MLL as max(LL)
      R_MPL_tbl$R <- R_MPL_tbl$LL - MLL ## Calculate R
      
      tc_hat_mpl <- n + which(R_MPL_tbl$LL == MLL)
      
      tc_range <- (n + 1):(n + fh)
      ## (R_MPL_tbl$R && R_MPL_tbl$R < 0 bør altid gælde!)
      #LI_IDs <- which(log_cutoff < R_MPL_tbl$R && R_MPL_tbl$R < 0) 
      
      LI = list(c(0,0), c(0,0), c(0,0)) ## Initialize list of LI's for each of three cutoffs
      for(i in 1:3) {
        LI_IDs <- which(log_cutoff[i] < R_MPL_tbl$R)
        LI_vect <- tc_range[LI_IDs] ## Values of tc that match criteria
        
        ## If LI_vect is empty, range() will return c(-Inf, Inf).
        ## So check that LI_vect is not empty.
        if(length(LI_vect != 0)) {
          LI[[i]] <- range(LI_vect)
        } else {
          warning("No values of tc in LI.\n")
          LI[[i]] = NA
        }
        
      }
      
      if(fb) {print("Done calculating R")}
      
      list(LI = LI, R = R_MPL_tbl$R, LL = R_MPL_tbl$LL, MLL = MLL, tc_hat_mpl = tc_hat_mpl)
    }
    
    mpl_output <- LI_MPL(
      cutoff = mpl_cutoff, 
      fit = fit,
      log_p = log_p,
      t = t, 
      fh = fh, 
      fb = fb)

    
    #   df_plot <- data.frame(tc = (n + 1):(n - hold_out + fh), log_mpl = li$LL)
    #   df_fit <- data.frame(tc_hat = mpl_fit_ARGARCH$fit[[1]][[1]])
    #   df_ML <- data.frame(tc_mpl = li$LI_vect[which(li$LL[li$LI_IDs] == max(li$LL[li$LI_IDs]))])
    #   df_LI_l <- data.frame(LI_l = li$LI[1])
    #   df_LI_u <- data.frame(LI_u = li$LI[2])
    #   mpl_plot <- ggplot(data = df_plot, aes(x = tc, y = log_mpl)) +
    #     geom_point() +
    #     geom_vline(data = df_fit, mapping = aes(xintercept = tc_hat), color = "red") +
    #     geom_vline(data = df_ML, mapping = aes(xintercept = tc_mpl), color = "red", linetype = "dashed") +
    #     geom_vline(data = df_LI_l, mapping = aes(xintercept = LI_l), color = "green3", linetype = "dashed") +
    #     geom_vline(data = df_LI_u, mapping = aes(xintercept = LI_u), color = "green3", linetype = "dashed")
    # }
    
    if(fb){print(paste0("Likelihood intervals calculated"))}
    if(fb){print(paste0("Generating MPL plot"))}
    tc_seq <- (n(t) + 1):(n(t) + fh)
    df_plot <- data.frame(tc = tc_seq, log_mpl = mpl_output$LL)
    tc_hat = fit[[1]]$tc
    #tc_mpl = li$LI_vect[which(li$LL[li$LI_IDs] == max(li$LL[li$LI_IDs]))]
    tc_mpl =  tc_seq[which(mpl_output$LL == max(mpl_output$LL))]
    li1_l = mpl_output$LI[[1]][1]
    li1_u = mpl_output$LI[[1]][2]
    li2_l = mpl_output$LI[[2]][1]
    li2_u = mpl_output$LI[[2]][2]
    li3_l = mpl_output$LI[[3]][1]
    li3_u = mpl_output$LI[[3]][2]
    mpl_plot <- ggplot(data = df_plot, aes(x = tc, y = log_mpl)) +
      geom_point() +
      geom_vline(aes(xintercept = tc_hat), color = "blue") +
      geom_vline(aes(xintercept = tc_mpl), color = "blue", linetype = "dashed") +
      geom_vline(aes(xintercept = li1_l), color = "green3", linetype = "dashed") +
      geom_vline(aes(xintercept = li1_u), color = "green3", linetype = "dashed") +
      geom_vline(aes(xintercept = li2_l), color = "orange", linetype = "dotted") +
      geom_vline(aes(xintercept = li2_u), color = "orange", linetype = "dotted") +
      geom_vline(aes(xintercept = li3_l), color = "red", linetype = "dotdash") +
      geom_vline(aes(xintercept = li3_u), color = "red", linetype = "dotdash")
  }
  if(fb){print(paste0("MPL plot generated"))}
  
  ## Contour data ====
  ## Generate SSE data for contour plot of SSE wrt m and omega.
  contour_data <- NULL
  contour_plot <- NULL
  surface_plot <- NULL
  fit_plot <- NULL
  param_plot <- NULL
  matrix_plot <- NULL
  if(cp == 1 || sp == 1) {
    if(fb){print(paste0("Generating contour plot data..."))}
    tc_val <- fit[[1]]$tc
    x_contour = seq(lower[1], upper[1], length.out = 101) ## m. lower and upper are filter range.
    y_contour =  seq(lower[2], upper[2], length.out = 101) ## omega
    sse_contour <- function(x, y) {
      SSE2(par = c(x, y), tc = tc_val, log_p = log_p, t = t)
    }
    z_contour <- t(outer(x_contour, y_contour, Vectorize(sse_contour))) ## z matrix must be transposed for correct contour plot! (Don't know why)
    contour_data <- list(x = x_contour, y = y_contour, z = z_contour)
    if(fb){print(paste0("Fit contour data generated"))}
  }
  
  ## Contour plot ====
  if(cp == 1) {
    if(fb){print(paste0("Generating contour plot..."))}
    point <- c(fit[[1]]$m, fit[[1]]$omega)
    z_coord <- SSE2(par = c(point[1], point[2]), tc = tc_val, log_p = log_p, t = t)
    
    contour_plot <- plot_ly(
      x = contour_data$x, 
      y = contour_data$y, 
      z = contour_data$z,
      type = "contour",
      contours = list(
        showlabels = FALSE,
        exponentformat = "e"),
      ncontours = 10,
      colorbar = list(exponentformat = "e", title = "SSE")
    )
    contour_plot <- contour_plot %>% layout(
      title = paste0("tc = ", tc_val),
      xaxis = list(title = "m"),
      yaxis = list(title = "omega")
    ) %>% add_trace(
      x = point[1], y = point[2], type = "scatter", mode = "markers", showlegend = FALSE) %>%
      add_trace(x = c(point[1], point[1], lower[1]), y = c(lower[2], point[2], point[2]), type = "scatter", mode = "lines", line = list(dash = "dash", width = 1, color = "orange"), showlegend = FALSE) %>%
      add_annotations(
        xref = "x",
        yref = "y",
        x = point[1],
        y = point[2],
        #xanchor = "middle",
        xanchor = "left",
        yanchor = "bottom",
        text = paste0("(", point[1], ", ", point[2], ", ", z_coord,")"),
        font = list(color = 'white', size = 10),
        showarrow = FALSE
      )
    
    
    if(fb){print(paste0("Contour plot generated"))}
  }
  ## Note:  
  ##  Hvis exponentformat = "e": scientific notation  
  ## exponentformat = "b":
  ## p for pico: $1e-12$  
  ## f for femto: $1e-15$ 
  ## osv (se https://en.wikipedia.org/wiki/Names_of_small_numbers )
  
  ## Surface plot ====
  if(sp == 1) {
    if(fb){print(paste0("Generating surface plot..."))}
    surface_plot <- plot_ly(
      x = contour_data$x, 
      y = contour_data$y, 
      z = contour_data$z
    ) %>% add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorbar = list(exponentformat = "e", title = "SSE")
    )
    surface_plot <- surface_plot %>% layout(
      title = paste0("tc = ", tc_val),
      scene = list(
        camera=list(
          eye = list(x=1.1, y=-1.3, z=0.1)
        ),
        xaxis = list(title = "m"),
        yaxis = list(title = "omega"),
        zaxis = list(title = "SSE")
      )
    )
    if(fb){print(paste0("Surface plot generated"))}
  }
  
  
  ## Trace plot ====
  if(tp[1] == 1 || tp[2] == 1 || tp[3] == 1) { ## boolean vector: c(1 on/off, 2, on/off, 3 on/off)
    if(mode == "F2" || mode == "MPL") {
      tc_val <- fit[[1]]$tc

      ## Calculate trace 

      set.seed(fit2_best$ID) ## Set same seed as was used for best fit.
      s <- seq_len(opt2_counts[[fit2_best$ID]])
      opt <- sapply(
        s, 
        # function(i) { 
        #   set.seed(1234)
        #   optim(par = c(runif(1, lower[1], upper[1]), runif(1, lower[2], upper[2])), ## Random init values for m and omega
        #         fn = SSE2, 
        #         tc = tc_val,
        #         log_p = log_p,
        #         t = t,
        #         method = "Nelder-Mead",
        #         control = list(maxit = i))$par
        # }
        
        function(i) { 
          set.seed(fit2_best$ID)
          opt_i <- optim(par = c(runif(1, lower[1], upper[1]), runif(1, lower[2], upper[2])), ## Random init values for m and omega
                         fn = SSE2, 
                         tc = tc_val,
                         log_p = log_p,
                         t = t,
                         lower = c(lower[1], lower[2]), 
                         upper = c(upper[1], upper[2]),
                         method = "L-BFGS-B",
                         control = list(maxit = i))$par
          ## Calculate B and add to vector if selector = 2 or 3
          if(tp[2] == 1 || tp[3] == 1){
            opt_i <- c(opt_i, beta_calculator(log_p, t, tc_val, opt_i[1], opt_i[2])[2])
          }
          opt_i
        }
      )
      if(tp[2] == 1 || tp[3] == 1){ ## selector 2 or 3: three rows (m, omega, B)
        opt <- matrix(c(c(m_init, o_init, beta_calculator(log_p, t, tc_val, m_init, o_init)[2]), opt), byrow = FALSE, nrow = 3)
      } else { ## selector 1: only two rows (m, omega)
        opt <- matrix(c(c(m_init, o_init,), opt), byrow = FALSE, nrow = 2)
      }
      tp_gen <- function(selector) { ## 1: (m, omega), 2: (B, m), 3: (B, omega)
        if(fb){print(paste0("Generating trace plot ", selector, "..."))}
        lattice_dim <- 100

        ## Generate lattice
        if(selector == 1) {
          lattice <- as.matrix(expand.grid(m = seq(min(lower[1], min(opt[1,])), max(upper[1], max(opt[1,])), length.out = lattice_dim), omega = seq(min(lower[2], min(opt[2,])), max(upper[2], max(opt[2,])), length.out = lattice_dim)))
          
          lattice_vals <- apply(lattice, 1, SSE2, tc = tc_val,
                                log_p = log_p,
                                t = t)
        } 
        if(selector == 2) {
          lattice <- as.matrix(expand.grid(m = seq(min(lower[1], min(opt[1,])), max(upper[1], max(opt[1,])), length.out = lattice_dim), B = seq(min(lower[3], min(opt[3,])), max(upper[3], max(opt[3,])), length.out = lattice_dim)))
          
          beta_tp <- beta_calculator(log_p, t, tc_val, opt[1, ][1], opt[2, ][1])

          SSE_tp <- function(par, a, c1, c2, tc, omega, log_p, t) {SSE(par = list(A = a, B = par[2], C1 = c1, C2 = c2, tc = tc, m = par[1], omega = omega), t = t, log_p = log_p)}
          
          if(tp_id < 1 || tp_id > ncol(opt)) {
            tp_id <- 1
            warning("tp_id for (B,m) trace plot is <1 or too high. Replaced with tp_id[1]=1\n")
          }
          lattice_vals <- apply(lattice, 1, SSE_tp,
            a = beta_tp[1],
            c1 = beta_tp[3],
            c2 = beta_tp[4],
            tc = tc_val,
            omega = opt[2, ][tp_id],
            log_p = log_p,
            t = t)
        } 
        if(selector == 3) {
          lattice <- as.matrix(expand.grid(omega = seq(min(lower[2], min(opt[2,])),max(upper[2], max(opt[2,])), length.out = lattice_dim), B = seq(min(lower[3], min(opt[3,])), max(upper[3], max(opt[3,])), length.out = lattice_dim)))
          
          beta_tp <- beta_calculator(log_p, t, tc_val, opt[1, ][1], opt[2, ][1])
          
          SSE_tp <- function(par, a, c1, c2, tc, m, log_p, t) {SSE(par = list(A = a, B = par[2], C1 = c1, C2 = c2, tc = tc, m = m, omega = par[1]), t = t, log_p = log_p)}
          
          if(tp_id > ncol(opt)) {
            tp_id <- 1
            warning("tp_id for (B,omega) trace plot too high. Replaced with tp_id[2]=1\n")
          }
          lattice_vals <- apply(lattice, 1, SSE_tp,
                                a = beta_tp[1],
                                c1 = beta_tp[3],
                                c2 = beta_tp[4],
                                tc = tc_val,
                                m = opt[1, ][tp_id],
                                log_p = log_p,
                                t = t)
        } 
        
        lattice <- as.data.frame(lattice)
        lattice$vals <- lattice_vals
        
        if(selector == 1) {
          opt_df <- data.frame(m = opt[1,], omega = opt[2,])

          opt_df$step <- 1:nrow(opt_df)
          
          trace_plot <- ggplot(lattice, aes(m, omega)) + 
            geom_raster(aes(fill = vals)) +
            geom_contour(aes(z = vals), col = 'white') +
            geom_path(data = opt_df, col = 'orange') +
            geom_point(data = opt_df, col = 'orange') +
            geom_point(data = head(opt_df, 1), shape = 21, fill = "red", color = "orange", size = 3) +
            geom_point(data = tail(opt_df, 1), shape = 23, fill = "green", color = "orange", size = 3) +
            scale_fill_viridis_c() +
            coord_cartesian(expand = FALSE)
        }
        if(selector == 2) {
          opt_df <- data.frame(m = opt[1,], B = opt[3,])  ## Get fitted values from optim
          
          opt_df$step <- 1:nrow(opt_df)
          
          trace_plot <- ggplot(lattice, aes(m, B)) + 
            geom_raster(aes(fill = vals)) +
            geom_contour(aes(z = vals), col = 'white') +
            geom_path(data = opt_df, col = 'orange') +
            geom_point(data = opt_df, col = 'orange') +
            geom_point(data = head(opt_df, 1), shape = 21, fill = "red", color = "orange", size = 3) +
            geom_point(data = tail(opt_df, 1), shape = 23, fill = "green", color = "orange", size = 3) +
            scale_fill_viridis_c() +
            coord_cartesian(expand = FALSE)
        }
        if(selector == 3) {
          opt_df <- data.frame(omega = opt[2,], B = opt[3,])
          
          opt_df$step <- 1:nrow(opt_df)
          
          trace_plot <- ggplot(lattice, aes(omega, B)) + 
            geom_raster(aes(fill = vals)) +
            geom_contour(aes(z = vals), col = 'white') +
            geom_path(data = opt_df, col = 'orange') +
            geom_point(data = opt_df, col = 'orange') +
            geom_point(data = head(opt_df, 1), shape = 21, fill = "red", color = "orange", size = 3) +
            geom_point(data = tail(opt_df, 1), shape = 23, fill = "green", color = "orange", size = 3) +
            scale_fill_viridis_c() +
            coord_cartesian(expand = FALSE)
        }
        if(fb){print(paste0("Trace plot generated"))}
        trace_plot
      }
    } else {warning(paste0("Trace plot only generated when mode = F2\n"))}
  }
  ## Generate trace plots
  if(tp[1] == 1) {trace_plot_mo <- tp_gen(1)} else {trace_plot_mo <- NULL} 
  if(tp[2] == 1) {trace_plot_bm <- tp_gen(2)} else {trace_plot_bm <- NULL} 
  if(tp[3] == 1) {trace_plot_bo <- tp_gen(3)} else {trace_plot_bo <- NULL} 
  
  
  ## Fit plot ====
  ## Red vertical line: T2
  ## Green vertical line: tc_hat
  if(fp == 1) {
    if(fb){print(paste0("Generating fit plot..."))}
    
    ## Plotting full data set in [T1, T2 + hold_out]
    plot_data <- data.frame(ID = time_ID, log_p = log_price)
    
    fit_plot <- ggplot(plot_data, aes(x = ID, y = log_p)) + 
      geom_line(color = "royalblue1") + 
      geom_function(fun = LPPLS, n = fit[[1]]$tc - 1, args = list(
        A = fit[[1]]$A, 
        B = fit[[1]]$B, 
        C1 = fit[[1]]$C1, 
        C2 = fit[[1]]$C2, 
        tc = fit[[1]]$tc, 
        m = fit[[1]]$m, 
        omega = fit[[1]]$omega
      ), color = "red") +
      geom_vline(aes(xintercept = length(t)), color = "green", linetype = "dashed") +
      geom_vline(aes(xintercept = fit[[1]]$tc), color = "red", linetype = "dashed") +
      labs(x = "tidsindeks", y = "log-pris")
    
    if(fb){print(paste0("Fit plot generated"))}
  }
  
  ## Parameter plot ====
  if(pp == 1){
    #if(mode == "F2") {
    if(fb){print(paste0("Generating parameter plot..."))}
    plot_df_merged <- gather(data = fit[[2]], key = param, value = estimate, -c(ID, tc))
    plot_df_lines <- tibble(
      param = c("m", "omega", "A", "B", "C1", "C2", "D", "value"), 
      estimate = c(
        fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc, fit[[1]]$tc)
      ) %>% group_by(param, estimate)
    
    param_plot <- ggplot(plot_df_merged, aes(tc, estimate)) +
      geom_point(size = 0.5) +
      geom_vline(data = plot_df_lines, mapping = aes(xintercept = estimate), color = "red") +
      facet_wrap(~param, scales = "free_y", ncol = 1)
    # } else {warning(paste0("Paramater plot only generated, when mode = F2\n"))}
    if(fb){print(paste0("Parameter plot generated"))}
  }
  
  ## Matrix plot ====
  if(mp == 1) {
    if(mode == "F2" || mode == "MPL") {
      if(fb){print(paste0("Generating matrix plot data..."))}
      tc_seq <- n(t) + 1:fh
      m_seq <- c(lower[1], (lower[1]+upper[1])/2, upper[1])
      omega_seq <- c(lower[2], (lower[2]+upper[2])/2, upper[2])
      
      tbl <- tibble(tc = 1, m = 1, omega = 1, sse = 1)
      
      row <- 1
      for(i in seq_along(m_seq)) {
        for(j in seq_along(omega_seq)) {
          for(k in seq_along(tc_seq)) {
            tbl[row, "sse"] <- SSE2(c(m_seq[i], omega_seq[j]), tc_seq[k], log_p, t)
            tbl[row, "m"] <- m_seq[i]
            tbl[row, "omega"] <- omega_seq[j]
            tbl[row, "tc"] <- tc_seq[k]
            row <- row + 1
          }
        }
      }
      if(fb){print(paste0("Matrix plot data generated"))}
      if(fb){print(paste0("Generating matrix plot..."))}
      
      matrix_plot <- ggplot(tbl, aes(tc, sse)) +
        geom_line(mapping = aes(x = tc, y = sse)) +
        facet_wrap(omega ~ m, as.table = FALSE, labeller = "label_both")
    } else {warning(sprintf("Matrix plot only generated when mode = F2 og mode = MPL\n"))}
    if(fb){print(paste0("Matrix plot generated"))}
  }
  
  ## Check filtering of preliminary fits
  num_estimates <- ((mode == "F2") + (mode == "MPL")) * fh * num_searches + (mode == "F1") * num_searches
  if(length(out_of_range_tracker$B) > 0) {warning(paste0(length(out_of_range_tracker$B), " of ", num_estimates, " estimates of B were out of range.\n", "See out_of_range_tracker output.\n"))}
  if(length(out_of_range_tracker$D) > 0) {warning(paste0(length(out_of_range_tracker$D), " of ", num_estimates, " estimates of B were out of range.\n", "See out_of_range_tracker output.\n"))}

  ## Check filtering of best fit
  if(beta_vals[[2]] >= upper[3]) {warning(paste0("tc_hat (final fit): B out of range.\n"))}
  if(damp <= lower[4]) {warning(paste0("tc_hat (final fit): D out of range.\n"))}
  
  end_time <- Sys.time()
  print(end_time - start_time)
  list(fit = fit, mpl_output = mpl_output, mpl_plot = mpl_plot, fit_plot = fit_plot, contour_data = contour_data, contour_plot = contour_plot, surface_plot = surface_plot, trace_plot_mo = trace_plot_mo, trace_plot_bm = trace_plot_bm, trace_plot_bo = trace_plot_bo, param_plot = param_plot, matrix_plot = matrix_plot, out_of_range_tracker = out_of_range_tracker)
}

## LPPLS function
LPPLS <- function(t, A, B, C1, C2, tc, m, omega, mode = 0, T1 = 500, T2 = 1990, omega2 = 0, omega3 = 0) {
  if (mode == 0) {
    ## Filimonov
    d <- omega * log(tc - t)
    A + (tc - t)^m * (B + C1 * cos(d) + C2 * sin(d))
  } else if (mode == 1) {
    ## First order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale Invariance), equation (19).
    ## See also Johansen et al. (2000), equation (18).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## C2 instead of phi
    ## m instead of beta
    
    A + (tc - t)^m * (B + C1 * cos(omega * log(tc - t) + C2))
    
  } else if (mode == 2) {
    ## Second order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale Invariance), equation (20).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## A instead of A_2
    ## B instead of B_2
    ## C1 instead of C_2
    ## C2 instead of phi_2
    ## m instead of beta
    ## omega instead of omega_1
    
    tau <- tc - t
    A + tau^m / sqrt(1 + (tau / T1)^(2*m)) * 
      (B + C1 * cos(omega * log(tau) + (omega2 / (2 * m)) * log(1 + (tau / T1)^(2*m)) + C2))
    
  } else if (mode == 3) {
    ## Third order expansion.
    ## See Johansen 1999 (Predicting Financial Crashes Using Discrete Scale Invariance), equation (22).
    ## I am repurposing the variable names from the Filimonov implementation:
    ## A instead of A_3
    ## B instead of B_3
    ## C1 instead of C_3
    ## C2 instead of phi_3
    ## m instead of beta
    ## omega instead of omega_1
    
    tau <- tc - t
    A + (
      tau^m / sqrt(1 + (tau / T1)^(2 * m) + (tau / T2)^(4 * m)) *
        (B + C1 * cos(
              omega + log(tau) + (omega2 / (2 * m)) *
                log(1 + (tau / T1)^(2 * m)) +
                (omega3 / (4 * m)) *
                log(1 + (tau / T2)^(4 * m)) +
                C2)))
    
  } else {stop("LPPLS mode must be 0, 1, 2 or 3.")}
}

## General SSE for custom contour plots
## par list is {tc, m, omega, A, B, C1, C2}
SSE <- function(par, log_p, t, mode = 0) {
  res <- log_p - LPPLS(t, par$A, par$B, par$C1, par$C2, par$tc, par$m, par$omega, mode)
  drop(sum(res^2, na.rm = TRUE))
}

## Generates contour plot for SSE wrt two parameters.
## Provide a par-list of fixed values for 5 of the following 7 parameters:
##      {tc, m, omega, A, B, C1, C2}
## IMPORTANT: Parameters should be in the order above.
##      eg:
##      par = list(tc = 2000, A = 4, B = -0.0015, C1 = -0.01, C2 = 0.1)
## Provide list of names of x and y variabes as character strings.
##      eg:
##      vars = list(x = "m", y = "omega")
## lower and upper are vectors of ranges for x and y axes:
##      lower = c(m_lower, B_lower)
##      upper = c(m_upper, B_upper)
## mode, LPPLS mode. mode = 0 is Filimonov.
contour_plot <- function(log_p, t, par, vars, lower = c(0.0, -0.001), upper = c(1.0, 0.001), cp = 1, sp = 0, fb = 0, point = c(NA,NA), mode = 0) {
  ## Generate SSE data for contour plot of SSE wrt m and B
  if(length(log_p) != length(t)){stop("Length of log_p and t must be same.")}
  if(fb){print(paste0("Generating contour plot data..."))}
  contour_plot <- NULL
  surface_plot <- NULL
  x_contour = seq(lower[1], upper[1], length.out = 101) ## m. lower and upper are filter range.
  y_contour =  seq(lower[2], upper[2], length.out = 101) ## B
  sse_contour <- function(x, y) {
    ## Make a list for x and y
    xy_list <- list(x, y)
    names(xy_list) <- vars[1:2]
    SSE(par = c(par, xy_list), log_p = log_p, t = t, mode = mode)
  }
  z_contour <- t(outer(x_contour, y_contour, Vectorize(sse_contour))) ## z matrix must be transposed for correct contour plot! (Don't know why)
  contour_data <- list(x = x_contour, y = y_contour, z = z_contour)
  if(fb){print(paste0("Fit contour data generated"))}
  
  ## Contour plot
  if(cp == 1) {
    if(fb){print(paste0("Generating contour plot..."))}
    contour_plot <- plot_ly(
      x = contour_data$x, 
      y = contour_data$y, 
      z = contour_data$z,
      type = "contour",
      contours = list(
        showlabels = FALSE,
        exponentformat = "e"),
      ncontours = 10,
      colorbar = list(exponentformat = "e", title = "SSE")
    )
    contour_plot <- contour_plot %>% layout(
      title = paste0("SSE wrt. ", vars[[1]], " and ", vars[[2]]),
      xaxis = list(title = vars[[1]]),
      yaxis = list(title = vars[[2]])
    )
    if(!is.na(sum(point))) {
      #df_point <- data.frame(x = point[1], y = point[2])
      z_coord <- SSE(par = c(par, list(m = point[1], omega = point[2])), log_p = log_p, t = t)
      contour_plot <- contour_plot %>% 
        add_trace(x = point[1], y = point[2], type = "scatter", mode = "markers", showlegend = FALSE) %>%
        add_trace(x = c(point[1], point[1], lower[1]), y = c(lower[2], point[2], point[2]), type = "scatter", mode = "lines", line = list(dash = "dash", width = 1, color = "orange"), showlegend = FALSE) %>%
        add_annotations(
          xref = "x",
          yref = "y",
          x = point[1],
          y = point[2],
          xanchor = "middle",
          yanchor = "bottom",
          text = paste0("(", point[1], ", ", point[2], ", ", z_coord,")"),
          font = list(color = 'white', size = 10),
          showarrow = FALSE
        )
    }
    if(fb){print(paste0("Contour plot generated"))}
    warning("(Don't worry about all these scatter warnings...)")
  }
  
  ## Surface plot
  if(sp == 1) {
    if(fb){print(paste0("Generating surface plot..."))}
    surface_plot <- plot_ly(
      x = contour_data$x, 
      y = contour_data$y, 
      z = contour_data$z
    ) %>% add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      ),
      colorbar = list(exponentformat = "e", title = "SSE")
    )
    surface_plot <- surface_plot %>% layout(
      title = paste0("SSE wrt. ", vars[[1]], " and ", vars[[2]]),
      scene = list(
        camera=list(
          eye = list(x=1.1, y=-1.3, z=0.1)
        ),
        xaxis = list(title = vars[[1]]),
        yaxis = list(title = vars[[2]]),
        zaxis = list(title = "SSE")
      )
    )
    if(fb){print(paste0("Surface plot generated"))}
  }
  list(contour_plot, surface_plot)
}

## Fit plot ====
## Input:
## 1) fit, list of parametervalues for A, B, C1, C2, tc, m and omega.
## 2) time_ID, vector of time IDs.
## 3) log_price, vector of log-prices.
## 4) mode, LPPLS mode. mode = 0 is Filimonov.
fit_plot <- function(fit, time_ID, log_price, mode) {
  ## Plotting full data set in [T1, T2 + hold_out]
  plot_data <- data.frame(ID = time_ID, log_p = log_price)
  
  fit_plot <- ggplot(plot_data, aes(x = ID, y = log_p)) + 
    geom_line(size = 0.1, color = "royalblue1") + 
    geom_function(fun = LPPLS, args = list(
      A = fit$A, 
      B = fit$B, 
      C1 = fit$C1, 
      C2 = fit$C2, 
      tc = fit$tc, 
      m = fit$m, 
      omega = fit$omega,
      mode = mode
    ), color = "red") +
    labs(x = "tidsindeks", y = "log-pris")
  fit_plot
}












