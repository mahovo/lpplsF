# lpplsF
Filimonov-variant of the LPPLS model


Implementation of the LPPLS model for analysing financial bubbles.  
By Martin Hoshi Vognsen.

This implementation is based on the following articles:  

1) A stable and robust calibration scheme of the log-periodic power law model  
Filimonov, V., Sornette, D.  
Physica A: Statistical Mechanics and its Applications, 2013  

2) Modified profile likelihood inference and interval forecast of the burst of financial bubbles  
Filimonov, V., Demos, G., Sornette, D.
Quantitative Finance, 2017  


Calculate complete model based on input time ID and price vectors.  

Inputs:  
1) time_ID, time index  
        vector, T1:T2 where T1 = 1  
2) log_price, price  
        vector, must be same length as t  
3) fh, length of forecast horizon. fh = s means forecast period is [T2+1, T2+s], s>1.  
     Unit is indexes.  
     Default is 90 time units.  
4) hold_out, number of time units for holdout data.   
        That is, the number of data points after T2.  
5) lower, lower limits of parameter filter  
     {m, omega, B, D}  
     Filimonov 2017: B < 0, D >= 0.8, 0.1 <= m <= 0.9, 6 <= omega <= 13   
     Damping, D = (m * abs(B) / (omega * sqrt(C1^2 + C2^2))  
     Defaults: lower = c(0.1, 6, -1e14, 0.8), upper = c(0.9, 13, -1e-14, 1e6)  
     NOTE: For trace_plot mode 2 and 3 (bm and bo), lower for B should be set to  
           the desired minimum value of the B-axis in the contour plot.  
           Example:  
           lower = c(0.1, 6, -0.03, 0.8)  
           If the optimization algorithm goes below the lower limits at any point  
           the lower limit for the B-axis will be extended automatically.  
6) upper, upper limits of parameter filter  
7) tc_init, initial tc value for optim().  
        tc_init should be in [T2+1, T2+s]  
8) m_init, initial m value for optim()  
        m_init should be in [lower[1], upper[1]]  
        The algorithm will search for initial values of m and omega that give B<0.  
        If such values are found, the given init values will be overwritten.  
9) o_init, initial omega value for optim()  
       Omega_init should be in [lower[2], upper[2]]  
       Default is 13, which increases the probability of B<0  
       (of course not for random initial values).  
       The algorithm will search for initial values of m and omega that give B<0.  
       If such values are found, the given init values will be overwritten.  
10) num_searches, number of times to repeat the optimization wrt tc, m and omega.  
         First iteration uses given initial values.  
         Any subsequent iterations use random initial   values.
11) mode (see also under output below),   
         mode = "F1", optimize wrt tc, m and omega simultaneously  
         mode = "F2", optimize wrt m and omega simultaneously, then wrt tc  
         mode = "MPL", Modified Profile Likelihood with likelihood intervals and plot.  
12) mpl_cutoff, cutoff level (c) for likelihood interval. See Filimonov2017 equation (39).  
              vector with 3 elements, for 3 levels.  
              Default is c = c(0.05, 0.1, 0.5).  
              c = 0.05 indicates that values of tc outside the likelihood interval  
              have a probability of 0.05.  
              Note that the intervals for level c = 0.05 can be extremely wide.  
12) mpl_plot, boolean for Modified Profile Likelihood plot with likelihood intervals.  
13) cp, boolean for contour plot.  
         tc value is fixed to the estimate with the lowest objective function value.  
14) sp, boolean for contour plot.  
         tc value is fixed to the estimate with the lowest objective function value.  
15) tp, boolean for vector trace plot selector.  
         Any combination of three plots are selected with a vector:  
         c(selector 1 on/off, selector 2 on/off, selector 3 on/off).  
         selector 1: m, omega  
         selector 2: B, m  
         selector 3: B, omega  
         Example  
         c(1, 0, 1) selects trace plot 1 and 3.  
         Traces convergence of parameters during optimization of F2 objective function.  
         Only works when mode = "F2".  
         The plotted trace path is the one for the value of tc that produced the  
         smallest value of SSE2.  
         NOTE: Contour plots in trace_plot 2 and 3 only apply to the first step in the trace.  
               In a trace plot wrt B and m, omega is not fixed. (And vice versa for m and omega.)  
16) tp_id, integer  
         Indicate which step of the trace to use for contour plot.  
         If the indicated id is too high, 1 will be chosen.  
17) fp, boolean for fit plot.  
18) mp, boolean for matrix plot.  
             x-axis is tc, y-axis is sse.  
             3 columns of plots for smallest and biggest value of m, and value between the two  
             3 columns of plots for smallest and biggest value of omega, and value between the two  
19) factr, factr argument for optim control list when using L-BFGS-B. See optim documentation.  
20) fb, boolean for printing feedback during execution.  



Outputs:  
1)  
mode = F1:  
   fit, list of fits with random initial values.  
        fit[[1]]:  
             list of best fit coefficients optimized wrt m and omega.  
        fit[[2]]:  
             Tibble of fits with random starting parameters.  
             Sorted by objective function value (best at top).  
             Output format:  
             [  
                  list(tc, m, omega, A, B, C1, C2, D, value_min, ID),  
                  list(tc, m, omega, A, B, C1, C2, D, value, ID),  
                  ...,  
                  list(tc, m, omega, A, B, C1, C2, D, value_max, ID)  
             ]  
        fit[[3]]:  
             If no fits passed the filter: Only unfiltered fits returned.  
   ID is an identifier indicating the order pre-sorting.  
   These are the indexes reported in "out of range" warnings.  
mode = F2:  
   fit, list with two or three elements:  
        fit[[1]]:  
             list of best fit coefficients optimized wrt m and omega.  
        fit[[2]]:  
             If at least one fit passed the filter:  
                  Tibble of best fits for objective function F2 for each value of tc.  
                  Columns are: {ID, value, tc, m, omega, A, B, C1, C2, damp}  
                  Sorted by objective function value (best at top).  
                  For each value of tc the best fit is picked from fit2_tmp.  
             If no fits passed the filter:   
                  Only unfiltered fits returned.  
        fit[[3]]:  
             If at least one fit passed the filter:  
                  Unfiltered fits returned.  
             If no fits passed the filter:   
                  Nothing is returned in fit[[3]]  
        fit[[4]]:  
             List of all tmp fits  
             For each tc value, there is a tibble with all fits using random starting points.  
2) mpl_output  
        mpl = list(LI, R, LL, MLL)  
             LI, vector of 3 elements, one for each value of mpl_cutoff.  
             R, number, relative likelihood.  
             LL, vector of MPL log-likelihoods for each tc value.  
             MLL, maximum of LL.  
             tc_hat_mpl, tc for which LL takes it's maximum.  
3) mpl_plot, Modified Profile Likelihood plot with likelihood intervals.  
4) contour_data, list of x, y and z data for contour plot of SSE wrt m and omega.  
5) fit_plot, object containing fit plot  
6) contour_plot, object containing contour plot  
7) surface_plot, object containing surface plot  
8) trace_plot_mo, object containing trace plot for m and omega  
9) trace_plot_bm, object containing trace plot for B and m  
10) trace_plot_bo, object containing trace plot for B and omega  
11) matrix_plot, object containing matrix plot  
12) out_of_range_tracker  
        F1 mode: list of random iteration IDs what resulted in out-of-range parameters.  
        F2 mode: list of random tc value IDs and iteration IDs what resulted in   
                 out-of-range parameters.  
                 tc value ID start with 1 for T2+1. That is, the first time step  
                 after the modelling period, or the first time step of the  
                 forecasting horizon.  
