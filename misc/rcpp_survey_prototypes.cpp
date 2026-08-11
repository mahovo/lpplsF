// Prototypes for the package-wide Rcpp survey (2026-08-10).
//
// Extends misc/bench_beta_prototypes.cpp, which covered the beta solver and the
// two MPL matrices, with the routines the survey identified as the remaining
// hot spots: the SSE2 objective that optim() drives, the grid scans that
// create_contour_data() / create_matrix_plot() / create_trace_plot() run over
// it, and the AR(1)-GARCH(1,1) path simulator in lr_trend_test().
//
// Plain Rcpp only (LinkingTo: Rcpp), no RcppArmadillo, with the same
// hand-rolled Gauss-Jordan solve, so every measurement here reflects the
// lightest possible dependency.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Shared: 4x4 Gauss-Jordan with partial pivoting (as in bench_beta_prototypes)
// ---------------------------------------------------------------------------
static bool solve4(double *A, double *b) {
  const int p = 4;
  for (int col = 0; col < p; ++col) {
    int piv = col; double best = std::fabs(A[col * p + col]);
    for (int r = col + 1; r < p; ++r) {
      double v = std::fabs(A[r * p + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (!(best > 1e-300)) return false;
    if (piv != col) {
      for (int c = 0; c < p; ++c) std::swap(A[col * p + c], A[piv * p + c]);
      std::swap(b[col], b[piv]);
    }
    const double d = A[col * p + col];
    for (int c = 0; c < p; ++c) A[col * p + c] /= d;
    b[col] /= d;
    for (int r = 0; r < p; ++r) {
      if (r == col) continue;
      const double f = A[r * p + col];
      if (f == 0.0) continue;
      for (int c = 0; c < p; ++c) A[r * p + c] -= f * A[col * p + c];
      b[r] -= f * b[col];
    }
  }
  return true;
}

// Accumulate the normal equations for the LPPLS design
// X = [1, tau^m, tau^m cos(w log tau), tau^m sin(w log tau)], tau = tc - t.
// Optionally writes the three non-constant basis columns into `basis`
// (layout f[0..n), g[0..n), h[0..n)) so a second pass can form residuals
// without repeating pow/log/cos/sin.
static bool normal_eq(const double *y, const double *t, int n,
                      double tc, double m, double omega,
                      double *A, double *b, double *basis) {
  double sf = 0, sg = 0, sh = 0, sff = 0, sgg = 0, shh = 0;
  double sfg = 0, sfh = 0, sgh = 0, sy = 0, syf = 0, syg = 0, syh = 0;
  for (int i = 0; i < n; ++i) {
    const double tau  = tc - t[i];
    const double taum = std::pow(tau, m);
    const double lt   = omega * std::log(tau);
    const double f = taum, g = taum * std::cos(lt), h = taum * std::sin(lt);
    const double yi = y[i];
    sf += f; sg += g; sh += h;
    sff += f * f; sgg += g * g; shh += h * h;
    sfg += f * g; sfh += f * h; sgh += g * h;
    sy += yi; syf += yi * f; syg += yi * g; syh += yi * h;
    if (basis) { basis[i] = f; basis[n + i] = g; basis[2 * n + i] = h; }
  }
  A[0] = (double)n; A[1] = sf;  A[2] = sg;  A[3] = sh;
  A[4] = sf;        A[5] = sff; A[6] = sfg; A[7] = sfh;
  A[8] = sg;        A[9] = sfg; A[10] = sgg; A[11] = sgh;
  A[12] = sh;       A[13] = sfh; A[14] = sgh; A[15] = shh;
  b[0] = sy; b[1] = syf; b[2] = syg; b[3] = syh;
  return solve4(A, b);
}

// Scratch buffer for the basis, grown once and reused across calls so the
// two-pass form costs no per-call allocation.
static std::vector<double> g_basis;
static double *basis_buf(int n) {
  if ((int)g_basis.size() < 3 * n) g_basis.resize(3 * n);
  return g_basis.data();
}

// ---------------------------------------------------------------------------
// The SSE2 objective, fused. This is what optim() calls fh x num_searches x
// (iterations) times per fit; in pure R it is beta_calculator() + eval_lppls()
// + sum((log_p - fitted)^2), which walks the series three times and rebuilds
// the same tau^m / cos / sin basis twice.
// ---------------------------------------------------------------------------

// Variant A: two passes, basis buffered. Residuals formed explicitly, so the
// arithmetic matches the R version term for term.
// [[Rcpp::export]]
double sse2_cpp(NumericVector par, double tc, NumericVector log_p,
                NumericVector t) {
  const int n = t.size();
  const double m = par[0], omega = par[1];
  double A[16], b[4];
  double *bs = basis_buf(n);
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, omega, A, b, bs))
    return R_PosInf;
  const double *f = bs, *g = bs + n, *h = bs + 2 * n;
  double sse = 0.0;
  for (int i = 0; i < n; ++i) {
    const double fit = b[0] + b[1] * f[i] + b[2] * g[i] + b[3] * h[i];
    const double r = log_p[i] - fit;
    if (R_finite(r)) sse += r * r;      // mirrors na.rm = TRUE
  }
  return sse;
}

// Variant B: two passes, basis recomputed (no buffer). Included to measure
// what the buffer is worth.
// [[Rcpp::export]]
double sse2_cpp_nobuf(NumericVector par, double tc, NumericVector log_p,
                      NumericVector t) {
  const int n = t.size();
  const double m = par[0], omega = par[1];
  double A[16], b[4];
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, omega, A, b, NULL))
    return R_PosInf;
  double sse = 0.0;
  for (int i = 0; i < n; ++i) {
    const double tau  = tc - t[i];
    const double taum = std::pow(tau, m);
    const double lt   = omega * std::log(tau);
    const double fit = b[0] + taum * (b[1] + b[2] * std::cos(lt) +
                                              b[3] * std::sin(lt));
    const double r = log_p[i] - fit;
    if (R_finite(r)) sse += r * r;
  }
  return sse;
}

// Variant C: single pass, SSE from the normal equations alone via
// SSE = y'y - beta'(X'y). Algebraically exact, but it differences two large
// nearly-equal quantities, so it loses relative precision. Measured for
// completeness; NOT proposed for use.
// [[Rcpp::export]]
double sse2_cpp_alg(NumericVector par, double tc, NumericVector log_p,
                    NumericVector t) {
  const int n = t.size();
  const double m = par[0], omega = par[1];
  double A[16], b[4], rhs[4];
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, omega, A, b, NULL))
    return R_PosInf;
  // normal_eq overwrites b with the solution, so rebuild X'y for the dot
  // product. Cheap relative to the pass, and keeps the variant self-contained.
  double sy = 0, syf = 0, syg = 0, syh = 0, yy = 0;
  for (int i = 0; i < n; ++i) {
    const double tau  = tc - t[i];
    const double taum = std::pow(tau, m);
    const double lt   = omega * std::log(tau);
    const double yi = log_p[i];
    sy += yi; syf += yi * taum;
    syg += yi * taum * std::cos(lt); syh += yi * taum * std::sin(lt);
    yy += yi * yi;
  }
  rhs[0] = sy; rhs[1] = syf; rhs[2] = syg; rhs[3] = syh;
  return yy - (b[0] * rhs[0] + b[1] * rhs[1] + b[2] * rhs[2] + b[3] * rhs[3]);
}

// The beta solver on its own, returning the coefficients (same contract as
// beta_cpp in bench_beta_prototypes.cpp; repeated here so this file stands
// alone and so a fit can use one .Call for the objective and one for beta).
// [[Rcpp::export]]
NumericVector beta_cpp2(NumericVector log_p, NumericVector t,
                        double tc, double m, double omega) {
  const int n = t.size();
  double A[16], b[4];
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, omega, A, b, NULL))
    stop("singular normal equations");
  return NumericVector::create(_["A"] = b[0], _["B"] = b[1],
                               _["C1"] = b[2], _["C2"] = b[3]);
}

// Objective and coefficients together: what fit_lppls needs right after the
// optim() call, without paying for a second pass over the series.
// [[Rcpp::export]]
NumericVector sse2_beta_cpp(NumericVector par, double tc, NumericVector log_p,
                            NumericVector t) {
  const int n = t.size();
  const double m = par[0], omega = par[1];
  double A[16], b[4];
  double *bs = basis_buf(n);
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, omega, A, b, bs))
    return NumericVector::create(_["value"] = R_PosInf, _["A"] = NA_REAL,
                                 _["B"] = NA_REAL, _["C1"] = NA_REAL,
                                 _["C2"] = NA_REAL);
  const double *f = bs, *g = bs + n, *h = bs + 2 * n;
  double sse = 0.0;
  for (int i = 0; i < n; ++i) {
    const double r = log_p[i] - (b[0] + b[1] * f[i] + b[2] * g[i] + b[3] * h[i]);
    if (R_finite(r)) sse += r * r;
  }
  return NumericVector::create(_["value"] = sse, _["A"] = b[0], _["B"] = b[1],
                               _["C1"] = b[2], _["C2"] = b[3]);
}

// ---------------------------------------------------------------------------
// Grid scan: the whole (m, omega) lattice in one .Call. This is what
// create_contour_data() (101x101), create_matrix_plot() (fh x 3 x 3) and
// create_trace_plot() (100x100) each drive one R-level cell at a time.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix sse2_grid_cpp(NumericVector m_seq, NumericVector omega_seq,
                            double tc, NumericVector log_p, NumericVector t) {
  const int nm = m_seq.size(), nw = omega_seq.size(), n = t.size();
  NumericMatrix out(nm, nw);
  double *bs = basis_buf(n);
  const double *y = log_p.begin(), *tt = t.begin();
  for (int a = 0; a < nm; ++a) {
    for (int c = 0; c < nw; ++c) {
      double A[16], b[4];
      if (!normal_eq(y, tt, n, tc, m_seq[a], omega_seq[c], A, b, bs)) {
        out(a, c) = R_PosInf; continue;
      }
      const double *f = bs, *g = bs + n, *h = bs + 2 * n;
      double sse = 0.0;
      for (int i = 0; i < n; ++i) {
        const double r = y[i] - (b[0] + b[1] * f[i] + b[2] * g[i] + b[3] * h[i]);
        if (R_finite(r)) sse += r * r;
      }
      out(a, c) = sse;
    }
  }
  return out;
}

// Profile SSE over a tc sequence at fixed (m, omega) -- the shape
// create_matrix_plot() actually scans (fh values of tc for each of 9 (m,omega)
// corners).
// [[Rcpp::export]]
NumericVector sse2_tc_seq_cpp(NumericVector tc_seq, double m, double omega,
                              NumericVector log_p, NumericVector t) {
  const int nt = tc_seq.size(), n = t.size();
  NumericVector out(nt);
  double *bs = basis_buf(n);
  const double *y = log_p.begin(), *tt = t.begin();
  for (int k = 0; k < nt; ++k) {
    double A[16], b[4];
    if (!normal_eq(y, tt, n, tc_seq[k], m, omega, A, b, bs)) {
      out[k] = R_PosInf; continue;
    }
    const double *f = bs, *g = bs + n, *h = bs + 2 * n;
    double sse = 0.0;
    for (int i = 0; i < n; ++i) {
      const double r = y[i] - (b[0] + b[1] * f[i] + b[2] * g[i] + b[3] * h[i]);
      if (R_finite(r)) sse += r * r;
    }
    out[k] = sse;
  }
  return out;
}

// ---------------------------------------------------------------------------
// eval_lppls(), mode 0 only. Mostly here so the fitted-value path can be
// compared like for like; the R version is already vectorised.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector eval_lppls0_cpp(NumericVector t, double A, double B, double C1,
                              double C2, double tc, double m, double omega) {
  const int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double tau = tc - t[i];
    const double d   = omega * std::log(tau);
    out[i] = A + std::pow(tau, m) * (B + C1 * std::cos(d) + C2 * std::sin(d));
  }
  return out;
}

// ---------------------------------------------------------------------------
// The MPL matrices, unchanged from bench_beta_prototypes.cpp (renamed so both
// files can be sourced into one session).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix compute_X_cpp2(NumericVector Psi, double tc, NumericVector t) {
  const int n = t.size();
  const double m = Psi[0], omega = Psi[1], b = Psi[3], c1 = Psi[4], c2 = Psi[5];
  NumericMatrix X(n, 6);
  for (int i = 0; i < n; ++i) {
    const double tau = tc - t[i];
    const double lt  = std::log(tau);
    const double tm  = std::pow(tau, m);
    const double ct  = std::cos(omega * lt), st = std::sin(omega * lt);
    X(i, 0) = tm * lt * (b + c1 * ct + c2 * st);
    X(i, 1) = tm * lt * (-c1 * st + c2 * ct);
    X(i, 2) = 1.0;
    X(i, 3) = tm;
    X(i, 4) = tm * ct;
    X(i, 5) = tm * st;
  }
  return X;
}

// [[Rcpp::export]]
NumericMatrix compute_H_cpp2(NumericVector Psi, double tc,
                             NumericVector log_p, NumericVector t) {
  const int n = t.size();
  const double m = Psi[0], omega = Psi[1], a = Psi[2], b = Psi[3],
               c1 = Psi[4], c2 = Psi[5];
  NumericMatrix H(6, 6);
  double h11 = 0, h12 = 0, h14 = 0, h15 = 0, h16 = 0, h22 = 0, h25 = 0, h26 = 0;
  for (int i = 0; i < n; ++i) {
    const double tau = tc - t[i];
    const double lt  = std::log(tau);
    const double tm  = std::pow(tau, m);
    const double ct  = std::cos(omega * lt), st = std::sin(omega * lt);
    const double res = log_p[i] - (a + tm * (b + c1 * ct + c2 * st));
    const double lt2 = lt * lt;
    h11 += res * tm * lt2 * (b + c1 * ct + c2 * st);
    h12 += res * tm * lt2 * (-c1 * st + c2 * ct);
    h14 += res * tm * lt;
    h15 += res * tm * lt * ct;
    h16 += res * tm * lt * st;
    h22 += res * (-1.0) * tm * lt2 * (c1 * ct + c2 * st);
    h25 += res * (-1.0) * tm * lt * st;
    h26 += res * tm * lt * ct;
  }
  H(0, 0) = h11;
  H(0, 1) = h12; H(1, 0) = h12;
  H(0, 3) = h14; H(3, 0) = h14;
  H(0, 4) = h15; H(4, 0) = h15;
  H(0, 5) = h16; H(5, 0) = h16;
  H(1, 1) = h22;
  H(1, 4) = h25; H(4, 1) = h25;
  H(1, 5) = h26; H(5, 1) = h26;
  return H;
}

// ---------------------------------------------------------------------------
// Fused MPL moments. compute_mpl_loglik() only ever uses the n x 6 gradient
// matrices through crossprod(), so the matrices themselves need never be
// materialised: one pass accumulates X'X, X_hat'X_tc and H directly. Returns
// the three 6x6 blocks plus s_tc and leaves det() to R, so the determinants
// stay on R's LAPACK path and the guard (det1 <= 0 -> NA) is unchanged.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List mpl_moments_cpp(double tc, double tc_hat, NumericVector Psi_hat_tc,
                     NumericVector Psi_hat, NumericVector log_p,
                     NumericVector t) {
  const int n = t.size(), p = 6;
  const double m = Psi_hat_tc[0], w = Psi_hat_tc[1], a = Psi_hat_tc[2],
               b = Psi_hat_tc[3], c1 = Psi_hat_tc[4], c2 = Psi_hat_tc[5];
  const double mh = Psi_hat[0], wh = Psi_hat[1], bh = Psi_hat[3],
               c1h = Psi_hat[4], c2h = Psi_hat[5];

  // beta at (tc, m, w) for the fitted values -- as compute_mpl_loglik() does
  double A4[16], b4[4];
  if (!normal_eq(log_p.begin(), t.begin(), n, tc, m, w, A4, b4, NULL))
    return List::create(_["ok"] = false);

  NumericMatrix XtX(p, p), Xc(p, p), H(p, p);
  double sse = 0.0;
  double xt[6], xh[6];
  for (int i = 0; i < n; ++i) {
    const double tau = tc - t[i];
    const double lt  = std::log(tau);
    const double tm  = std::pow(tau, m);
    const double ct  = std::cos(w * lt), st = std::sin(w * lt);

    // fitted value from the freshly solved beta
    const double fit = b4[0] + tm * (b4[1] + b4[2] * ct + b4[3] * st);
    const double rf  = log_p[i] - fit;
    if (R_finite(rf)) sse += rf * rf;

    // gradient row at (tc, Psi_hat_tc), Filimonov (2017) eq. (36)
    xt[0] = tm * lt * (b + c1 * ct + c2 * st);
    xt[1] = tm * lt * (-c1 * st + c2 * ct);
    xt[2] = 1.0; xt[3] = tm; xt[4] = tm * ct; xt[5] = tm * st;

    // gradient row at (tc_hat, Psi_hat)
    const double tauh = tc_hat - t[i];
    const double lth  = std::log(tauh);
    const double tmh  = std::pow(tauh, mh);
    const double cth  = std::cos(wh * lth), sth = std::sin(wh * lth);
    xh[0] = tmh * lth * (bh + c1h * cth + c2h * sth);
    xh[1] = tmh * lth * (-c1h * sth + c2h * cth);
    xh[2] = 1.0; xh[3] = tmh; xh[4] = tmh * cth; xh[5] = tmh * sth;

    for (int r = 0; r < p; ++r) {
      const double xr = xt[r], hr = xh[r];
      for (int c = 0; c < p; ++c) {
        XtX(r, c) += xr * xt[c];      // crossprod(X_hat_tc, X_hat_tc)
        Xc(r, c)  += hr * xt[c];      // crossprod(X_hat,    X_hat_tc)
      }
    }

    // residual-weighted second derivatives, eq. (37)
    const double res = log_p[i] - (a + tm * (b + c1 * ct + c2 * st));
    const double lt2 = lt * lt;
    H(0, 0) += res * tm * lt2 * (b + c1 * ct + c2 * st);
    H(0, 1) += res * tm * lt2 * (-c1 * st + c2 * ct);
    H(0, 3) += res * tm * lt;
    H(0, 4) += res * tm * lt * ct;
    H(0, 5) += res * tm * lt * st;
    H(1, 1) += res * (-1.0) * tm * lt2 * (c1 * ct + c2 * st);
    H(1, 4) += res * (-1.0) * tm * lt * st;
    H(1, 5) += res * tm * lt * ct;
  }
  H(1, 0) = H(0, 1);
  H(3, 0) = H(0, 3);
  H(4, 0) = H(0, 4);
  H(5, 0) = H(0, 5);
  H(4, 1) = H(1, 4);
  H(5, 1) = H(1, 5);

  return List::create(_["ok"] = true, _["XtX"] = XtX, _["Xcross"] = Xc,
                      _["H"] = H, _["s_tc"] = sse / (double)n);
}

// ---------------------------------------------------------------------------
// .sim_ar_garch(): two serial recursions over n + burn, so no R vectorisation
// is available. Uses R's RNG so a given seed reproduces the R version exactly.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector sim_ar_garch_cpp(int n, double ar1, double om0, double alpha1,
                               double beta1, int burn = 500) {
  const int N = n + burn;
  NumericVector z = rnorm(N);          // same stream as stats::rnorm(n + burn)
  std::vector<double> s2(N), ee(N), y(N, 0.0);
  double denom = 1.0 - alpha1 - beta1;
  if (denom < 1e-6) denom = 1e-6;
  s2[0] = om0 / denom;
  ee[0] = std::sqrt(s2[0]) * z[0];
  for (int i = 1; i < N; ++i) {
    s2[i] = om0 + alpha1 * ee[i - 1] * ee[i - 1] + beta1 * s2[i - 1];
    ee[i] = std::sqrt(s2[i] > 0 ? s2[i] : 0.0) * z[i];
  }
  for (int i = 1; i < N; ++i) y[i] = ar1 * y[i - 1] + ee[i];
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = y[burn + i];
  return out;
}
