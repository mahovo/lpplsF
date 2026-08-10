// Prototypes for the Rcpp question. Plain Rcpp only (LinkingTo: Rcpp), no
// Armadillo, so the measured numbers reflect the lightest dependency option.
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Gauss-Jordan with partial pivoting on a small dense system.
static bool solve_small(double *A, double *b, int p) {
  for (int col = 0; col < p; ++col) {
    int piv = col; double best = std::fabs(A[col * p + col]);
    for (int r = col + 1; r < p; ++r) {
      double v = std::fabs(A[r * p + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < 1e-300) return false;
    if (piv != col) {
      for (int c = 0; c < p; ++c) std::swap(A[col * p + c], A[piv * p + c]);
      std::swap(b[col], b[piv]);
    }
    double d = A[col * p + col];
    for (int c = 0; c < p; ++c) A[col * p + c] /= d;
    b[col] /= d;
    for (int r = 0; r < p; ++r) {
      if (r == col) continue;
      double f = A[r * p + col];
      if (f == 0.0) continue;
      for (int c = 0; c < p; ++c) A[r * p + c] -= f * A[col * p + c];
      b[r] -= f * b[col];
    }
  }
  return true;
}

// Fused: one pass accumulating X'X and X'y directly, so the n x 4 design
// matrix and every length-n temporary are never materialised.
// [[Rcpp::export]]
NumericVector beta_cpp(NumericVector log_p, NumericVector t,
                       double tc, double m, double omega) {
  const int n = t.size();
  double sf = 0, sg = 0, sh = 0, sff = 0, sgg = 0, shh = 0;
  double sfg = 0, sfh = 0, sgh = 0, sy = 0, syf = 0, syg = 0, syh = 0;
  for (int i = 0; i < n; ++i) {
    const double tau  = tc - t[i];
    const double taum = std::pow(tau, m);
    const double lt   = omega * std::log(tau);
    const double f = taum, g = taum * std::cos(lt), h = taum * std::sin(lt);
    const double y = log_p[i];
    sf += f; sg += g; sh += h;
    sff += f * f; sgg += g * g; shh += h * h;
    sfg += f * g; sfh += f * h; sgh += g * h;
    sy += y; syf += y * f; syg += y * g; syh += y * h;
  }
  double A[16] = { (double)n, sf,  sg,  sh,
                   sf,        sff, sfg, sfh,
                   sg,        sfg, sgg, sgh,
                   sh,        sfh, sgh, shh };
  double b[4] = { sy, syf, syg, syh };
  if (!solve_small(A, b, 4)) stop("singular normal equations");
  return NumericVector::create(_["A"] = b[0], _["B"] = b[1],
                               _["C1"] = b[2], _["C2"] = b[3]);
}

// Gradient matrix, Filimonov (2017) eq. (36). Replaces an explicit R loop.
// [[Rcpp::export]]
NumericMatrix compute_X_cpp(NumericVector Psi, double tc, NumericVector t) {
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

// Residual-weighted second derivatives, Filimonov (2017) eq. (37). Symmetric.
// [[Rcpp::export]]
NumericMatrix compute_H_cpp(NumericVector Psi, double tc,
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
