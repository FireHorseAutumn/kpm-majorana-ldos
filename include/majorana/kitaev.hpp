#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace majorana {

/// CSR sparse matrix, FP64. Same layout as kpm-tightbinding's
/// kpm::CsrMatrix; the integration decision is DESIGN.md section 6.
struct CsrMatrix {
  int n = 0;                  // matrix dimension
  std::vector<int> row_ptr;   // size n + 1
  std::vector<int> col_idx;   // size nnz
  std::vector<double> values; // size nnz

  int nnz() const { return static_cast<int>(values.size()); }
};

/// Kitaev chain: L open sites, hopping t, imposed p-wave pairing delta,
/// chemical potential mu. The BdG Hamiltonian is 2L x 2L; row i is the
/// electron component of site i, row L+i the hole component. Defaults
/// sit at the sweet spot (delta = t, mu = 0).
struct KitaevParams {
  int L = 64;
  double t = 1.0;
  double delta = 1.0;
  double mu = 0.0;
};

/// Local-moment parameters. a = 0 means unset: the caller must supply a
/// spectral bound — sqrt((2t + |mu|)^2 + 4*delta^2) covers the chain.
struct LdosParams {
  int n_cheb = 2048;
  double a = 0.0;
  double b = 0.0;
};

/// R1: build the chain as CSR. BdG blocks [[h, Delta], [-Delta, -h]]:
/// h tridiagonal with -mu on-site and -t hopping, Delta antisymmetric
/// p-wave with Delta_{i,i+1} = -delta (the pairing sign is a gauge
/// choice; the spectrum and end weights carry |delta| only).
inline CsrMatrix kitaev_chain(const KitaevParams& p) {
  if (p.L < 1)
    throw std::invalid_argument("kitaev_chain: L must be >= 1");
  const int L = p.L;
  CsrMatrix H;
  H.n = 2 * L;
  H.row_ptr.reserve(2 * L + 1);
  H.row_ptr.push_back(0);
  H.col_idx.reserve(10 * L - 8);
  H.values.reserve(10 * L - 8);
  auto emit = [&H](int col, double val) {
    H.col_idx.push_back(col);
    H.values.push_back(val);
  };
  for (int i = 0; i < L; ++i) { // electron row of site i
    if (i > 0)
      emit(i - 1, -p.t);
    emit(i, -p.mu);
    if (i < L - 1)
      emit(i + 1, -p.t);
    if (i > 0)
      emit(L + i - 1, p.delta);
    if (i < L - 1)
      emit(L + i + 1, -p.delta);
    H.row_ptr.push_back(H.nnz());
  }
  for (int i = 0; i < L; ++i) { // hole row of site i
    if (i > 0)
      emit(i - 1, -p.delta);
    if (i < L - 1)
      emit(i + 1, p.delta);
    if (i > 0)
      emit(L + i - 1, p.t);
    emit(L + i, p.mu);
    if (i < L - 1)
      emit(L + i + 1, p.t);
    H.row_ptr.push_back(H.nnz());
  }
  return H;
}

/// Dense ED eigensystem: values ascending, column k of vectors the
/// eigenvector of values[k]. Eigen backs the oracle only.
struct Eigenpairs {
  std::vector<double> values;
  Eigen::MatrixXd vectors;
};

/// R1: dense exact-diagonalization oracle. Throws if the eigensolve
/// fails.
inline Eigenpairs exact_eigenpairs(const CsrMatrix& A) {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(A.n, A.n);
  for (int row = 0; row < A.n; ++row)
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      H(row, A.col_idx[p]) = A.values[p];
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("exact_eigenpairs: eigensolve failed");
  Eigenpairs out;
  const auto& ev = solver.eigenvalues();
  out.values.assign(ev.data(), ev.data() + ev.size());
  out.vectors = solver.eigenvectors();
  return out;
}

/// R1: eigenvalues only, ascending.
inline std::vector<double> exact_eigenvalues(const CsrMatrix& A) {
  return exact_eigenpairs(A).values;
}

/// R2: local Chebyshev moments mu_n(i) = <e_i|T_n(H')|e_i>, with i a
/// row of the 2L Nambu basis and H' = (H - b*I)/a the rescaled chain.
/// Deterministic; no stochastic trace. Vector recurrence t_0 = e_i,
/// t_1 = H' e_i, t_{n+1} = 2 H' t_n - t_{n-1}, so mu_n = t_n[i].
/// Returns mu_0..mu_n_cheb (size n_cheb + 1).
///
/// The recurrence is bounded exactly when the spectrum of H lies in
/// [b - a, b + a]: |T_n| <= 1 on [-1, 1], so ||t_n|| <= ||t_0|| = 1.
/// If a under-covers the spectrum, T_n(H') grows exponentially and the
/// moments are silently wrong, so every recurrence vector is norm-checked
/// (NaN fails the same check, catching a NaN a or b). Roundoff keeps
/// ||t_n||^2 at 1 + O(n * eps); the 1e-6 tolerance is five decades above
/// that and decades below any real violation.
inline std::vector<double> local_moments(const CsrMatrix& A, int site,
                                         const LdosParams& p) {
  if (p.n_cheb < 1)
    throw std::invalid_argument("local_moments: n_cheb must be >= 1");
  if (p.a <= 0.0)
    throw std::invalid_argument(
        "local_moments: a must be a positive spectral bound (a = 0 is unset)");
  if (site < 0 || site >= A.n)
    throw std::invalid_argument("local_moments: site out of range");
  auto spmv_scaled = [&A, a = p.a, b = p.b](const std::vector<double>& x,
                                            std::vector<double>& y) {
    for (int row = 0; row < A.n; ++row) {
      double acc = 0.0;
      for (int q = A.row_ptr[row]; q < A.row_ptr[row + 1]; ++q)
        acc += A.values[q] * x[A.col_idx[q]];
      y[row] = (acc - b * x[row]) / a;
    }
  };
  const auto check_bounded = [](double norm2) {
    if (!(norm2 <= 1.0 + 1e-6))
      throw std::runtime_error(
          "local_moments: recurrence vector norm exceeded 1; the spectrum "
          "of H is not covered by [b - a, b + a] (increase a)");
  };
  std::vector<double> mu(p.n_cheb + 1);
  std::vector<double> t_prev(A.n, 0.0), t_cur(A.n), t_next(A.n);
  t_prev[site] = 1.0;
  mu[0] = 1.0;
  spmv_scaled(t_prev, t_cur);
  double norm2_1 = 0.0;
  for (double x : t_cur)
    norm2_1 += x * x;
  check_bounded(norm2_1);
  mu[1] = t_cur[site];
  for (int n = 1; n < p.n_cheb; ++n) {
    spmv_scaled(t_cur, t_next);
    double norm2 = 0.0;
    for (int i = 0; i < A.n; ++i) {
      t_next[i] = 2.0 * t_next[i] - t_prev[i];
      norm2 += t_next[i] * t_next[i];
    }
    check_bounded(norm2);
    mu[n + 1] = t_next[site];
    std::swap(t_prev, t_cur);
    std::swap(t_cur, t_next);
  }
  return mu;
}

/// R3: Jackson damping factors g_0..g_n_cheb (Weisse, Wellein,
/// Alvermann, Fehske, Rev. Mod. Phys. 78, 275 (2006), Eq. 71).
/// g_0 = 1, g_n -> 0 as n -> n_cheb; the damping that kills the Gibbs
/// oscillations of the truncated expansion.
inline std::vector<double> jackson_kernel(int n_cheb) {
  if (n_cheb < 1)
    throw std::invalid_argument("jackson_kernel: n_cheb must be >= 1");
  std::vector<double> g(n_cheb + 1);
  const double denom = static_cast<double>(n_cheb) + 1.0;
  const double cot = 1.0 / std::tan(std::numbers::pi / denom);
  for (int n = 0; n <= n_cheb; ++n) {
    const double phi = std::numbers::pi * n / denom;
    g[n] =
        ((n_cheb - n + 1.0) * std::cos(phi) + std::sin(phi) * cot) / denom;
  }
  return g;
}

/// R3: local DOS on an energy grid from damped local moments:
///   rho'(eps) = [g_0 mu_0 + 2 sum_{n=1}^{N-1} g_n mu_n T_n(eps)]
///               / (pi sqrt(1 - eps^2)),   eps in (-1, 1),
/// with N + 1 = moments.size(); the grid is n_grid points uniform in
/// theta = arccos(eps). Unscaled: rho(E) = rho'((E - b)/a) / a with
/// E = a*eps + b. T_n by the recurrence T_{n+1} = 2 eps T_n - T_{n-1}
/// (never cos(n arccos eps)). Returns (E, rho) pairs sorted ascending
/// in E.
inline std::vector<std::pair<double, double>>
ldos(const std::vector<double>& moments, const std::vector<double>& jackson,
     double a, double b, int n_grid) {
  if (moments.size() < 2 || jackson.size() < moments.size())
    throw std::invalid_argument("ldos: moments/jackson size mismatch");
  if (a <= 0.0)
    throw std::invalid_argument("ldos: a must be positive");
  if (n_grid < 2)
    throw std::invalid_argument("ldos: n_grid must be >= 2");
  const int n_cheb = static_cast<int>(moments.size()) - 1;
  std::vector<double> w(n_cheb);
  for (int n = 1; n <= n_cheb - 1; ++n)
    w[n] = 2.0 * jackson[n] * moments[n];
  std::vector<std::pair<double, double>> out;
  out.reserve(n_grid);
  for (int k = 0; k < n_grid; ++k) {
    const double theta = std::numbers::pi * (k + 0.5) / n_grid;
    const double eps = std::cos(theta);
    double sum = jackson[0] * moments[0];
    double t_prev = 1.0, t_cur = eps;
    for (int n = 1; n <= n_cheb - 1; ++n) {
      sum += w[n] * t_cur;
      const double t_next = 2.0 * eps * t_cur - t_prev;
      t_prev = t_cur;
      t_cur = t_next;
    }
    const double rho_prime = sum / (std::numbers::pi * std::sin(theta));
    out.emplace_back(a * eps + b, rho_prime / a);
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace majorana
