#pragma once

#include <Eigen/Dense>

#include <cstdint>
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
/// row of the 2L Nambu basis. Deterministic; no stochastic trace.
inline std::vector<double> local_moments(const CsrMatrix&, int,
                                         const LdosParams&) {
  throw std::logic_error("local_moments: not implemented (DESIGN.md R2)");
}

/// R3: Jackson damping factors g_0..g_n_cheb.
inline std::vector<double> jackson_kernel(int) {
  throw std::logic_error("jackson_kernel: not implemented (DESIGN.md R3)");
}

/// R3: local DOS on an energy grid from damped local moments.
inline std::vector<std::pair<double, double>>
ldos(const std::vector<double>&, const std::vector<double>&, double, double,
     int) {
  throw std::logic_error("ldos: not implemented (DESIGN.md R3)");
}

} // namespace majorana
