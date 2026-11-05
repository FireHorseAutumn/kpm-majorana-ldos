// R1 gates (DESIGN.md section 3) at L = 32, t = 1: particle-hole
// symmetry of the ED spectrum, the exact sweet-spot zero mode, the
// trivial-phase gap against the derived bulk bound, and Tr H = 0.
// Every check prints its measured value; any failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

double phs_max(const std::vector<double>& ev) {
  double m = 0.0;
  const int n = static_cast<int>(ev.size());
  for (int k = 0; k < n; ++k)
    m = std::max(m, std::fabs(ev[k] + ev[n - 1 - k]));
  return m;
}

double trace(const majorana::CsrMatrix& A) {
  double tr = 0.0;
  for (int r = 0; r < A.n; ++r)
    for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
      if (A.col_idx[p] == r)
        tr += A.values[p];
  return tr;
}

} // namespace

int main() {
  using namespace majorana;
  char buf[256];
  const int L = 32;
  const int n = 2 * L;

  const KitaevParams sweet{L, 1.0, 1.0, 0.0};
  const KitaevParams triv{L, 1.0, 1.0, 3.0};
  const CsrMatrix Hs = kitaev_chain(sweet);
  const CsrMatrix Ht = kitaev_chain(triv);

  std::snprintf(buf, sizeof buf, "n=%d, nnz=%d (expect %d)", Hs.n, Hs.nnz(),
                10 * L - 8);
  check(Hs.n == n && Hs.nnz() == 10 * L - 8, "builder dimensions", buf);

  // The BdG matrix is symmetric by construction; checked bitwise.
  std::vector<double> dense(n * n, 0.0);
  for (int r = 0; r < n; ++r)
    for (int p = Hs.row_ptr[r]; p < Hs.row_ptr[r + 1]; ++p)
      dense[r * n + Hs.col_idx[p]] = Hs.values[p];
  bool symmetric = true;
  for (int r = 0; r < n && symmetric; ++r)
    for (int c = 0; c < n; ++c)
      if (dense[r * n + c] != dense[c * n + r]) {
        symmetric = false;
        break;
      }
  check(symmetric, "H symmetric (bitwise)", "densified CSR compared");

  const double tr_s = trace(Hs);
  const double tr_t = trace(Ht);
  std::snprintf(buf, sizeof buf, "sweet %.3e, mu=3 %.3e (bound 1e-13)", tr_s,
                tr_t);
  check(std::fabs(tr_s) < 1e-13 && std::fabs(tr_t) < 1e-13, "Tr H = 0", buf);

  const Eigenpairs es = exact_eigenpairs(Hs);
  const Eigenpairs et = exact_eigenpairs(Ht);

  const double phs_s = phs_max(es.values);
  const double phs_t = phs_max(et.values);
  std::snprintf(buf, sizeof buf, "sweet %.3e, mu=3 %.3e (bound 1e-12)", phs_s,
                phs_t);
  check(phs_s < 1e-12 && phs_t < 1e-12, "PHS max|E_k + E_2L-1-k|", buf);

  // Sweet spot: exact zero modes, one per end; the null space is
  // supported on the end sites only, so any returned eigenbasis of it
  // has its full weight there. End-site rows: 0, L-1 (electron) and
  // L, 2L-1 (hole).
  int n_zero = 0;
  double e_zero = 0.0;
  double w_min = 1.0;
  for (int k = 0; k < n; ++k) {
    if (std::fabs(es.values[k]) >= 1e-12)
      continue;
    ++n_zero;
    e_zero = std::max(e_zero, std::fabs(es.values[k]));
    const auto& V = es.vectors;
    const double w = V(0, k) * V(0, k) + V(L - 1, k) * V(L - 1, k) +
                     V(L, k) * V(L, k) + V(2 * L - 1, k) * V(2 * L - 1, k);
    w_min = std::min(w_min, w);
  }
  std::snprintf(buf, sizeof buf,
                "%d zero modes, max |E| %.3e, min end weight %.17g (bound 1 "
                "- 1e-12)",
                n_zero, e_zero, w_min);
  check(n_zero >= 1 && w_min >= 1.0 - 1e-12,
        "sweet-spot zero mode on the two end sites", buf);

  // Same derivation as the gap bound below, at mu = 0: E^2(k) = 4t^2,
  // the bulk band is flat at |E| = 2t. Pins the hopping scale, which
  // PHS and the trace cannot see.
  double flat = 0.0;
  for (int k = 0; k < n; ++k) {
    if (std::fabs(es.values[k]) < 1e-12)
      continue;
    flat = std::max(flat, std::fabs(std::fabs(es.values[k]) - 2.0));
  }
  std::snprintf(buf, sizeof buf, "max ||E| - 2t| %.3e (bound 1e-12)", flat);
  check(flat < 1e-12, "sweet-spot bulk modes at |E| = 2t", buf);

  // Bulk bound derived from the built chain's Bloch Hamiltonian
  // [[h(k), D(k)], [D*(k), -h(k)]] with h(k) = -mu - 2t cos k and
  // D(k) = -2i delta sin k: E(k) = ±sqrt((2t cos k + mu)^2 +
  // 4 delta^2 sin^2 k). At delta = t this is E^2(k) = 4t^2 + mu^2 +
  // 4t mu cos k, monotonic in cos k; for mu > 0 the minimum sits at
  // k = pi, so E_gap = |mu - 2t| = 1 at mu = 3. The open chain samples
  // discrete momenta and misses k = pi, so its min |E| sits above the
  // bound; the 1e-12 slack covers rounding only. Derived from the
  // chain, not quoted from a source (DESIGN.md section 7).
  double e_min_t = 0.0;
  for (int k = 0; k < n; ++k)
    e_min_t = (k == 0) ? std::fabs(et.values[k])
                       : std::min(e_min_t, std::fabs(et.values[k]));
  const double e_gap = std::fabs(triv.mu - 2.0 * triv.t);
  std::snprintf(buf, sizeof buf, "min|E| %.17g vs derived bulk gap %.17g",
                e_min_t, e_gap);
  check(e_min_t >= e_gap - 1e-12,
        "trivial phase: min|E| bounded by the analytic bulk gap", buf);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
