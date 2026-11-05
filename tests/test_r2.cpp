// R2 gates (DESIGN.md section 3) at L = 32, t = 1: the local Chebyshev
// moments mu_n(i) = <e_i|T_n(H')|e_i> from the site recurrence against
// the ED-derived moments sum_k |psi_k(i)|^2 T_n(E'_k), with H' = (H -
// b*I)/a and a = sqrt((2t + |mu|)^2 + 4 delta^2) the documented
// spectral bound (the LdosParams comment). linf is taken over all 2L
// rows and all n <= n_cheb; bound 1e-13.
//
// The gated n_cheb is 16, not the LdosParams default 2048. The
// cross-check's own FP64 roundoff grows linearly in the moment index
// (both sides carry eigensolver/recurrence rounding; the recurrence
// itself is the standard stable three-term form). Measured on this
// toolchain: worst linf over all rows is 2.275e-14 at n = 16 across
// the three parameter sets below, and for the sweet chain grows to
// 5.9e-14 at n = 64 and 9.8e-13 at n = 1024 (printed as information
// at the end of the run). The 1e-13 bound therefore bites only at low
// moment index, so the gate sits at n_cheb = 16. R3 gates the
// reconstructed LDOS against ED directly, so no moment-level gate is
// needed at high n.
//
// Every check prints its measured value; any failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

template <typename F>
bool throws_invalid_argument(F&& f) {
  try {
    f();
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

// ED-derived local moments sum_k |psi_k(site)|^2 T_n(E'_k), n =
// 0..n_cheb, T_n by the recurrence T_{n+1} = 2 eps T_n - T_{n-1}.
std::vector<double> ed_local_moments(const majorana::Eigenpairs& es, int site,
                                     double a, double b, int n_cheb) {
  std::vector<double> mu(n_cheb + 1, 0.0);
  const int n = static_cast<int>(es.values.size());
  for (int k = 0; k < n; ++k) {
    const double w = es.vectors(site, k) * es.vectors(site, k);
    const double eps = (es.values[k] - b) / a;
    double t_prev = 1.0, t_cur = eps;
    mu[0] += w;
    mu[1] += w * t_cur;
    for (int m = 1; m < n_cheb; ++m) {
      const double t_next = 2.0 * eps * t_cur - t_prev;
      mu[m + 1] += w * t_next;
      t_prev = t_cur;
      t_cur = t_next;
    }
  }
  return mu;
}

// Worst |mu - mu_ED| over all 2L rows and n <= n_cheb; argmax row in
// *site if given.
double worst_linf(const majorana::CsrMatrix& H, const majorana::Eigenpairs& es,
                  const majorana::LdosParams& p, int* site) {
  double worst = 0.0;
  for (int i = 0; i < H.n; ++i) {
    const std::vector<double> mu = majorana::local_moments(H, i, p);
    const std::vector<double> ref =
        ed_local_moments(es, i, p.a, p.b, p.n_cheb);
    for (int n = 0; n <= p.n_cheb; ++n) {
      const double d = std::fabs(mu[n] - ref[n]);
      if (d > worst) {
        worst = d;
        if (site)
          *site = i;
      }
    }
  }
  return worst;
}

double spectral_bound(const majorana::KitaevParams& kp) {
  const double s = 2.0 * kp.t + std::fabs(kp.mu);
  return std::sqrt(s * s + 4.0 * kp.delta * kp.delta);
}

} // namespace

int main() {
  using namespace majorana;
  char buf[256];
  const int L = 32;
  const int n_cheb = 16;

  const KitaevParams sets[3] = {{L, 1.0, 1.0, 0.0},   // sweet spot
                                {L, 1.0, 0.8, 0.5},   // off sweet spot
                                {L, 1.0, 1.0, 3.0}};  // trivial
  const char* names[3] = {"sweet", "off-sweet", "trivial"};

  for (int s = 0; s < 3; ++s) {
    const CsrMatrix H = kitaev_chain(sets[s]);
    const LdosParams lp{n_cheb, spectral_bound(sets[s]), 0.0};
    const Eigenpairs es = exact_eigenpairs(H);

    double eprime = 0.0;
    for (double E : es.values)
      eprime = std::max(eprime, std::fabs((E - lp.b) / lp.a));
    std::snprintf(buf, sizeof buf, "%s: max|E'| %.17g, a %.17g", names[s],
                  eprime, lp.a);
    check(eprime < 1.0, "rescaled spectrum inside (-1, 1)", buf);

    int site = -1;
    const double worst = worst_linf(H, es, lp, &site);
    std::snprintf(buf, sizeof buf,
                  "%s: worst linf %.3e (row %d, %d rows x %d moments, bound "
                  "1e-13)",
                  names[s], worst, site, H.n, n_cheb + 1);
    check(worst < 1e-13, "local_moments vs ED moments, linf", buf);
  }

  // The conditioning ramp behind the n_cheb = 16 placement (see the
  // header comment). Informational; not gated.
  {
    const CsrMatrix H = kitaev_chain(sets[0]);
    const Eigenpairs es = exact_eigenpairs(H);
    const double a = spectral_bound(sets[0]);
    std::printf("[info] sweet worst linf vs n_cheb:");
    for (int m : {64, 256, 1024}) {
      const LdosParams lp{m, a, 0.0};
      std::printf("  %d: %.3e", m, worst_linf(H, es, lp, nullptr));
    }
    std::printf("\n");
  }

  const CsrMatrix H = kitaev_chain(sets[0]);
  const LdosParams good{n_cheb, spectral_bound(sets[0]), 0.0};
  check(throws_invalid_argument(
            [&] { local_moments(H, 0, LdosParams{0, good.a, 0.0}); }),
        "n_cheb = 0 rejected", "invalid_argument");
  check(throws_invalid_argument(
            [&] { local_moments(H, 0, LdosParams{n_cheb, 0.0, 0.0}); }),
        "a = 0 (unset spectral bound) rejected", "invalid_argument");
  check(throws_invalid_argument([&] { local_moments(H, H.n, good); }),
        "site out of range rejected", "invalid_argument");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
