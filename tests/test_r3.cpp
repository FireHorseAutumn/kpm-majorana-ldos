// R3 gates (DESIGN.md section 3): the zero-mode figure. L = 1000,
// t = delta = 1, mu = 0 (the sweet spot); n_cheb = 4000; the rescaling
// a = sqrt(8), b = 0 puts the bulk bands at eps = +-1/sqrt(2) and the
// two exact zero modes at eps = 0. End-site row: 0 (electron component
// of the left end); bulk row: L/2. The window is |eps| <= w = 0.2
// (|E| <= 0.566), deep inside the gap between the bulk bands.
//
// Every bound is pinned from the Jackson kernel arithmetic and the ED
// oracle BEFORE the KPM moments are computed; no bound is derived from
// a KPM output. Derivation: integrating the damped reconstruction of a
// unit delta at eps_k over |eps| <= w is exact arithmetic on the
// moments mu_n = T_n(eps_k). With th0 = arccos(w),
//
//   K(eps_k) = (pi - 2 th0)/pi
//              - (4/pi) sum_{n even, 2 <= n <= N-1} g_n T_n(eps_k)
//                        sin(n th0)/n,
//
// the odd n vanishing under the symmetric window; the truncation
// N - 1 matches ldos(). Then W_end = sum_k w_k K(eps_k) exactly: the
// zero modes contribute w_ED * K(0), the bulk bands the rest.
//
//   * Peak weight (the DESIGN gate): |W_end - w_ED| <=
//     w_ED*|1 - K(0)| + (1 - w_ED)*K_leak + 1e-10, with K_leak = max
//     over bulk eigenvalues of |K(eps_k)|. The first term is the known
//     broadening — the weight the Jackson peak leaves outside the
//     window; the second the bulk-band leakage in; the last is
//     roundoff, pinned at 10x the estimate below. The bound is
//     saturated by the known-broadening term by construction: the gate
//     certifies there is no weight beyond the known shape.
//   * Bulk leakage: W_bulk <= K_leak + K(0)*wz_bulk + 1e-10, wz_bulk
//     the zero-mode weight on the bulk row (asserted < 1e-20).
//   * Peak position: |E_peak| <= 1.5*dE0 with dE0 = a*pi/n_grid the
//     grid spacing at E = 0 — the peak may land in the bin containing
//     0 or in either neighbor ("within one grid bin").
//   * Machinery: |W - W_pred| <= 1e-9 per row with W_pred = sum_k w_k
//     K(eps_k) the full ED+kernel prediction; 1e-9 is 50x the combined
//     roundoff estimate.
//
// The R2 conditioning ramp is respected by construction: no gate at
// this rung compares single moments at n = 4000, where the R2
// cross-check measured a linear roundoff ramp reaching ~1e-12 per
// moment. The gated window weights are sums of moments whose
// coefficient l1 norm is <= (4/pi) sum_m 1/(2m) ~ 4.4, so the moment
// noise enters the weights at <= ~5e-12 — twenty times below the 1e-10
// slack, far below the 1e-9 machinery bound. LDOS grid values carry
// <= (2N/pi)/sin(th) * 1e-12 ~ 3e-9, immaterial against O(10..100)
// peak heights and the one-bin position gate.
//
// Heights are resolution-limited — they grow linearly in n_cheb/a —
// so they are reported, never gated (DESIGN.md section 3).
//
// The carried R2 gate (LDOS summed over all 2L rows vs the ED
// histogram, kpm-tightbinding L1 discipline) runs at its native
// L = 32, off-sweet and trivial: L1(KPM summed, hist) <= L1(ED-moments
// reconstruction, hist) + 1e-6, the bound pinned from the oracle and
// the kernel — same reconstruction, same bins — before the KPM sum is
// formed.
//
// Every check prints its measured value; any failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
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

// K(eps_k), the Jackson-broadened weight of a unit delta at eps_k
// inside |eps| <= w, th0 = arccos(w); derivation in the header. T_n by
// the recurrence, never cos(n arccos eps).
double window_weight(double eps_k, double th0,
                     const std::vector<double>& g) {
  const int n_cheb = static_cast<int>(g.size()) - 1;
  double k = (std::numbers::pi - 2.0 * th0) / std::numbers::pi;
  double t_nm2 = 1.0;   // T_{n-2}
  double t_nm1 = eps_k; // T_{n-1}
  for (int n = 2; n <= n_cheb - 1; n += 2) {
    const double t_n = 2.0 * eps_k * t_nm1 - t_nm2;
    const double t_np1 = 2.0 * eps_k * t_n - t_nm1;
    k -= (4.0 / std::numbers::pi) * g[n] * t_n * std::sin(n * th0) / n;
    t_nm2 = t_n;
    t_nm1 = t_np1;
  }
  return k;
}

// The same window weight assembled from a computed moment sequence:
// W = mu_0 (pi - 2 th0)/pi - (4/pi) sum g_n mu_n sin(n th0)/n over
// even 2 <= n <= N-1.
double kpm_window_weight(const std::vector<double>& mu, double th0,
                         const std::vector<double>& g) {
  const int n_cheb = static_cast<int>(mu.size()) - 1;
  double w = mu[0] * (std::numbers::pi - 2.0 * th0) / std::numbers::pi;
  for (int n = 2; n <= n_cheb - 1; n += 2)
    w -= (4.0 / std::numbers::pi) * g[n] * mu[n] * std::sin(n * th0) / n;
  return w;
}

// Normalized DOS histogram, bins of width dE centered at -emax,
// -emax + dE, ..., integrating to 1.
std::vector<std::pair<double, double>>
histogram(const std::vector<double>& evals, double dE, double emax) {
  const int n_bins = static_cast<int>(std::llround(2.0 * emax / dE)) + 1;
  std::vector<std::pair<double, double>> hist(n_bins);
  for (int k = 0; k < n_bins; ++k)
    hist[k] = {-emax + k * dE, 0.0};
  for (double E : evals) {
    const int k = static_cast<int>(std::floor((E + emax + 0.5 * dE) / dE));
    if (k < 0 || k >= n_bins)
      throw std::out_of_range("histogram: eigenvalue outside range");
    hist[k].second += 1.0;
  }
  for (auto& bin : hist)
    bin.second /= static_cast<double>(evals.size()) * dE;
  return hist;
}

// Mean rho of the fine KPM grid inside the histogram bin centered at
// center, width dE.
double bin_average(const std::vector<std::pair<double, double>>& dos,
                   double center, double dE) {
  double acc = 0.0;
  int cnt = 0;
  for (const auto& p : dos)
    if (std::fabs(p.first - center) <= 0.5 * dE) {
      acc += p.second;
      ++cnt;
    }
  return cnt > 0 ? acc / cnt : 0.0;
}

// Window weight by the midpoint rule in theta = arccos(E/a) of an
// ldos() grid over |E| <= e_win: an independent path from moments to
// window weight that shares nothing but the moments. Every other weight
// gate consumes the analytic window formula on both sides, so an error
// in that formula — or in ldos()'s normalization (the 2/pi, the 1/a) —
// cancels out of them; this check crosses the two paths. Measured
// quadrature error at this file's settings (L = 500..1000, n_cheb =
// 4000, w = 0.2): 6.8e-11 at n_grid = 2048, ~1e-13 at n_grid = 16384;
// the 1e-8 pin carries 100x margin over the coarser-grid measurement.
double trapz_window(const std::vector<std::pair<double, double>>& dos,
                    double a, double e_win, int n_grid) {
  const double dth = std::numbers::pi / n_grid;
  double acc = 0.0;
  for (const auto& p : dos) {
    if (std::fabs(p.first) > e_win)
      continue;
    acc += p.second * a * std::sin(std::acos(p.first / a)) * dth;
  }
  return acc;
}

double l1_vs_hist(const std::vector<std::pair<double, double>>& dos,
                  const std::vector<std::pair<double, double>>& hist,
                  double dE) {
  double l1 = 0.0;
  for (const auto& bin : hist)
    l1 += std::fabs(bin_average(dos, bin.first, dE) - bin.second) * dE;
  return l1;
}

double spectral_bound(const majorana::KitaevParams& kp) {
  const double s = 2.0 * kp.t + std::fabs(kp.mu);
  return std::sqrt(s * s + 4.0 * kp.delta * kp.delta);
}

// Jackson kernel formula gates, independent of any KPM use. Every
// weight gate in this file and test_r4.cpp consumes g through
// window_weight on both sides, so a formula transcription error cancels
// there (verified: a denominator off-by-one N + 1 -> N passes all of
// them). These checks pin the formula itself against Weisse et al.,
// Rev. Mod. Phys. 78, 275 (2006): g_0 = 1; g_1 = cos(pi/(N+1)) exactly
// (Eq. 68); positivity on n < N with g_N = 0 and monotone decrease (the
// Fejer-square construction, Eq. 67); and the delta-variance identity
// (1 - g_2)/2 = N (1 - cos(2 pi/(N+1))) / (2 (N+1)) (Eq. 75 at a = 0,
// T_2(0) = -1). Measured on this toolchain: |g_1 - cos| = 0 bitwise at
// every tested N; the sigma^2 identity matches to 7.6e-10 relative at
// N = 16384, exactly at small N (roundoff in 1 - g_2); pins carry
// 100x margin over the measured values.
void jackson_formula_gates(int n_cheb) {
  char buf[256];
  const auto g = majorana::jackson_kernel(n_cheb);
  const double d = n_cheb + 1.0;

  std::snprintf(buf, sizeof buf, "N = %d, g_0 = %.17g (exact 1)", n_cheb,
                g[0]);
  check(g[0] == 1.0, "Jackson g_0 == 1 (bitwise)", buf);

  const double g1_err = std::fabs(g[1] - std::cos(std::numbers::pi / d));
  std::snprintf(buf, sizeof buf, "N = %d, |g_1 - cos(pi/(N+1))| = %.3e "
                "(RMP Eq. 68; bound 1e-15)", n_cheb, g1_err);
  check(g1_err <= 1e-15, "Jackson g_1 = cos(pi/(N+1))", buf);

  double min_g = 2.0;
  bool monotone = true;
  for (int n = 1; n < n_cheb; ++n) {
    min_g = std::min(min_g, g[n]);
    if (!(g[n] < g[n - 1]))
      monotone = false;
  }
  std::snprintf(buf, sizeof buf, "N = %d, min g[1..N-1] = %.3e, monotone "
                "decrease: %s", n_cheb, min_g, monotone ? "yes" : "NO");
  check(min_g > 0.0 && monotone,
        "Jackson kernel positive and decreasing on n < N", buf);

  std::snprintf(buf, sizeof buf, "N = %d, |g_N| = %.3e (bound 1e-14)",
                n_cheb, std::fabs(g[n_cheb]));
  check(std::fabs(g[n_cheb]) <= 1e-14, "Jackson g_N = 0", buf);

  const double s2 = (g[0] - g[2]) / 2.0;
  const double s2_rmp =
      n_cheb * (1.0 - std::cos(2.0 * std::numbers::pi / d)) / (2.0 * d);
  const double s2_rel = std::fabs(s2 - s2_rmp) / s2_rmp;
  std::snprintf(buf, sizeof buf,
                "N = %d, rel err of (1 - g_2)/2 vs RMP Eq. 75 = %.3e "
                "(bound 1e-8)",
                n_cheb, s2_rel);
  check(s2_rel <= 1e-8, "Jackson delta variance at 0 (RMP Eq. 75)", buf);
}

// The carried R2 gate at L = 32: LDOS summed over all 2L rows vs the
// ED histogram. The bound is pinned first: the same reconstruction
// built from the ED-derived total moments against the same histogram,
// plus 1e-6 of roundoff slack.
void carried_r2_gate(const majorana::KitaevParams& kp, const char* name) {
  using namespace majorana;
  char buf[256];
  const int n_cheb = 512, n_grid = 2048;
  const double dE = 0.1;
  const CsrMatrix H = kitaev_chain(kp);
  const double a = spectral_bound(kp);
  const Eigenpairs es = exact_eigenpairs(H);
  const std::vector<double> g = jackson_kernel(n_cheb);
  const LdosParams lp{n_cheb, a, 0.0};
  const auto hist = histogram(es.values, dE, a);

  // The pin: oracle + kernel only, no KPM output.
  std::vector<double> mu_ed(n_cheb + 1, 0.0);
  for (double E : es.values) {
    const double eps = E / a;
    double t_prev = 1.0, t_cur = eps;
    mu_ed[0] += 1.0;
    mu_ed[1] += t_cur;
    for (int n = 1; n < n_cheb; ++n) {
      const double t_next = 2.0 * eps * t_cur - t_prev;
      mu_ed[n + 1] += t_next;
      t_prev = t_cur;
      t_cur = t_next;
    }
  }
  for (double& m : mu_ed)
    m /= H.n;
  const auto dos_ed = ldos(mu_ed, g, a, 0.0, n_grid);
  const double l1_ed = l1_vs_hist(dos_ed, hist, dE);
  const double bound = l1_ed + 1e-6;

  // The KPM side: per-row LDOS, summed over all 2L rows, normalized.
  std::vector<double> summed(n_grid, 0.0);
  for (int row = 0; row < H.n; ++row) {
    const auto dos = ldos(local_moments(H, row, lp), g, a, 0.0, n_grid);
    for (int k = 0; k < n_grid; ++k)
      summed[k] += dos[k].second;
  }
  std::vector<std::pair<double, double>> dos_kpm(n_grid);
  for (int k = 0; k < n_grid; ++k)
    dos_kpm[k] = {dos_ed[k].first, summed[k] / H.n};
  const double l1_kpm = l1_vs_hist(dos_kpm, hist, dE);

  std::snprintf(buf, sizeof buf,
                "%s: L1 %.6f (pin %.6f = ED-moments L1 %.6f + 1e-6)",
                name, l1_kpm, bound, l1_ed);
  check(l1_kpm <= bound, "R2 carried: row-summed LDOS vs ED histogram", buf);
}

} // namespace

int main() {
  using namespace majorana;
  char buf[256];
  const int L = 1000;
  const int n = 2 * L;
  const int n_cheb = 4000;
  const int n_grid = 16384;
  const double w = 0.2;
  const KitaevParams sweet{L, 1.0, 1.0, 0.0};
  const double a = std::sqrt(8.0);
  const CsrMatrix H = kitaev_chain(sweet);
  const Eigenpairs es = exact_eigenpairs(H);

  // Structural anchors, R1 lineage at the R3 scale.
  int n_zero = 0;
  double e_zero = 0.0, w_ed = 0.0, wz_bulk = 0.0, bulk_dev = 0.0;
  for (int k = 0; k < n; ++k) {
    const double E = es.values[k];
    if (std::fabs(E) < 1e-12) {
      ++n_zero;
      e_zero = std::max(e_zero, std::fabs(E));
      w_ed += es.vectors(0, k) * es.vectors(0, k);
      wz_bulk += es.vectors(L / 2, k) * es.vectors(L / 2, k);
    } else {
      bulk_dev = std::max(bulk_dev, std::fabs(std::fabs(E) - 2.0));
    }
  }
  std::snprintf(buf, sizeof buf, "%d modes, max |E| %.3e (bound 1e-12)",
                n_zero, e_zero);
  check(n_zero == 2 && e_zero < 1e-12, "two exact zero modes at L = 1000",
        buf);
  std::snprintf(buf, sizeof buf, "w_ED = %.17g (bound: 1/2 +- 1e-12)", w_ed);
  check(std::fabs(w_ed - 0.5) < 1e-12,
        "zero-mode end weight on row 0 (R1: exact end Majoranas)", buf);
  std::snprintf(buf, sizeof buf, "max ||E| - 2t| %.3e (bound 1e-11)",
                bulk_dev);
  check(bulk_dev < 1e-11, "bulk bands at |E| = 2t", buf);
  std::snprintf(buf, sizeof buf, "%.3e (bound 1e-20)", wz_bulk);
  check(wz_bulk < 1e-20, "zero-mode weight on the bulk row", buf);

  // ---- Pins: kernel arithmetic and the ED oracle only. ----
  const std::vector<double> g = jackson_kernel(n_cheb);
  const double th0 = std::acos(w);
  const double k0 = window_weight(0.0, th0, g);
  double k_leak = 0.0;
  double w_pred_end = 0.0, w_pred_bulk = 0.0;
  for (int k = 0; k < n; ++k) {
    const double kw = window_weight(es.values[k] / a, th0, g);
    const double w_end = es.vectors(0, k) * es.vectors(0, k);
    const double w_blk = es.vectors(L / 2, k) * es.vectors(L / 2, k);
    w_pred_end += w_end * kw;
    w_pred_bulk += w_blk * kw;
    if (std::fabs(es.values[k]) >= 1e-12)
      k_leak = std::max(k_leak, std::fabs(kw));
  }
  const double bound_weight =
      w_ed * std::fabs(1.0 - k0) + (1.0 - w_ed) * k_leak + 1e-10;
  const double bound_bulk = k_leak + k0 * wz_bulk + 1e-10;
  std::printf("[info] pins: K(0) = %.12f (1 - K(0) = %.3e), K_leak = %.3e, "
              "W_pred end %.12f, bulk %.6e\n",
              k0, 1.0 - k0, k_leak, w_pred_end, w_pred_bulk);

  // ---- The KPM side. ----
  const LdosParams lp{n_cheb, a, 0.0};
  const std::vector<double> mu_end = local_moments(H, 0, lp);
  const std::vector<double> mu_bulk = local_moments(H, L / 2, lp);
  const double w_end = kpm_window_weight(mu_end, th0, g);
  const double w_bulk = kpm_window_weight(mu_bulk, th0, g);

  std::snprintf(buf, sizeof buf, "|W_end - w_ED| = %.6e (pin %.6e)",
                std::fabs(w_end - w_ed), bound_weight);
  check(std::fabs(w_end - w_ed) <= bound_weight,
        "end-site peak weight vs ED zero-mode end weight", buf);
  std::snprintf(buf, sizeof buf, "|W_end - W_pred| = %.3e (bound 1e-9)",
                std::fabs(w_end - w_pred_end));
  check(std::fabs(w_end - w_pred_end) <= 1e-9,
        "end-site window weight vs ED+kernel prediction", buf);
  std::snprintf(buf, sizeof buf, "W_bulk = %.6e (pin %.6e)", w_bulk,
                bound_bulk);
  check(w_bulk <= bound_bulk,
        "bulk-site in-gap weight below the kernel leakage bound", buf);
  std::snprintf(buf, sizeof buf, "|W_bulk - W_pred| = %.3e (bound 1e-9)",
                std::fabs(w_bulk - w_pred_bulk));
  check(std::fabs(w_bulk - w_pred_bulk) <= 1e-9,
        "bulk-site window weight vs ED+kernel prediction", buf);

  // Peak position on the reconstructed grid; heights reported, never
  // gated.
  const auto dos_end = ldos(mu_end, g, a, 0.0, n_grid);
  const auto dos_bulk = ldos(mu_bulk, g, a, 0.0, n_grid);
  const double dE0 = a * std::numbers::pi / n_grid;
  int k_peak = 0;
  for (int k = 1; k < n_grid; ++k)
    if (dos_end[k].second > dos_end[k_peak].second)
      k_peak = k;
  std::snprintf(buf, sizeof buf, "|E_peak| = %.6e (bound %.6e = 1.5 bins)",
                std::fabs(dos_end[k_peak].first), 1.5 * dE0);
  check(std::fabs(dos_end[k_peak].first) <= 1.5 * dE0,
        "end-site LDOS peak within one grid bin of E = 0", buf);

  // Quadrature cross-gate: the analytic window weight against the
  // direct trapezoid of the ldos() reconstruction (helper comment).
  const double q_end = trapz_window(dos_end, a, w * a, n_grid);
  const double q_blk = trapz_window(dos_bulk, a, w * a, n_grid);
  std::snprintf(buf, sizeof buf,
                "end %.3e, bulk %.3e (bound 1e-8; quadrature error at this "
                "grid ~1e-13)",
                std::fabs(q_end - w_end), std::fabs(q_blk - w_bulk));
  check(std::fabs(q_end - w_end) <= 1e-8 &&
            std::fabs(q_blk - w_bulk) <= 1e-8,
        "ldos() trapezoid vs analytic window weight", buf);

  double h_zero = 0.0, h_band_end = 0.0, h_zero_bulk = 0.0,
         h_band_bulk = 0.0;
  for (int k = 0; k < n_grid; ++k) {
    const double E = dos_end[k].first;
    if (std::fabs(E) < 0.1) {
      h_zero = std::max(h_zero, dos_end[k].second);
      h_zero_bulk = std::max(h_zero_bulk, dos_bulk[k].second);
    }
    if (std::fabs(std::fabs(E) - 2.0) < 0.1) {
      h_band_end = std::max(h_band_end, dos_end[k].second);
      h_band_bulk = std::max(h_band_bulk, dos_bulk[k].second);
    }
  }
  std::printf("[info] resolution-limited heights (n_cheb = %d, a = sqrt(8)): "
              "end-site zero peak %.4f, end-site bulk band %.4f; bulk-site "
              "zero region %.6f, bulk-site bulk band %.4f\n",
              n_cheb, h_zero, h_band_end, h_zero_bulk, h_band_bulk);
  std::printf("[info] bulk-site argmax E = %.6f (the bulk band; no zero "
              "peak)\n",
              std::max_element(dos_bulk.begin(), dos_bulk.end(),
                               [](const auto& p, const auto& q) {
                                 return p.second < q.second;
                               })
                  ->first);

  // ---- The carried R2 gate at its native L = 32. ----
  carried_r2_gate(KitaevParams{32, 1.0, 0.8, 0.5}, "off-sweet");
  carried_r2_gate(KitaevParams{32, 1.0, 1.0, 3.0}, "trivial  ");

  // ---- Jackson formula gates (independent of the KPM paths). ----
  jackson_formula_gates(7);
  jackson_formula_gates(n_cheb);

  // Contract checks.
  check(throws_invalid_argument([&] { jackson_kernel(0); }),
        "jackson_kernel(0) rejected", "invalid_argument");
  check(throws_invalid_argument(
            [&] { ldos({1.0}, jackson_kernel(4), a, 0.0, 16); }),
        "ldos with a single moment rejected", "invalid_argument");
  check(throws_invalid_argument(
            [&] { ldos({1.0, 0.0}, jackson_kernel(4), 0.0, 0.0, 16); }),
        "ldos with a = 0 rejected", "invalid_argument");
  check(throws_invalid_argument(
            [&] { ldos({1.0, 0.0}, jackson_kernel(4), a, 0.0, 1); }),
        "ldos with n_grid < 2 rejected", "invalid_argument");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
