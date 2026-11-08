// R4 gates (DESIGN.md section 3): the transition sweep. t = delta = 1
// fixed, mu swept across 2t = 2: a coarse grid 0.0:0.1:3.0 (31 points)
// plus a fine bracket 1.82:0.02:2.18 (19 points) around the transition.
// L = 256, bulk row L/2 = 128, n_cheb = 16384, window |E| <= E_win =
// 0.05 with th0 = acos(E_win/a) per point; a = sqrt((2t + mu)^2 +
// 4 delta^2) is the documented spectral bound (the LdosParams comment),
// b = 0. mu >= 0 only: the staggered gauge transform c_j -> (-1)^j c_j
// sends mu -> -mu exactly, so the sweep covers |mu| < 2t.
//
// Resolution budget, stated per Kill 1 (DESIGN.md section 4):
//   * E_win = 0.05 is the zero-mode window: a split pair below E_win
//     reads as a zero mode. From the built chain at delta = t the
//     zero-mode Majorana recurrence is u_{i+1} = -(mu/2t) u_i, so the
//     pair splitting sits at ~r^L with r = mu/2t and the row-0 pair
//     weight at (1 - r^2)/2 + O(r^{2L}) — the R3 sweet-spot value 1/2
//     at r = 0; the ED in-window end weight tracks (1 - r^2)/2 to the
//     printed digits at every swept point. At the deepest gated
//     present point mu = 1.8 the splitting is ~0.9^256 ~ 2e-12, ten
//     decades under E_win; measured splittings are printed per point.
//     The present/absent gates are window-weight statements,
//     insensitive to any splitting inside the window.
//   * KPM broadening dE_KPM = pi a/n_cheb runs 5.4e-4 (mu = 0) to
//     1.04e-3 (mu = 3): E_win/dE_KPM >= 48 at every point, so window
//     captures are exact to the K(eps_k) arithmetic below.
//   * mu resolution: the crossover bracket half-width is bounded at
//     0.10, hard from the two crossing arithmetics. Up side: the
//     closed-gap bulk-leakage band ends where the bulk gap reopens
//     past the window — E_1 >= mu - 2t > E_win for mu >= 2.06 (the
//     open chain's lowest level sits above the periodic gap, R1), so
//     mu_up <= 2.06. Down side: the exact row-0 weight (1 - r^2)/2
//     crosses w_min = 0.05 at mu = 2t sqrt(1 - 2 w_min) = 1.897, one
//     fine step gives mu_down >= 1.88. Worst bracket [1.88, 2.06]:
//     half-width 0.09. The down-side offset is the smooth weight decay
//     of the locator, not a KPM resolution effect. At mu = 2t the open
//     chain's lowest bulk level E_1 (measured 0.0122, the finite-size
//     floor at the pi t/L scale) sits inside the window — the signal
//     of the closed gap.
//
// Pins are kernel arithmetic and the ED oracle, computed per point
// BEFORE the KPM moments run. K(eps_k) is the Jackson-broadened weight
// of a unit delta inside |eps| <= E_win/a; derivation in the
// test_r3.cpp header.
//   * Cross-gate (Kill 2): |W - W_pred| <= 1e-9 on both rows at every
//     point, W_pred = sum_k w_k K(eps_k). The R2 conditioning ramp at
//     n = 16384 extrapolates to ~1.6e-11 per moment; the window weights
//     see <= ~7e-11 of it (coefficient l1 norm <= 4.4), 14x under the
//     bound.
//   * Quadrature cross-gate: |trapz(ldos) - W| <= 1e-7 at one point per
//     regime (mu = 1.0 end row, 2.0 bulk row, 2.6 end row). R4 otherwise
//     never calls ldos(): every other gate is a window-weight identity
//     that consumes the analytic K formula on both sides, so a formula
//     or ldos() normalization error would cancel. Measured quadrature
//     error at n_grid = 16384: <= 3.3e-9.
//   * Present (coarse, mu <= 1.8): W_end >= w_min = 0.05 against the
//     exact row-0 weight (1 - r^2)/2 = 0.095 at mu = 1.8 — 1.9x of
//     margin.
//   * Absent (coarse, mu >= 2.2): W_end <= sum_k w_k|K(eps_k)| + 1e-9,
//     the exact leakage ceiling; no state sits inside the window (bulk
//     gap >= 0.2 = 4 E_win). Printed per point: decades under w_min.
//   * Crossover located: mu_down = max{mu: W_end >= w_min}, mu_up =
//     min{mu: W_end <= 1e-6}, the ceiling pinned at w_min/5e4, between
//     the absent-side ceilings (printed, <= 1.3e-8 over the sweep) and
//     the zone floor (bulk leakage at the closed gap, measured >=
//     7e-4). The bracket must contain 2t with half-width <= 0.10 per
//     the budget above; measured [1.88, 2.06], half-width 0.09.
//   * Bulk gap closing, same window: W_bulk stays under the same exact
//     leakage ceiling on the open side (|mu - 2| >= 0.15) and jumps at
//     mu = 2t — gated W_bulk(2t) >= 1e-6 and >= 100x the open-side
//     maximum (expected ~0.016, the in-window bulk levels at ~2/L per
//     level), with the argmax mu* of W_bulk within 0.05 of 2t. The
//     window is the E_win of the zero-mode gates: the same resolution.
//
// Record of the run: the first sweep failed the crossover half-width
// gate as first written (bound 0.075, measured 0.090). Cause: a
// factor-2 error in this header's end-weight expectation — it quoted
// 1 - r^2 where the split pair carries (1 - r^2)/2 on row 0 (each of
// psi± = (gamma_L ± i gamma_R)/sqrt(2) carries A^2/4), which the ED
// in-window weight confirms to the printed digits at every point and
// which R3's w_ED = 1/2 already encoded. The down crossing therefore
// sits at 1.897, not 1.949. No physics or machinery discrepancy: the
// Kill 2 cross-gate passed at 7.6e-15 on the same run, and no bound
// that compares KPM against the oracle was moved. The resolution
// statement above carries the corrected, hard-bounded arithmetic.
//
// Every check prints its measured value; any failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

// K(eps_k), the Jackson-broadened weight of a unit delta at eps_k
// inside |eps| <= w, th0 = acos(w); derivation in the test_r3.cpp
// header. T_n by the recurrence, never cos(n acos eps).
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

// The same window weight assembled from a computed moment sequence.
double kpm_window_weight(const std::vector<double>& mu, double th0,
                         const std::vector<double>& g) {
  const int n_cheb = static_cast<int>(mu.size()) - 1;
  double w = mu[0] * (std::numbers::pi - 2.0 * th0) / std::numbers::pi;
  for (int n = 2; n <= n_cheb - 1; n += 2)
    w -= (4.0 / std::numbers::pi) * g[n] * mu[n] * std::sin(n * th0) / n;
  return w;
}

// Window weight by the midpoint rule in theta = arccos(E/a) of an
// ldos() grid over |E| <= e_win; the independent path from moments to
// window weight (derivation and motivation in the test_r3.cpp helper).
// Measured quadrature error at this file's settings (n_cheb = n_grid =
// 16384): 3.3e-9 worst over mu in {1.0, 2.0, 2.6} and both rows; the
// 1e-7 pin carries 30x margin.
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

double spectral_bound(const majorana::KitaevParams& kp) {
  const double s = 2.0 * kp.t + std::fabs(kp.mu);
  return std::sqrt(s * s + 4.0 * kp.delta * kp.delta);
}

struct SweepPoint {
  double mu, a, th0;
  double w_end = 0.0, w_bulk = 0.0;     // KPM window weights
  double pred_end = 0.0, pred_bulk = 0.0; // ED+kernel predictions (pins)
  double leak_end = 0.0, leak_bulk = 0.0; // exact ceilings sum w|K| + 1e-9
  double e1 = 0.0, e2 = 0.0;            // two smallest |E| (the pair)
  double w_pair_end = 0.0;              // ED in-window end weight
  double e_gap_ed = -1.0;               // ED bulk-row support edge
  bool coarse = false;
};

} // namespace

int main() {
  using namespace majorana;
  char buf[256];
  const int L = 256;
  const int bulk_row = L / 2;
  const int n_cheb = 16384;
  const double e_win = 0.05;
  const double w_min = 0.05;
  const double cross_ceil = 1e-6;

  std::vector<SweepPoint> pts;
  for (int i = 0; i <= 30; ++i) {
    SweepPoint sp;
    sp.mu = i / 10.0;
    sp.coarse = true;
    pts.push_back(sp);
  }
  for (int k = 0; k <= 18; ++k) {
    SweepPoint sp;
    sp.mu = 1.82 + 0.02 * k;
    pts.push_back(sp);
  }
  std::sort(pts.begin(), pts.end(),
            [](const SweepPoint& p, const SweepPoint& q) {
              return p.mu < q.mu;
            });

  const std::vector<double> g = jackson_kernel(n_cheb);
  std::printf("[info] sweep: %zu points, L = %d, n_cheb = %d, E_win = %.3f; "
              "dE_KPM = pi a/n_cheb in [%.3e, %.3e]\n",
              pts.size(), L, n_cheb, e_win,
              std::numbers::pi * spectral_bound(KitaevParams{L, 1, 1, 0}) /
                  n_cheb,
              std::numbers::pi * spectral_bound(KitaevParams{L, 1, 1, 3}) /
                  n_cheb);
  std::printf("[info]   mu      W_end      W_pred_end   |E_pair| min,max   "
              "w_pair_end(ED)  W_bulk     E_gap_ED\n");

  double xgate_end = 0.0, xgate_bulk = 0.0;
  double xquad = 0.0;
  int n_quad = 0;
  for (auto& sp : pts) {
    const KitaevParams kp{L, 1.0, 1.0, sp.mu};
    const CsrMatrix H = kitaev_chain(kp);
    sp.a = spectral_bound(kp);
    sp.th0 = std::acos(e_win / sp.a);
    const Eigenpairs es = exact_eigenpairs(H);
    const int n = 2 * L;

    // ---- Pins: kernel arithmetic and the ED oracle only. ----
    double e1 = 1e30, e2 = 1e30, e_gap = 1e30;
    for (int k = 0; k < n; ++k) {
      const double E = es.values[k];
      const double kw = window_weight(E / sp.a, sp.th0, g);
      const double we = es.vectors(0, k) * es.vectors(0, k);
      const double wb = es.vectors(bulk_row, k) * es.vectors(bulk_row, k);
      sp.pred_end += we * kw;
      sp.pred_bulk += wb * kw;
      sp.leak_end += we * std::fabs(kw);
      sp.leak_bulk += wb * std::fabs(kw);
      if (std::fabs(E) <= e_win)
        sp.w_pair_end += we;
      if (wb >= 1e-6)
        e_gap = std::min(e_gap, std::fabs(E));
      if (std::fabs(E) < e1) {
        e2 = e1;
        e1 = std::fabs(E);
      } else if (std::fabs(E) < e2) {
        e2 = std::fabs(E);
      }
    }
    sp.leak_end += 1e-9;
    sp.leak_bulk += 1e-9;
    sp.e1 = e1;
    sp.e2 = e2;
    sp.e_gap_ed = e_gap;

    // ---- The KPM side, after the pins. ----
    const LdosParams lp{n_cheb, sp.a, 0.0};
    const std::vector<double> mu_end = local_moments(H, 0, lp);
    const std::vector<double> mu_blk = local_moments(H, bulk_row, lp);
    sp.w_end = kpm_window_weight(mu_end, sp.th0, g);
    sp.w_bulk = kpm_window_weight(mu_blk, sp.th0, g);
    xgate_end = std::max(xgate_end, std::fabs(sp.w_end - sp.pred_end));
    xgate_bulk = std::max(xgate_bulk, std::fabs(sp.w_bulk - sp.pred_bulk));

    // Quadrature cross-gate at one coarse point per sweep regime:
    // topological (mu = 1.0, end row), closed gap (mu = 2.0, bulk row),
    // trivial (mu = 2.6, end row). The coarse flag excludes the fine
    // bracket's own mu = 2.0.
    if (sp.coarse && (sp.mu == 1.0 || sp.mu == 2.0 || sp.mu == 2.6)) {
      const int n_grid = 16384;
      const bool end_row = (sp.mu != 2.0);
      const auto& moments = end_row ? mu_end : mu_blk;
      const double w = end_row ? sp.w_end : sp.w_bulk;
      const auto dos = ldos(moments, g, sp.a, 0.0, n_grid);
      const double dq =
          std::fabs(trapz_window(dos, sp.a, e_win, n_grid) - w);
      xquad = std::max(xquad, dq);
      ++n_quad;
      std::printf("[info] quadrature at mu = %.1f, %s row: |trapz - "
                  "analytic| = %.3e\n",
                  sp.mu, end_row ? "end " : "bulk", dq);
    }

    std::printf("[info] %5.2f  %.4e  %.4e    %.3e %.3e   %.4e      %.4e  "
                "%.4f\n",
                sp.mu, sp.w_end, sp.pred_end, sp.e1, sp.e2, sp.w_pair_end,
                sp.w_bulk, sp.e_gap_ed);
  }

  std::snprintf(buf, sizeof buf,
                "worst |W - W_pred| over %zu points: end %.3e, bulk %.3e "
                "(bound 1e-9)",
                pts.size(), xgate_end, xgate_bulk);
  check(xgate_end <= 1e-9 && xgate_bulk <= 1e-9,
        "Kill 2 cross-gate: window weights vs ED+kernel prediction", buf);

  std::snprintf(buf, sizeof buf,
                "worst |trapz - analytic| over %d points: %.3e (bound 1e-7)",
                n_quad, xquad);
  check(n_quad == 3 && xquad <= 1e-7,
        "ldos() trapezoid vs analytic window weight", buf);

  // Present below (coarse, mu <= 1.8) and absent above (mu >= 2.2).
  double worst_present = 1e30;
  double mu_present = -1.0;
  double worst_absent = -1e30, worst_absent_pin = 0.0, mu_absent = -1.0;
  int n_present = 0, n_absent = 0;
  for (const auto& sp : pts) {
    if (!sp.coarse)
      continue;
    if (sp.mu <= 1.8 + 1e-9) {
      ++n_present;
      if (sp.w_end < worst_present) {
        worst_present = sp.w_end;
        mu_present = sp.mu;
      }
    }
    if (sp.mu >= 2.2 - 1e-9) {
      ++n_absent;
      if (sp.w_end - sp.leak_end > worst_absent) {
        worst_absent = sp.w_end - sp.leak_end;
        worst_absent_pin = sp.leak_end;
        mu_absent = sp.mu;
      }
    }
  }
  std::snprintf(buf, sizeof buf,
                "min W_end %.6f at mu = %.2f over %d points (pin %.2f; exact "
                "row-0 weight (1 - r^2)/2 = 0.095 at mu = 1.8)",
                worst_present, mu_present, n_present, w_min);
  check(n_present == 19 && worst_present >= w_min,
        "zero-mode peak present for mu < 2t", buf);
  std::snprintf(buf, sizeof buf,
                "worst W_end - ceiling %.3e at mu = %.2f over %d points "
                "(ceiling there %.3e)",
                worst_absent, mu_absent, n_absent, worst_absent_pin);
  check(n_absent == 9 && worst_absent <= 0.0,
        "zero-mode peak absent for mu > 2t", buf);

  // Crossover bracket from the full scan.
  double mu_down = -1.0, mu_up = -1.0;
  for (const auto& sp : pts) {
    if (sp.w_end >= w_min)
      mu_down = sp.mu;
    if (mu_up < 0.0 && sp.w_end <= cross_ceil)
      mu_up = sp.mu;
  }
  std::snprintf(buf, sizeof buf,
                "bracket [%.2f, %.2f], half-width %.3f (bound 0.10; "
                "derivation in the header)",
                mu_down, mu_up, 0.5 * (mu_up - mu_down));
  check(mu_down > 0.0 && mu_up > 0.0 && mu_down < 2.0 && mu_up > 2.0 &&
            0.5 * (mu_up - mu_down) <= 0.10,
        "crossover bracket contains 2t within the stated resolution", buf);

  // Bulk gap closing in the same window.
  double open_max = 0.0, w_bulk_at_2 = -1.0, mu_star = -1.0, w_star = -1.0;
  double worst_open = -1e30, open_pin = 0.0, mu_open = -1.0;
  for (const auto& sp : pts) {
    if (std::fabs(sp.mu - 2.0) < 1e-9)
      w_bulk_at_2 = sp.w_bulk;
    if (sp.w_bulk > w_star) {
      w_star = sp.w_bulk;
      mu_star = sp.mu;
    }
    if (std::fabs(sp.mu - 2.0) >= 0.15 - 1e-9) {
      open_max = std::max(open_max, sp.w_bulk);
      if (sp.w_bulk - sp.leak_bulk > worst_open) {
        worst_open = sp.w_bulk - sp.leak_bulk;
        open_pin = sp.leak_bulk;
        mu_open = sp.mu;
      }
    }
  }
  std::snprintf(buf, sizeof buf,
                "worst W_bulk - ceiling %.3e at mu = %.2f (ceiling there "
                "%.3e)",
                worst_open, mu_open, open_pin);
  check(worst_open <= 0.0,
        "open side (|mu - 2t| >= 0.15): bulk window weight under the exact "
        "leakage ceiling",
        buf);
  std::snprintf(buf, sizeof buf,
                "W_bulk(mu = 2t) %.6e vs open-side max %.3e (pins: >= 1e-6 "
                "and >= 100x open side)",
                w_bulk_at_2, open_max);
  check(w_bulk_at_2 >= 1e-6 && w_bulk_at_2 >= 100.0 * open_max,
        "bulk gap closed at mu = 2t within the window resolution", buf);
  std::snprintf(buf, sizeof buf,
                "argmax of W_bulk at mu* = %.2f (bound |mu* - 2t| <= 0.05)",
                mu_star);
  check(mu_star > 0.0 && std::fabs(mu_star - 2.0) <= 0.05,
        "gap-closing locator within the stated resolution", buf);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
