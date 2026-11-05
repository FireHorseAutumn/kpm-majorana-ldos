// Scaffold acceptance test: the interface carries the documented
// defaults. R1 (kitaev_chain, exact_eigenvalues), R2 (local_moments),
// and R3 (jackson_kernel, ldos) are implemented and gated in
// test_r1.cpp, test_r2.cpp, and test_r3.cpp. Every check prints its
// measured value; any failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <cstdio>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

} // namespace

int main() {
  char buf[256];

  const majorana::KitaevParams kp;
  std::snprintf(buf, sizeof buf, "L=%d t=%.1f delta=%.1f mu=%.1f", kp.L, kp.t,
                kp.delta, kp.mu);
  check(kp.L == 64 && kp.t == 1.0 && kp.delta == 1.0 && kp.mu == 0.0,
        "KitaevParams defaults (the sweet spot)", buf);

  const majorana::LdosParams lp;
  std::snprintf(buf, sizeof buf, "n_cheb=%d a=%.2f b=%.1f", lp.n_cheb, lp.a,
                lp.b);
  check(lp.n_cheb == 2048 && lp.a == 0.0 && lp.b == 0.0,
        "LdosParams defaults (a=0: unset, caller supplies the bound)", buf);

  const majorana::CsrMatrix empty;
  std::snprintf(buf, sizeof buf, "n=%d nnz=%d", empty.n, empty.nnz());
  check(empty.n == 0 && empty.nnz() == 0, "CsrMatrix empty state", buf);

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
