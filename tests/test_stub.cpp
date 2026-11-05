// Scaffold acceptance test: the interface stubs carry the documented
// defaults, and the R2/R3 entry points throw std::logic_error until
// their rungs land. R1 (kitaev_chain, exact_eigenvalues) is implemented
// and gated in test_r1.cpp. Every check prints its measured value; any
// failure exits nonzero.

#include "majorana/kitaev.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

template <typename F>
bool throws_not_implemented(F&& f) {
  try {
    f();
  } catch (const std::logic_error& e) {
    return std::string(e.what()).find("not implemented") != std::string::npos;
  }
  return false;
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

  check(throws_not_implemented([&] { majorana::local_moments(empty, 0, lp); }),
        "local_moments stub throws", "DESIGN.md R2");
  check(throws_not_implemented([&] { majorana::jackson_kernel(16); }),
        "jackson_kernel stub throws", "DESIGN.md R3");
  check(throws_not_implemented([&] { majorana::ldos({}, {}, 1.0, 0.0, 16); }),
        "ldos stub throws", "DESIGN.md R3");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
