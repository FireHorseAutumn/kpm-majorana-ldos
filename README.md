# kpm-majorana-ldos

Edge local density of states and the zero-mode peak of a Kitaev chain —
and later, optionally, a small 2D topological BdG island — computed with
the Kernel Polynomial Method, validated against exact diagonalization
and against the analytic end-mode signature.

Boundary, stated once: mean-field BdG with imposed pairing, like every
Kitaev-chain numerics exercise. No braiding, no fusion rules, no
topological-computation claim lives here.

## Relationship to kpm-tightbinding

Reuses the core of
[kpm-tightbinding](https://github.com/FireHorseAutumn/kpm-tightbinding):
the KPM recurrence, the Jackson machinery, the CSR contract. LDOS
replaces the stochastic trace with site start vectors — a small
extension of that core's pattern, not a new method. Nothing is copied,
vendored, or modified there. R1 took no dependency: the chain builder
and the ED oracle touch nothing in that core, so there was nothing to
integrate. R2 recorded the fallback: that repo's CMake forces CUDA and
exports nothing, so the CSR surface stays layout-identical and local —
the integration decision is deferred, not made (DESIGN.md section 6).

## Status

R4 done (2026-11-05). The transition sweep at L = 256, t = delta = 1,
mu across 2t (coarse 0.0:0.1:3.0 plus a fine bracket 1.82:0.02:2.18,
50 points; window |E| <= 0.05, n_cheb = 16384): the zero-mode peak is
present for mu <= 1.8 (min end-site window weight 0.095 — the exact
row-0 pair weight (1 - r^2)/2 — against a 0.05 pin) and absent for
mu >= 2.2 (under the exact kernel-leakage ceilings); the crossover is
bracketed at [1.88, 2.06], containing 2t with half-width 0.090 against
the hard-derived 0.10 resolution; the bulk gap closing is visible in
the same window — bulk-site window weight 7.8e-3 at mu = 2t against
5.7e-10 worst on the open side, the closing argmax at mu* = 1.98. The
Kill 2 cross-gate against the ED+kernel prediction passes at 7.6e-15
over both rows and all 50 points (bound 1e-9). The first run failed
the crossover half-width as first stated (0.075): the header's
end-weight expectation was off by a factor 2, confirmed by the ED
column and re-derived hard — no oracle-comparison bound moved (record
in the test_r4.cpp header). R3 (2026-11-05). The full pipeline is implemented: `local_moments`
(site-start Chebyshev recurrence, deterministic), `jackson_kernel`,
`ldos` (Jackson-damped reconstruction on an arccos-uniform grid). The
zero-mode figure gates pass at L = 1000, sweet spot, n_cheb = 4000: the
end-site peak sits within one grid bin of E = 0; the integrated peak
weight matches the ED zero-mode end weight (1/2, the R1 exact end
Majoranas) to 4.0e-9 against a pin of 4.2e-9 derived from the kernel's
known broadening before the KPM run; the bulk-site in-gap weight is
1.5e-10 against the pinned kernel-leakage bound 2.5e-10. Heights are
resolution-limited and reported, never gated: zero peak 90.61 vs bulk
band 64.08 on the end site; the bulk site shows no zero peak. R2
(2026-11-05): moments vs the ED-derived sum_k |psi_k(i)|^2 T_n(E'_k),
worst linf 2.275e-14 (bound 1e-13) at L = 32. R1 (2026-11-04):
`kitaev_chain` builds the 2L x 2L Nambu CSR; `exact_eigenpairs` /
`exact_eigenvalues` are the Eigen dense ED oracle; PHS worst 1.2e-14
(bound 1e-12); two exact sweet-spot zero modes, |E| <= 4.4e-17, with
end-site weight 1 - 4.4e-16; min|E| = 1.0239 at mu = 3 against the
derived bulk gap |mu - 2t| = 1; Tr H = 0 exactly. CTest 4/4, zero
warnings. The design record — chain scope, the ED and analytic end-mode
oracles, the gated rung decomposition, kill criteria, and the
primary-source reading list — is [`DESIGN.md`](DESIGN.md).
