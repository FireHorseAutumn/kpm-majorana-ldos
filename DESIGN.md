# DESIGN.md — kpm-majorana-ldos

Design record for Kitaev-chain edge LDOS and the zero-mode peak by KPM.
Status: R1 done (2026-11-04) — the chain builder and the Eigen ED
oracle are implemented; the R2/R3 interfaces remain throwing stubs. The
R1 gates pass at L = 32, t = 1: PHS worst 1.2e-14 (bound 1e-12); two
exact zero modes at the sweet spot, |E| <= 4.4e-17, end-site weight
1 - 4.4e-16 (bound 1e-12); min|E| = 1.0239 at mu = 3 against the bulk
gap; Tr H = 0 exactly (bound 1e-13). The gap bound is derived from the
built chain's Bloch Hamiltonian — E(k) = ±sqrt((2t cos k + mu)^2 +
4 delta^2 sin^2 k), which at delta = t is monotonic in cos k with the
minimum at k = pi, giving E_gap = |mu - 2t| = 1 — not quoted from the
section 7 sources, which remain unread. The open chain misses k = pi
and sits 2.4% above the bound. This note fixes the scope, the oracles,
the rung decomposition with pre-committed gates, and the kill criteria,
so each rung starts in five minutes.

Ordering: third of the three KPM extensions — the for-fun one, and the
cheapest physics in the portfolio (section 5).

## 1. Physics scope

Kitaev chain: spinless fermions on L open sites, hopping t, imposed
p-wave pairing Delta, chemical potential mu. Mean-field BdG, a 2L x 2L
Nambu matrix, open boundaries. The topological phase carries unpaired
Majorana end modes; the trivial phase does not. The computed objects are
local densities of states: at the end sites and at bulk sites, by KPM
with site start vectors — local moments mu_n(i) = <e_i|T_n(H')|e_i>,
deterministic, no stochastic trace.

The figure the rungs build toward: end-site LDOS in the topological
phase showing a zero-energy peak that bulk-site LDOS does not.

Exact anchors, to verify from Kitaev's paper (section 7) before any
number quotes them: topological for |mu| < 2t and trivial above; the
bulk gap closes at |mu| = 2t; the sweet spot Delta = t, mu = 0, where
the zero mode sits exactly on the end sites; away from the sweet spot
the end modes split by ~e^{-L/xi} with localization length xi.

## 2. Oracles

1. **Exact diagonalization.** Dense ED (Eigen) of the 2L matrix. L = 512
   is a 1024 x 1024 dense solve — cheap on host. The oracle for spectra
   and eigenvectors at every gated size.
2. **Analytic end-mode signature.** At the sweet spot the zero mode is
   exact: E = 0, weight fully on the end sites. Away from it:
   localization length xi and splitting ~e^{-L/xi}. The bulk gap closing
   at |mu| = 2t. All to verify from section 7 before use.
3. **Particle-hole symmetry.** The BdG spectrum is symmetric about 0
   exactly — the cheapest gate, run first.

## 3. Rungs and gates

Each rung is one sitting. Gates are pre-committed; tolerances are not
adjusted after a run.

- **R1 — chain builder + ED oracle (host).** `kitaev_chain` and
  `exact_eigenvalues`. Gates at L = 32, t = 1: PHS
  max|E_k + E_{2L-1-k}| < 1e-12; at the sweet spot (Delta = 1, mu = 0) a
  zero mode with |E| < 1e-12 whose weight on the two end sites is
  >= 1 - 1e-12; in the trivial phase (mu = 3) the minimum |E| bounded
  away from 0 by the analytic bulk gap, the bound pinned from the source
  at R1; Tr H = 0 to 1e-13.
- **R2 — local-moment machinery (host).** `local_moments` from the start
  vector e_i. Gates at L = 32: moments against the ED-derived local
  moments sum_k |psi_k(i)|^2 T_n(E'_k), linf < 1e-13; LDOS summed over
  all 2L rows equals the total DOS, gated against the ED histogram under
  the kpm-tightbinding L1 regression discipline.
- **R3 — the zero-mode figure (host).** L ~ 1000, n_cheb of a few
  thousand — 1D, seconds. Gates: the end-site LDOS peak sits within one
  grid bin of E = 0; the integrated peak weight matches the ED zero-mode
  end weight within a pre-committed bound — the Jackson kernel broadens
  the delta into a known shape, so weight is gated and height is
  reported as resolution-limited, never gated; bulk-site LDOS weight
  inside the gap stays below the leakage bound computed from the kernel
  at the chosen n_cheb.
- **R4 — the transition sweep (host).** mu swept across 2t at fixed
  t = Delta. Gates: the zero-mode peak is present for |mu| < 2t and
  absent above, the crossover located within the stated resolution; the
  bulk gap closing is visible within the same resolution.
- **R5 — optional, later, unscheduled:** a small 2D topological BdG
  island with edge modes. Chains never need the GPU; 2D might. Decided
  if opened.

## 4. Kill criteria

- **Kill 1 — resolution vs splitting.** Away from the sweet spot the end
  modes split by ~e^{-L/xi}; below the KPM resolution a "zero mode" is
  indistinguishable from a split pair. If the interesting mu regime puts
  the splitting below resolution at achievable (L, n_cheb), the claim
  closes to "zero mode resolved at the sweet spot," recorded plainly
  with the resolution budget.
- **Kill 2 — oracle disagreement.** ED and KPM disagree beyond the
  pre-committed tolerance and no bug is found: stop. No tolerance
  tuning.
- **Scope guard.** Braiding, fusion rules, Majorana-operator algebra: a
  different repo, if ever. This one computes LDOS.

## 5. Cost notes

The cheapest physics in the portfolio. 1D, 2L dimension, a CSR with a
handful of diagonals: L = 1000 with 4000 moments is seconds of host
time, one start vector, no stochastic average. R1–R4 are one sitting
each. No GPU rung is planned for 1D.

## 6. Relationship to kpm-tightbinding

Reuses its core: the KPM recurrence, the Jackson machinery, the CSR
contract. LDOS replaces the stochastic trace with site start vectors — a
small extension of that core's pattern, not a new method. Same stance as
the sibling designs: nothing copied, vendored, or modified there.
`majorana::CsrMatrix` carries the same layout. R1 decision
(2026-11-04): no integration yet — R1 is the chain builder plus the
Eigen ED oracle and touches nothing in the KPM core, so there is no
surface to depend on or reimplement. The choice binds at R2, when
`local_moments` needs the recurrence: default a CMake build dependency
against a pinned kpm-tightbinding commit, fallback a local
reimplementation of the small CSR surface under cross-gate. That repo is
never modified by this work.

## 7. Primary sources needed before physics claims

To fetch and verify before any number derived from them is quoted:

- Kitaev, *Unpaired Majorana fermions in quantum wires*, Phys. Usp. 44,
  131 (2001) — the phase boundary, the sweet spot, the localization
  length. The section 1 anchors are placeholders until this is on hand
  and read.
- Weisse, Wellein, Alvermann, Fehske, *The kernel polynomial method*,
  Rev. Mod. Phys. 78, 275 (2006) — the local-moment/LDOS conventions
  and the kernel-resolution bounds.
- A standard BdG reference for the Nambu conventions (shared with the
  kpm-bdg reading list).

## 8. Honest limitations

Mean-field, noninteracting, imposed pairing. No braiding, no fusion, no
topological-computation claim. Open chains; ED gates at moderate L; the
zero-mode claim is exactly as strong as the resolution budget of
section 4.
