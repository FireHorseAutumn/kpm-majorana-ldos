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
integrate. The decision binds at R2, when the local moments need the
KPM recurrence; the default is a CMake build dependency against a pinned
commit.

## Status

R1 done (2026-11-04). `kitaev_chain` builds the 2L x 2L Nambu CSR;
`exact_eigenpairs` / `exact_eigenvalues` are the Eigen dense ED oracle.
The DESIGN.md section-3 gates pass at L = 32, t = 1: PHS worst 1.2e-14
(bound 1e-12); two exact sweet-spot zero modes, |E| <= 4.4e-17, with
end-site weight 1 - 4.4e-16; min|E| = 1.0239 at mu = 3 against the
derived bulk gap |mu - 2t| = 1; Tr H = 0 exactly. CTest 2/2. The R2/R3
interfaces (local moments, Jackson, LDOS) remain throwing stubs; the
stub CTest checks they still throw. The design record — chain scope,
the ED and analytic end-mode oracles, the gated rung decomposition,
kill criteria, and the primary-source reading list — is
[`DESIGN.md`](DESIGN.md).
