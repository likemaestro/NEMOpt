# Stability (buckling/bifurcation) diagnostics for NEMOpt: implementation, verification, architecture and cost

Scope: how to add microscopic (Bloch–Floquet), cell-periodic (k = 0) and macroscopic (strong-ellipticity) stability diagnostics to the planned NEMOpt package without disturbing the FEAD2026 reproduction profile. The notes build on the refactoring report (`reports/NEMOpt refactoring strategy literature.md`, NDR-0001 to NDR-0044) and propose material for NDR-0045 onward.

Evidence labels used below:
- **[PR]** peer-reviewed paper (bibliographic data checked against Crossref on 2026-10-12 unless marked otherwise).
- **[Doc]** official software documentation or source code.
- **[Exp]** a small experiment run in this study (Julia 1.11.7, SparseArrays CHOLMOD `BUILD_VERSION` 5.2.0, one BLAS thread, Linux container). The scripts (`exp1.jl`, `exp2.jl`, `exp3.jl`, `asm.jl`, `lanczos.jl`; stdlib only) are preserved in [`evidence/stability_experiments/`](evidence/stability_experiments/), with the recorded output of `exp1.jl` in `exp1.out`. Run each from that folder with Julia 1.11.7. Each experiment is also described here well enough to be re-created.
- **[Slides]** a non-peer-reviewed presentation.
- **[Preprint]** an arXiv preprint.
- **Inference** marks reasoning from cited facts; **Unverified** marks claims recalled but not checked in this session.

Code facts (verified in the working tree; `git diff 30d0842 HEAD` shows no change to `Full_Domain/` or `Domain_Reduction/`):
- The full-cell solver fixes the last node, the center node, by dropping its two DOFs (`Kff = K[1:end-2, 1:end-2]`) ([NEMOpt.jl L342](../../Full_Domain/NEMOpt.jl#L342)).
- The reduced tangent is formed as `vcat(hcat(Kmm + Kms*D1 + D1'*Ksm + D1'*Kss*D1, …), …)` and solved with `\` ([NEMOpt.jl L357–367](../../Full_Domain/NEMOpt.jl#L357)).
- D1 has unit entries mapping slaves to masters. D2 rows encode each slave's lattice offset (L1, L2 or L1 + L2) ([NEMOpt.jl L195–203](../../Full_Domain/NEMOpt.jl#L195)).
- At 80×80 the reduced system has 2·(159 + 6240) = 12,798 DOFs; at 100×100 it has 2·(199 + 9800) = 19,998 (inference from the numbering code).
- The half cell has 6 prescribed nodes and a D1 with −1 entries, i.e., point-symmetric (odd-fluctuation) constraints ([NEMOpt-DR.jl L242–299](../../Domain_Reduction/NEMOpt-DR.jl#L242)).
- The constitutive law is P = μF − (μ − λ ln J)F⁻ᵀ, i.e., W = μ/2(I₁ − 2) − μ ln J + λ/2 (ln J)². The tangent is 𝔸_ijkl = μδ_ikδ_jl + λF⁻ᵀ_ijF⁻ᵀ_kl + (μ − λ ln J)F⁻ᵀ_ilF⁻ᵀ_kj. Low-density elements blend to the small-strain linear law ℂ ([NEMOpt.jl L281–299](../../Full_Domain/NEMOpt.jl#L281)).
- The article itself names "instabilities that might arise particularly under compressive macroscopic loading" as open future work (article p. 16, [Güven & Özdemir 2026](https://doi.org/10.1016/j.finel.2026.104599)).
- Every published case is tensile or shear: §5.2 uses F_M = 1.075·I and F_M = [1.025 0.025; 0.025 1.025] (refactoring report, "same results" section). The legacy Newton scale `Rep_MU = |max(F_M − I)|·Lx` is 0 for purely compressive F_M and then loops forever ([NEMOpt.jl L375, L388](../../Full_Domain/NEMOpt.jl#L375)). **Compressive benchmarks therefore cannot run in the FEAD2026 profile at all**; they need the NDR-0031/0032 Newton-control options.

---

## Q1. Inertia and positive-definiteness checks with sparse direct solvers

### Takeaway
CHOLMOD's supernodal Cholesky failure (`issuccess == false`) is a reliable binary "not positive definite" flag, but it says nothing about how many eigenvalues are negative. A singular positive-semidefinite matrix can still pass. CHOLMOD's LDLᵀ is simplicial and **unpivoted**, so its negative-pivot count gives the inertia by Sylvester's law only when the factorization does not break down or grow; it carries no accuracy guarantee for indefinite matrices. It is also 10–21× slower than Cholesky at NEMOpt sizes. The low-risk design is:
1. Cholesky (or a Cholesky attempt) on every converged tangent, as a free flag.
2. Shifted-Cholesky bisection or shift-invert Lanczos for margins.
3. An LDLᵀ negative-pivot count only on failure, cross-checked by a pivoted solver (dense Bunch–Kaufman on small cases; MUMPS with sequential root, or MA57/Pardiso, through an optional extension).

### Cited Findings
**CHOLMOD semantics [Doc]**
- CHOLMOD "includes a non-supernodal LDLᵀ factorization method that can factorize symmetric indefinite matrices if all of their leading submatrices are well-conditioned (D is diagonal)" — [CHOLMOD User Guide, abstract (SuiteSparse dev branch, fetched 2026-10-12)](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CHOLMOD/Doc/CHOLMOD_UserGuide.tex).
- "CHOLMOD does not include a supernodal LDLᵀ factorization, so if a supernodal factorization is selected, it will be in the form LLᵀ. The LDLᵀ method can be used to factorize positive definite matrices and indefinite matrices whose leading minors are well-conditioned (2-by-2 pivoting is not supported)." To factorize a large indefinite matrix one must force the simplicial LDLᵀ, which "will be significantly slower than a supernodal LLᵀ factorization" — [CHOLMOD User Guide, `cholmod_analyze`/`cholmod_factorize` sections](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CHOLMOD/Doc/CHOLMOD_UserGuide.tex).
- For `cholmod_rowfac` (the simplicial kernel): "No 2-by-2 pivoting (or any other pivoting) is performed." The LDL code comments add that for indefinite A "no guarantee of accuracy is provided, since no numeric pivoting is performed" — [CHOLMOD User Guide, `cholmod_rowfac`](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CHOLMOD/Doc/CHOLMOD_UserGuide.tex); [LDL user guide mirror](https://svn.r-project.org/R-packages/branches/Matrix-mer2/inst/doc/UFsparse/ldl_userguide.pdf). This **resolves the [UNVERIFIED] item** in `efficient_fe_adjoint.md`: CHOLMOD's LDLᵀ is unpivoted.
- On failure: "If the matrix is not positive definite the routine returns TRUE, but Common->status is set to CHOLMOD_NOT_POSDEF and L->minor is set to the column at which the failure occurred. Columns L->minor to L->n-1 are set to zero." The factorization routines accept a `beta*I` shift ("only the real part, beta[0], is used"). Complex matrices are supported, with "most methods treat[ing] the matrix as Hermitian" — [CHOLMOD User Guide](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CHOLMOD/Doc/CHOLMOD_UserGuide.tex).
- `cholmod_dbound` can bound the magnitude of D entries ("ensures that the magnitude of the entries of D are greater than or equal to Common->dbound"). A non-zero `dbound` would silently alter the pivots, so it must stay at its default for inertia counting (inference from the same guide).

**Julia wrappers [Doc]**
- Julia's `ldlt` for sparse matrices forces the simplicial path with `@cholmod_param final_ll = false` and `@cholmod_param supernodal = 0` ("Really make sure it's an LDLt by avoiding supernodal factorization") — [SparseArrays `cholmod.jl` L1661–1663 in the Julia 1.11.7 stdlib](https://github.com/JuliaSparse/SparseArrays.jl/blob/main/src/solvers/cholmod.jl).
- The `ldlt` methods available in Julia 1.11.7 are CHOLMOD sparse (real and complex) and `SymTridiagonal` only. Dense Bunch–Kaufman is `bunchkaufman`, and `LinearAlgebra.inertia(::BunchKaufman; atol, rtol)` exists, with tolerance-based zero classification — [Exp, `methods(ldlt)`, `methods(LinearAlgebra.inertia)` in 1.11.7]. The CHOLMOD `D` component is extracted with `diag(sparse(F.LD))`; `sparse(F.D)` throws `CanonicalIndexError` [Exp].
- `cholesky(A; check=true)` on an indefinite matrix throws `PosDefException` with the message "matrix is not positive definite; Factorization failed." and no column index [Exp].

**Experiments [Exp]**
The test matrix is the plane-strain Q4 stiffness of a periodic n×n cell, E = 100, ν = 0.3, assembled with the same constraint structure as NEMOpt (periodic identification plus a pinned center node).
- `exp2.jl`, heterogeneous densities 0.02–1.0, n = 14, 390 DOFs, shifts σ placed in 28 well-separated spectral gaps (1 to 370 eigenvalues below σ):
  - `cholesky(...; shift=-σ, check=false)` failed on all 28 indefinite shifts.
  - The `ldlt` negative-pivot count equaled the exact count in **28/28** cases, and still did with σ within 1e−13 relative of an eigenvalue.
  - Dense Bunch–Kaufman `inertia` reported (385, 4, 1) at the 1e−13 gap because of its zero tolerance.
- `exp1.jl`:
  - `ldlt([0 1; 1 0])` returns `issuccess == false`, a breakdown on a nonsingular indefinite matrix, which is the expected failure mode of unpivoted LDLᵀ.
  - The **unpinned** periodic tangent, which is singular (two translation null vectors with eigenvalues −2.7e−13 and −7e−15), passed both `cholesky` and `ldlt` (min |D| = 4.8e−12). Cholesky success therefore means "no detectably negative eigenvalue", not "positive definite with margin".
- Timing at the legacy reduced sizes (minimum of 3 runs; same sparsity as the NEMOpt reduced tangent):

| Mesh (DOFs) | nnz(K) | nnz(L) | UMFPACK LU | CHOLMOD Cholesky | `cholesky!` numeric refactor | `ldlt` (simplicial) | 1 solve | 4-RHS solve |
|---|---|---|---|---|---|---|---|---|
| 80×80 (12,798) | 179,146 | 1,138,379 | 0.091 s | 0.029 s | 0.018 s | 0.311 s | 1.9 ms | 4.0 ms |
| 100×100 (19,998) | 279,948 | 2,515,224 | 0.181 s | 0.050 s | 0.049 s | 1.047 s | 2.4 ms | 5.4 ms |

**Pivoted alternatives [Doc]**
- MUMPS 5.9.1 user guide (2026-07-20):
  - INFO(12)/INFOG(12) give "number of negative pivots … if SYM=1 or 2". With the default ICNTL(13)=0 this "excludes pivots from the parallel root node treated by ScaLAPACK"; the user "should set ICNTL(13)=1 or use a single processor in order to get the exact number". Null pivots detected under ICNTL(24) are excluded from the count.
  - With SYM=1, "numerical pivoting is switched off", and "the successful factorization of a symmetric matrix with SYM=1 is thus not an indication that the matrix provided was symmetric positive definite … check that the number of negative pivots or inertia (INFOG(12)) is 0".
  - "For the complex version, the value SYM=1 is currently treated as SYM=2. We do not have a version for Hermitian matrices in this release", and "for complex symmetric matrices (SYM=1 or 2), INFO(12) will be 0".
  - Source: [MUMPS 5.9.1 Users' guide](https://mumps-solver.org/doc/userguide_5.9.1.pdf).
- MUMPS.jl: latest v1.7.0, released 2026-10-02; v1.6.0 on 2026-01-19 ([MUMPS.jl releases](https://github.com/JuliaSmoothOptimizers/MUMPS.jl/releases)). Its README describes LDLᵀ/LU support but does not document an inertia accessor ([MUMPS.jl README](https://github.com/JuliaSmoothOptimizers/MUMPS.jl)).
- Pardiso.jl: latest v1.1.2 (2025-12-10). It requires "either having a valid license for Panua Pardiso or having the MKL library installed" ([Pardiso.jl releases](https://github.com/JuliaSparse/Pardiso.jl/releases); [README](https://github.com/JuliaSparse/Pardiso.jl)).
- QDLDL.jl is a factorization "for quasi-definite linear systems" ([QDLDL.jl README](https://github.com/osqp/QDLDL.jl)).
- Latest registered versions in the local General registry snapshot, which is at least as recent as 2026-10-02 because it contains MUMPS 1.7.0: HSL.jl 0.5.6, QDLDL.jl 0.4.1, LDLFactorizations.jl 0.10.2 [Exp, registry tarball].

**Theory [PR]**
- Bunch & Kaufman introduced "stable methods for calculating inertia and solving symmetric linear systems" with 1×1/2×2 pivoting — [Bunch & Kaufman 1977, Math. Comp. 31(137):163–179](https://doi.org/10.1090/S0025-5718-1977-0428694-0).
- Bounded Bunch–Kaufman variants are analyzed by [Ashcraft, Grimes & Lewis 1998, SIAM J. Matrix Anal. Appl. 20:513–561](https://doi.org/10.1137/S0895479896296921).
- Haynsworth inertia additivity, In(K) = In(K₁₁) + In(K/K₁₁) for the Schur complement, is the basis of the condensation in Q3 — [Haynsworth 1968, Linear Algebra Appl. 1:73–81](https://doi.org/10.1016/0024-3795(68)90050-5).
- The CHOLMOD algorithm reference is [Chen, Davis, Hager & Rajamanickam 2008, ACM TOMS 35(3), Algorithm 887](https://doi.org/10.1145/1391989.1391995).

### Inferences
- **Sylvester's law (spectrum slicing).** If K − σI = LDLᵀ exists, the number of negative D entries equals the number of eigenvalues of K below σ (Golub & Van Loan, *Matrix Computations*, 4th ed., 2013, and Parlett, *The Symmetric Eigenvalue Problem*, SIAM 1998; standard results whose bibliographic details were not re-verified here). Three practical consequences for NEMOpt:
  - (a) "Cholesky of K − σI succeeds" ⇔ λ_min(K) > σ, up to round-off, so bisection on σ with `cholesky!` (symbolic analysis reused) brackets λ_min at about 0.02–0.05 s per step at 80×80–100×100.
  - (b) The sign of λ_min is ambiguous within roughly n·u·‖K‖ (`exp1`: a PSD-singular matrix passes). Diagnostics must report a **three-state verdict** (stable / indeterminate / unstable) with an explicit tolerance, not a boolean.
  - (c) A failed Cholesky during Newton iterations at a non-equilibrium iterate is not a stability statement. Only failures at converged equilibria (and in χ/adjoint solves at the converged U) carry physical meaning, so the log must tag the state.
- **The unpivoted LDLᵀ count is a diagnostic, not a proof.** It agreed in all 28 tests, but CHOLMOD documents no accuracy guarantee, and breakdown (as in `ldlt([0 1;1 0])`) or growth is possible. Use it only after a Cholesky failure, require `issuccess`, report min|D| and max|D| as a growth indicator, and cross-check:
  - dense `bunchkaufman` + `inertia` when n ≲ 3,000 (unit and regression tests);
  - MUMPS with ICNTL(13)=1, SYM=2, CNTL(1) > 0, or MA57/Pardiso, in an optional extension, when an exact count matters for production.
- **Complex Hermitian Bloch matrices.** CHOLMOD handles them directly (Q3). MUMPS cannot report inertia for complex matrices. A pivoted real solver can still be used on the real embedding [A −B; B A] of H = A + iB, whose eigenvalues are those of H each counted twice, so the inertia doubles (standard identity).
- **The cost ordering at NEMOpt sizes is decisive**: Cholesky (0.03–0.05 s) ≪ LU (0.09–0.18 s) ≪ simplicial LDLᵀ (0.31–1.05 s). LDLᵀ should therefore never be the routine fallback for the Newton solve: an indefinite converged tangent should fall back to UMFPACK LU for the solve, while LDLᵀ runs only as a logged diagnostic.

### Gaps
- MKL PARDISO reports positive and negative eigenvalue counts in iparm(22)/iparm(23) for symmetric indefinite matrices, and its pivot perturbation may affect those counts (Unverified; recalled, not re-read).
- HSL MA57/MA97 return inertia and need an HSL licence (Unverified).
- Whether MUMPS.jl exposes `INFOG(12)` through a documented accessor was not checked in its source.
- No error bound for the inertia computed by an unpivoted sparse LDLᵀ was found; the Higham (2002, Ch. 10) Cholesky success condition was not re-read.

---

## Q2. Eigenvalue computation in Julia: shift-invert, packages and cost

### Takeaway
For the few smallest eigenvalues of the reduced (real or complex Hermitian) tangent, use shift-invert Lanczos with one reused CHOLMOD factor. KrylovKit.jl is the actively maintained default (v0.10.4, June 2026; Lanczos for Hermitian maps and `geneigsolve` for Hermitian-definite pencils). ArnoldiMethod.jl (v0.4.0, Feb 2024, Krylov–Schur, documented shift-and-invert recipe) is the fallback. Arpack.jl (last GitHub release 2021) should not be a dependency. Measured cost for 6 eigenpairs is about 26–29 solves, i.e., 3–5 factorization-equivalents. The center-node pin injects two spurious near-zero "pinned-translation" eigenvalues, so magnitudes need translation deflation rather than pinning.

### Cited Findings
**Package status [Doc]**

| Package | Latest version (release date) | Relevant capability |
|---|---|---|
| KrylovKit.jl | v0.10.4 (2026-06-23); v0.10.3 (2026-04-29); v0.10.0 (2025-08-02) | `eigsolve` with `:LM/:SR/…` selectors and `EigSorter`; for "a self-adjoint linear map (hermitian or, in the real case, symmetric), the eigenvalues are real"; `geneigsolve` "for symmetric or hermitian generalized eigenvalue problems where the matrix or linear map `B`" is positive definite |
| ArnoldiMethod.jl | v0.4.0 (2024-02-22); the registry snapshot (≥ 2026-10-02) has nothing newer | Matrix-free Krylov–Schur for "general" matrices, "symmetric or non-symmetric and either real or complex"; for Hermitian matrices "the partial Schur form coincides with the partial eigendecomposition"; documents "Shift-and-invert with LinearMaps.jl" and a `ShiftAndInvert(factorize(A), B, …)` example for generalized problems |
| Arpack.jl | Latest GitHub release v0.5.3 (2021-07-01); latest registered v0.5.4 (date not found) | README: "Users running into issues with this package may want to try KrylovKit.jl or ArnoldiMethod.jl" |

Sources: [KrylovKit releases feed](https://github.com/Jutho/KrylovKit.jl/releases); [KrylovKit eig docs](https://github.com/Jutho/KrylovKit.jl/blob/master/docs/src/man/eig.md); [ArnoldiMethod releases feed](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl/releases); [ArnoldiMethod docs](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl/blob/master/docs/src/index.md); [Arpack.jl releases/tags](https://github.com/JuliaLinearAlgebra/Arpack.jl/tags); [Arpack.jl README](https://github.com/JuliaLinearAlgebra/Arpack.jl).

**Algorithms [PR]**
- Krylov–Schur restarting — [Stewart 2002, SIAM J. Matrix Anal. Appl. 23(3):601–614](https://doi.org/10.1137/S0895479800371529).
- Spectral-transformation Lanczos — [Ericsson & Ruhe 1980, Math. Comp. 35:1251](https://doi.org/10.2307/2006390).
- Shifted block Lanczos for sparse symmetric generalized problems (the NASTRAN/industrial standard) — [Grimes, Lewis & Simon 1994, SIAM J. Matrix Anal. Appl. 15(1):228–272](https://doi.org/10.1137/S0895479888151111).
- ARPACK — [Lehoucq, Sorensen & Yang 1998, SIAM](https://doi.org/10.1137/1.9780898719628).

**Cost measurements [Exp]** (`exp1.jl`/`exp3.jl`; naive shift-invert Lanczos with full reorthogonalization and a fixed start vector, 6 eigenpairs to a relative Ritz residual of 1e−8, warm runs)

| Mesh | Real k = 0 | Complex Bloch θ = (0.7, 1.3) | Factorization for comparison |
|---|---|---|---|
| 80×80 | 29 solves, 0.102 s | 26 solves, 0.265 s | real Cholesky 0.029 s; complex Hermitian Cholesky 0.052 s; complex solve 2.5 ms |
| 100×100 | 29 solves, 0.140 s | 26 solves, 0.348 s | real Cholesky 0.050 s; complex Hermitian Cholesky 0.121 s; complex solve 4.5 ms |

- The first call cost 4.2 s, which is JIT compilation; in a long-running process this is a one-time cost.

**Spectrum at k = 0 with the center pin [Exp]** (homogeneous cell, identity metric)
- The six smallest eigenvalues at 80×80 are [0.01172, 0.01172, 0.2371, 0.2371, 0.2615, 0.2615]; at 100×100 [0.00718, 0.00718, 0.1518, 0.1518, 0.1667, 0.1667].
- Multiplied by n², the third and fourth (1517.6 at 80×80, 1517.9 at 100×100) match the continuum value μ(2π)² = 1518.4 of the lowest cell-periodic shear wave.
- The first pair (×n² ≈ 75.0 → 71.8, falling slowly like 1/ln n) are **pinned-translation pseudo-modes**: nearly rigid translations with a point singularity at the pinned node.

### Inferences
- **Recommended solver call.** Use `F = cholesky(Symmetric(K))` (or `Hermitian` for Bloch), then `eigsolve(x -> F \ x, x0, nev, :LM, Lanczos(; krylovdim = 30, tol = 1e-10))` and invert the Ritz values.
  - When K is indefinite, factor K − σI with σ < 0 below the expected λ_min (Cholesky succeeds), or use LU as the operator (Arnoldi).
  - Fix x0 from a private `StableRNG` so the diagnostic never touches the global RNG.
  - This fits the planned extension: KrylovKit as a weak dependency (Q7).
- **Cost model.** One margin computation costs about one factorization plus 26–30 solves, i.e., about 4–6 factorization-equivalents at these sizes. On the coordinator's machine (LU 0.15 s, Cholesky + 4 RHS 0.08 s, i.e., 1.6–2.4× slower than this container), that is about 0.2–0.3 s at 80×80 and 0.3–0.4 s at 100×100 per state.
- **Metric and generalized problems.**
  - Inertia, and hence the stability verdict, is independent of the SPD metric M: K x = λ M x has the same sign pattern as K by congruence.
  - Magnitudes need a mesh-independent normalization. On NEMOpt's uniform periodic grid every reduced node carries the same tributary area h², so a geometric (unit-density) lumped mass is exactly h²·I. The identity metric is therefore already the lumped-mass metric, and λ_min(K)/h² is a mesh-converging energy-per-L² margin; the 1517.6 → 1517.9 → 1518.4 sequence above confirms this.
  - Avoid density-weighted mass: void elements create spurious low-mass modes.
  - No generalized solver is needed. If a non-uniform metric is wanted, scale by M^{−1/2} (diagonal) or use KrylovKit's `geneigsolve`.
- **Do not use the pinned matrix for margins.** For eigenvalue magnitudes at k = 0, use the unpinned periodic tangent with the two translation vectors t₁, t₂ (exact null vectors, since K t = 0) deflated by projection: operate on P K P with P = I − T(TᵀT)⁻¹Tᵀ, or factor K + αTTᵀ. The pinned matrix is still exact for inertia (Q3), but shift-invert would converge first to the two pseudo-modes and waste two eigenpair slots. Under refinement those slots also tend to zero, which would masquerade as a "softening" design.
- **A linearized critical-load estimate is possible with existing kernels.** Q3 recommends direct tangent testing along the load path. As an option, the planned third derivative Z = ∂𝔸/∂F gives the directional derivative of the tangent along the load path, dK/dt = Σ_e ∫ Bᵀ(Z·∇(dU/dt))B. Solving (K + t* dK/dt)φ = 0 by shift-invert then estimates the next critical load factor t* from one converged state, the finite-strain analogue of the linearized buckling of Thomsen et al. 2018. This is a cheap per-evaluation predictor whose accuracy must be validated against stepping and bisection.

### Gaps
- The publication date of Arpack.jl v0.5.4 and the commit activity of ArnoldiMethod.jl after February 2024 could not be retrieved, because the GitHub API was blocked and the registry snapshot has no dates.
- KrylovKit and ArnoldiMethod were not benchmarked (packages not installed). The cost numbers come from a hand-written Lanczos and should be within a small factor (inference).

---

## Q3. Bloch–Floquet (microscopic) stability on the unit cell

### Takeaway
Three tests together are complete in the sense of Geymonat–Müller–Triantafyllidis:
- **Cell-periodic, k = 0.** NEMOpt's existing reduced tangent, with the pinned center, is exactly the k = 0 Bloch stability matrix for inertia.
- **Microscopic, k ≠ 0.** Bloch phases on the existing D1/D2 map, **without** the center pin, give a complex Hermitian tangent.
- **Macroscopic, k → 0.** Strong ellipticity of 𝔸_M (Q4).

The implementation is a small change at scatter time: multiply each slave row of D1 by e^{i(n₁θ₁ + n₂θ₂)}, with (n₁, n₂) read from the D2 row. Time-reversal symmetry halves the zone. The four real "TRIM" points θ ∈ {0, π}² catch all period-doubling modes, including the classic 2×2 hole-array pattern, with real arithmetic. One wave vector costs one complex Hermitian Cholesky (0.05 s at 80×80, 0.12 s at 100×100 in this container).

### Cited Findings
**Theory [PR]**
- Microscopic bifurcation and macroscopic loss of rank-one convexity are the two stability notions for periodic hyperelastic media. Microscopic stability is assessed by Bloch waves over all wave vectors; the long-wavelength limit connects to rank-one convexity of the homogenized tangent — [Geymonat, Müller & Triantafyllidis 1993, Arch. Ration. Mech. Anal. 122(3):231–290](https://doi.org/10.1007/BF00380256).
- Onset of failure is "the first occurrence of a bifurcation in the fundamental solution, for which all cells deform identically", computed by "a Bloch wave analysis on the deformed unit cell". The paper describes "a numerical algorithm which reduces the size of stability matrices by an order of magnitude" — [Triantafyllidis, Nestorović & Schraad 2006, J. Appl. Mech. 73(3):505–515](https://doi.org/10.1115/1.2126695) (abstract via Crossref; issued online 2005).
- In layered hyperelastic composites "the macroscopic stable region consistently envelops the microscopic stable region", and the loss-of-ellipticity condition coincides with long-wavelength buckling — [Triantafyllidis & Maker 1985, J. Appl. Mech. 52(4):794–800](https://doi.org/10.1115/1.3169148), as summarized in [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371) [Preprint]. See also [Triantafyllidis & Schnaidt 1993, JMPS 41(9):1533–1565](https://doi.org/10.1016/0022-5096(93)90039-I).
- Real arithmetic with standard FE codes: Bloch conditions can be imposed by duplicating the model into real and imaginary parts coupled through the phase relations — [Åberg & Gudmundson 1997, J. Acoust. Soc. Am. 102(4):2007–2013](https://doi.org/10.1121/1.419652). Bertoldi's group used the same two-copy coupling in ABAQUS ("one for the real part and one for the imaginary part, coupled by the Bloch type displacement") — [Bertoldi, MSM workshop slides, Univ. Twente, 22 Sept 2008](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) [Slides].
- Topology optimization with Bloch–Floquet linear buckling:
  - Thomsen, Wang & Sigmund use "homogenization, a linearized stability criterion, and Floquet-Bloch theory" to maximize the most critical buckling strength of square, triangular and hexagonal cells under biaxial, uniaxial and shear loads — [Thomsen, Wang & Sigmund 2018, CMAME 339:115–136](https://doi.org/10.1016/j.cma.2018.04.031) (abstract via the [DTU Orbit record](https://orbit.dtu.dk/en/publications/buckling-strength-topology-optimization-of-2d-periodic-materials-/); full text not accessible here).
  - In the follow-up, "both short and long wavelength buckling is captured by employing the Floquet-Bloch boundary conditions". The critical mode is captured "by sweeping k-vectors along the boundaries of the irreducible Brillouin zone", justified only by citation. Spurious void modes are suppressed by interpolating the stress stiffness as ρ̄^p E₁ without the E₀ offset used in K₀ (p = 3, E₀ = 1e−5 E₁). The eigensolver is SLEPc — [Wang & Sigmund, arXiv:2210.00003](https://arxiv.org/abs/2210.00003) [Preprint].
- **The IBZ-contour shortcut is not guaranteed.**
  - Spectral extrema "can occur inside the Brillouin zone" — [Harrison, Kuchment, Sobolev & Winn 2007, J. Phys. A 40:7597–7618](https://doi.org/10.1088/1751-8113/40/27/011).
  - The probability that a band-gap extremum lies on the IBZ contour has been quantified for the 17 plane crystallographic groups — [Maurin, Claeys, Deckers & Desmet 2018, IJSS 135:26–36](https://doi.org/10.1016/j.ijsolstr.2017.11.006).
  - Cell-periodic-only buckling in microstructure topology optimization goes back to [Neves, Sigmund & Bendsøe 2002, IJNME 54:809–834](https://doi.org/10.1002/nme.449).

**Practice of Zhang, Feng & Khandelwal [Preprint]** ([arXiv:2010.02371](https://arxiv.org/abs/2010.02371))
- For a parallelogram cell the wave-vector space "can always be chosen as k₁, k₂ ∈ [0,1)".
- "For buckling modes that are larger than one unit cell, i.e. k ≠ 0, the rigid-body motion is implicitly suppressed by the Bloch boundary conditions … so there is no need for adding extra constraints. However … when k = 0, extra constraints are needed for constraining the rigid-body translation, e.g. fixing one arbitrary point."
- They sample "a 100×100 uniform mesh over [0,1)×[0,1)" plus three refined 100×100 zones near the axes and origin, then bisect in load to refine the first bifurcation point.
- At the k → 0 limit "the infimum … is not reachable and can only be computed as the limit k → 0, and the microscale stability surface is discontinuous at origin (singular point)".
- The condensation of Triantafyllidis et al. needs "the inversion of K_NN … only once at each loading step" and "can be more efficient … especially when the sampling grid in the k-space is large".

**Hole-array benchmark result [Slides]**
- Under uniaxial compression of the square array of circular holes, "a critical instability is detected at an engineering strain of −0.03 with k₀₁ = 1/(2l₁) and k₀₂ = 1/(2l₂)". "Y-periodicity is broken and an enlarged cell pY with p = (2,2) is found." Only then, "after the 1st Eigen mode at an engineering strain of −0.055", is loss of ellipticity found — [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf).

**Experiments [Exp]** (`exp1.jl`, `exp3.jl`; Bloch phase applied at element scatter, K_red[I,J] += conj(φ_a)·Ke[a,b]·φ_b)
- Homogeneous periodic cell, n = 12, unpinned:
  - θ = (1e−6, 0) gives three smallest eigenvalues [9.3e−13, 1.3e−12, 10.31]. The translations are recovered continuously as θ → 0.
  - At θ = (π, 0) the unpinned spectrum starts [2.6211, 2.6211, 9.17]. **Pinning the center node there changes it to [2.6211, 4.3652, 9.17]**: the pin is a spurious constraint at k ≠ 0 and raises eigenvalues, i.e., it is non-conservative.
- Boundary condensation (Haynsworth), with a heterogeneous cell (n = 16, 512 DOFs), a 12×12 θ grid, and a shift chosen so that 97 of 144 wave vectors are unstable: the inertia of K(θ) − σI equaled In(K_ii − σI) + In(S(θ)) at **all 144** points, with K_ii k-independent and S the dense Schur complement on the master DOFs.
- Condensation timing:
  - 80×80: 318 master DOFs. K_ii Cholesky 0.049 s; one 318-RHS complex solve 1.55 s; dense S formation 0.028 s; dense `eigvals` 0.045 s; dense Bunch–Kaufman 0.073 s (first calls, including compile).
  - 100×100: 398 master DOFs. 0.044 s, 2.50 s, 0.042 s, 0.043 s and 0.027 s respectively.
  - Direct complex Hermitian Cholesky of the full K(θ): 0.052 s at 80×80 and 0.121 s at 100×100 (warm, `exp1`).

### Inferences
- **The k = 0 meaning of NEMOpt's existing tangent is exact for inertia.**
  - In the unpinned periodic space V, translations t satisfy K t = 0, and V = V_c ⊕ T with V_c the pinned subspace. Hence q(v + t) = q(v) and In(K_full) = In(K_pinned) + (0, 0, 2).
  - A failed Cholesky of the converged reduced Newton tangent therefore means instability against **cell-periodic (1×1) perturbations at fixed F_M**, which is exactly the k = 0 Bloch test.
  - It does not detect period-multiplying modes, such as the 2×2 pattern of hole arrays, or macroscopic modes.
  - The same argument holds for the half cell's 6 prescribed nodes within its symmetry class (Q5).
- **Bloch map on NEMOpt's numbering.**
  - Write θ_j = k·L_j ∈ (−π, π]. For each slave DOF row r, read the integer offset (n₁, n₂) ∈ {0,1}² from the D2 row: column block 1:2 means an L1 jump and 3:4 an L2 jump. Corner slaves carry L1, L1+L2 and L2; edge-3 slaves carry L2 and edge-4 slaves L1 ([NEMOpt.jl L200–203](../../Full_Domain/NEMOpt.jl#L200)).
  - Set D1(θ) = diag(e^{i(n₁θ₁+n₂θ₂)})·D1 and form K(θ) = T(θ)ᴴ K T(θ) in reduced numbering. With the planned one-triangle assembler (NDR-0021), apply the phases at scatter time.
  - **Do not drop the center DOFs for θ ≠ 0. At θ = 0, either pin, which is exact for inertia, or deflate translations, which gives exact magnitudes.**
  - The affine D2(F_M − I)L term plays no role, because perturbations are homogeneous in the jump.
- **Cheap structure.**
  - (i) K(−θ) = conj K(θ) (time reversal), so only half the zone is needed for any design, even one without symmetry.
  - (ii) At the four TRIM points θ ∈ {0, π}², all phases are ±1, so K(θ) is **real** and real Cholesky/LDLᵀ apply. These points cover the 2×1, 1×2 and 2×2 supercell modes.
  - (iii) K(θ) = Σ_{m∈{−1,0,1}²} e^{i m·θ} K_m with real sparse K_m precomputed once at the converged state. Per-k assembly becomes a sparse linear combination, and per-k cost is dominated by the complex Cholesky.
  - (iv) Exact metamorphic identity (Bloch theorem): the θ = 0 periodic spectrum of a p₁×p₂ supercell equals the union of the 1×1 spectra at θ_j ∈ {2πm/p_j}. A 2×2 supercell test validates the phase map against TRIM points to round-off.
- **The neighborhood of k = 0 needs separate treatment.** Along a direction n, λ_min(K(θ)) → 0 like |θ|² for the acoustic branch even in stable cells, so grids must exclude a small ball around θ = 0 and delegate that region to the macroscopic test (Q4).
  - Consistency check linking the two modules: with the identity metric on the uniform grid, lim_{θ→0} λ_min(K(θ))·N_el²/|θ|² along n should approach min eig Q(n) of the computed 𝔸_M. N_el is the number of elements per side, since the Rayleigh quotient of a long-wave Bloch mode is h²|k|²·aᴴQa/|a|². The prediction matches the homogeneous `exp1` values (about 1e−12 at θ = 1e−6, n = 12).
- **Sampling recommendation.** Use the following levels:
  - Level 1: TRIM points, real arithmetic, 3 extra factorizations.
  - Level 2: the IBZ contour of the cell's actual point group plus TRIM, for example 15 points per segment. For a general p2 (C2-only) design this is a larger half-zone contour.
  - Level 3, at final designs: a uniform 16×16 or 32×32 grid on the half zone (with ZFK-style refinement near the axes), the macroscopic test near 0, and bisection in load factor for the first critical point.
  - Contour-only scans must be flagged as not guaranteed (Harrison et al.; Maurin et al.).
- **Condensation is only worth it for dense grids at these sizes.** Per-k direct complex Cholesky (0.05–0.12 s) is already comparable to the condensed per-k cost measured with a naive implementation (0.03–0.07 s). With the S(θ) = Σ e^{i m·θ} G_m decomposition (four offset classes, about 1,300–1,600 real solves of setup, roughly 3–8 s at coordinator speed), per-k cost would fall to an O(n_b³) dense factorization of about 318–398 unknowns, a few ms. That pays off only above a few hundred k-points per state. In-production use requires K_ii to be positive definite (checked once by Cholesky), otherwise Haynsworth's correction term must be tracked.
- **Spurious modes.** NEMOpt's voids use the linear law εℂ with no stress-dependent term, so classical void pseudo-modes are less likely than in SIMP-with-geometric-stiffness schemes. However, the nonlinear "soft-element band" (ρ̃ about 0.37–0.46 at final p and β, per the refactoring report) can host spurious local modes, so stability verdicts on final designs should use thresholded re-analysis, as the report already proposes for comparisons. If optimization ever uses stability constraints, Wang–Sigmund-style separate interpolation of the stress stiffness is the documented remedy.

### Gaps
- The full text of Thomsen et al. 2018 (k-sampling details, verification examples, cost) and of Michel et al. 2007 and Triantafyllidis et al. 2006 (numerical algorithms and tabulated critical values) could not be retrieved: DTU Orbit returned HTTP 403 and the polytechnique.fr server returned 503 or had a TLS chain error.
- The exact condensation algorithm of Triantafyllidis et al. 2006 was not read; the Haynsworth-based version above is this study's construction.
- Whether Zhang, Feng & Khandelwal's preprint has a peer-reviewed version (believed to be in IJNME, 2021) is Unverified.

---

## Q4. Macroscopic check: strong ellipticity (rank-one convexity) of 𝔸_M

### Takeaway
The acoustic tensor test is effectively free: microseconds for an exact quartic test, under 1 ms for a 360-angle scan. It needs no sampling compromise, because in 2D det Q(n) is a homogeneous quartic in n and can be tested exactly by polynomial roots; a 0.25° scan is the published practice. For the FEAD2026 constituent itself there is a closed-form acoustic tensor: **the solid never loses ellipticity in compression**. Any macroscopic loss of ellipticity in NEMOpt is therefore purely structural, and the closed form gives an exact unit test.

### Cited Findings
- Strong ellipticity "implies that Q_ij u_i u_j > 0 for all vectors u and n such that n ⊗ u ≠ 0", with Q_ij = 𝔸_piqj n_p n_q; macroscopic instability occurs "at the first point along the loading path" where Q_ij u_i u_j = 0 — [Slesarenko & Rudykh 2017, JMPS 99:471–482](https://doi.org/10.1016/j.jmps.2016.11.002) (Eqs. 21–22).
- An ellipticity indicator B(λ) = min over unit m, M of (m⊗M):𝔸:(m⊗M) is evaluated "at every π/720 radian increment in both m and M" in [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371) [Preprint].
- Homogenized porous solids can lose strong ellipticity even when the matrix material preserves it — [Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134](https://doi.org/10.1016/0020-7683(81)90092-5), as cited in Zhang et al.
- Bertoldi et al. monitor "the loss of ellipticity of the homogenized tangent modulus" alongside Bloch analysis — [Bertoldi, Reis, Willshaw & Mullin 2010, Adv. Mater. 22(3):361–366](https://doi.org/10.1002/adma.200901956) ([open copy](https://dash.harvard.edu/server/api/core/bitstreams/7312037c-e969-6bd4-e053-0100007fdf3b/content)).
- **Closed form for the FEAD2026 energy [Exp, derived and verified].**
  - From the code's tangent, Q(n) = μI + c·m⊗m with m = F⁻ᵀn and c = λ + μ − λ ln J. Hence min_n λ_min(Q) = min(μ, μ + c/σ_min², μ + c/σ_max²), with σ the singular values of F.
  - `exp2.jl` compared this with a 3,600-angle scan for six states and found agreement to all printed digits (7 significant):
    - equibiaxial 0.5: 38.46154;
    - uniaxial 0.3: 38.46154;
    - equibiaxial 2.0: 38.46154;
    - uniaxial 12: −8.744614;
    - simple shear 0.8: 38.46154;
    - biaxial 6×3: 30.61728.
  - With E = 100 and ν = 0.3 (λ = 57.69, μ = 38.46), ellipticity is lost only when λ(ln J − 1) − μ ≥ μσ_min². That gives **uniaxial tension at stretch exp(1 + 2μ/λ) = 10.3123** (scan at 1.001× that value: −5.8e−2), and **never under any compression**, since c > 0 whenever J ≤ e^{1+μ/λ}.
  - Cost: a naive allocating 360-angle scan takes 6.8e−4 s per state.
- NEMOpt stores 𝔸_M in the Voigt order mapI = [1,2,1,2], mapJ = [1,2,2,1] (11, 22, 12, 21) ([NEMOpt.jl L44](../../Full_Domain/NEMOpt.jl#L44)).

### Inferences
- **Exact test.** In 2D write n = (cos φ, sin φ). Then Q₁₁(n) is a quadratic form and det Q(n) = a₀n₁⁴ + a₁n₁³n₂ + a₂n₁²n₂² + a₃n₁n₂³ + a₄n₂⁴. Strong ellipticity ⇔ [Q₁₁(n) > 0 ∀n] ∧ [det Q(n) > 0 ∀n]. That is: a 2×2 positive-definiteness check, plus a₀ > 0, a₄ > 0, and no real root of the quartic in t = tan φ (companion-matrix eigenvalues).
  - Return the margin min_φ λ_min(Q) by a 1° scan refined with Brent, cross-checked against the root test. The ZFK 0.25° scan is the documented fallback.
  - Use the full 2×2×2×2 𝔸_M, which has no minor symmetry under prestress, never a symmetrized version.
- **Positive definiteness of 𝔸_M versus ellipticity.**
  - Also log the eigenvalues of the 4×4 major-symmetric Voigt matrix 𝔸_M. Its positive definiteness is the relevant condition under mixed or traction control (strain-driven cells are not governed by it).
  - The FEAD2026 objective pushes ν toward −1, i.e., κ → 0⁺, which approaches a macroscopic dilatational mechanism. That is a near-singular 𝔸_M, not loss of ellipticity: in 2D isotropy, ellipticity needs only μ > 0 and κ + μ > 0.
  - Reporting both separates "auxetic near-mechanism" from "macroscopically unstable".
- **Where it runs.** After every homogenization at negligible cost, including inside the FEAD2026 profile as a read-only observer; it touches no state.
- **Combined verdict.** Complete stability at a converged state requires: k = 0 tangent positive definite, Bloch PD for θ outside a small ball, and Q(n) PD for all n. This follows the Geymonat et al. definitions; the region near θ = 0 is covered by the ellipticity test.

### Gaps
- No published NEMOpt-like (compressible neo-Hookean, Bonet–Wood form) acoustic-tensor closed form was found. The derivation here is this study's and needs an independent check: the scan agreement is numerical evidence, not proof.
- For constant-tangent periodic Q4 meshes, whether the discrete FE tangent loses positivity exactly at continuum loss of ellipticity was not established. The test in Q6 should establish it numerically under refinement.

---

## Q5. Half-cell (point-symmetric) reduction and stability

### Takeaway
The half cell represents only the C2-symmetric subspace: odd fluctuations, w(−y) = −w(y). A successful Cholesky on HalfCell tangents therefore proves nothing about antisymmetric modes or any k ≠ 0 mode. Recommended policy:
- Reconstruct the full-cell state from the converged half-cell solution, which is an exact full-cell equilibrium by Palais' principle of symmetric criticality.
- Re-assemble the full-cell tangent once, and run the full diagnostic set on it.
- Use the parity decomposition at k = 0 and at TRIM points only as a cheap optional pre-check and as a metamorphic test.

Group-theoretic block-diagonalization exists in the literature, but gives only a 2-block split at four wave vectors for C2, so it saves little.

### Cited Findings
- The half-domain code prescribes 6 nodes and uses D1 entries of −1 and +1, i.e., point-symmetric constraints ([NEMOpt-DR.jl L242–299](../../Domain_Reduction/NEMOpt-DR.jl#L242)). The refactoring report states that "full-versus-half agreement ends at the first symmetry-breaking bifurcation" and that the half-domain design equals the full-domain design translated by half a period ([refactoring report](../../reports/NEMOpt%20refactoring%20strategy%20literature.md)).
- Principle of symmetric criticality: for a group-invariant functional, under conditions that hold for compact groups such as C2, a critical point of the functional restricted to symmetric fields is a critical point of the full functional — [Palais 1979, Commun. Math. Phys. 69:19–30](https://doi.org/10.1007/BF01941322).
- Group-theoretic computational bifurcation, i.e., block-diagonalization by isotypic decomposition — [Healey 1988, CMAME 67:257–295](https://doi.org/10.1016/0045-7825(88)90049-7) (note: the correct DOI ends in −7, not −8).
- Number of unit cells needed in an RVE for microstructural bifurcation — [Saiki, Terada, Ikeda & Hori 2002, CMAME 191(23–24):2561–2585](https://doi.org/10.1016/S0045-7825(01)00413-3). That this paper uses group theory is Unverified; only the title and a secondary summary were seen.
- Symmetric bifurcation in cellular solids: the work posits that at onset "flipping [the microscopic velocity's] sign does not change macroscopic states" and derives "the onset conditions for symmetric bifurcation" — [Ohno, Okumura & Noguchi 2002, JMPS 50(5):1125–1153](https://doi.org/10.1016/S0022-5096(01)00106-5), content via the [KAKEN 13650084 grant summary](https://kaken.nii.ac.jp/grant/KAKENHI-PROJECT-13650084) (secondary).
- Honeycomb post-bifurcation and stability with symmetry-based analyses — [Combescure, Henry & Elliott 2016, IJSS 88–89:296–318](https://doi.org/10.1016/j.ijsolstr.2016.02.016); [Combescure, Elliott & Triantafyllidis 2020, JMPS 142:103976](https://doi.org/10.1016/j.jmps.2020.103976). The group-theoretic content is Unverified from abstracts.
- Bertoldi's group found that RVEs of increasing size built from "even" primitive cells share the same first eigenmode. They called this enlarged-RVE eigen analysis "not rigorous … Rigorous confirmation performing Bloch wave analyses", citing Saiki, Terada et al. CMAME 2002 — [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) [Slides].

### Inferences
- **Representation theory for C2 (Inference, standard).** The 180° rotation acts on displacement fields as (Rw)(y) = −w(−y).
  - At k = 0, and at the other TRIM points where −k ≡ k, R commutes with the Bloch tangent of a C2-symmetric design at a C2-symmetric state. The space then splits into R = +1 (odd fluctuations, the half-cell space) and R = −1 (even fluctuations, which contain the translations).
  - In(K_full(θ_TRIM)) = In(K_sym) + In(K_anti).
  - At general k, R maps k to −k. Combined with time reversal this only reproduces K(−θ) = conj K(θ); **there is no block split**.
  - The half-cell's 6 prescribed nodes are the C2 centers and their images, where odd fields vanish identically; they are not extra constraints within the symmetric subspace.
- **Recommended HalfCell stability procedure.**
  1. Map half-cell nodal displacements to the full mesh by the C2 extension, applying the half-period design offset.
  2. Assert full-cell residual ‖R_full‖ ≈ 0 at the reconstructed state as a check of the mapping, by Palais.
  3. Assemble the full-cell tangent once.
  4. Run the k = 0 check (pinned), TRIM, contour or grid, and the 𝔸_M ellipticity test, which needs nothing new because 𝔸_M is the same.

  Cost: one full assembly plus diagnostics on 2× the half-cell DOFs, i.e., the full-cell numbers in Q7.
- **Optional parity pre-check.** Build the antisymmetric half-cell tangent by flipping the parity map: D1 sign pattern inverted, C2-center nodes freed, and one node pinned to remove translations. Its Cholesky plus the existing HalfCell Cholesky gives the full k = 0 verdict at half-cell cost. The same split at the three non-zero TRIM points gives period-doubling checks.
  - This is worth implementing first as a **metamorphic test**: In(full, k = 0) = In(sym) + In(anti). Only if profiling shows the full-cell reassembly matters should it become a production path.
- **Interpretation for FEAD2026 results.** A HalfCell optimum is guaranteed to be an equilibrium of the full cell but not a stable one. Full-versus-half agreement in L2 tests (refactoring report) is expected exactly as long as the full-cell diagnostics are stable. A full-cell instability at a published HalfCell optimum would mean the published design sits on an unstable symmetric branch: worth reporting, not a reproduction failure.

### Gaps
- No paper was found that treats Bloch stability of point-symmetric half-cell reductions in finite-strain topology optimization specifically.
- The group-theoretic content of Saiki et al. 2002 and Combescure et al. 2016/2020 was not verified from full text.
- Ikeda & Murota's *Imperfect Bifurcation in Structures and Materials* (Springer) is the standard monograph, but its edition and DOI were not checked.

---

## Q6. Published numerical benchmarks for verifying an implementation

### Takeaway
The strongest benchmarks are exact or semi-analytic and can be run inside NEMOpt's own material model:
1. the FEAD2026 acoustic-tensor closed form (Q4);
2. Sylvester and Bloch identities (Q1, Q3);
3. **laminates**, where finite-strain homogenization is exact for piecewise-homogeneous fields and the incompressible neo-Hookean limit has the closed form λ_cr = (1 − μ̆/μ̄)^{1/4};
4. the square array of circular holes, for qualitative and semi-quantitative agreement: first instability at k = (½, ½), i.e., a 2×2 pattern, for porosity above about 0.34, at an engineering strain of about −0.03, with loss of ellipticity at about −0.055 for Bertoldi's material and geometry.

Exact hole-array numbers require the original papers' material laws (Rivlin, Gent, or near-incompressible neo-Hookean), which differ from NEMOpt's compressible law. Hole-array checks are therefore cross-code tests, not reproductions.

### Cited Findings
**Laminates [PR]**
- For plane-strain compression along the layers of an incompressible neo-Hookean laminate, the macroscopic critical stretch is λ_macr = (1 − μ̆/μ̄)^{1/4}. Here μ̄ = v_f μ_f + v_m μ_m and, for layers, μ̆ is the harmonic mean (v_f/μ_f + v_m/μ_m)⁻¹ — Eqs. 29–30 of [Slesarenko & Rudykh 2017, JMPS 99:471–482](https://doi.org/10.1016/j.jmps.2016.11.002), attributing the 2D result to [Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148). The 3D analogue has exponent 1/3 (Eq. 28). Eq. 30 was garbled in text extraction and is read here as the harmonic mean, which is the standard laminate shear modulus (Inference).
- Example values computed from that formula (Inference, arithmetic):

| v_f | μ_f/μ_m | λ_cr |
|---|---|---|
| 0.5 | 10 | 0.9045 |
| 0.2 | 10 | 0.8668 |
| 0.1 | 15 | 0.8574 |
| 0.5 | 100 | 0.9900 |

- In 3D fiber composites at small v_f, the microscopic (finite-wavelength) curves "branch out" from the macroscopic ones, i.e., micro instability precedes macro — [Slesarenko & Rudykh 2017](https://doi.org/10.1016/j.jmps.2016.11.002). For layered 2D composites, Triantafyllidis & Maker found that the macroscopic stable region envelops the microscopic one ([via Zhang et al.](https://arxiv.org/abs/2010.02371)). Combined normal and shear loading of layered composites is treated by [Nestorović & Triantafyllidis 2004, JMPS 52(4):941–974](https://doi.org/10.1016/j.jmps.2003.06.001).

**Square array of circular holes**
- Uniaxial compression gives a 2×2 pattern of "mutually orthogonal ellipses" from "an elastic instability above a critical value of the applied strain".
  - Numerical void fractions studied: 0.34, 0.39, 0.49, 0.59, 0.69. "For F < 0.34 (where localization occurs prior to the pattern transformation) the Poisson's ratio is always positive during loading." The lowest asymptotic ν was −0.904 at F = 0.70.
  - Material: "a two-term I1-based Rivlin model, which was modified to include compressibility with a high bulk modulus"; CPE6H elements.
  - Experimental specimen: holes of "9.95 ± 0.02 mm" from 10 mm pins at "spacing of 12.0 mm", i.e., porosity π·4.975²/12² = 0.540 (Inference, arithmetic).
  - ε_c versus F appears only in a figure (Fig. 3a).
  - Source: [Bertoldi, Reis, Willshaw & Mullin 2010, Adv. Mater. 22:361–366](https://doi.org/10.1002/adma.200901956) ([open copy](https://www.epfl.ch/labs/flexlab/wp-content/uploads/2017/12/22_2009_Bertoldi_AdvMat_NegativePoissonRatioElasticInstability.pdf)).
- For the infinite square array under uniaxial compression: critical Bloch instability at engineering strain −0.03 with k = (1/(2l₁), 1/(2l₂)), giving a 2×2 enlarged cell; loss of ellipticity at −0.055. The PSM-4 rubber was fitted by a two-term Rivlin model with "C1 = 0.55 MPa, C2 = 0.3 MPa", and "for the two-terms Rivlin and the Neo-Hookean material models almost identical stress-strain behavior" — [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) [Slides]. These numbers correspond to the JMPS paper below but were read only from the slides.
- Primary papers whose bibliographic data were verified but whose full text was not retrieved:
  - [Bertoldi, Boyce, Deschanel, Prange & Mullin 2008, JMPS 56(8):2642–2668](https://doi.org/10.1016/j.jmps.2008.03.006) (square and oblique lattices of circular voids, rectangular lattices of elliptical voids; Bloch and loss-of-ellipticity analyses);
  - [Mullin, Deschanel, Bertoldi & Boyce 2007, PRL 99:084301](https://doi.org/10.1103/PhysRevLett.99.084301);
  - [Bertoldi & Boyce 2008, PRB 77:052105](https://doi.org/10.1103/PhysRevB.77.052105);
  - [Overvelde & Bertoldi 2014, JMPS 64:351–366](https://doi.org/10.1016/j.jmps.2013.11.014);
  - [Michel, Lopez-Pamies, Ponte Castañeda & Triantafyllidis 2007, JMPS 55(5):900–938](https://doi.org/10.1016/j.jmps.2006.11.006);
  - [Triantafyllidis, Nestorović & Schraad 2006](https://doi.org/10.1115/1.2126695) (failure surfaces of porous neo-Hookean solids with circular holes in strain space).
- Zhang et al. report first bifurcation points in the same notation, e.g., λ_c = 0.9270 at (k₁, k₂) = (0.5, 0.5) for a rectangular RVE and λ_c = 0.9934 with the minimum at the origin for another case. These are usable as cross-code targets only if their geometry and material (in the preprint) are replicated — [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371) [Preprint].

**Wallin & Tortorelli golden case**
- It contains a negative tangent component 𝔸₁₂₁₂ = −0.9861 ([refactoring report](../../reports/NEMOpt%20refactoring%20strategy%20literature.md); [Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324)).

### Inferences
- **Benchmark B-S1, laminate oracle (highest value).**
  - Use a two-density layered design aligned with the mesh: filter bypassed (the NDR-0043 pure-SIMP test mode), both phases in the nonlinear branch, i.e., ρ^p well above the 0.05 switch threshold, stiffness ratio g(ρ₂).
  - The exact finite-strain solution is piecewise homogeneous with a rank-one jump a⊗N at the interfaces. Q4 elements represent it exactly when interfaces lie on element edges, so FE 𝔸_M must equal the semi-analytic laminate tangent (a 2-unknown nonlinear solve for a) to round-off, at any ν.
  - That makes three checks: (i) 𝔸_M correctness at finite strain, an independent oracle for NDR-0020/0021; (ii) the acoustic-tensor module, via loss of ellipticity under compression along the layers; (iii) the incompressible limit (ν → 0.5, e.g., 0.49, 0.499) against λ_cr = (1 − μ̆/μ̄)^{1/4}.
  - Bloch scans on the same cell should find no instability before macroscopic loss of ellipticity, if the 2D Triantafyllidis–Maker envelope result carries over to the compressible case (to be confirmed).
  - Running it needs a compressive Newton (NDR-0031/0032 options), so it lives outside FEAD2026.
- **Benchmark B-S2, hole array (semi-quantitative).**
  - Use a circular hole at porosity 0.5–0.6 on 80×80 and 100×100, with voids as ρ = 0 (linear law with ε = 5e−4·E) under uniaxial plane-strain compression.
  - Expectations: (a) the first instability at the M point θ = (π, π), real arithmetic; (b) the critical strain in the low-percent range (about −0.03 for Bertoldi's Rivlin rubber at porosity of about 0.5–0.6); (c) macroscopic loss of ellipticity later; (d) for porosity below about 0.34, macroscopic first.
  - Quantitative acceptance should be against an independent body-fitted computation with NEMOpt's exact material (e.g., a Ferrite/Gridap or commercial-code script with Bloch conditions), because staircase holes, residual void stiffness and the compressible Bonet–Wood law shift ε_c by an amount not known a priori. Expected magnitude: a few percent (Inference, unquantified).
- **Benchmark B-S3, homogeneous solid.**
  - Under any compressive F_M, Q(n) ≥ μ (closed form), so 𝔸_M of a uniform solid cell (= 𝔸(F_M), an existing L2 oracle) must pass the ellipticity test.
  - Under uniaxial tension, failure must occur at s = 10.3123 to the bisection tolerance.
  - The Bloch K(θ) of the uniform cell should stay PD at all θ under modest compression. This is expected from the continuum (constant-coefficient ellipticity ⇒ stability against periodic perturbations), but **for the discrete Q4 tangent, which is not PD as a 4×4 Voigt matrix under compression because of the stress term, it must be confirmed under mesh refinement** (Inference).
- **Benchmark B-S4, regression vectors.** Store the diagnostics outputs of the Wallin & Tortorelli §4.1 case and of the §5.2 published optima (inertia at k = 0, TRIM verdicts, min Q eigenvalue, 𝔸_M eigenvalues) as hash-stamped golden values once accepted. They protect the diagnostics themselves against regressions.

### Gaps
- Exact published critical strains or stretches with full parameter sets for neo-Hookean hole arrays could not be extracted. Michel et al. 2007, Triantafyllidis et al. 2006, Bertoldi et al. 2008, Bertoldi & Boyce 2008 and Overvelde & Bertoldi 2014 were not accessible (403/503, TLS) or report values only in figures.
- Recalled but **Unverified**: Mullin et al. 2007 hole diameter 8.67 mm with 9.97 mm spacing (porosity about 0.59); the Overvelde & Bertoldi 2014 circular-pore critical strain.
- Geymonat et al. 1993 give general theory and a laminate example; their numerical values were not checked.
- The compressible neo-Hookean laminate critical stretch has no closed form found here; it must come from the semi-analytic laminate solve.

---

## Q7. Architecture hooks, logging, provenance, tests and run-time budget, without changing FEAD2026

### Takeaway
Implement diagnostics as **read-only observers** on a converged-state object, outside the optimizer's call graph. Keep the stdlib-only parts (Cholesky flag, shifted-Cholesky bisection, LDLᵀ count, TRIM and Bloch via CHOLMOD complex, acoustic tensor) in the core. Put KrylovKit, MUMPS and Pardiso in package extensions. Log every factorization with its state tag. Accept the whole feature by one Regime-B test: all FEAD2026 golden masters and trajectories bitwise identical with diagnostics on and off.

Cost at coordinator speed:
- Per-evaluation k = 0, TRIM and ellipticity checks: under 0.5 s at 80×80 and under 0.8 s at 100×100, i.e., under 1% of a legacy evaluation.
- Full Bloch scans on final designs: about 6–50 s at 80×80 and 15–120 s at 100×100, depending on path versus grid.

### Cited Findings
- Package extensions are "a module in a file … automatically loaded when some other set of packages are loaded" and require Julia 1.9+ — [Pkg docs, "Conditional loading of code in packages (Extensions)"](https://pkgdocs.julialang.org/v1/creating-packages/). The refactoring report already plans `ext/NEMOptPlotsExt.jl`.
- NonconvexMMA has no per-iteration callback. The FEAD2026 profile must preserve its call order (objective before constraint at the same x) — refactoring report, NDR-0026, citing [NonconvexCore](https://github.com/JuliaNonconvex/NonconvexCore.jl/blob/master/src/functions/value_jacobian.jl).
- The legacy reduced tangent is asymmetric at about 1.4e−14 and is routed to UMFPACK LU; 59 factorizations per evaluation; 49 s per 40×40 evaluation; the thesis reports 156.02 s per iteration at 80×80 with domain reduction — [refactoring report](../../reports/NEMOpt%20refactoring%20strategy%20literature.md); [coordinator benchmark](../NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl).
- Coordinator measurements (Julia 1.11.7, about 13k-DOF reduced tangent): one LU about 0.15 s; one Cholesky with 4 RHS about 0.08 s (task statement). This container measured 0.091 s and 0.033 s respectively for the same sizes and pattern, so the coordinator machine is about 1.6–2.4× slower [Exp].
- Homogenized results "are only valid before the bifurcation occurs" — [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371), as recorded in `numerical_method_choices.md`.

### Inferences

**Where the diagnostics plug in**
- `ConvergedState` (immutable view): design hash, F_M, load-step index, U, cached Gauss-point 𝔸/P (NDR-0017), the reduced-tangent factor if one exists (NDR-0016/0022), 𝔸_M, CellSymmetry.
- `diagnose(state, ::StabilityConfig) -> StabilityReport` is pure. It never mutates caches or factors (it copies or refactors), and draws random start vectors only from its own `StableRNG(seed)`. It scopes any CHOLMOD parameter changes (as Julia's `@cholmod_param` does) and never calls `BLAS.set_num_threads`.
- Trigger points, configured per profile:
  1. **after each load step's Newton convergence** (cheap k = 0 flag only; tag the state);
  2. **after each objective evaluation** (k = 0 verdict, 𝔸_M ellipticity and Voigt eigenvalues; optionally TRIM);
  3. **every N outer iterations and at continuation-stage ends** (TRIM plus contour);
  4. **on final designs**, raw and thresholded (Level-3 grid, margins, critical-load bisection).

  In FEAD2026 the trigger runs after the objective value and gradient are computed and returned, never between `Objective` and `Objective_Grad`. Default: off, except triggers 2 and 4 as opt-in.
- Fallback policy for the Default profile's Cholesky solver: Cholesky, then (if it fails at a converged state) the logged LDLᵀ inertia diagnostic, then UMFPACK LU for the actual solve. The FEAD2026 profile keeps LU for the solve and runs the Cholesky attempt **only on a separate copy**, `Symmetric(K, :L)`, logging ‖K − Kᵀ‖/‖K‖.

**Logging schema (one record per factorization)**
- Run ID, evaluation and iteration counters, load step, Newton iteration, state tag (`:iterate | :converged | :chi | :adjoint | :diagnostic`), θ.
- Matrix size and nnz; solver path (`:cholmod_ll | :cholmod_ldlt | :umfpack | :mumps`); success; CHOLMOD failure column (minor) if obtainable; negative-pivot count; min|D| and max|D|; shift; symmetrization residual; wall time.
- Aggregates per evaluation: Cholesky attempts, failures by tag, fallback counts.

**Provenance fields (extending NDR-0041)**
- Julia version; SparseArrays and `CHOLMOD.BUILD_VERSION` (5.2.0 in 1.11.7); BLAS vendor and threads.
- Eigen package and version; algorithm parameters (krylovdim, tol, nev, maxiter, convergence flags and residuals); RNG seed; metric (`:identity_lumped`), deflation mode (`:pin | :translation_projection`).
- k-sampling spec (type, N per segment or grid, time-reversal and point-group reduction, excluded radius around 0); acoustic-tensor method (`:quartic_roots | :scan(Δφ)`); tolerances defining the three-state verdict.
- Design state (raw or thresholded); CellSymmetry and how a full-cell state was reconstructed from a HalfCell.

**Packaging**
- Core: `stability/` with `inertia.jl` (Cholesky flag, shifted bisection, LDLᵀ count), `bloch.jl` (phase map from D1/D2, TRIM, contour, grid, optional condensation), `ellipticity.jl` (acoustic tensor, quartic test) and `report.jl`.
- Extensions:
  - `NEMOptKrylovKitExt` (margins, modes);
  - `NEMOptArnoldiMethodExt` (fallback);
  - `NEMOptMUMPSExt` (exact real inertia with ICNTL(13)=1);
  - optionally `NEMOptPardisoExt`.
- Arpack.jl is not used.

**Run-time budget**
Estimated at the coordinator's speed (this container's measurements × 2; plausible range × 1.6–2.4):

| Diagnostic | 80×80 (12,798 DOF) | 100×100 (19,998 DOF) | When |
|---|---|---|---|
| D0: k = 0 Cholesky flag on converged tangent | ≈ 0 extra in Default (factor reused); ≈ 0.06 s in FEAD2026 (copy) | ≈ 0 / ≈ 0.10 s | every evaluation |
| D1: macroscopic ellipticity + 𝔸_M eigenvalues | < 1 ms | < 1 ms | every evaluation |
| D2: TRIM points (π,0), (0,π), (π,π), real Cholesky each + phase assembly | ≈ 0.2–0.5 s | ≈ 0.3–0.7 s | every evaluation (Default) or every N iterations |
| D3: inertia count on failure (simplicial LDLᵀ) | ≈ 0.6 s | ≈ 2.1 s | only after a failure |
| D4: λ_min margin (factor + about 29-solve shift-invert Lanczos, 6 eigs, deflated) | ≈ 0.26 s real; ≈ 0.6 s complex per θ | ≈ 0.38 s real; ≈ 0.9 s complex per θ | final designs; on demand |
| D5: Bloch contour, about 60 θ (complex Cholesky each) | ≈ 6 s | ≈ 15 s | every N iterations; final designs |
| D6: Bloch half-zone grid 16×16 (≈ 129 θ) / 32×32 (≈ 513 θ) | ≈ 13 s / ≈ 53 s | ≈ 31 s / ≈ 123 s | final designs only |
| D7: critical load by bisection (about 10–15 load levels × D2 or D5) | ≈ 3–90 s plus Newton re-solves | ≈ 5–220 s plus Newton re-solves | final designs, compressive studies |
| HalfCell: full-cell reconstruction + one full assembly | one full-cell assembly (legacy triplet assembly 1.36 s at 80×80; much less with the Tensors.jl kernel) | proportionally more | whenever a HalfCell result is diagnosed |

- Against the legacy evaluation (49 s at 40×40; 156 s per iteration at 80×80 with domain reduction), D0–D2 cost under 1% per evaluation. After Phase 3–4 speedups bring an evaluation to a few seconds (unquantified), D2 every iteration could approach 10–20%, hence "every N iterations" as the default.

**Tests (IDs extend the report's scheme)**
- K-S1: closed-form FEAD2026 acoustic tensor versus scan; 6 states; relative error ≤ 1e−12.
- K-S2: loss-of-ellipticity threshold in uniaxial tension at s = exp(1 + 2μ/λ) to bisection tolerance; no loss under a dense set of compressive F.
- K-S3: quartic-root test versus 0.25° scan on random major-symmetric tensors (metamorphic).
- E-S1: Sylvester tests. Shifted Cholesky succeeds iff σ < λ_min away from ties; the LDLᵀ negative count equals the dense eigenvalue count on small heterogeneous cells (as in `exp2`: 28/28); a three-state verdict near ties.
- E-S2: In(pinned) = In(unpinned) − 2 zeros at θ = 0; the pinned spectrum contains exactly two extra pseudo-modes that drift down under refinement.
- E-S3: Bloch metamorphic tests. Hermitian to round-off; K(−θ) = conj K(θ); 2π periodicity; TRIM matrices real; the 2×2 supercell spectrum at θ = 0 equals the union of the 1×1 TRIM spectra; cell-symmetry permutations of θ.
- E-S4: long-wave consistency, λ_min(K(θ))·N_el²/|θ|² → min eig Q(n) of 𝔸_M as θ → 0.
- E-S5: HalfCell parity, In(full, k = 0) = In(sym) + In(anti); the C2 reconstruction gives full-cell residual ≈ 0.
- E-S6: condensed-versus-direct inertia over a θ grid (`exp3`: 144/144).
- B-S1 to B-S4: as in Q6.
- **R-S0 (gate): every FEAD2026 golden master and trajectory bitwise identical with all diagnostics enabled versus disabled.**

**Proposed decision records (all opened as "proposed")**

| ID | Proposed decision | Class | FEAD2026 effect | Acceptance evidence |
|---|---|---|---|---|
| NDR-0045 | Stability diagnostics as a pure read-only observer on `ConvergedState`, called outside the optimizer call graph; private RNG; scoped CHOLMOD parameters; three-state verdicts | B | none; off by default | R-S0 bitwise on all golden masters and trajectories with diagnostics on |
| NDR-0046 | Factorization log with state tags, failure column, fallback counters and symmetrization residual; FEAD2026 runs the Cholesky attempt on a copy only | B | none | R-S0; log present for every solve |
| NDR-0047 | k = 0 stability = Cholesky flag on the converged reduced tangent; on failure, simplicial LDLᵀ negative-pivot count, cross-checked by dense Bunch–Kaufman in tests and by MUMPS (ICNTL(13)=1) in an extension; margins by translation-deflated shift-invert, never by the pinned spectrum | B (new output) | none | E-S1, E-S2 |
| NDR-0048 | Macroscopic strong-ellipticity test of 𝔸_M (exact quartic test plus scan margin) and 𝔸_M Voigt eigenvalues, logged per evaluation | B (new output) | none | K-S1 to K-S3, B-S3 |
| NDR-0049 | Bloch–Floquet module: phases from D1/D2 offsets, no pin for θ ≠ 0, time-reversal halving, TRIM real checks, contour or grid options flagged "not guaranteed" for contour-only; complex Hermitian CHOLMOD; optional Haynsworth condensation | B (new capability) | none | E-S3, E-S4, E-S6, B-S1, B-S2 |
| NDR-0050 | HalfCell stability policy: C2 reconstruction plus full-cell diagnostics by default; parity-split half-cell tangents as test and pre-check only | B | none | E-S5; residual check |
| NDR-0051 | Eigen-solver and packaging: KrylovKit (Lanczos) as an extension, ArnoldiMethod as fallback, no Arpack.jl; MUMPS and Pardiso as optional extensions; stdlib-only core | P | none | dependency audit; extension load tests |
| NDR-0052 | Stability verification suite: analytic (K-S*), metamorphic (E-S*) and benchmark (B-S*) tests; hole-array acceptance against an independent body-fitted computation with NEMOpt's material, not against figure-read literature values | P | none | suite green; documented tolerances |
| NDR-0053 | Diagnostics provenance fields and report schema (extends NDR-0041) | P | none | provenance audit |
| NDR-0054 | Compressive and critical-load studies (bisection in load, optional Z-based linearized predictor) require the NDR-0031/0032 Newton options; never part of FEAD2026 | A (new studies) | n/a: legacy Newton hangs for compressive F_M | B-S1/B-S2 run in the Default profile |
| NDR-0055 | (Deferred research) stability-constrained optimization (Bloch or ellipticity constraints, stress-stiffness interpolation against pseudo-modes) | A | n/a | separate study |

### Gaps
- The runtime of KrylovKit, ArnoldiMethod, MUMPS and Pardiso on these matrices was not measured; the packages were not installed, so the cost table uses a hand-written Lanczos and stdlib CHOLMOD.
- The coordinator-machine scaling factor (1.6–2.4×) comes from only two paired measurements.
- Per-k phase-assembly cost for the planned assembler is estimated, not measured.
- Whether Julia's CHOLMOD wrapper exposes `L->minor` through a public API was not checked; the `PosDefException` message does not contain it.
- Whether any published FEAD2026 optimum is unstable is unknown until NDR-0047/0049 are run on recovered designs.

## Coordinator verification note (added after re-running the evidence scripts)

The coordinator re-ran `exp1.jl` and `exp2.jl` from `evidence/stability_experiments/` (Julia 1.11.7, CHOLMOD 5.2.0) and reproduced the recorded output. One caveat applies to `exp1.jl`.

- **Its shifted-inertia lines (m = 1, 3 and 7) are not a valid inertia test.** The script places the shift σ midway between `ev[m]` and `ev[m+1]` and assumes this leaves exactly m eigenvalues below σ. The square-symmetric test cell has eigenvalues in equal pairs, and all three shifts land on a repeated pair, so K − σI is singular in every case. *(Corrected after the addendum writer's re-computation, `evidence/synthesis_checks/exp1_degeneracy_check.py`: strictly negative / zero counts are 0/2, 2/2 and 6/2 for m = 1, 3 and 7. The first version of this note named only m = 3 and m = 7.)*
- **What the two methods report there.** Unpivoted CHOLMOD LDLᵀ gives 0, 4 and 6 negative pivots: it counted both near-zero eigenvalues as negative only at m = 3; at m = 1 and m = 7 its count equals the strictly negative count. Tolerance-based dense Bunch–Kaufman `inertia` reports (286, 0, 0), (282, 2, 2) and (279, 6, 1): it flags two, one or none of the paired eigenvalues as zero depending on round-off. Both behaviors are round-off-dependent at an exactly singular shift. *(Corrected: the first version of this note said LDLᵀ counted the near-zero pivots as negative at both m = 3 and m = 7.)*
- **The reliability evidence is `exp2.jl`.** It uses a heterogeneous cell, skips nearly equal eigenvalues, and gives 28 shifts with 0 mismatches. Its near-coincident test also shows that at a relative distance of 1e−13 from an eigenvalue, the tolerance-based Bunch–Kaufman count changes from (385, 5, 0) to (385, 4, 1).
- **Practical consequence.** Inertia counts are unreliable when a shift or load level lies within round-off of an eigenvalue. That supports the three-state verdict (stable / indeterminate / unstable) with a margin, rather than a raw negative-pivot count.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. The items below correct it after the supplied papers and open-access sources were read in full. Each item gives the corrected statement and the verification note that supports it. All notes are in [NEMOpt supplied papers verification](../NEMOpt%20supplied%20papers%20verification/). Line numbers refer to this file as of 2026-10-12.

- **L146 ("the finite-strain analogue of the linearized buckling of Thomsen et al. 2018").** Corrected: "a finite-strain analogue of linear buckling (Thomsen et al. 2018, Eq. 10), closest to the 'one point approach' of Dalklint et al. 2021, K ≈ K_o + π_jK_g (Eqs. 18 and 21)". Thomsen et al.'s geometric stiffness is a linear prestress stiffness at small strain. Source: [optimization note](../NEMOpt%20supplied%20papers%20verification/open_access_optimization_papers.md), §1 edit 9.
- **L157 ("Three tests together are complete in the sense of Geymonat–Müller–Triantafyllidis").** Corrected: "Three tests together cover every mode class that Geymonat–Müller–Triantafyllidis distinguish; near θ = 0 the coverage holds only in the limit sense of their Cor. 3.5 and along a load path that starts stable". Source: [GMT note](../NEMOpt%20supplied%20papers%20verification/geymonat_muller_triantafyllidis_1993.md), E17.
- **L171 (Thomsen et al. quote and source).** Corrected quote: "homogenization theory, a linearized stability criterion and Floquet-Bloch theory". Corrected source: the DTU Orbit accepted manuscript, read in full, replaces the abstract via the DTU Orbit record. Source: optimization note, §1 edit 10.
- **L225 (full texts of Thomsen et al. 2018, Michel et al. 2007 and Triantafyllidis et al. 2006 not retrieved).** Corrected: all three were later read in full.
  - Thomsen et al. 2018, from the DTU Orbit accepted manuscript. It reports no CPU times.
  - Michel et al. 2007, from the HAL deposit (submitted version). The PDF supplied under that name was in fact Lopez-Pamies & Ponte Castañeda 2007.
  - Triantafyllidis et al. 2006, from a user-supplied PDF.
  - Michel et al. and Triantafyllidis et al. report critical values only in figures, with no tables; digitized values are in the Michel and TNS notes.
  - Sources: optimization note, §1 edit 11; [Michel note](../NEMOpt%20supplied%20papers%20verification/michel_et_al_2007.md), E12; [TNS note](../NEMOpt%20supplied%20papers%20verification/triantafyllidis_nestorovic_schraad_2006.md), §5 research-note edits.
- **L226 ("The exact condensation algorithm of Triantafyllidis et al. 2006 was not read").** Corrected: Triantafyllidis et al. 2006 condense onto the left, bottom and corner DOFs (their Eq. 3.12). They sign-test the Schur complement by LDU without stating that K_ii must be positive definite. The Haynsworth-based version above generalizes it. Source: TNS note, §5 research-note edits.
- **L239 ([Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134](https://doi.org/10.1016/0020-7683(81)90092-5), "as cited in Zhang et al.").** Wrong paper. Corrected: "[Abeyaratne & Triantafyllidis 1984, J. Appl. Mech. 51(3):481–486](https://doi.org/10.1115/1.3167661), as cited in Zhang et al. and in Geymonat et al. 1993, p. 233". Triantafyllidis et al. 2006 also cite only the 1984 paper for porous solids. Sources: GMT note E18; TNS note, §2.2.
- **L263 ("This follows the Geymonat et al. definitions; the region near θ = 0 is covered by the ellipticity test.").** Corrected: "These conditions are necessary by Geymonat et al.'s definitions. The ellipticity test covers the region near θ = 0 only in the sense that, at a state with Λ ≥ 0, α_𝕄 = 0 if and only if Λ₄ = Λ₅ = 0 (Thm 3.3(i), Cor. 3.5). They give no radius below which α_𝕄 > 0 guarantees Bloch positivity, so the ball radius is an engineering choice." Source: GMT note E19.
- **L353 (Triantafyllidis, Nestorović & Schraad 2006: "failure surfaces of porous neo-Hookean solids with circular holes in strain space").** Corrected: "failure surfaces in log-strain space for compressible neo-Hookean two-phase solids, κ/μ = 98 and 9.8, with circular inclusions at volume fraction π/16 and μ_f/μ_m from 0.02 ("voids") to 50". Their "voids" are soft inclusions, not empty holes. The Michel note's interim wording (matrix law unverified) is superseded by the TNS full text. Sources: TNS note, §5 research-note edits (claim N11); Michel note E12.
- **L368, expectation (d) ("for porosity below about 0.34, macroscopic first").** Corrected: no macroscopic-first expectation. At porosity π/16 with κ/μ = 10, Michel et al. 2007 find the 2×2 mode first under uniaxial-strain and equibiaxial compression. In the addendum, the hole-array benchmark is now split into B-S2a to B-S2d. Source: Michel note E12 (claim 15).
- **L377 (published critical strains for neo-Hookean hole arrays could not be extracted).** Corrected:
  - Michel et al. 2007 give full parameter sets for compressible neo-Hookean hole arrays (κ/μ = 10), with critical strains in figures only. They are digitized in the Michel note from the HAL submitted version and must be checked against the journal version before tolerances are fixed.
  - Triantafyllidis et al. 2006 report values only in figures; digitized values are in the TNS note.
  - Bertoldi et al. 2008, Bertoldi & Boyce 2008 and Overvelde & Bertoldi 2014 remain unread.
  - Sources: Michel note E12; TNS note, §5 research-note edits.
- **L379 ("Geymonat et al. 1993 give general theory and a laminate example; their numerical values were not checked").** Corrected: "… their quoted loss-of-ellipticity stretch λ = 0.927 (stiff fraction 0.5, p. 288) was reproduced as 0.9267 from their Eqs. 6.6, 6.55 and 6.56". Source: GMT note E20; scripts in `evidence/gmt_laminate/` of the verification folder.
