# Prove equivalence first, then accelerate NEMOpt

Seven specialists are best placed to answer this question: a finite-strain homogenization mechanician, a density-based topology-optimization methodologist, a design-sensitivity (adjoint) analyst, a sparse-linear-algebra and floating-point specialist, a scientific-software verification and validation (V&V) engineer, a Julia performance and package architect, and a research-reproducibility specialist. They converge on one strategy. **Freeze and characterize the legacy code. Restructure it into a typed Julia package that replays the legacy arithmetic bit for bit. Take the speedups that do not change a single bit. Only then change the arithmetic (verified to tolerance) and the formulation (behind named legacy switches), and record every step as a numbered decision record opened as "proposed".**

The central finding is that "same results as the article" cannot be a single test. The repository cannot reproduce any published optimization result exactly. The initial designs and seeds, the Table 2 densities, the package versions, the per-case optimizer variant and the linear-response code are all unrecorded, and neither shipped driver matches a published case ([Güven & Özdemir 2026](https://doi.org/10.1016/j.finel.2026.104599); [driver settings](../Full_Domain/NEMOpt_main.jl#L10)). What can be proved is a chain of four links:

| Link | Claim | Kind of proof |
|---|---|---|
| E1 | The refactored code equals the legacy code bitwise in a pinned environment | Hash equality |
| E2 | Modernized numerics stay inside the legacy code's own measured round-off noise | Statistical |
| E3 | The legacy code reproduces the article's verification claim (adjoint-versus-finite-difference differences of at most 5.6e−8) and its Poisson's-ratio and bulk-modulus values within pre-registered bands | Evidential |
| E4 | Every deliberate deviation is re-run on the article's cases and reported | Documented re-runs |

The full text of Wallin & Tortorelli (2020), the formalism the article says it follows, strengthens both ends of the chain. Most of the code's apparent quirks are that paper's reference formulation, not NEMOpt inventions: the clipped SIMP law, the blend evaluated on ρ̃^p, the symmetric sensitivity and the fixed center node. Its fully specified 4×4 verification case gives the refactor an independent oracle that needs no unpublished input. That case was reproduced to every printed digit in this study and independently re-run by the coordinator. It also shows that **the code is correct and the article's printed Eq. 20 is a misprint** ([Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324)).

The measured baseline justifies the effort. A 40×40 full-domain evaluation takes **49 s**. Roughly 22 s of that is a constitutive call that takes 40 µs and allocates 21.9 KB, executed 87 times per Gauss point; a Tensors.jl kernel does the same work in 0.08 µs. A 1.4e−14 round-off asymmetry also routes all **59 sparse factorizations** per evaluation to LU ([coordinator benchmark](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl)).

Seven verified code behaviors change results if altered, so they must survive as switches. Three claims from the earlier review are corrected here:
- H_steep(ρ̃^p) matches the 2023 thesis and Wallin & Tortorelli, so it is not a defect.
- Only Table 3, not Table 4, uses ω as the steepness.
- The extra Newton correction changes results rather than only wasting work.

## Seven specialists converge on freeze, prove, then accelerate

The question sits where finite-strain mechanics, optimization methodology, numerical analysis and software engineering meet, so no single discipline can answer it. The panel below is the smallest group whose combined expertise covers every decision the refactor forces. Each row states what that specialist would insist on and the primary evidence behind it.

Three evidence labels are used throughout:
- **Verified** means checked in this study against one of four sources: the article; Güven's thesis (Güven, M. (2023), M.Sc. thesis, İzmir Institute of Technology, [https://gcris.iyte.edu.tr/handle/11147/13904](https://gcris.iyte.edu.tr/handle/11147/13904)); the code at commit `30d0842`; or the coordinator's instrumented measurement.
- **Inference** means derived from verified facts but not tested.
- **Unverified** means taken from a source that could not be checked.

The full text of Wallin & Tortorelli (2020) was read for this report and compared with the code item by item. The article names it as its primary formalism, and the paper is open access under CC BY-NC-ND. The comparison note's scripts were re-run independently by the coordinator ([Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324); [comparison note](../research_notes/NEMOpt%20refactoring%20strategy%20literature/wallin_tortorelli_2020.md); [evidence scripts](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/wallin_tortorelli_2020/)).

| Specialist | What they would say | Anchoring evidence |
|---|---|---|
| Finite-strain computational-homogenization mechanician | Preserve the mechanics exactly: Bloch-type total-displacement boundary conditions, the fixed center node (full domain) or six prescribed nodes (half domain), the symmetric total-field χ-form of 𝔸_M, and the parity rules of the point-symmetric half cell (odd fluctuations, χ fields and adjoints; even densities). A 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear load case, but only on the symmetric branch (Ohno et al. 2001 assume a point-symmetric internal state; the restriction to the symmetric branch is derived, not stated there). Full-versus-half agreement ends at the first symmetry-breaking bifurcation. | Article §3, Eqs. 41–50 ([Güven & Özdemir 2026](https://doi.org/10.1016/j.finel.2026.104599)); Wallin & Tortorelli Eqs. 25, 32–35 and footnote 3 ([2020](https://doi.org/10.1016/j.mechmat.2020.103324)); [Ohno et al. 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) |
| Density-based topology-optimization methodologist | The result is a function of the optimization path. Rounding the objective, counting continuation per function evaluation, and switching between MMA87 and GCMMA each change the path, so none of them may change silently. Count continuation in design-loop iterations, as every source does, and compare designs after thresholding with a common evaluator. | [Stolpe & Svanberg 2001](https://doi.org/10.1007/s001580050177); [Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w); [Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2) (abstract only) |
| Design-sensitivity (adjoint) analyst | Keep the symmetric sensitivity; the article's printed one-sided Eq. 20 is wrong. Wallin & Tortorelli and the legacy code both solve one adjoint per tensor component. One adjoint per scalar response is a proposed improvement that is exact by linearity, but it must be shown to agree with the per-component path to round-off. Verify gradients on unrounded functions with Taylor-remainder slopes and the complex step. | Wallin & Tortorelli Eqs. 50–51 ([2020](https://doi.org/10.1016/j.mechmat.2020.103324)); [Tortorelli & Michaleris 1994](https://doi.org/10.1080/174159794088027573); [Farrell et al. 2013](https://doi.org/10.1137/120873558); [Martins et al. 2003](https://doi.org/10.1145/838250.838251) |
| Sparse-linear-algebra and floating-point specialist | Make the reduced tangent exactly symmetric by construction, choose Cholesky explicitly, and reuse factors. Never let Julia's `\` pick the solver from an exact symmetry test. Expect last-bit differences whenever operation order changes, and set tolerances from conditioning or measured noise, not from √eps. | [SparseArrays `\` dispatch source](https://github.com/JuliaSparse/SparseArrays.jl); [CHOLMOD, Chen et al. 2008](https://doi.org/10.1145/1391989.1391995); [UMFPACK, Davis 2004](https://doi.org/10.1145/992200.992205); [Higham 2002](https://doi.org/10.1137/1.9780898718027) (bibliographic details not re-verified) |
| V&V and scientific-software-testing engineer | Characterize before changing anything, and never refactor and change behavior in the same commit. Golden masters show "same as before"; analytic oracles, metamorphic relations and published benchmarks show "correct". When bitwise identity is lost, accept a change against the spread of round-off-perturbed trusted runs. | [Kanewala & Bieman 2014](https://doi.org/10.1016/j.infsof.2014.05.006); [Kanewala & Chen 2019](https://doi.org/10.1109/MCSE.2018.2875368); [Rosinski & Williamson 1997](https://doi.org/10.1137/S1064827594275534); [Baker et al. 2015](https://doi.org/10.5194/gmd-8-2829-2015) |
| Julia performance and package architect | The 500× kernel gap is architectural: untyped globals, closures and heap-allocated 2×2 tensors. Pass state through typed structs, put kernels on stack-allocated Tensors.jl types, replace two forked scripts with dispatch on cell symmetry, and enforce type stability and zero allocation in continuous integration (CI). | [Pelenitsyn et al. 2021](https://doi.org/10.1145/3485527); [Bezanson et al. 2017](https://doi.org/10.1137/141000671); [Carlsson & Ekre 2019](https://doi.org/10.5334/jors.182); [Verdugo & Badia 2022](https://doi.org/10.1016/j.cpc.2022.108341) |
| Research-reproducibility and publication-integrity specialist | Pin and archive the environment and give each release a DOI. Record each decision before running its evidence, and make every results-changing decision a row in the next paper's deviation table. Define "reproduced" explicitly: bitwise identity is the expectation only for identical computations. | [NASEM 2019](https://doi.org/10.17226/25303); [Sandve et al. 2013](https://doi.org/10.1371/journal.pcbi.1003285); [Smith et al. 2016](https://doi.org/10.7717/peerj-cs.86); [Barker et al. 2022](https://doi.org/10.1038/s41597-022-01710-x); [Haftka et al. 2019](https://doi.org/10.1007/s00158-019-02298-4) |

### The consolidated recommendation

The seven positions fit together because they act on different layers. The mechanician and the optimization methodologist define **what must not change**: the formulation and the path-defining semantics. The numerical analyst and the Julia architect define **how the arithmetic may change, and what each change costs in bitwise identity**. The V&V and reproducibility specialists define **how each change is proven and recorded**.

The consolidated strategy therefore has three layers:

| Layer | Content | Lifetime |
|---|---|---|
| **Semantic profile** named FEAD2026 | Every results-defining legacy behavior, bundled as a named switch | Permanent |
| **Numerics layer** | First replays the legacy arithmetic bit for bit; later moves to modern kernels under tolerance-based evidence | Legacy arithmetic retired after Phase 4 |
| **Decision trail** | Numbered records, each opened as "proposed" with acceptance criteria written down before any evidence is run | Permanent |

The order of work matters as much as the content. It follows Wilson et al.'s rule to optimize only after the code works correctly ([Wilson et al. 2014](https://doi.org/10.1371/journal.pbio.1001745)), and the separation of refactoring from behavior change that the characterization-testing literature recommends ([Kanewala & Bieman 2014](https://doi.org/10.1016/j.infsof.2014.05.006)). The sequence is:
1. Freeze and recover.
2. Characterize.
3. Restructure bitwise.
4. Take bitwise-preserving speedups.
5. Modernize the numerics to tolerance.
6. Introduce deliberate deviations.
7. Release.

The panel rejects four tempting shortcuts for the reproduction track:
- **A rewrite from the article.** The article has at least ten internal inconsistencies, and the code is the only executable specification.
- **A migration to Ferrite.jl.** Any library change renumbers degrees of freedom and reorders sums. Released Ferrite (v1.7.0 and earlier) also rejects constraints in which a degree of freedom is both master and slave ([Ferrite CHANGELOG](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md)).
- **A swap to NLopt's MMA.** It is an independent CCSA implementation and therefore a different algorithm ([NLopt algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/)).
- **Warm-started or inexact Newton as a default.** The published evidence reports only "very similar" optima, not identical ones ([Senne et al. 2022](https://doi.org/10.1007/s11081-022-09769-6)).

## Corrections the notes require before anything is preserved

The six research notes largely agree, but several statements in them, and in the coordinator's earlier review, are wrong, overstated or superseded by the Wallin & Tortorelli full text. The decision to preserve or change a behavior depends on these points, so the corrected positions are stated once here and used in every later section.

| Earlier statement | Corrected position | Evidence (status) |
|---|---|---|
| Coordinator item 7: evaluating the linear/nonlinear blend as H_steep(ρ̃^p), on the filtered, unprojected density raised to p, is a defect. | **Not a defect.** The thesis writes H_{β̂,ω̂}(ρ̃^p) with β̂ = 200 and ω̂ = 0.05. Wallin & Tortorelli's Eq. 73, the source of article Eq. 30, uses the same argument (filtered volume fraction raised to p) and the same threshold 0.05. Only the steepness differs: 500 there, 200 in the article, thesis and code, with no stated rationale. Changing the argument alters results (class A). | [Thesis p. 40, Eqs. 5.23–5.24](https://gcris.iyte.edu.tr/handle/11147/13904); [Wallin & Tortorelli 2020, p. 9](https://doi.org/10.1016/j.mechmat.2020.103324); [NEMOpt.jl L297–299](../Full_Domain/NEMOpt.jl#L297) (verified) |
| Coordinator item 8: Tables 3 and 4 both use ω as the projection steepness. | **Only Table 3** ramps ω from 1 to 32 as the steepness, with β = 0.5 as the threshold. Table 4 ramps β from 1 to 32 as the steepness, with ω = 0.5 as the threshold. That matches the code, the thesis and Wallin & Tortorelli (β steepness, ω threshold). The article also states that Table 3's r is in element sizes. | [Article pp. 12–13, Tables 3–4](https://doi.org/10.1016/j.finel.2026.104599); [NEMOpt.jl L53](../Full_Domain/NEMOpt.jl#L53) (verified) |
| The extra Newton correction is only a wasted factorization. | **It changes results.** `solve` returns the norm of the residual *before* the correction, and the correction is then applied unconditionally. So each load step ends with one more correction, reassembly and factorization than thesis Algorithm 5.2, which stops on the latest residual. The delivered U is more converged than TOL implies. Removing the extra correction is class A, and it plausibly degrades the ~1e−8 finite-difference agreement of Table 2 (inference). Wallin & Tortorelli do not describe their Newton solver, so the code is the only specification of the tolerance and the load stepping. | [NEMOpt.jl L365–367, L387–395](../Full_Domain/NEMOpt.jl#L365); [thesis p. 36, Algorithm 5.2](https://gcris.iyte.edu.tr/handle/11147/13904) (verified; the Table 2 link is inference) |
| Verification note: the "published-case configuration" is F_M = 1.025·I, 80×80, N_ITER = 120, V = 0.5, MMA87. | That is the **as-shipped full-domain driver, which matches no published case**. It reproduces the thesis Chapter 7 schedule. The half-domain driver is the closest to §5.2/Fig. 9 (same V = 0.4, same schedule, 300 iterations), but its F_M, mesh and filter radius differ. | [NEMOpt_main.jl L10–26](../Full_Domain/NEMOpt_main.jl#L10); [NEMOpt-DR_main.jl L8–24](../Domain_Reduction/NEMOpt-DR_main.jl#L8) (verified) |
| Verification note: the published runs used PNG initial designs (`ITER_1.png`). | **Unrecorded.** Both drivers default to `ITER_1.png`, but neither the article nor the repository says which image or seed produced any published result. Article §5.1 speaks of "same random initial densities". | [NEMOpt_main.jl L28–29](../Full_Domain/NEMOpt_main.jl#L28); [article §5.1](https://doi.org/10.1016/j.finel.2026.104599) (verified) |
| Verification note: the half-domain `K_Upper` applies the constraint only on the column side, which suggests a reduced system that is non-symmetric by construction. | **`K_Upper` is computed and never used.** The half-domain code solves `vcat(K_Middle, K_Lower)`, which has the same form as the full-domain system and is symmetric in exact arithmetic. The asymmetry is round-off (1.4e−14 measured for the full domain), which is enough to send `\` to UMFPACK LU. | [NEMOpt-DR.jl L572–575, 664–667, 804–807](../Domain_Reduction/NEMOpt-DR.jl#L572) (verified) |
| Architecture note: the Tensors.jl kernel gave "identical results", so expect bitwise or near-bitwise equality. | The benchmark checked agreement only to **rtol 1e−10**, for P and the tangent, not bitwise. Z = ∂𝔸/∂F was not cross-checked. Treat the kernel swap as a tolerance-level change that needs its own Z test. | [Benchmark L176–177](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl) (verified) |
| Verification-note gap: whether MMA87 calls the objective during inner iterations. | **Resolved.** NonconvexMMA evaluates value and Jacobian once at setup and once per inner iteration. MMA87 does two inner iterations in the first outer iteration and one thereafter; MMA02 repeats until its approximation is conservative or `maxinner` = 40 is reached. The article's "every 5 iterations" therefore corresponds in code to "every 5 objective evaluations, counting the setup evaluation". Under MMA87 that is offset by about two evaluations; under MMA02 it depends on how many iterations are rejected. | [NonconvexMMA `mma_algorithm.jl`](https://github.com/JuliaNonconvex/NonconvexMMA.jl/blob/master/src/mma_algorithm.jl); [`mma_approx.jl`](https://github.com/JuliaNonconvex/NonconvexMMA.jl/blob/master/src/mma_approx.jl) (verified by source reading) |
| Coordinator measurement: "4-RHS solve 0.61 s (LU each) vs 0.08 s (one Cholesky)". | The comparison **conflates two effects**: switching from LU to Cholesky, and reusing one factor for four right-hand sides. The 7.6× ratio is the combined gain, not a Cholesky-versus-LU ratio. | [Benchmark L137–138](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl) (verified) |
| Verification note: in legacy mode the refactored code "reproduces the FEAD 2026 runs bit-for-bit". | Bitwise identity can be proved **only against the legacy code at `30d0842` in a pinned environment**. The published runs used unrecorded inputs and environments. Commit `30d0842` dates from 2026-05-07, after the 6 May 2026 revision. The §5.1 and §5.3 schedules appear in no executable line, which suggests some figures came from earlier script versions (inference). | [GitHub repository](https://github.com/likemaestro/NEMOpt); [article history](https://doi.org/10.1016/j.finel.2026.104599) (dates verified; provenance is inference) |
| Several notes: the Wallin & Tortorelli full text is inaccessible, so its switch variable, adjoints and settings are unverified. | **Superseded.** The paper is open access (CC BY-NC-ND). Its full text was read and compared item by item. | [Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324); [comparison note](../research_notes/NEMOpt%20refactoring%20strategy%20literature/wallin_tortorelli_2020.md) (verified) |
| Numerical-choices note: the clip `max(H^p, ε)` is "neither classic nor modified SIMP", which makes it a NEMOpt-specific choice, to be switched to modified SIMP by default. | It is **Wallin & Tortorelli's ersatz-material law (Eq. 71) verbatim**: the reference formulation, not a NEMOpt invention. Modified SIMP is therefore an *optional deviation* behind a switch, not a bug fix. Only the floor differs: ε = 1e−5 there and 5e−4 in the code, a value stated in neither the article nor the thesis. The larger floor moves the clip boundary in H from 0.0215 to 0.0794 at p = 3 and enlarges the zero-gradient zone. | [Wallin & Tortorelli 2020, pp. 9, 13](https://doi.org/10.1016/j.mechmat.2020.103324); [NEMOpt.jl L59–60](../Full_Domain/NEMOpt.jl#L59) (verified; the clip values are computed) |
| Numerical-choices note: make a physical-density void switch the new default, because ρ̃^p departs from Wang et al. (2014), and treat the band of soft hyperelastic elements as a NEMOpt problem. | **Re-attributed.** The departure from Wang et al. is Wallin & Tortorelli's own choice: their Eq. 73 uses the filtered ν̃^p. The "band of soft hyperelastic elements" exists in their scheme too: ν̃ between 0.368 and 0.452 at their final p = 3 and β = 40, against 0.368 to 0.462 in NEMOpt (computed). It is not a NEMOpt defect. Changing the switch would depart from both sources, so it stays the default and the alternative is offered only as an option. | [Wallin & Tortorelli 2020, pp. 2, 9](https://doi.org/10.1016/j.mechmat.2020.103324); [Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4) |
| By analogy with Eq. 16, the one-sided article Eq. 20 might be taken as equivalent to the code's symmetric sensitivity. | **It is not; the code is correct and Eq. 20 is a misprint.** The one-sided 𝔸_M of Eq. 16 does equal the symmetric form at χ-equilibrium (Wallin & Tortorelli Eq. 24). That identity does not carry over to the sensitivity. On Wallin & Tortorelli's benchmark, Eq. 20 gives ∂A₁₁₁₁/∂ν₁ = 10.146857 instead of the printed 26.548821. The symmetric form used by the code, by thesis Eq. 4.25 and by Wallin & Tortorelli Eq. 51 reproduces every Table 1 entry, as do central finite differences; the coordinator confirmed this by an independent re-run. The article's Table 2 therefore validates the code, not the printed equation. Whether to issue a corrigendum is the authors' decision. | [Thesis p. 29](https://gcris.iyte.edu.tr/handle/11147/13904); [Wallin & Tortorelli 2020, pp. 5–7, 12](https://doi.org/10.1016/j.mechmat.2020.103324); [evidence scripts](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/wallin_tortorelli_2020/) (verified against published values) |

## "Same results" is a four-link chain, not a single test

The constraint that the refactored code "provably produce the same results as the article" can be met only by splitting it into links, each with its own kind of proof.

**Link E1** states that the refactored code, running in the FEAD2026 profile with the legacy arithmetic, is bitwise identical to commit `30d0842` in the pinned reference environment. That is a hard proof by hash equality of every characterization output and every iterate.

**Link E2** states that the refactored code with modern arithmetic but FEAD2026 semantics stays within the measured round-off noise of the legacy code. This holds per quantity at the kernel, evaluation and sensitivity levels, and inside an ensemble of perturbed legacy runs at the trajectory level. That is a statistical proof, the same logic climate modelers use to validate ports ([Rosinski & Williamson 1997](https://doi.org/10.1137/S1064827594275534); [Baker et al. 2015](https://doi.org/10.5194/gmd-8-2829-2015)).

**Link E3** states that the legacy code itself reproduces the article: the Table 2 accuracy claim in kind, and the reported ν and κ values within pre-registered bands. This link is evidential, because the article's inputs are partly unknown.

**Link E4** covers every deliberate deviation: each one is re-run on the article's cases, and its change against the published values is reported.

The National Academies explicitly endorse relaxing bitwise identity to "reproducible results within an accepted range of variation" once the numerics are understood (Finding 4-3 of [NASEM 2019](https://doi.org/10.17226/25303)). E1 and E2 are genuine proofs of equivalence to the code. E3 is the strongest claim available about the article, and Phase 0's recovery of the authors' original inputs is the only way to harden it.

Two tolerance regimes follow from this chain.

**Regime B (bitwise)** applies to pure refactoring in the reference environment: same Julia, SuiteSparse and BLAS versions, `BLAS.set_num_threads(1)`, no `@fastmath`, and no new `muladd` reassociation. Only hash equality passes.

**Regime T (tolerance)** applies to any deliberate change of operation order, solver, kernel or library. Its bands are set empirically, by running the legacy code under round-off-sized perturbations (±1 ulp on the design, LU versus Cholesky, BLAS thread count). They are cross-checked against the a priori bound that relative forward error scales roughly as κ(K)·u.

Julia's default `isapprox` tolerance (√eps ≈ 1.5e−8) cannot be used for either regime. It is too loose to detect a broken pure refactor, and its norm-wise test can hide large relative errors in near-zero shear entries of 𝔸 ([Julia `isapprox`](https://docs.julialang.org/en/v1/base/math/)). Each band must therefore be stated per quantity. A change in Newton iteration count near the 1e−7 threshold is an expected Regime-T event to log, not a failure in itself.

| Level | Quantities compared | FEAD2026 profile, legacy arithmetic (Regime B) | Modern arithmetic (Regime T) | Independent oracles that stay valid after changes | Basis |
|---|---|---|---|---|---|
| L0 environment | Julia version, Manifest hash, BLAS vendor and threads, SuiteSparse version, solver path taken | Identical | Recorded on every run | — | [Sandve et al. 2013](https://doi.org/10.1371/journal.pcbi.1003285); [NASEM 2019](https://doi.org/10.17226/25303) |
| L1 kernels | P, 𝔸, Z, H, H′, g, g′, Ĥ, filter element matrices, node numbering, D1/D2 | Hash equality | Per-entry error within a small multiple of n·u; P and 𝔸 against AD to the measured 1e−10 | 𝔸(F = I) = ℂ; Z against central differences of 𝔸; objectivity; the thesis' printed 8×8 numberings | [Higham 2002](https://doi.org/10.1137/1.9780898718027); [Kanewala & Bieman 2014](https://doi.org/10.1016/j.infsof.2014.05.006) |
| L2 single evaluation | U, four χ fields, P_M, 𝔸_M, ν, κ, Vf, Newton iteration counts and residuals | Hash equality | Within the measured legacy noise band and the κ-based bound | Uniform solid gives the exact 𝔸_m(F_M); uniform void gives εℂ; power-of-two scaling of E; periodic shift; 90° rotation; full ↔ half domain modulo a half-period shift; Wallin & Tortorelli §4.1 cases A and B | [Rosinski & Williamson 1997](https://doi.org/10.1137/S1064827594275534); [Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324) |
| L3 sensitivities | ∂𝔸_M/∂ρ (all 16 components), ∇J, ∇g | Hash equality | L2 band propagated; aggregated versus per-component adjoint to about 1e−12 | Wallin & Tortorelli Table 1 (printed digits); article Table 1/2 replication; Taylor-remainder slope 2; complex step; dot-product tests | [Farrell et al. 2013](https://doi.org/10.1137/120873558); [Martins et al. 2003](https://doi.org/10.1145/838250.838251) |
| L4 full runs | J, Vf, ν, κ histories; iterates; MMA counters; final design | Bitwise trajectory replay | Inside an ensemble of round-off-perturbed legacy runs, on pre-registered metrics | Article values at printed precision; thresholded re-analysis, measure of non-discreteness (M_nd), intersection-over-union (IoU); Voigt-type energy bound at F_M; Hashin–Shtrikman bulk bound only on the linear-limit tangent 𝔸_M(F_M = I) of thresholded designs (NDR-0064) | [NASEM 2019](https://doi.org/10.17226/25303); [Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2) (abstract only); [Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x) |

At the kernel and evaluation levels, several analytic oracles are exact:
- **Uniform solid (ρ ≡ 1).** The code's own definitions give g = H(1)^p = 1 and Ĥ(1) = 1. The exact periodic solution is the homogeneous deformation, so 𝔸_M must equal the neo-Hookean first elasticity tensor at F_M, and P_M the first Piola–Kirchhoff stress (article Eqs. 10, 16, 24–25).
- **Uniform void (ρ ≡ 0).** Both tanh projections vanish exactly, because tanh(βω) + tanh(−βω) = 0, so 𝔸_M = εℂ to round-off.
- **Doubling E.** This is exact in binary floating point and must leave ν bit-identical. It also exposes the legacy stopping test, whose scale `Rep_MU = |max(F_M − I)|·Lx` is a displacement while the residual is a force ([NEMOpt.jl L375, L388](../Full_Domain/NEMOpt.jl#L375)).
- **Full versus half domain.** The half-domain physical design equals the full-domain design translated by half a period (Ly/2); this was verified by emulating both `Generate_Densities` routines. Tests must therefore compare 𝔸_M, P_M, ν, Vf and gradients, not raw arrays. They must also use a shear F_M, so that mistakes in the odd/even parity rules cannot pass silently ([NEMOpt.jl L822–834](../Full_Domain/NEMOpt.jl#L822); [NEMOpt-DR.jl L178–184](../Domain_Reduction/NEMOpt-DR.jl#L178)).

The sensitivity level has two published anchors of very different strength.

**The article's own check** fixes these parameters:

| Parameter | Value |
|---|---|
| F_M | [1.05 0.05; 0.05 1.05] |
| L₁, L₂ | 3, 3 |
| E, ν | 100, 0.3 |
| r | 2 |
| p | 3 |
| (β, ω) | (32, 0.5) |
| (β̂, ω̂) | (200, 0.05) |
| Mesh | 4×4, with four design variables |
| Finite difference | Forward, step 1e−6 |

Table 2 reports adjoint-minus-finite-difference differences no larger than **5.59858e−8** (ρ₃, component 1111). The thesis prints the identical numbers but labels them "δA(FDM) − δA(FEM)", which flips the sign convention, so acceptance should compare magnitudes ([article §4.1](https://doi.org/10.1016/j.finel.2026.104599); [thesis pp. 42–43, Tables 5.1–5.2](https://gcris.iyte.edu.tr/handle/11147/13904)). Neither source gives ρ₁…ρ₄, and the unit of Table 1's r (elements or absolute length) is not stated. The replication must therefore:
- use documented densities;
- test both interpretations of r;
- call `Stiffness_Sensitivity_Unfiltered` for all 16 components, because `Objective_Grad` computes only δ𝔸₁₁₁₁, δ𝔸₂₂₂₂ and δ𝔸₁₁₂₂ ([NEMOpt_main.jl L109](../Full_Domain/NEMOpt_main.jl#L109)).

**Wallin & Tortorelli's §4.1 case** is the stronger, fully reproducible test and should become a golden test. Its specification:
- a 4×4 periodic cell on [−1, 1]²;
- element groups with the given volume fractions 0.3, 0.5, 0.7 and 0.9;
- a pure SIMP law S = ν³∇Ψ̃, with no filter, projection or blend;
- F_M = diag(1.3, 0.8) (case A), plus a uniaxial case B;
- their 2D energy;
- central differences with step 1e−5 and a reported agreement of six significant figures.

The scripts in the evidence folder reproduce the published values, and the coordinator re-ran them independently: S_M; 𝔸_M; the case-B values F₂₂ = 0.872199, S₁₁ = 1.43773, E₁ = 3.28191 and ν₁₂ = 0.356115; and every Table 1 sensitivity, by both central differences and the symmetric adjoint. The reproduction requires Q4 elements, 2×2 Gauss points, a fixed center node, and B = μ/8 exactly as printed. That value is a typo relative to the consistent μ/2, which their own A = (λ + μ)/4 implies, but the published numbers reproduce only with μ/8, so the benchmark must use it uncorrected ([Wallin & Tortorelli 2020, pp. 10–12](https://doi.org/10.1016/j.mechmat.2020.103324); [evidence scripts](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/wallin_tortorelli_2020/)).

The case has three useful properties. It contains a negative tangent component, 𝔸₁₂₁₂ = −0.9861, which exercises sign handling in Z. Its pure power law satisfies a degree-3 Euler identity, Σνᵢ∂(·)/∂νᵢ = 3(·), which is a cheap invariant test. It needs only a pluggable energy, a pure-SIMP mode with a filter bypass, centroid-based group assignment, and the index map from Wallin & Tortorelli's (11, 21, 12, 22) ordering to NEMOpt's Voigt positions (1, 4, 3, 2).

Independent verification also needs pure, unrounded functions ρ ↦ 𝔸_M and ρ ↦ J. The legacy objective is quantized at 1e−3, which makes finite differences below that step meaningless. The `max(H^p, ε)` clip and the rounding are also exactly the non-analytic constructs that complex-step differentiation must avoid ([Martins et al. 2003](https://doi.org/10.1145/838250.838251)).

At the trajectory level, link E3 needs explicit acceptance targets:

| Case | Load | Reported value |
|---|---|---|
| §5.2 case 1 | F_M = 1.075·I; V ≤ 0.4; 100×100 | ν = −0.9008 |
| §5.2 case 2 | F_M = [1.025 0.025; 0.025 1.025]; V ≤ 0.4; 100×100 | ν = −0.9132 |
| §5.3, three cases | V ≤ 0.6 | κ = 28.4368, 27.9679 and 24.1839 |

These come from [article §5.2–5.3](https://doi.org/10.1016/j.finel.2026.104599). From Eq. 55, the two ν targets imply objective values of about 0.984 and 0.753 (inference). One plausibility check supports the code's bulk definition: κ = 28.4368 is 98.6% of the two-dimensional Hashin–Shtrikman upper bound of 28.85 for a porous solid at 60% solid ([Bendsøe & Sigmund 1999, Eq. 6](https://link.springer.com/doi/10.1007/s004190050248); inference, from λ + μ = 96.15). This comparison is meaningful only because §5.3 case 1 (F_M = 1.0001·I) is essentially linear. Hashin–Shtrikman bounds are theorems for linear elastic moduli, so no such bound applies to the finite-strain tangent values of the other cases. Three further points apply:
- the ε = 5e−4 floor raises the bound to 28.905 (computed);
- p = 3 SIMP itself exceeds the bulk bound by up to 1.6% in gray elements under plane strain (computed). This is expected, because p = 3 is below the 2D threshold p\* = 3.5 of Bendsøe & Sigmund (1999, Eq. 11) for the in-plane Poisson's ratio 0.4286;
- the bound applies to the hydrostatic-strain modulus κ = ¼(𝔸₁₁₁₁ + 𝔸₂₂₂₂ + 𝔸₁₁₂₂ + 𝔸₂₂₁₁) of any anisotropic two-phase cell, square-symmetric or not. Milton (2002, §23.4, Eq. 23.30, p. 466 of the author's posted copy) gives the optimal 2D upper energy bound for any geometry of two well-ordered isotropic phases, and at hydrostatic strain it reduces to the 2D Hashin–Shtrikman bulk bound; it reproduces 28.85 and 28.905, with the void as the limit of the ε-material ([Milton, *The Theory of Composites*](https://gwmilton.com/TheoryCompositesNOPRINT.pdf); verified). Bendsøe & Sigmund's (1999) strain-energy bound for void/solid composites, Eq. 21, gives the same result at hydrostatic strain (derived).

NDR-0064 records how bounds are used from now on.

Starting bands of ±0.02 on ν and 1–2% on κ are experience values. They must be replaced by the measured spread of the legacy ensemble before any comparison is run, because the thesis itself notes that a different initial design can lead to a different optimum topology ([thesis §7.1.2, pp. 54–55](https://gcris.iyte.edu.tr/handle/11147/13904)). Sigmund notes that, without practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors" ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), abstract; full text not read). This report therefore recommends re-analyzing final designs with one frozen evaluator, thresholding them, and comparing them on M_nd ([Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x)) and IoU.

Four items cannot be turned into tests:
- The linear-response value ν = −0.9218 needs the absent LEMOpt code.
- The Fig. 7 timings are unlabeled bars on unstated hardware.
- The thesis results form a third parameter family (r = 1.5, V = 0.5) and must never be mixed into the article's golden masters.
- Wallin & Tortorelli's design examples publish figures only. They also differ in objective, Poisson definition, energy and kinematics, so they are not quantitative targets.

## Seven verified behaviors change results and must become switches

The reproduction-baseline audit classified fifteen code findings. Seven of them change numerical results if "fixed" ([NEMOpt_main.jl](../Full_Domain/NEMOpt_main.jl#L44); [NEMOpt.jl](../Full_Domain/NEMOpt.jl#L297)):
1. the output rounding fed to MMA;
2. continuation counted per function evaluation;
3. the Newton stopping logic;
4. MMA87 versus MMA02;
5. image ingestion;
6. the blend argument;
7. the MMA options.

The Wallin & Tortorelli comparison sorts these behaviors into two groups. Some are **inherited reference formulation**: the clipped SIMP law, the ρ̃^p blend and its soft-element band, and plain MMA in the full domain. Others **remain likely defects** that no source supports:
- the rounding of J and ν before MMA;
- p and β changed inside the objective before the gradient is computed;
- continuation counted per function evaluation;
- the unused `options` local;
- the full-domain schedule's β overshoot to 33.

The rule that follows is the same for both groups: **every class-A behavior is preserved verbatim in the FEAD2026 profile, and any new default is introduced only through a decision record that names its switch.** The column "Effect if changed" separates real results changes from logging or maintainability issues.

| Behavior (location, `30d0842`) | Status and origin | Effect if changed | FEAD2026 profile | Proposed option or default |
|---|---|---|---|---|
| J and ν rounded to 3 decimals; the rounded J is returned to MMA and the rounded ν enters ∇J = 200(ν − ν₀)∇ν ([L79–80, L104, L117](../Full_Domain/NEMOpt_main.jl#L79)) | Verified likely defect; in no source | **A.** The returned value carries ±5e−4 noise, far above the slack epsimin = 1e−7 that the GPL Python port of Svanberg's MATLAB code adds in its GCMMA conservativeness check (also the port's subproblem tolerance; Svanberg's 2007 notes not re-read); NonconvexMMA's check has no tolerance at all ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python), `concheck` and `examples/gcmma_beam.py`). The gradient factor error is up to 5e−4/\|ν+1\|, about 0.5% at ν = −0.90. With kkt = 1e−8, KKT termination becomes practically unattainable (inference). | Keep | Unrounded by default; round only in logs |
| p/β continuation advanced per objective evaluation, inside `Objective`, before `Objective_Grad` runs in the same rrule. The full-domain schedule (`β < 32`, then +2 from odd values) overshoots to β = 33 ([L45, L90–102, L213–216](../Full_Domain/NEMOpt_main.jl#L90)). | Verified likely defects; the article, thesis and Wallin & Tortorelli all count design iterations and state terminal values | **A.** At each step, J uses (p_old, β_old) but the gradient uses (p_new, β_new). This happens on about 50 of 300 half-domain evaluations (computed). The schedule also depends on how often the optimizer rejects iterations. | Keep, including the overshoot in the as-shipped preset | Per MMA outer iteration, applied between iterations; schedules stored as `{istart, max, every, delta}` data with explicit terminal values ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)) |
| Newton control: residual tested before the update (one extra correction per load step); `Rep_MU = \|max(F_M − I)\|·Lx`, which is 0 for any compressive F_M and then loops forever; no iteration cap or line search; 10 fixed load steps from U = 0 ([NEMOpt.jl L372–396](../Full_Domain/NEMOpt.jl#L372)) | Verified; Wallin & Tortorelli describe no solver; differs from thesis Algorithm 5.2 | **A** for moving the check, changing the scale or warm-starting. **B** for guards that fire only where legacy would hang or diverge. For F_M = 1.0001·I the absolute target is 2e−11, a plausible but unverified hang risk. | Keep, plus an explicit error instead of a hang | Check after the update; a scale that is never zero; tighter tolerance; iteration cap; bisection fallback on failure |
| Optimizer: MMA87 (full domain) versus MMA02 (half domain); `s_init` 0.25 and `maxinner` 40, against library defaults of 0.5 and 10 ([L235–237](../Full_Domain/NEMOpt_main.jl#L235); [DR L234](../Domain_Reduction/NEMOpt-DR_main.jl#L234)) | Verified. MMA87 matches Wallin & Tortorelli's plain MMA with "the numerical parameter values provided in Svanberg (1987)", wording the article echoes. MMA02 follows the thesis. No source documents `s_init` or `maxinner`. | **A.** The two shipped drivers cannot produce the "same optimum" of Fig. 8 under different optimizers. | Keep per preset | Explicit configuration; identical optimizer for any full-versus-half comparison |
| Stiffness interpolation `max(H^p, ε)` with zero derivative below ε = 5e−4 ([L59–60](../Full_Domain/NEMOpt.jl#L59)) | Verified reference formulation: Wallin & Tortorelli Eq. 71 (ε = 1e−5 there). The code's ε is stated in neither the article nor the thesis. Classic SIMP instead floors the density variable, 0 < ρ_min ≤ ρ ([Bendsøe & Sigmund 1999, p. 638](https://link.springer.com/doi/10.1007/s004190050248); verified). | **A.** The clip creates a non-C¹ kink and a zero-sensitivity zone, while GCMMA's convergence theory assumes continuously differentiable functions ([Svanberg 2002](https://doi.org/10.1137/S1052623499362822); assumption not re-verified). Stiffness differs from modified SIMP by at most ε. | Keep, with ε recorded | **Keep the clip as the default.** Offer modified SIMP ε + (1 − ε)H^p ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)) only as an opt-in deviation behind the `stiffness_interp` switch, recorded as a deviation from both Wallin & Tortorelli and FEAD 2026, not as a bug fix. |
| Void switch Ĥ(ρ̃^p) with β̂ = 200, ω̂ = 0.05 ([L297–299](../Full_Domain/NEMOpt.jl#L297)) | Verified reference formulation. Matches the thesis and Wallin & Tortorelli Eq. 73 in argument and threshold; steepness is 200 against their 500, with no rationale anywhere. | **A** | Keep | Keep as default. Offer β̂ = 500 and a physical-density switch, Ĥ(H(ρ̃)^p) with (500, 0.01) ([Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4); [Dalklint et al. 2021, Eq. 9](https://doi.org/10.1007/s00158-021-03012-z), from Wallin and Tortorelli's group), as exploratory options. Add distortion diagnostics. |
| Image ingestion: crop rows `floor(n/2):n` give one extra row; `imresize` interpolation; rescaling by `min(V_MAX/Vf, 1)` ([L780–798](../Full_Domain/NEMOpt.jl#L780)) | Verified | **A** with respect to the image and package versions | Keep, and store the ingested design vectors as hashed fixtures | Corrected loader as an option |
| Unused local `options = MMAOptions(s_init = 0.5/(β+1), …)` inside `Objective` ([L100](../Full_Domain/NEMOpt_main.jl#L100)) | Verified likely defect (code reading, standard Julia scoping): the intended β-dependent update never reaches the optimizer | None today; actually implementing it would be A | Preserve the no-op | Delete; any β-scaled move limits ([Guest et al. 2011](https://doi.org/10.1007/s00158-011-0676-1)) need a new record |
| `Volume_Constraint_Grad` reads globals set by `Objective`; `Volume_Constraint` repeats the filter solve ([L128–164](../Full_Domain/NEMOpt_main.jl#L128)) | Verified | **B** today, because NonconvexCore evaluates the objective before the constraint at the same x ([NonconvexCore `value_jacobian.jl`](https://github.com/JuliaNonconvex/NonconvexCore.jl/blob/master/src/functions/value_jacobian.jl)). Latent A if the call order changes. | Preserve the call order | Explicit state |
| ρ̃^p with fractional p throws `DomainError` if ρ̃ < 0 ([L297, L618](../Full_Domain/NEMOpt.jl#L297)) | Verified in code; no filter undershoot observed in tests | Latent crash | Detect and raise | Detect and raise; clamp only by decision |
| Logged "G": A₁₂₁₂ + A₂₁₂₁ in the full domain versus A₁₂₁₂ in the half domain ([L71](../Full_Domain/NEMOpt.jl#L71); [DR L67](../Domain_Reduction/NEMOpt-DR.jl#L67)) | Verified. The half-domain value equals Wallin & Tortorelli's μ₁₂; because 𝔸 has no minor symmetry, the two differ under prestress. | **B** (logging only) | Rename both | Rename both |
| Square-mesh index assumptions; transposed Jacobian in the half-domain filter ([L692](../Full_Domain/NEMOpt.jl#L692); [DR L388](../Domain_Reduction/NEMOpt-DR.jl#L388)) | Verified | **B** for square, axis-aligned meshes, which covers every published case | Assert | Generalize |
| Eqs. 53 and 56 commented out; the Eq. 56 gradient references δA₂₂₁₁, which is never computed ([L70–72, L109–121](../Full_Domain/NEMOpt_main.jl#L70)) | Verified | Needed for the §5.1 and §5.3 presets | Not active | Objective types; δ𝔸₂₂₁₁ := δ𝔸₁₁₂₂ by major symmetry (inference; verify numerically) |
| No Project/Manifest; `cd` into the output folder; 𝜌 (U+1D70C) used in place of ρ (U+03C1); dead code ([.gitignore L24](../.gitignore#L24)) | Verified | **B** | Fix | Fix |

The article, thesis, code and Wallin & Tortorelli disagree in ways that decide which value a reproduction uses. The governing rule is: **for formulas, code first, thesis second and article third; Wallin & Tortorelli where all three are silent; for case parameters, the article's own table or text.** Every article-internal error is listed as an erratum; none of them changes what the code computes.

| Topic | Article | Thesis | Code `30d0842` | Wallin & Tortorelli 2020 | Treat as |
|---|---|---|---|---|---|
| Strain energy | Eq. 22 prints (1/μ)(I₁ − 2) | μ/2 (Eq. 5.4, p. 32), "following Bonet & Wood" | μ/2, implied by P | Split neo-Hookean (Eq. 97); 2D test prints B = μ/8, a typo for μ/2 that their published numbers nonetheless use | μ/2; erratum. Their energy is a deliberate different choice. |
| Projection names | α/β (Eq. 28, with an α(α − β) typo; Table 1); **Table 3 alone** uses ω as steepness | β steepness, ω threshold (Eq. 5.21) | Same as the thesis | Same (Eq. 72) | Code/thesis names; values from each case's table |
| SIMP law and floor | Eq. 29 `max`; ε not stated | Eq. 5.20; ε not stated | `max`, ε = 5e−4 | Eq. 71, ε = 1e−5 | `max` form as the reference default; ε = 5e−4 recorded as a FEAD2026 value; modified SIMP opt-in only |
| Blend argument and steepness | Eq. 30 omits the argument; (200, 0.05) | ρ̃^p, (200, 0.05) (p. 40) | ρ̃^p, (200, 0.05) | ν̃^p, (500, 0.05) (Eq. 73) | ρ̃^p with (200, 0.05); not a defect |
| Sensitivity form | One-sided Eq. 20 | Symmetric (Eq. 4.25, p. 29) | Symmetric | Symmetric (Eq. 51) | Symmetric. Eq. 20 is a misprint with an O(1) error; a corrigendum is the authors' decision. |
| Displacement adjoints | Eq. 19, one per component | Eq. 4.21 | One per component | One per component (Eq. 50), no aggregation | Per component in FEAD2026; aggregation is a proposed T-class improvement |
| Newton residual sign | R = F_ext − F_int, but RHS −(R_I + D₁ᵀR_D) | Eq. 5.15 K·dU = R is correct; the reduced Eq. 5.16 repeats the extra minus (p. 37) | R = F_int, which gives the correct Newton step | Not described | Code; erratum |
| Newton control | — | Algorithm 5.2 stops on the latest residual (p. 36) | One extra correction per load step | Not described | Code, for reproduction (class A) |
| Volume normalization | Eq. 33 and A.1 omit 1/Ω; Eqs. 53, 55, 57 include it | 1/\|Ω\| | 1/(Lx·Ly) | 1/\|Ω\| (Eq. 96) | Normalized |
| Table 2 label | "Adj − FD" | "FDM − FEM" for the same numbers (p. 42) | — | — | Compare magnitudes |
| Filter radius | Tables 3–4: r = 1 element; §5.2–5.3 text: 2.0 × Lx/100 | r = 1.5 | r = 1 element | r = 2 elements, the same as the article's text value | Unresolved; test r = 1 and r = 2 elements |
| Continuation | Per "optimization iteration": §5.1 p +0.2/5; §5.2 p +0.1/5 then β +1/5; §5.3 p +0.5/5 then β +4/5 | p +0.5/5, β +2/10 (p. 49) | Per evaluation; full domain p +0.5/5, β +2/10 (overshoots to 33); half domain p +0.1/5, β +1/5 | p +0.05/5; β doubled every 20 to 40; per design iteration | The article's text per section, with per-evaluation counting in FEAD2026 |
| Optimizer | MMA with "default parameters" [32] | GCMMA, Svanberg defaults | MMA87 / MMA02, `s_init` 0.25, `maxinner` 40 | Plain MMA, Svanberg (1987) parameters; nothing else stated | Per preset; record `s_init` and `maxinner` as undocumented values |
| Poisson's ratio | Eq. 54, an objective toward −1; the article credits Kazemi et al. (2020), whose Eq. 24 (p. 8) is a different, 3D linear formula for cubic lattices, ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)], so the attribution is incorrect ([Kazemi, Vaziri & Norato 2020](https://doi.org/10.1016/j.cma.2020.112895); verified) | Eq. 7.5, with the same incorrect attribution | Eq. 54 | Tangent uniaxial-stress −S₂₂₁₁E₁ (Eqs. 69, 96), as a constraint ν ≤ ν₀ | Eq. 54, a deliberate and documented deviation. On their uniaxial state, Eq. 54 gives 0.5596 against their 0.3561 (computed), so the two must never be compared numerically. Whether to correct the attribution in an erratum is the authors' decision. |
| Objective rounding; box constraint | Not mentioned; ρ_min ≤ ρ ≤ 1 | — | 3 decimals; [0, 1] | —; lower bound via ε | Code |
| Linear-response comparison (ν = −0.9218) | Reported | LEMOpt.jl named | Absent | — | Not reproducible from the repository |

## Most of the 49-second evaluation is overhead the literature already removes

The coordinator's instrumented run of the unmodified kernels (Julia 1.11.7, four cores) fixes the baseline:

| Measurement | Value |
|---|---|
| 40×40 full-domain evaluation | **41.8 s objective + 7.4 s gradient** |
| Constitutive calls per Gauss point | **87**: 60 Newton, 4 χ, 17 homogenization, 6 gradient |
| Sparse factorizations per evaluation | **59**: 50 Newton over 10 load steps, 4 χ, 1 filter, 4 gradient; all UMFPACK LU, because the reduced matrix is asymmetric at 1.4e−14 |
| One constitutive call | About 40 µs and 21.9 KB allocated, against 0.08 µs for a Tensors.jl automatic-differentiation equivalent |
| Tangent assembly at 80×80 | 4.57 s with incremental insertion, against 1.36 s with triplets |

These numbers come from the [coordinator benchmark](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl). Multiplying out (87 calls × 6,400 Gauss points × 40 µs) puts **about 22 s, close to half the evaluation, in the constitutive kernel alone** (inference).

The literature diagnoses each remaining cost:
- Preallocation and vectorized triplet assembly gave top88 its hundredfold speedup. The same paper symmetrizes K because round-off asymmetry otherwise forces MATLAB from Cholesky onto slower LU ([Andreassen et al. 2011](https://doi.org/10.1007/s00158-010-0594-7)).
- Assembling one triangle cuts assembly time and memory by about half ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)).
- A preallocated assembler inserts an element matrix in 83.4 ns, against 606.4 µs for naive `K[dofs,dofs] += Ke` ([Ferrite assembly docs](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/assembly/)).
- In geometrically nonlinear topology optimization, the nonlinear solve and objective account for 74.5–77.1% of runtime ([Senne et al. 2022](https://doi.org/10.1007/s11081-022-09769-6)).

The changes fall into three classes with very different proof burdens. Class B is the most valuable for this project: those changes remove work without changing any arithmetic operation that produces a result, so they can be proven by hash equality. The bitwise claims in the table are inferences from code reading. Each must be confirmed by a hash test, and any change that fails is reclassified as tolerance-level.

| Change | Literature support | Baseline evidence | Class | Expected effect |
|---|---|---|---|---|
| Assemble and factor the Helmholtz filter once per run; share it across objective, constraint and gradients | top88's PDE variant factors once outside the loop ([Andreassen et al. 2011](https://doi.org/10.1007/s00158-010-0594-7)); [Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072) (filter definition; full text not read) | Reassembled on every call, and repeated inside `Volume_Constraint` ([L128–130](../Full_Domain/NEMOpt_main.jl#L128)) | B (same factorization algorithm) | Removes every per-evaluation filter assembly and factorization |
| Reuse the final Newton assembly, which is already the tangent at the converged U, and one factorization for the four χ solves and all adjoint solves | Multiple right-hand sides against one factorization ([Behrou et al. 2021](https://doi.org/10.1016/j.compstruct.2021.113729)) | `χ_Solver` and `FindTrialU` reassemble the same matrix at the same U ([L401–462](../Full_Domain/NEMOpt.jl#L401), [L520–596](../Full_Domain/NEMOpt.jl#L520)) | B (same matrix, same UMFPACK factors; inference) | Removes 7 assemblies and 6 factorizations per evaluation |
| Cache the converged-state Gauss-point constitutive results for homogenization and sensitivities | Deterministic recomputation | 17 + 6 calls per Gauss point ([benchmark](../research_notes/NEMOpt%20refactoring%20strategy%20literature/evidence/bench_legacy_kernels.jl)) | B | Constitutive calls drop from 87 to about 61 per Gauss point |
| Move the PNG and text output out of the objective | [Wilson et al. 2014](https://doi.org/10.1371/journal.pbio.1001745) | A plot on every evaluation ([L75](../Full_Domain/NEMOpt_main.jl#L75)) | B | Wall time |
| Triplet assembly instead of incremental CSC insertion | [Andreassen et al. 2011](https://doi.org/10.1007/s00158-010-0594-7); [Engblom & Lukarski 2016](https://doi.org/10.1016/j.parco.2016.04.001) | 4.57 s → 1.36 s at 80×80 | B if duplicates are summed in element order (inference), otherwise T | About 3.4× faster assembly |
| Tensors.jl constitutive kernel (AD for P and 𝔸; hand-coded Z tested against AD and finite differences) | [Carlsson & Ekre 2019](https://doi.org/10.5334/jors.182); [Pelenitsyn et al. 2021](https://doi.org/10.1145/3485527) | 40 µs → 0.08 µs, agreement to 1e−10 | T | About 500× per call |
| Preallocated one-triangle assembly in reduced master–slave numbering, wrapped as `Symmetric` | [Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w); [Ferrite assembler](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/assembly/) | 1.4e−14 asymmetry from the Tᵀ K T products | T | Exact symmetry by construction; no sparse matrix products |
| Explicit CHOLMOD Cholesky with symbolic reuse (`cholesky!`) and a logged LDLᵀ/LU fallback | [Chen et al. 2008](https://doi.org/10.1145/1391989.1391995); [SparseArrays source](https://github.com/JuliaSparse/SparseArrays.jl) | Every solve uses LU today | T | 7.6× combined with reuse at 4 RHS (measured); Cholesky alone not yet measured |
| Aggregated adjoint: one state adjoint and one filter back-substitution per scalar function, with weights ∂J/∂𝔸 | One adjoint per response function ([Tortorelli & Michaleris 1994](https://doi.org/10.1080/174159794088027573)). Wallin & Tortorelli and the code both solve one per component (Eq. 50), so this is a proposed improvement, exact by linearity of the adjoint in its right-hand side. | Three adjoint solves for ν, two for Eq. 53, four for Eq. 56 | T (summation order; must agree with the per-component path to round-off, about 1e−12) | Fewer back-substitutions; simpler code |
| Precomputed reference geometry; deterministic colored threading | [Träff 2023](https://backend.orbit.dtu.dk/ws/files/338184592/Erik_Albert_Tr_ff.pdf); [Ferrite threading how-to](https://ferrite-fem.github.io/Ferrite.jl/stable/howto/threaded_assembly/) | — | T | Modest at 80×80 |
| Remove the extra Newton correction | [Thesis p. 36, Algorithm 5.2](https://gcris.iyte.edu.tr/handle/11147/13904) | 10 of 50 Newton factorizations (inference) | A | About 20% fewer Newton factorizations, but a less converged U |
| Warm start, adaptive bisection, modified or inexact Newton, reanalysis | [Senne et al. 2022](https://doi.org/10.1007/s11081-022-09769-6) (30–69% less objective time, "very similar" optima); [Amir & Sigmund 2011](https://doi.org/10.1007/s00158-010-0586-7); [Thillaithevan et al. 2024](https://doi.org/10.1007/s00158-024-03761-7); [Zhang et al. 2023](https://doi.org/10.1016/j.cma.2023.116119) | — | A | Large, but unquantified for homogenization; opt-in only |
| Iterative or multigrid solvers, matrix-free operators | [Amir et al. 2014](https://doi.org/10.1007/s00158-013-1015-5); [Aage et al. 2015](https://doi.org/10.1007/s00158-014-1157-0); nonlinear matrix-free application costs 1075.4 ms against 1.316 ms for the linear case ([Träff 2023](https://backend.orbit.dtu.dk/ws/files/338184592/Erik_Albert_Tr_ff.pdf)) | 3·10³–1.3·10⁴ DOFs | Not justified | — |

Once the class-B and class-T changes are in, about 50 Newton factorizations per evaluation remain and become the dominant cost. Ferrari & Sigmund observe the same: 85–90% of 2D runtime goes to the equilibrium solve once everything else is optimized ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)). Any further large gain therefore requires class-A changes to Newton itself, which is why those changes belong in decision records and not in performance commits.

Explicit Cholesky also gives a free diagnostic. A failed Cholesky of a converged tangent signals a loss of positive definiteness, which has physical meaning as a microstructural instability (inference). Logging the fallback count turns a solver detail into a stability indicator.

If the next paper reports cost, it should follow Hoefler and Belli's rules ([Hoefler & Belli 2015](https://doi.org/10.1145/2807591.2807644); [Chen & Revels 2016](https://arxiv.org/abs/1608.04295)):
- separate per-evaluation speedup from time to solution;
- report absolute times for both codes on the same machine;
- use BenchmarkTools minima with distributions.

## A typed package with two profiles, not a rewrite

The target is one Julia package that replaces the two forked solvers (856 and 1,000 lines) and their two nearly identical drivers (249 and 255 lines). Dispatch happens only at the true variation points:
- **`CellSymmetry`**, with `FullCell` and `HalfCell` methods, supplies the constraint maps D1/D2, the map from design vector to elements, and the homogenization domain weights. `HalfCell` keeps the paired scale factors (2 outside the sensitivity sum, 2 inside the volume gradient) and the half-period design offset ([NEMOpt-DR.jl L863, L934](../Domain_Reduction/NEMOpt-DR.jl#L863)).
- **Material types** are pluggable. The FEAD2026 energy (Bonet & Wood form) provides stress and tangent through Tensors.jl automatic differentiation, with a hand-coded Z fast path tested against AD. Wallin & Tortorelli's 2D energy, with B as a parameter, is used for the benchmark.
- **Interpolation types** carry their parameters explicitly: `ClipSIMP` (the default), `ModifiedSIMP` (opt-in), and a `PureSIMP` mode with a filter bypass for the benchmark. Projection, void-switch and filter types do the same.
- **Objective types** (`PoissonTarget`, `ShearEq53`, `BulkEq56`) replace the commented-out lines, implementing Wilson et al.'s rule not to control behavior by commenting code out ([Wilson et al. 2017](https://doi.org/10.1371/journal.pcbi.1005510)).
- **`Continuation` and `NewtonControl` structs** carry their counting and stopping semantics.
- **A `Profile` struct** bundles every class-A switch: `FEAD2026()` reproduces the legacy semantics, and `Default()` holds the proposed new defaults.

All state is passed through immutable parameter structs read from TOML. That removes the global-variable pathology which the Julia manual and the type-stability literature identify as the main obstacle to compiler optimization ([Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/); [Pelenitsyn et al. 2021](https://doi.org/10.1145/3485527)). Dispatch on concrete singleton types is resolved at compile time, so this generality costs nothing at runtime, as Gridap's peer-reviewed design demonstrates for finite-element assembly ([Verdugo & Badia 2022](https://doi.org/10.1016/j.cpc.2022.108341)).

```text
NEMOpt/
  Project.toml, Manifest.toml        # committed; reproduce/FEAD2026/ carries the reference Manifest
  src/  mesh  symmetry  material  interpolation  projection  filter
        assembly  newton  homogenization  sensitivity  objectives
        continuation  optimizer  profiles  provenance
  ext/  NEMOptPlotsExt.jl             # plotting and video as a package extension
  config/  fead2026_s5_2_case1.toml ... as_shipped_fd.toml  as_shipped_dr.toml  wt2020_s4_1.toml
  test/  kernels  evaluation  sensitivities  trajectories  benchmarks/wt2020  golden/ (hash-stamped JLD2)
  docs/decisions/  NDR-0001.md ...    # decision records
  reproduce/FEAD2026/                 # presets, environment, acceptance scripts
```

Two layer choices deserve explicit justification.

**Finite elements: keep the in-house structured-grid assembly,** with the master–slave map applied at scatter time in reduced numbering, rather than adopt Ferrite.jl now. Four facts support this:
- Released Ferrite (v1.7.0 and earlier) refuses affine constraints in which a degree of freedom is both master and slave. The fix was merged on 2026-09-28 but was unreleased on 2026-10-12.
- Ferrite v1.6.0 fixed a condensation-transpose bug that "gave silently wrong results" on non-symmetric matrices ([Ferrite CHANGELOG](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md)).
- Any library adoption changes numbering and summation order, so it can never be bitwise.
- The alternatives fit worse: Gridap is still 0.x and evolving rapidly, GridapTopOpt is level-set based, and TopOpt.jl supports linear elasticity only ([Wegert et al. 2025](https://doi.org/10.1007/s00158-024-03927-3); [TopOpt.jl](https://github.com/JuliaTopOpt/TopOpt.jl)).

**Optimization: pin NonconvexMMA in the FEAD2026 profile and preserve its call order.** The package was revived as v1.1.2 in August 2026 after about 27 months without commits. It depends on Zygote and Optim and has no per-iteration callback; its documented hook for hand-written gradients is `CustomGradFunction` ([Nonconvex docs](https://julianonconvex.github.io/Nonconvex.jl/stable/gradients/user_defined/); [NonconvexMMA](https://github.com/JuliaNonconvex/NonconvexMMA.jl)). For new studies, a plain `(x, grad)` interface admits two options:
- a clean-room MMA/GCMMA implementation written from [Svanberg 1987](https://doi.org/10.1002/nme.1620240207) and [Svanberg 2002](https://doi.org/10.1137/S1052623499362822), which avoids the GPLv3 terms of Svanberg's MATLAB code ([Alexandersen 2022](https://arxiv.org/abs/2207.13695));
- NLopt's CCSAQ.

Either is a different algorithm and therefore needs its own record.

| Layer | FEAD2026 reference | Refactored target (status on 2026-10-12) |
|---|---|---|
| Language | Pinned capture version (the baseline used 1.11.7; the original is unknown) | Julia 1.13.x (stable 1.13.1, released 2026-09-25); LTS 1.10.12 only if needed; the upgrade is its own T-class record ([Julia downloads](https://julialang.org/downloads/manual-downloads/)) |
| Constitutive | Legacy kernel (Tullio, closures) | Tensors.jl v1.17.1 (ForwardDiff duals, stack-allocated) ([Carlsson & Ekre 2019](https://doi.org/10.5334/jors.182)) |
| FE and constraints | Legacy assembly, Tᵀ K T | In-house reduced-numbering assembler; Ferrite deferred |
| Linear algebra | UMFPACK LU via `\` | CHOLMOD Cholesky, explicit, with symbolic reuse ([Chen et al. 2008](https://doi.org/10.1145/1391989.1391995)) |
| Optimizer | NonconvexMMA, pinned, MMA87/MMA02 | Same for FEAD2026; clean-room MMA/GCMMA or NLopt CCSAQ by record |
| Testing | Hash-stamped golden masters | Test + TestItemRunner v1.3.2; JET v0.12.3 (`@test_opt`); AllocCheck v0.2.7; Aqua v0.8.18; BenchmarkTools v1.8.0; StableRNGs ([JET](https://github.com/aviatesk/JET.jl); [AllocCheck](https://github.com/JuliaLang/AllocCheck.jl)) |
| Golden files | — | Thin in-house JLD2 helper with explicit rtol, atol and norm, and SHA-256 for bitwise checks, because ReferenceTests.jl has had little activity since its 2024-11-08 release ([ReferenceTests.jl](https://github.com/JuliaTesting/ReferenceTests.jl)) |
| Provenance | Timestamped folder, no commit recorded | DrWatson-style `tagsave` (git commit and patch), Manifest hash, `versioninfo()`, solver path, profile flags ([Datseris et al. 2020](https://doi.org/10.21105/joss.02673)) |
| Docs, format, latency, CI | — | Documenter v1.19.0; Runic v1.13.0; PrecompileTools v1.3.4; GitHub Actions `setup-julia@v3` |

One small hazard needs a rule of its own. The code writes 𝜌 (U+1D70C, mathematical italic small rho) on 63 lines of the full-domain solver and 67 lines of the half-domain solver. Julia normalizes identifiers only with NFC plus a few hard-coded look-alike mappings, so 𝜌 and the ρ that `\rho` produces are **different variables** ([Julia manual, allowed variable names](https://docs.julialang.org/en/v1/manual/variables/#man-allowed-variable-names)). The rename to ASCII names in public APIs must be done atomically, and a lint test must forbid characters from the Mathematical Alphanumeric Symbols block.

Threading, when added, uses only Ferrite-style deterministic coloring with fixed, recorded BLAS threads. Atomic and lock-based accumulation make the summation order depend on scheduling, which would break every regression band ([Ferrite threading how-to](https://ferrite-fem.github.io/Ferrite.jl/stable/howto/threaded_assembly/)).

## Six evidence-gated phases from freeze to paper

The phases below order the work so that every phase exits on evidence the next one depends on. One sequencing insight cuts the cost of validation substantially. The ensembles of round-off-perturbed runs needed for link E3 are expensive at legacy speed: the thesis reports 156.02 s per iteration at 80×80 even with domain reduction ([thesis p. 48, Fig. 6.4](https://gcris.iyte.edu.tr/handle/11147/13904)). **The Phase 3 bridge release is proven bitwise identical to the legacy code, so the expensive ensembles can run on it instead, without weakening the proof.** Phase 1 therefore captures only short golden masters, and the long acceptance runs move to the end of Phase 3.

| Phase | Goal | Exit gate (evidence) | Records |
|---|---|---|---|
| 0. Freeze and recover | Tag `30d0842` and archive it with a DOI. Pin the capture environment. From the authors' own records, recover the per-figure initial designs and seeds, the Table 2 densities, the filter radius and MMA variant per case, the script versions, the Julia and package versions, and LEMOpt. Document every FEAD2026 value no source states. | Legacy code runs in the pinned environment; an inventory lists every recovered and every missing input | NDR-0001 to NDR-0005, NDR-0044 |
| 1. Characterize | Capture full-precision golden masters GM-0 to GM-10 by wrapping the unmodified source; measure the legacy noise floor; replicate Tables 1–2; check full ↔ half domain equivalence; run short trajectories under both MMA variants | Golden masters hashed; noise bands per quantity; Table 2-type agreement at or below 1e−7 | NDR-0006 to NDR-0010 |
| 2. Restructure (Regime B) | Package skeleton, profiles, typed state, `CellSymmetry` dispatch, ASCII identifiers, I/O moved out of the objective; add the pluggable-energy and pure-SIMP test modes needed for the Wallin & Tortorelli golden test | Every golden master and every as-shipped trajectory bitwise identical, on every commit | NDR-0011 to NDR-0014, NDR-0026, NDR-0043 |
| 3. Bitwise-preserving speedups | Filter factored once; converged tangent reused; constitutive caching; triplet assembly; guards that never fire on valid inputs | Still bitwise; benchmark recorded; bridge release archived with a DOI; acceptance ensembles for the §5.2 and §5.3 presets run on the bridge release | NDR-0015 to NDR-0019 |
| 4. Modern numerics (Regime T) | Tensors.jl kernel, reduced-numbering one-triangle assembly, explicit Cholesky, aggregated adjoint, precomputed geometry, colored threads, Julia upgrade | L1–L3 within the measured bands; aggregated adjoint agrees with the per-component path to round-off; Wallin & Tortorelli §4.1 reproduced to printed digits; L4 inside the legacy ensemble; analytic and metamorphic oracles; Taylor slope 2; complex step | NDR-0020 to NDR-0025 |
| 5. Deliberate deviations (class A) | New defaults behind switches for the likely defects (unrounded objective, consistent gradient, outer-iteration continuation with exact terminal values, Newton control, optimizer options); opt-in options for the reference-formulation choices (modified SIMP, void-switch variants) | Per record: article presets re-run; change against the published values reported; a deviation-register row written | NDR-0027 to NDR-0036, NDR-0038 |
| 6. Release and paper support | Verification report with a DOI, `CITATION.cff`, deviation register, provenance audit | Every figure and table value of the new paper traces to a tagged result file | NDR-0037, NDR-0039 to NDR-0042 |

Phases 2 and 3 must never mix refactoring with behavior change. A commit that fails the hash test is either reverted or reclassified to Regime T, with a new record. After Phase 4 is accepted, the legacy arithmetic survives only in the archived bridge release and its CI job. The FEAD2026 semantic switches stay in the live code permanently, because they define the formulation that the next paper deviates from.

## Forty-four proposed decisions map onto the next paper

The decision trail uses lightweight, immutable, sequentially numbered records in the form Nygard introduced ([Nygard 2011](https://www.cognitect.com/blog/2011/11/15/documenting-architecture-decisions)). They add the "Considered Options" and "Confirmation" sections of MADR 4.0.0 ([MADR](https://adr.github.io/madr/)), plus numerics and paper-traceability fields.

Two empirical findings shape the process. In open-source projects, about half of the repositories that use decision records hold only 1–5 of them ([Buchgeher et al. 2023](https://research.jku.at/en/publications/using-architecture-decision-records-in-open-source-projects-an-ms/)). An ICSA 2026 study of about 5,800 records found that about 63% were created directly as "accepted", skipping the deliberation the format exists to capture ([ICSA 2026](https://conf.researchr.org/details/icsa-2026/icsa-2026-papers/34/Architecture-Decision-Records-Adoption-Impact-and-Developer-Engagement-in-Open-Sou)). Each NEMOpt record therefore opens as **proposed**, with its acceptance criteria written before any evidence is run, and moves to accepted only with that evidence attached.

The records map directly onto what the next topology-optimization paper needs:
- **A methods paragraph** that cites FEAD 2026 and the release DOI, with a parameter table per result. The table should include the FEAD2026 values no source documents: ε = 5e−4, β̂ = 200, `s_init` and `maxinner`.
- **A verification statement** that cites NASEM Finding 4-3 for the tolerance-based definition, and the Wallin & Tortorelli golden test as independent evidence.
- **A deviation register**, "Differences from Güven & Özdemir (2026)", with one row per class-A record. Where a change also departs from Wallin & Tortorelli, the row says so.

That set serves:
- the "Replication of results" section that *Structural and Multidisciplinary Optimization* requires ([Haftka et al. 2019](https://doi.org/10.1007/s00158-019-02298-4), p. 405). The editorial fixes no content but recommends source code with a tagged version and DOI, all tuning and convergence parameters, plot data and intermediate results, and repeated runs to show scatter (pp. 406–409);
- the version-specific software citation principle ([Smith et al. 2016](https://doi.org/10.7717/peerj-cs.86));
- the FAIR principles for research software ([Barker et al. 2022](https://doi.org/10.1038/s41597-022-01710-x));
- Elsevier's research-data options for FEAD or CMAME ([Elsevier data guidelines](https://www.elsevier.com/authors/author-resources/research-data/data-guidelines)). The specific option of each journal could not be retrieved and must be checked before submission.

```markdown
---
id: NDR-0029                       # sequential, never reused (Nygard 2011)
title: Trigger p/β continuation per MMA outer iteration, between iterations
status: proposed                   # proposed | accepted | rejected | superseded-by NDR-xxxx
date: 2026-10-12
deciders: [ ... ]
class: A                           # P process | B bitwise-preserving | T tolerance-level | A results-changing
supersedes: []
related: [NDR-0006, NDR-0010, NDR-0027]
---
## Context
Legacy behavior with file:line and commit (e.g., Full_Domain/NEMOpt_main.jl L45, L90–102 @30d0842).
Article, thesis and Wallin & Tortorelli statements (equation, table, page) and any inconsistency between them.
Origin: reference formulation (W&T 2020) | FEAD 2026 choice | likely defect.
Decision drivers: fidelity to FEAD 2026, correctness, cost, maintainability.
## Considered options
Option 1 (FEAD2026 legacy) / Option 2 / Option 3, each with pros and cons.
## Decision
Chosen option; scope (FullCell, HalfCell, both); profiles affected.
## Rationale and primary sources
Citations with DOIs; each supporting claim labeled verified | inference | unverified.
## Effect on numerical results
affects-results: no | bitwise-only | round-off (Regime T) | yes (class A)
FEAD2026 switch value: profile.continuation = (counter = :objective_evaluations, timing = :inside_objective)
proposed default:      (counter = :outer_iterations, timing = :between_iterations)
deviation-from-FEAD2026-formulation: yes -> deviation-register row text: "..."
deviation-from-W&T-2020-formulation: yes | no
## Pre-registered acceptance criteria (written before evidence is run)
Levels L0–L4 affected; quantities; tolerance bands and their derivation (noise-floor run IDs,
condition estimate); ensemble size, perturbation types and metrics for L4; article values to be
compared at printed precision.
## Verification evidence (filled at acceptance)
Test IDs (GM-x, K-x, E-x, S-x, R-x, WT-x); result-file IDs with git describe and Manifest hash;
outcome per criterion; FEAD 2026 table/figure: reproduced to N digits | changed (old -> new).
## Performance evidence (only if cost will be reported)
Hardware, OS, Julia/BLAS/SuiteSparse versions, threads; BenchmarkTools minimum and distribution;
absolute legacy and new times (Hoefler & Belli 2015, rules 1, 5, 9).
## Paper traceability
Future-paper section/figure/table; methods sentence; release tag and DOI.
## Consequences
Positive, negative, follow-up records, risks.
## Links
Commits/PRs, CHANGELOG entry, issues.
```

The initial log follows. The classes are P (process or documentation), B (bitwise-preserving), T (tolerance-level) and A (results-changing). Each class-A row becomes one row of the next paper's deviation register when its option is used. **Every decision is proposed, and none is accepted until its evidence exists.**

| ID | Proposed decision | Class | FEAD2026 value → proposed default or option | Acceptance evidence | Status |
|---|---|---|---|---|---|
| NDR-0001 | Adopt numbered, immutable decision records using the template above; open each as proposed with pre-registered criteria; one umbrella record per batch of bitwise-verified refactoring commits | P | — | Records exist before their evidence runs | proposed |
| NDR-0002 | Tag `30d0842` as `legacy-30d0842` and archive it with a DOI; call it "the FEAD 2026 code" only if Phase 0 confirms that | P | — | Archived tag and DOI | proposed |
| NDR-0003 | Recover from the authors' records the per-figure initial designs and seeds, Table 2 densities, filter radius and MMA variant per case, script versions, environment and LEMOpt; record each unrecoverable item as a gap | P | — | Inventory table | proposed |
| NDR-0004 | Commit Project/Manifest for the capture environment (Julia 1.11.7 unless the original is recovered); BLAS threads = 1; log the SuiteSparse version and solver path | P | — | GM-0 environment fixture | proposed |
| NDR-0005 | Ground truth: code, then thesis, then article for formulas, with Wallin & Tortorelli where all three are silent; the article's tables and text for case parameters. Publish the errata list: Eq. 20 (misprint with an O(1) error; the code is correct), Eq. 22, Eq. 28, Eqs. 51–52 and Eq. 33. Note that a corrigendum is the authors' decision. | P | — | Errata table reviewed | proposed |
| NDR-0006 | Define "same results" as chain E1–E4 with Regime B and Regime T; bands come from measured noise floors, never from default `isapprox` | P | — | Noise-floor runs | proposed |
| NDR-0007 | Capture characterization data at full precision by wrapping the unmodified legacy source (GM-0 to GM-10) | P | — | Hash-stamped JLD2 fixtures | proposed |
| NDR-0008 | Keep presets separate: as-shipped full and half domain; article §5.1, §5.2 and §5.3 with r = 1 and r = 2 elements; thesis Chapter 7. Never mix them. | P | — | Preset files and tests | proposed |
| NDR-0009 | Article Table 1/2 protocol: documented ρ₁…ρ₄; r as elements and as absolute length; forward and central differences; all 16 components | P | — | \|adjoint − FD\| ≤ 1e−7 (published maximum 5.59858e−8) | proposed |
| NDR-0010 | Full-run acceptance: a legacy ensemble (size, perturbations, and metrics J, ν, κ, Vf, M_nd, thresholded ν, IoU), following SMO's advice to make repeated, slightly perturbed runs to show scatter (Haftka et al. 2019, p. 407), plus the article values at printed precision | P | — | Pre-registered bands | proposed |
| NDR-0011 | Package skeleton; all state in typed immutable structs read from TOML; no globals, no `cd` | B | — | Golden masters and as-shipped trajectories bitwise | proposed |
| NDR-0012 | Merge full and half domain via `CellSymmetry` dispatch, preserving parity rules, the paired half-domain scale factors and the half-period design offset | B | — | Bitwise per variant; equivalence modulo the shift under a shear F_M | proposed |
| NDR-0013 | ASCII public identifiers; lint test against mathematical-alphanumeric and Greek confusables | B | — | Atomic rename, bitwise | proposed |
| NDR-0014 | Move plotting and file output from the objective into a logger; write full precision; round only for display | B | — | Bitwise | proposed |
| NDR-0015 | Assemble and factor the filter once per run; share it across objective, constraint and gradients | B | — | Bitwise | proposed |
| NDR-0016 | Reuse the final Newton assembly and one factorization for the four χ solves and all adjoint solves | B | — | Hash equality, otherwise reclassify as T | proposed |
| NDR-0017 | Cache converged-state Gauss-point constitutive results for homogenization and sensitivities | B | — | Bitwise | proposed |
| NDR-0018 | Guards that never fire on terminating cases: error when Rep_MU = 0, a Newton iteration cap that raises, detection of ρ̃ < 0 before fractional powers, a square-mesh assertion | B | hang → explicit error | Bitwise on all golden masters; GM-10 now raises | proposed |
| NDR-0019 | Triplet assembly instead of incremental CSC insertion | B or T | — | Hash test, otherwise a Regime-T band | proposed |
| NDR-0020 | Tensors.jl constitutive kernel (AD for P and 𝔸; hand-coded Z tested against AD and central differences); zero allocations | T | legacy kernel kept in the bridge release | L1 within the measured band; L2–L4 bands | proposed |
| NDR-0021 | Preallocated pattern in reduced master–slave numbering; one-triangle assembly; `Symmetric` wrapper | T | — | L2 bands; exact symmetry test | proposed |
| NDR-0022 | Explicit CHOLMOD Cholesky with symbolic reuse; logged LDLᵀ/LU fallback; never rely on `\` dispatch | T | UMFPACK LU → Cholesky | L2 bands; fallback counts logged | proposed |
| NDR-0023 | Aggregated adjoint: one state adjoint and one filter back-substitution per scalar function, with weights ∂J/∂𝔸. This is a proposed algorithmic improvement: neither Wallin & Tortorelli (Eq. 50) nor the code aggregates, and the step is exact by linearity. | T | per component; weights use the rounded ν under the FEAD2026 value of NDR-0028 | Agreement with the per-component path to round-off (about 1e−12); Taylor slope 2 | proposed |
| NDR-0024 | Precomputed reference geometry; deterministic colored threading only; BLAS threads recorded | T | serial | Results invariant to thread count | proposed |
| NDR-0025 | Upgrade from the reference Julia version to 1.13.x as a separate step | T | reference version → 1.13.x | Noise-floor comparison | proposed |
| NDR-0026 | Replace rrules and globals with explicit `(x, grad)` callbacks while preserving NonconvexMMA's call order; pin NonconvexMMA for FEAD2026 | B | — | Bitwise trajectories | proposed |
| NDR-0027 | Objective value returned to the optimizer (likely defect) | A | rounded to 3 decimals → unrounded | Article presets re-run; deviation row | proposed |
| NDR-0028 | ν used in ∇J (likely defect) | A | rounded ν → exact ν | Taylor test passes for the default | proposed |
| NDR-0029 | Continuation counter, timing and terminal values (likely defects: per-evaluation counting, p/β changed before the gradient, full-domain β overshoot to 33). Schedules are stored as `{istart, max, every, delta}`; the overshoot is kept only in the as-shipped preset. The article, thesis and Wallin & Tortorelli all count design iterations. | A | per evaluation, inside the objective → per outer iteration, between iterations, with exact terminal values | Article presets; ensemble comparison | proposed |
| NDR-0030 | Unused β-dependent `s_init` local (likely defect: the intended update never reaches the optimizer) | A if implemented | no effect → deleted | Bitwise for deletion | proposed |
| NDR-0031 | Newton stopping control (no source other than the code specifies tolerance or scaling) | A | check before update, extra correction, `Rep_MU = \|max(F_M − I)\|·Lx`, TOL = 1e−7 → check after update (thesis Algorithm 5.2), a scale from a norm of F_M − I that is never zero, and a tolerance chosen so that Table 2-type agreement does not degrade | Article Table 1/2 protocol; Wallin & Tortorelli golden test; L2–L4 | proposed |
| NDR-0032 | Load stepping and Newton variants | A | 10 fixed steps from zero with full Newton → unchanged default; bisection only as a failure fallback; warm start and inexact Newton opt-in | Branch-consistency study | proposed |
| NDR-0033 | Stiffness interpolation. The clip `max(H^p, ε)` is Wallin & Tortorelli's reference law (Eq. 71, ε = 1e−5 there). It equals classic SIMP evaluated at max(H, ε^{1/p}); classic SIMP, however, bounds the design variable itself, 0 < ρ_min ≤ ρ, "to avoid a singular FEM problem" (Bendsøe & Sigmund 1999, p. 638), whereas the clip lets H and the volume fall below ε^{1/p} at constant stiffness. Modified SIMP ε + (1 − ε)H^p has the form of Bendsøe & Sigmund's two-material SIMP (Eq. 23), which "will never satisfy the Hashin-Shtrikman bounds for all densities" (p. 647). At p = 3 under plane strain neither law is bound-compatible (both exceed the bulk bound by about 1.6% near H = 0.93), so bound compatibility does not decide between them. | A (option only) | `max(H^p, ε)` with ε = 5e−4 → **unchanged default**; modified SIMP ε + (1 − ε)H^p as an opt-in deviation behind `stiffness_interp`, recorded as departing from both sources | Article presets re-run when the option is used; deviation row | proposed |
| NDR-0034 | Void (linear/nonlinear) switch, which matches the thesis and Wallin & Tortorelli Eq. 73 in argument and threshold. Its soft-element band is inherited from that scheme. | A (option only) | Ĥ(ρ̃^p) with (200, 0.05) → unchanged default; exploratory options β̂ = 500 and Ĥ(H(ρ̃)^p) with (500, 0.01) (precedents: Kumar et al. 2021; Dalklint et al. 2021, Eq. 9); distortion diagnostics | Controlled study before any default change | proposed |
| NDR-0035 | Optimizer variant and options per preset. Full-domain MMA87 matches Wallin & Tortorelli's plain MMA; half-domain MMA02 follows the thesis. | A | MMA87 / MMA02, `s_init` 0.25, `maxinner` 40, kkt 1e−8 → explicit; identical optimizer in any full-versus-half comparison; new MMA implementations only by new record | Trajectory comparison | proposed |
| NDR-0036 | Initial-design ingestion and RNG | B/A | legacy crop, resize and rescale → hashed design-vector fixtures (B); corrected loader optional (A); StableRNGs for random starts | Fixture hashes | proposed |
| NDR-0037 | Effective-property definitions: keep ν (Eq. 54, a deliberate deviation from Wallin & Tortorelli's uniaxial-stress ν; for a square-symmetric linear tensor it equals the 2D uniaxial-stress ratio C₁₁₂₂/C₁₁₁₁; it does not appear in Kazemi et al. 2020, whose Eq. 24 (p. 8) is ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)] for 3D cubic lattices, so the documentation must not cite Kazemi et al. as its source; the attribution in the article and in thesis Eq. 7.5 is incorrect, and whether to issue an erratum is the authors' decision), κ and the shear sum. Rename the logged G to G₁₂ = 𝔸₁₂₁₂ (Wallin & Tortorelli's μ₁₂) and G_Σ = 𝔸₁₂₁₂ + 𝔸₂₁₂₁ (Eq. 53). Document them as tangent measures with ν in (−1, 1). Add an optional uniaxial-stress ν₁₂ diagnostic; never compare it numerically with Eq. 54. | P | — | Documentation and tests | proposed |
| NDR-0038 | Re-enable Eqs. 53 and 56 as objective types, with δ𝔸₂₂₁₁ = δ𝔸₁₁₂₂ by major symmetry | A (new capability) | commented out → selectable | Numerical symmetry check; gradient check | proposed |
| NDR-0039 | Filter length as a physical length with an explicit convention (Helmholtz r versus cone R ≈ 2√3·r, a second-moment match on an unbounded domain: [Lazarov, Wang & Sigmund 2016, Eq. 10](https://backend.orbit.dtu.dk/ws/files/140201797/main.pdf), attributing it to Lazarov & Sigmund 2011); Wallin & Tortorelli's r = 2 elements equals the article's text value | P | r·Lx/Nelx with r = 1 element | Preset tests | proposed |
| NDR-0040 | Box constraint [0, 1] as coded; document the article's ρ_min | P | [0, 1] | Documentation | proposed |
| NDR-0041 | Provenance on every result file (git describe and patch, Manifest hash, `versioninfo()`, BLAS, solver path, profile flags, parameters, design hash); `CITATION.cff`; Zenodo DOIs for the legacy, bridge and refactored releases | P | — | Audit of the figure-to-file map | proposed |
| NDR-0042 | Keep in-house structured-grid assembly; defer Ferrite until a release untangles master–slave chains and a T-class record justifies the move | P | — | — | proposed |
| NDR-0043 | Add Wallin & Tortorelli §4.1 as an independent golden test. Requires a pluggable energy (their 2D form with the printed B = μ/8 typo kept, since the published numbers depend on it), a pure-SIMP mode S = ν³P with filter and projection bypassed, element-group densities assigned by centroid, and the index map (1, 2, 3, 4) → (1, 4, 3, 2). Port the evidence scripts as reference implementations. | B (new test modes; FEAD2026 unaffected) | — | Case A and case B S_M and 𝔸_M (including 𝔸₁₂₁₂ = −0.9861, F₂₂ = 0.872199, ν₁₂ = 0.356115) to printed digits; Table 1 sensitivities by symmetric adjoint and central differences (1e−5) to about 1e−6 relative; degree-3 Euler identity | proposed |
| NDR-0044 | Record in the profile and in the next paper's methods the FEAD2026 values that no source documents or justifies: ε = 5e−4 (against 1e−5), β̂ = 200 (against 500), `s_init` 0.25, `maxinner` 40, TOL 1e−7 scaled by Rep_MU, and 10 load steps. Optionally add a "Wallin & Tortorelli formulation" preset (β̂ = 500, ε = 1e−5, r = 2 elements, doubling β schedule) for formulation-sensitivity studies only, without claiming to reproduce their numbers, since energy and kinematics differ (SMO: "detailed settings for results need to be documented", Haftka et al. 2019, p. 406). | P (preset: A) | — | Documentation; preset tests | proposed |

## Conclusion

The most consequential finding is that the hardest requirement, provable sameness with a published article, is limited less by software than by missing records. Equivalence to the legacy code can be proven bitwise. Equivalence to the article can only be supported by evidence, because the inputs behind the published numbers were never committed: the initial designs, the Table 2 densities and the environment. **The cheapest way to strengthen that proof is the authors' own archives, not more code.** Each input recovered in Phase 0 turns a statistical acceptance band into an exact replay. Wallin & Tortorelli's fully specified verification case fills the other gap: it gives the refactor a correctness oracle that does not depend on any unpublished input, and it shows that the code is right where the article's printed Eq. 20 is not.

Two further points matter for planning. First, NEMOpt is more faithful to its source formulation than the article's wording suggests. The clipped SIMP law, the ρ̃^p blend with its soft-element band, the symmetric sensitivity and the fixed center node are all inherited reference formulation, so changing them is a modeling choice for future studies, not a repair. What is genuinely NEMOpt-specific is ε = 5e−4, β̂ = 200, the continuation schedules, the MMA options and the Newton control, plus five likely defects: rounding, p/β updates before the gradient, per-evaluation counting, the unused `options` local and the β overshoot. That is exactly the set of switches the FEAD2026 profile must name. Second, reproduction and efficiency conflict less than they appear to. A substantial tranche of the speedup preserves every bit, and the resulting bridge release makes the expensive ensemble validation affordable without weakening any proof. Once named and recorded, the legacy behaviors become a citable methods baseline that the next paper can deviate from deliberately, with one decision record and one re-run number per change.
## Addendum: buckling and elastic stability

A separate addendum extends this report to buckling and elastic stability: [NEMOpt buckling and stability addendum](NEMOpt%20buckling%20and%20stability%20addendum.md). It covers:
- the stability checks NEMOpt needs: one-cell k = 0 inertia, macroscopic ellipticity of 𝔸_M, and Bloch–Floquet scans;
- when one-cell and half-cell results stop being valid;
- the risk to the published load cases;
- a read-only diagnostics architecture that leaves the FEAD2026 profile bitwise unchanged;
- verification benchmarks;
- a constrain-first formulation for a future stability-aware paper.

The decision log continues there as **NDR-0045 to NDR-0063**, all proposed. NDR-0064 and NDR-0065 are recorded in the next section, and the roadmap records NDR-0066 to NDR-0082 follow it. Two items in the addendum change records here:
- **The panel-table sentence** "Full-versus-half agreement ends at the first symmetry-breaking bifurcation" is amended. The *guarantee* of agreement ends there, but agreement itself may continue on an unstable symmetric branch. Agreement between the full and half cell is therefore never evidence of stability (addendum NDR-0051).
- **NDR-0022's "logged LDLᵀ/LU fallback"** is refined. UMFPACK LU serves the solve after a Cholesky failure, and LDLᵀ runs only as a logged diagnostic, because it is unpivoted and 10–21× slower than Cholesky (addendum NDR-0047).

## Decisions recorded after the addendum

### NDR-0064: elastic bounds are checks only where they hold

**Hashin–Shtrikman bounds are theorems for linear elastic moduli** ([Hashin & Shtrikman 1963, JMPS 11:127](https://garfield.library.upenn.edu/classics1988/A1988Q027300001.pdf)). They therefore cannot validate the finite-strain tangent values (𝔸_M, ν, κ at F_M ≠ I) that NEMOpt reports.

**Nonlinear extensions exist but do not supply a tangent check:**
- Talbot & Willis (1985) extended the Hashin–Shtrikman variational principles to nonlinear media, in a geometrically linear setting ([Milton, *Theory of Composites*](https://resolve.cambridge.org/core/books/theory-of-composites/variational-principles-and-inequalities/F55B22AA23653A981F56C1D0C9F22D78), p. 271 and §13.7, pp. 282–286 of the [author-hosted copy](https://gwmilton.com/TheoryCompositesNOPRINT.pdf); verified).
- Ponte Castañeda (1991, JMPS 39:45) gave Hashin–Shtrikman-type bounds through a linear comparison composite, and Ponte Castañeda (1989, Proc. R. Soc. A 422:147) extended the approach to finite elasticity ([summary](https://apps.dtic.mil/sti/pdfs/ADA231637.pdf); verified only through secondary sources: the 1991 comparison bounds in Milton 2002, p. 271 and pp. 285–286; the 1989 bounds on the effective energy as cited by [Geymonat et al. 1993, p. 237](https://doi.org/10.1007/BF00380256) and, with full bibliographic data, by [Michel et al. 2007](https://doi.org/10.1016/j.jmps.2006.11.006)).
- These bound the effective *energy*, not tangent moduli.
- They are posed for geometrically linear problems and rely on growth conditions, and on convexity for the dual principles (Milton 2002, §13.7, pp. 283–286); finite-strain hyperelastic energies are not convex in F. Without convexity the homogenized energy is an infimum over fluctuations periodic on k×k supercells for all k, which can lie strictly below the one-cell value ([Müller 1987, ARMA 99:189–212](https://doi.org/10.1007/BF00284506), as restated in [Geymonat et al. 1993, Thm 2.3 and Eqs. 2.6–2.7](https://doi.org/10.1007/BF00380256), in [Michel et al. 2007, Eq. 2.22](https://doi.org/10.1016/j.jmps.2006.11.006) and in Milton 2002, §1.7, p. 14, which describes computing it "from solutions in increasingly large blocks of unit cells"; verified as restated there, Müller's own text not read). The one-cell value is exact for convex energies, and from a state where it is exact it stays exact while the cell is stable against Bloch waves of every wave vector, under a continuity hypothesis on the cell minimizers ([Geymonat et al. 1993, Prop. 5.3](https://doi.org/10.1007/BF00380256); verified).
- The second-order finite-strain results for porous rubbers use a Hashin–Shtrikman-type *estimate* for the linear comparison composite, "not to be confused with the Hashin–Shtrikman bounds", so they are estimates, not bounds ([Michel et al. 2007, §3.2.1](https://doi.org/10.1016/j.jmps.2006.11.006); verified from the HAL submitted version, hal-00111501). Derived for random microstructures ([Lopez-Pamies & Ponte Castañeda 2004, J. Elasticity 76:247–287](https://dialnet.unirioja.es/ejemplar/108150), bibliographic data verified; [2007, JMPS 55:1677–1701](https://doi.org/10.1016/j.jmps.2007.01.007), read in full), they have also been specialized to square and hexagonal arrays of circular holes. There they estimate only macroscopic loss of ellipticity, and against finite elements they ranged from "remarkably good" (square, porosity ≈ 0.2) to "an unsatisfactory overestimation" (hexagonal, porosity 0.5). Random-microstructure estimates lose ellipticity later than periodic arrays of equal porosity ([Michel et al. 2007, §§4.3, 5](https://doi.org/10.1016/j.jmps.2006.11.006)). They are therefore not conservative for periodic cells, and no such estimate exists for NEMOpt's optimized topologies (inference).

**Computed against the FEAD2026 constants** (plane strain, E = 100, ν = 0.3, so κ_s = 96.15 and μ_s = 38.46; [script](../research_notes/NEMOpt%20elastic%20bounds%20check/hs_simp_check.py)):
- **2D Hashin–Shtrikman upper bulk bound for voids:** 15.38, 21.37 and 28.85 at solid fraction 0.4, 0.5 and 0.6. Modeling the void as the ε = 5e−4 material raises these to 15.44, 21.43 and 28.905. The void values equal Bendsøe & Sigmund's void-limit formula (1999, Eq. 6), and all six equal Milton's optimal 2D energy bound at hydrostatic strain (2002, Eq. 23.30), which holds for any anisotropic two-phase cell ([bound-check scripts](../research_notes/NEMOpt%20supplied%20papers%20verification/evidence/milton_bounds/); verified against both texts).
- **2D shear bounds for voids:** the isotropic Hashin–Shtrikman shear bound is 7.40, 10.12 and 13.42 at solid fraction 0.4, 0.5 and 0.6. These values come from the corrected scripts (commit 5aa11db; 2D form ζ = κμ/(κ + 2μ)). The earlier shear output, 7.96, 10.82 and 14.23, was 6–9% too high; no report quoted it, and the bulk values were unaffected. The isotropic shear bound needs an isotropic cell. For a single shear modulus of an anisotropic cell, Milton's Eq. 23.30 at pure-shear strain gives 12.41, 16.03 and 19.89 (19.905 at 0.6 with the ε-material) (computed).
- **p = 3 SIMP exceeds the bulk bound** by up to 1.6% (1.0163 times the bound at projected density H = 14/15 ≈ 0.93). The excess vanishes at p = 3.5. This is Bendsøe & Sigmund's 2D condition p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} (Eq. 11, p. 639), whose bulk term gives 3.5 here and whose shear term gives 2.8 ([Bendsøe & Sigmund 1999](https://link.springer.com/doi/10.1007/s004190050248); verified in the full text). The paper derives it for a generic 2D isotropic material from the 2D Hashin–Shtrikman bounds (Eq. 6) and κ⁰ = E⁰/(2(1 − ν⁰)), μ⁰ = E⁰/(2(1 + ν⁰)) (Eq. 9). So ν⁰ is the in-plane Poisson's ratio, which under plane strain is ν/(1 − ν) = 0.4286, giving p\* = 3.5; the 3D ν = 0.3 would give the plane-stress value p\* = 3.08, set by the shear term. In 2D, p\* = 3 "holds only for ν⁰ = 1/3". Under plane strain, p = 3 SIMP stays below the shear bound (ratio 0.99998).
- **Consequence:** a design that still contains gray elements can exceed the bound without any bug. Bendsøe & Sigmund (1999, pp. 636, 651) likewise treat bound compatibility of gray as optional when designs end black-and-white.

**The one rigorous and simple finite-strain bound is Voigt-type.** Periodic homogenization defines the effective energy as an infimum over fluctuations periodic on k×k supercells for all k, and NEMOpt computes the one-cell (k = 1) value, which bounds it from above ([Geymonat et al. 1993, Eqs. 2.6–2.7](https://doi.org/10.1007/BF00380256)). Zero fluctuation is admissible for every k, so W̄(F̄) ≤ Ŵ¹(F̄) ≤ (1/|Ω|)∫ W_ρ(x, F̄) dV (standard derivation). The computed one-cell energy obeys the bound only if the computed state is the one-cell energy minimizer, so a violation flags either a bug or a non-minimizing, possibly unstable, branch. In linear elasticity, exceeding the Voigt bound is likewise tied to loss of stability ([Kochmann & Milton 2014](https://arxiv.org/pdf/1401.4142)).

| ID | Proposed decision | Class | FEAD2026 value → proposed default | Acceptance evidence | Status |
|---|---|---|---|---|---|
| NDR-0064 | Use elastic bounds only where they hold. **(a) Linear limit:** apply the 2D Hashin–Shtrikman upper bound on the hydrostatic-strain bulk modulus κ = ¼(𝔸₁₁₁₁ + 𝔸₂₂₂₂ + 𝔸₁₁₂₂ + 𝔸₂₂₁₁) only to 𝔸_M(F_M = I). Apply it to thresholded two-phase designs, with the void modeled as the ε-material; there it holds for any cell geometry and any effective anisotropy (Milton 2002, Eq. 23.30) and can serve as a hard gate. Never apply it to raw gray designs or to finite-strain tangents. Gray elements can exceed the bound by up to 1.6%, because p = 3 is below the plane-strain threshold p\* = 3.5 of Bendsøe & Sigmund 1999, Eq. 11. Do not apply the isotropic shear bound (13.42 at solid fraction 0.6 with voids) to anisotropic cells, because it needs isotropy. For a single shear modulus (G₁₂ = 𝔸₁₂₁₂ at F_M = I) use Milton's Eq. 23.30 at pure-shear strain: 12.41, 16.03 and 19.89 at solid fraction 0.4, 0.5 and 0.6 with voids, and 19.905 at 0.6 with the ε-material (computed). Bendsøe & Sigmund's (1999) strain-energy bound (Eq. 21) gives the same limit, 𝔸₁₂₁₂ + 𝔸₂₁₂₁ ≤ 39.8 at 0.6 (derived). For square-symmetric cells, the 2D trace bound on the pair of shear moduli (derived from Milton's Eqs. 23.33 and 23.44) is optional. **(b) Finite strain:** add a homogenized-energy output (new output; FEAD2026 arithmetic untouched) and check the Voigt-type bound W̄(F̄) ≤ (1/\|Ω\|)∫W_ρ(x, F̄) dV. Treat a violation as a bug or a non-minimizing branch, cross-referenced with the stability records NDR-0045 to NDR-0050. **(c)** Cite nonlinear Hashin–Shtrikman-type bounds and second-order estimates for comparison only, never as acceptance tests. The second-order results are estimates, not bounds; they predict only macroscopic loss of ellipticity, which a periodic cell's microscopic bifurcation precedes; and for periodic arrays they overestimated finite-element onset strains by up to about a factor of two (Michel et al. 2007, hexagonal array, porosity 0.5; digitized), so they never certify stability. **(d)** This record amends the L4 oracle list and qualifies the §5.3 κ plausibility statement. | P (checks); B (energy output) | No bound checks; κ compared with the isotropic bound in prose → (a) and (b) as automated tests | Bound formulas reproduce the computed values (bulk 15.38 / 21.37 / 28.85, equal to Bendsøe & Sigmund 1999 Eq. 6 and to Milton 2002 Eq. 23.30 at hydrostatic strain; 28.905 with ε; isotropic shear 7.40 / 10.12 / 13.42 and single-modulus shear 12.41 / 16.03 / 19.89 at solid fraction 0.4 / 0.5 / 0.6 with voids, after the scripts' 2D shear formula was corrected to ζ = κμ/(κ + 2μ) in commit 5aa11db); F_M = I evaluations of thresholded regenerated designs stay below the bound within the measured noise band; the Voigt-type energy check passes on all golden-master states; the bulk bound's validity for the hydrostatic-strain modulus of anisotropic cells is confirmed (Milton 2002, §23.4, Eq. 23.30, p. 466 of the author's posted copy: 2D, two isotropic phases with κ₁ > κ₂ > 0 and μ₁ > μ₂ > 0, the void as a limit; Bendsøe & Sigmund's 1999 Eq. 21 gives the same result for void/solid cells, derived), so (a) can be a hard gate for thresholded two-phase designs at F_M = I | proposed |

### NDR-0065: optional 2D kinematics modes

FEAD2026 and the 2023 thesis use **plane strain**. Wallin & Tortorelli did not: their §5 cell is one layer of 3D bricks with the out-of-plane macro stretch left free, a **generalized plane strain** state with zero mean out-of-plane stress ([comparison note, Q8](../research_notes/NEMOpt%20refactoring%20strategy%20literature/wallin_tortorelli_2020.md)). Thin-sheet metamaterials call for **plane stress**.

**At finite strain, plane stress is not a change of constants.** Each Gauss point needs a local Newton solve for the out-of-plane stretch so that S₃₃ = 0, and the tangent must be condensed accordingly. This is the approach of Klinkel & Govindjee (2002), cited from a later paper's description ([summary](https://arxiv.org/pdf/2007.05877)). A general formulation for compressible hyperelastic laws with open-source code is given by Ahmadi et al. ([arXiv:2410.22562](https://arxiv.org/pdf/2410.22562)).

**Plane stress cannot represent out-of-plane deformation.** Thin perforated sheets often deform out of plane, which makes their Poisson's ratio much less negative than in-plane theory predicts ([Mizzi et al. 2020](https://eng.ox.ac.uk/media/9677/mizzi-ijms-2020.pdf)). Plane-stress stability verdicts are therefore in-plane only and optimistic for thin sheets.

The choice also changes the linear constants, and with them the NDR-0064 bounds ([script](../research_notes/NEMOpt%20elastic%20bounds%20check/plane_stress_comparison.py); computed for E = 100, ν = 0.3; shear ratios from the corrected 2D shear formula, commit 5aa11db):

| | Plane strain (FEAD2026) | Plane stress |
|---|---|---|
| 2D bulk modulus of the solid | 96.15 | 71.43 |
| In-plane Poisson's ratio of the solid | 0.4286 | 0.30 |
| Hashin–Shtrikman bulk bound at 40% solid | 15.38 | 13.51 |
| p = 3 SIMP against the Hashin–Shtrikman bulk and shear bounds | exceeds the bulk bound by up to 1.6%; stays below the shear bound (ratio 0.99998) | stays at or just below the bulk bound (ratio 0.9999) but exceeds the shear bound by up to 0.05% (ratio 1.00047 at density 0.988) |
| SIMP threshold p\* (Bendsøe & Sigmund 1999, Eq. 11) | 3.5 (bulk-governed) | 3.08 (shear-governed) |

| ID | Proposed decision | Class | FEAD2026 value → proposed default | Acceptance evidence | Status |
|---|---|---|---|---|---|
| NDR-0065 | Add a `Kinematics2D` option with three modes. **PlaneStrain** is the default and reproduces FEAD2026. **PlaneStress** enforces S₃₃ = 0 at each Gauss point by a local Newton solve for the out-of-plane stretch, with a condensed tangent; Z is obtained by differentiating through the local solve and checked against automatic differentiation. **GeneralizedPlaneStrain** uses one uniform out-of-plane stretch per cell, with zero mean out-of-plane stress (one extra scalar unknown); it is the Wallin & Tortorelli §5 setting and enables formulation comparisons with that paper. Requirements: label plane-stress stability verdicts "in-plane modes only"; re-derive the closed-form ellipticity test (NDR-0049) for each mode; recompute the NDR-0064 bound values per mode (the plane-strain bound stays a valid, possibly loose, upper bound in GeneralizedPlaneStrain, because the cell's in-plane moduli obey plane stress ≤ generalized plane strain ≤ plane strain, whereas the plane-stress bound is not guaranteed to bound generalized-plane-strain cells; inference). Implement after Phase 4; the option never alters the FEAD2026 profile. | A (option) | Plane strain only → three modes, plane strain default | FEAD2026 golden masters bitwise unchanged; homogeneous uniaxial and biaxial tests against the analytical 3D solution with σ₃₃ = 0 (plane stress) or zero mean σ₃₃ (generalized plane strain); in the linear limit, 𝔸_M(F_M = I) equals linear homogenization with the corresponding 2D constants; Z against automatic differentiation and central differences; NDR-0064 bounds recomputed per mode | proposed |

## Roadmap toward a general, high-performance tool (NDR-0066 to NDR-0082)

The authors want NEMOpt to become a general topology-optimization tool with high performance. These records turn that goal into proposed decisions. Every record is **proposed**; none is accepted until its evidence exists, and none may change the FEAD2026 profile.

**Niche.** A general tool competes with large established codes. NEMOpt can lead in **finite-strain (nonlinear) topology optimization of microstructures and structures, with verified sensitivities and reproducible results**. Every other capability enters through an extension point around that core.

**Why the order matters.** Records NDR-0069 to NDR-0073 are architecture choices that are cheap to reserve in the Phase 2 interfaces and expensive to retrofit later. Their implementation still waits until Phase 4 is accepted.

**Ecosystem context:**
- **Matrix-free CPU code:** a MATLAB framework handles 128 million hexahedra on a 64 GB PC (Wang, Aage, Wu, Sigmund & Westermann, SMO 2025; cited via [a 2026 preprint](https://engrxiv.org/preprint/download/8045/13267/11508), not read).
- **GPU code:** 65.5 million elements in about 2 hours on one GPU, for linear compliance ([Träff et al. 2023, CMAME 410:116043](https://hgpu.org/?p=28174)).
- **Julia GPU code:** [GPUTopOpt (2026, GPL-3.0, doi:10.21950/R5MOML)](https://portalcientifico.uned.es/documentos/6a356e0375030556601f7a67) already offers matrix-free multigrid PCG, MMA and STL input/output.
- **Caveat for nonlinear problems:** a nonlinear matrix-free operator application cost 1075.4 ms against 1.316 ms for the linear case ([Träff 2023, thesis](https://backend.orbit.dtu.dk/ws/files/338184592/Erik_Albert_Tr_ff.pdf)). GPU gains for finite-strain problems must therefore be measured, not assumed.

| ID | Proposed decision | Class | FEAD2026 value → proposed default | Acceptance evidence | Status |
|---|---|---|---|---|---|
| NDR-0066 | **Scope and extension-point principle.** State the niche above. Make generality come from dispatch interfaces, not from built-in cases: physics, material, element, kinematics (NDR-0065), cell or lattice geometry, design representation, filter, projection, response, solver, optimizer and output. Every new feature arrives as interface + tests + decision record. Track performance regressions in CI. | P | Two forked scripts → one package with documented extension points | Interface documentation; each later record names its extension point; FEAD2026 golden masters bitwise unchanged after every feature merge; performance-regression job in CI | proposed |
| NDR-0067 | **Software license (decision owner: the authors).** The repository currently has no license, so nobody may legally reuse it. Options: MIT or BSD-3 (permissive, the Julia-ecosystem norm, maximizes adoption); MPL-2.0 (file-level copyleft); LGPL or GPL-3.0 (strong copyleft, compatible with GPL code such as GPUTopOpt but limits adoption). A permissive license means GPL code (GPUTopOpt; Svanberg's MATLAB MMA) may not be copied, so MMA stays a clean-room implementation from the papers (NDR-0035). Check institutional and funder (TÜBİTAK 122M134) requirements first. | P | No license → license chosen by the authors; `LICENSE` and `CITATION.cff` (NDR-0041) committed | License file present; dependency-license audit; no copied GPL code if permissive | proposed (awaiting authors) |
| NDR-0068 | **Documentation and example gallery** with Documenter.jl: tutorials that reproduce the FEAD2026 cases, the Wallin & Tortorelli §4.1 benchmark (NDR-0043) and the stability diagnostics; API reference; a "how to add a material, response or element" guide. | P | README only → versioned documentation site | Docs build in CI; every tutorial runs as a test | proposed |
| NDR-0069 | **3D cells.** Add hexahedral elements (H8, 2×2×2 Gauss) and 3D periodic constraints. Extend the domain reduction to 3D symmetries; the article's conclusions name 3D orthotropic microstructures as the next step. Reserve the dimension as a type parameter in Phase 2; implement after Phase 4. | B (new capability) | 2D only → `dim` type parameter, 2D default | 3D patch tests; 3D linear homogenization against published benchmark tensors; a one-layer 3D cell with free out-of-plane stretch matching the generalized-plane-strain mode of NDR-0065; gradient Taylor tests | proposed |
| NDR-0070 | **Iterative solver path** for 3D and large 2D problems: PCG with geometric multigrid on structured grids (matrix-free option); algebraic multigrid on unstructured meshes; Krylov solves inside Newton, with tolerances tied to the sensitivity accuracy. Direct CHOLMOD solves stay the 2D default (NDR-0022). Sources: [Amir, Aage & Lazarov 2014](https://doi.org/10.1007/s00158-013-1015-5); [Aage et al. 2015](https://doi.org/10.1007/s00158-014-1157-0). | T | Direct solver only → solver chosen by problem size | Same results as the direct solver within the measured bands; documented scaling study; adjoint consistency with inexact solves | proposed |
| NDR-0071 | **Macroscale structural problems:** non-periodic boundary conditions (supports, loads) on the same finite-strain Newton and adjoint machinery; responses such as end compliance and compliant-mechanism output displacement ([Buhl, Pedersen & Sigmund 2000](https://doi.org/10.1007/s001580050089); [Bruns & Tortorelli 2001](https://doi.org/10.1016/S0045-7825(00)00278-4)). | B (new capability) | Unit cells only → `ProblemType`: `PeriodicCell` (default) or `Structure` | Reproduction of a published finite-strain structural benchmark within stated tolerance; Taylor tests; linear-limit agreement with a linear code | proposed |
| NDR-0072 | **General response layer.** Objectives and constraints become functions of solver outputs (𝔸_M, P_M, homogenized energy W̄, u, eigenvalues). Their partial derivatives come from automatic differentiation (ForwardDiff/Enzyme), feeding one combined adjoint (NDR-0023). The hand-coded FEAD2026 objectives remain in the FEAD2026 profile. | T | Hard-coded objective types → `Response` interface with AD partials | Automatic-differentiation partials match hand-coded ones for Eqs. 53, 55 and 56 to round-off; Taylor slope 2 for every shipped response | proposed |
| NDR-0073 | **General lattices:** periodic cells with arbitrary lattice vectors (rectangular, oblique, hexagonal, triangular), through skewed structured grids or Ferrite meshes (revisits NDR-0042). The Bloch module takes the lattice vectors (NDR-0050). | B (new capability) | Rectangular cell only → `Lattice` type, rectangular default | Linear homogenization of a hexagonal cell reproduces isotropy where symmetry requires it; Bloch identities (E-S3) on non-rectangular lattices | proposed |
| NDR-0074 | **Load-path objectives:** targets defined along the load path F(t), such as a stress–strain curve or a constant ν over strain levels ([Wang, Sigmund & Jensen 2014](https://doi.org/10.1016/j.jmps.2014.05.003); the article's outlook). Sensitivities are summed over the retained load steps. | A (option) | Single-state objectives → optional path objectives | Taylor tests across load steps; an optimized design meets its target curve within a stated tolerance on re-analysis | proposed |
| NDR-0075 | **Material library and near-incompressibility:** Mooney–Rivlin, Gent, Ogden, and Wallin & Tortorelli's split neo-Hookean law (also needed for NDR-0043); nearly incompressible materials through u/p mixed elements ([Zhang, Alberdi & Khandelwal 2018](https://doi.org/10.1002/nme.5834), the article's ref. [40]). | A (option) | One neo-Hookean law → `Material` interface, FEAD2026 law as default | Automatic-differentiation tangent and Z tests per material; homogeneous analytical solutions; no locking in u/p patch tests | proposed |
| NDR-0076 | **Multi-material designs** with more than two phases ([Kazemi, Vaziri & Norato 2020](https://doi.org/10.1016/j.cma.2020.112895), the article's ref. [39]); multiphase bound checks extend NDR-0064 where applicable. | A (option) | Solid–void only → phase-count parameter | Two-phase limit reproduces solid–void results; gradient checks; multiphase linear-limit checks | proposed |
| NDR-0077 | **Robust length-scale control:** eroded, intermediate and dilated projections ([Wang, Lazarov & Sigmund 2011](https://doi.org/10.1007/s00158-010-0602-y)). This is also the eroded design that stability constraints need (NDR-0059). | A (option) | Single projection → optional robust formulation | Minimum length scale measured on final designs; reproduction of a published robust 2D example within stated tolerance | proposed |
| NDR-0078 | **Parallel and GPU execution:** deterministic colored threaded assembly first (NDR-0024); GPU kernels through KernelAbstractions.jl for assembly and matrix-free operators, adopted only where measurements show a benefit for finite-strain problems; MPI distribution (e.g., PartitionedArrays.jl) only when 3D sizes require it. GPL code is never copied (NDR-0067). | T | Serial → threads by default; GPU and MPI optional | Results within bands on 1, N threads and GPU; measured speedups reported per Hoefler & Belli; FEAD2026 profile unchanged | proposed |
| NDR-0079 | **Checkpoint, restart and ensembles:** serialize the full optimizer state (design, MMA asymptotes and history, iteration counters, continuation stage, RNG state) so that a resumed run replays bitwise; batch runner for the round-off ensembles (NDR-0010) and parameter sweeps. | B | No restart → checkpoint every N iterations | An interrupted and resumed run is bitwise identical to an uninterrupted one; ensemble runner reproduces NDR-0010 statistics | proposed |
| NDR-0080 | **Additive-manufacturing overhang control** for 3D: the layer-wise filter of [Langelaar 2016](https://doi.org/10.1016/j.addma.2016.06.010) and [Langelaar 2017](https://doi.org/10.1007/s00158-016-1522-2) on structured grids; front propagation on unstructured meshes ([van de Ven et al. 2020](https://research.tudelft.nl/en/publications/overhang-control-based-on-front-propagation-in-3d-topology-optimi)). How to apply it to periodic cells is open (inference). | A (option) | None → optional filter | Printability check on final designs; reproduction of a published self-supporting 2D example | proposed |
| NDR-0081 | **Output and interoperability:** VTK field output (WriteVTK.jl); STL export of thresholded designs; body-fitted re-analysis export through Gmsh for verification (supports NDR-0054 and the stability verification); plotting as a Makie package extension. | B | PNG frames and a video → structured outputs; plots optional | Round-trip tests; body-fitted re-analysis of a regenerated FEAD2026 design within stated tolerance | proposed |
| NDR-0082 | **Alternative design representations** (level set, feature mapping or geometry projection) deferred; the design-representation interface is reserved in Phase 2. Level-set work can use GridapTopOpt ([Wegert et al. 2025](https://doi.org/10.1007/s00158-024-03927-3)). | P (deferred) | Density only → interface reserved | — | proposed (deferred) |

**Suggested order:**
1. NDR-0066 and NDR-0067 now.
2. The interfaces for NDR-0069 to NDR-0073 in Phase 2.
3. Implementation of NDR-0069 to NDR-0072 after Phase 4.
4. NDR-0073 to NDR-0081 by research need.

## Revisions after full-text verification (2026-10-12)

The supplied papers and the open-access sources were read in full and checked against this report. The verification notes are in [NEMOpt supplied papers verification](../research_notes/NEMOpt%20supplied%20papers%20verification/). Each note gives its recommended edits as "old text → new text". Short names below: "GMT note" (`geymonat_muller_triantafyllidis_1993.md`), "Michel note" (`michel_et_al_2007.md`), "B&S note" (`bendsoe_sigmund_1999.md`), "mechanics note" (`open_access_mechanics_sources.md`) and "optimization note" (`open_access_optimization_papers.md`). "Coordinator fact n" refers to the coordinator's instructions for this revision. Every decision record keeps its status: all remain proposed, NDR-0067 remains "proposed (awaiting authors)" and NDR-0082 remains "proposed (deferred)". Edits were located by text, because the notes' line numbers predate NDR-0064 to NDR-0082. Every "old" string matched the current wording unless a bullet says otherwise.

**Applied: panel, oracles and plausibility check**
- Panel table, mechanician row: the symmetric-branch restriction is now marked as derived; Ohno et al. 2001 only assume a point-symmetric internal state. Why: no open source on the reduction discusses stability. Source: mechanics note, Item 1 (main L36 edit).
- Panel table, methodologist row, and the L4 row of the tolerance table: Sigmund 2022 marked "(abstract only)". Why: no legal full text exists. Source: optimization note, §3 edit 2.
- κ plausibility paragraph: the 28.85 bound now cites Bendsøe & Sigmund 1999, Eq. 6. Source: B&S note §4.3 (optional edit).
- κ plausibility paragraph, list lead-in: "caveats" changed to "points", because the third item is now a confirmation. Second bullet: added that the 1.6% excess is expected, because p = 3 lies below the plane-strain threshold p\* = 3.5. Source: B&S note §4.3 (L173 edit).
- κ plausibility paragraph, third bullet: "an assumption still to be confirmed" replaced by the confirmation. Milton 2002, Eq. 23.30 bounds the hydrostatic-strain modulus of any anisotropic two-phase cell; Bendsøe & Sigmund Eq. 21 gives the same result (derived). Why: Milton's §23.4 states the bound for any geometry of two well-ordered isotropic phases. Sources: mechanics note, Item 3 (main L171 edit); B&S note §4.3 (L174 edit, applied in its confirmed form); coordinator fact 4.
- Sigmund paragraph after the acceptance-band sentence: now quotes the abstract with its condition ("without practical common benchmarks"). The frozen evaluator, thresholding, M_nd and IoU are now this report's own recommendation, with M_nd credited to Sigmund 2007. Why: the full text was not obtained. Source: optimization note, §3 edit 1.

**Applied: behavior, erratum and speedup tables**
- Rounding row: the 1e−7 tolerance is now attributed to the epsimin slack of the GPL Python port, not to "Svanberg's GCMMA code", and the port files are named. Why: Svanberg's 2007 notes could not be read. Source: optimization note, §5 edit 1.
- Clip row, "Status and origin" cell: added that classic SIMP floors the density variable (Bendsøe & Sigmund 1999, p. 638). Source: B&S note §4.4 (optional edit).
- Void-switch row and NDR-0034: Dalklint et al. 2021, Eq. 9 added as a second precedent for Ĥ(H(ρ̃)^p) with (500, 0.01). Source: optimization note, §2 edit 4.
- Disagreement table, Poisson's-ratio row: the article's credit to Kazemi et al. is now marked incorrect. Kazemi et al.'s Eq. 24 (p. 8) is ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)] for 3D cubic lattices. The thesis cell notes the same attribution in Eq. 7.5, and the "Treat as" cell adds that an erratum is the authors' decision. Sources: mechanics note, Item 2 (main L240 edit); coordinator fact 5. The thesis and "Treat as" additions go beyond the note's edit and follow coordinator fact 5.
- Speedup table, filter row: Lazarov & Sigmund 2011 marked "(filter definition; full text not read)". Source: optimization note, §6 edit 2.
- Paper-traceability list: the SMO "Replication of results" item now cites p. 405 and summarizes the editorial's recommendations (pp. 406–409). Why: the editorial was read in full, and it fixes no content. Source: optimization note, §4 edit 1.

**Applied: decision log NDR-0001 to NDR-0044**
- NDR-0010: added SMO's advice to make repeated, slightly perturbed runs (Haftka et al. 2019, p. 407). Source: optimization note, §4 edit 2.
- NDR-0033: added three facts, all quoting Bendsøe & Sigmund 1999. The clip equals classic SIMP evaluated at max(H, ε^{1/p}); modified SIMP has the form of their two-material SIMP (Eq. 23); bound compatibility does not separate the two laws. Decision and default unchanged. Source: B&S note §4.4.
- NDR-0037: added that Eq. 54 does not appear in Kazemi et al. 2020 and that the article's and thesis Eq. 7.5's attribution is incorrect; an erratum is the authors' decision. The decision to keep Eq. 54 is unchanged. Sources: mechanics note, Item 2 (NDR-0037 edit); coordinator fact 5.
- NDR-0039: the r = R/(2√3) convention now cites Lazarov, Wang & Sigmund 2016, Eq. 10 (a second-moment match on an unbounded domain). Source: optimization note, §6 edit 1.
- NDR-0044: appended the SMO editorial's "detailed settings for results need to be documented" (p. 406). Source: optimization note, §4 edit 3.

**Applied: NDR-0064 prose**
- Talbot & Willis bullet: now "Hashin–Shtrikman variational principles", "geometrically linear setting", and page numbers in Milton 2002 (p. 271; §13.7, pp. 282–286). The author-hosted copy is linked; the existing Cambridge link is kept. Source: mechanics note, Item 3 (Talbot bullet).
- Ponte Castañeda bullet: the secondary sources are now named: Milton 2002 for the 1991 bounds; Geymonat et al. 1993, p. 237 and Michel et al. 2007 for the 1989 paper. The label stays "secondary". Sources: mechanics note, Item 3 (optional Ponte Castañeda edit); GMT note, Q1 claim 19; Michel note, claim 3 (neither has an E-number).
- Convexity and multi-cell bullet:
  - The Müller (1987) link pointed to Duerinckx & Gloria (2016). It now points to doi:10.1007/BF00284506.
  - The convexity wording is refined: geometrically linear setting, growth conditions, and convexity for the dual principles.
  - The multi-cell infimum is stated with Geymonat et al.'s Thm 2.3 and Eqs. 2.6–2.7, Michel et al.'s Eq. 2.22 and Milton's §1.7.
  - Added that the one-cell value stays exact while the cell is Bloch-stable (Prop. 5.3).
  - Müller's evidence label changes from a reference-list citation to "verified as restated"; Müller's own text was not read.
  - Sources: GMT note E10 and its evidence-label list; Michel note E2; mechanics note, Item 3 (convexity bullet); coordinator fact 6. The three notes give overlapping replacements for this bullet. They were merged: the mechanics note's convexity wording, then GMT E10's content, with Michel's and Milton's restatements added.
- Porous-rubber bullet: replaced. The second-order results are Hashin–Shtrikman-type estimates, not bounds. They exist for periodic arrays too, are not conservative for periodic cells, and none exists for NEMOpt's topologies. Sources: Michel note E1; coordinator fact 7. Michel et al. 2007 is labeled "HAL submitted version", not the note's "author version", per coordinator fact 7. The supplied PDF, Lopez-Pamies & Ponte Castañeda 2007 (JMPS 55:1677–1701), is cited as read in full.
- "Computed" list, bulk bullet: added that the values equal Bendsøe & Sigmund Eq. 6 (void limit) and Milton Eq. 23.30 at hydrostatic strain, with a link to the bound-check scripts. Sources: B&S note, claim C2; mechanics note F3.1 and its evidence note.
- "Computed" list, new shear bullet: the corrected isotropic shear bounds (7.40 / 10.12 / 13.42) and the earlier, 6–9% too high values (7.96 / 10.82 / 14.23). It also gives the single-modulus bound from Milton Eq. 23.30 (12.41 / 16.03 / 19.89; 19.905 with the ε-material). Why: the scripts' 2D shear formula was wrong and was fixed in commit 5aa11db. Sources: coordinator facts 1 and 4; B&S note §1.4; mechanics note, side finding and F3.4.
- "Computed" list, p = 3 bullet: replaced. The hedges on the attribution and on the formula's exact form are removed. The condition is Bendsøe & Sigmund 1999, Eq. 11, p. 639, and ν⁰ is the in-plane Poisson's ratio: p\* = 3.5 under plane strain, 3.08 (set by shear) under plane stress. The bullet also gives the analytic peak (1.0163 at H = 14/15) and the plane-strain shear ratio 0.99998. Sources: B&S note §4.1 and §1.3; coordinator fact 3.
- "Consequence" bullet: appended Bendsøe & Sigmund's view that bound compatibility of gray is optional for black-and-white designs. Source: B&S note §4.1 (optional edit).
- Voigt-type paragraph: the effective energy is now an infimum over supercell-periodic fluctuations. NEMOpt's one-cell value bounds it from above, and the computed energy obeys the bound only at the one-cell minimizer. Source: GMT note E11.

**Applied: NDR-0064 and NDR-0065 records**
- NDR-0064 (a), scope: the bulk bound holds for any thresholded two-phase cell at F_M = I and can serve as a hard gate there. Source: coordinator fact 4 (mechanics note, Item 3).
- NDR-0064 (a), gray designs: the 1.6% excess is now explained by p = 3 < p\* = 3.5. Source: B&S note §4.2(a).
- NDR-0064 (a), shear: the isotropic shear bound needs isotropy. For a single shear modulus, use Milton Eq. 23.30 at pure shear (12.41 / 16.03 / 19.89; 19.905 with ε), which equals Bendsøe & Sigmund's Eq. 21 limit (39.8 on 𝔸₁₂₁₂ + 𝔸₂₁₂₁). The square-symmetric trace bound is optional. Sources: B&S note §4.2(b); mechanics note, Item 3 (NDR-0064 (a) edit); coordinator fact 4.
- NDR-0064 (c): second-order estimates are named. They predict only macroscopic onset and overestimated finite-element onset strains by up to about a factor of two, so they never certify stability. Source: Michel note E3.
- NDR-0064 acceptance cell:
  - Bound values tie to Bendsøe & Sigmund Eq. 6 and Milton Eq. 23.30, and the corrected shear values were added (commit 5aa11db).
  - "Confirm in the composites literature before a hard gate" is replaced by the Milton confirmation, so (a) can be a hard gate for thresholded two-phase designs at F_M = I.
  - Sources: B&S note §4.2(c) and (d); mechanics note, Item 3 (acceptance edit); coordinator facts 1 and 4.
- NDR-0065 table: the plane-stress cell said p = 3 SIMP "stays at or just below" the bound without qualification. It now says this holds for the bulk bound (ratio 0.9999), while p = 3 SIMP exceeds the shear bound by up to 0.05% (ratio 1.00047 at density 0.988). The row label now names both bounds, and the plane-strain cell adds the shear ratio 0.99998. Sources: coordinator fact 2; B&S note §4.7. No other prose in either report made the unqualified claim.
- NDR-0065 table, new row "SIMP threshold p\*": 3.5 (bulk-governed) under plane strain and 3.08 (shear-governed) under plane stress. Source: B&S note §4.7.
- NDR-0065 lead-in sentence: the shear ratios come from the corrected formula (commit 5aa11db). Source: coordinator fact 1.
- NDR-0065 record, requirements: added that the plane-strain bound stays valid, possibly loose, in generalized plane strain, while the plane-stress bound is not guaranteed there (inference). Source: B&S note §4.7 (generalized plane strain).

**Recommended edits not applied, or applied in modified form**
- B&S note §4.2(c) and §4.3 (L174): "cross-check Eq. 21 against Allaire & Kohn or Milton before (a) becomes a hard gate". Not applied as worded, because the cross-check is done: Milton 2002, Eq. 23.30 confirms the bound (coordinator fact 4).
- B&S note §4.6 (evidence-script edits): not applied here. Scripts were outside the files this revision may edit, and the fix is already in commit 5aa11db (coordinator fact 1).
- Mechanics note, Talbot bullet: the proposed text drops the Cambridge link. The link was kept and the author-hosted copy added, because no note shows the link to be wrong.
- GMT note E10, first sentence ("They require convexity that finite-strain hyperelastic energies lack"): replaced by the mechanics note's more precise wording, which rates that sentence "partly supported"; E10's substance is kept in full.
- Optimization note, §5 edit 5 (cite no Svanberg 2007 parameter value until a snapshot is read): a standing rule, not a text edit. The reports cite no such value.
- NDR-0005 errata list: not extended with the Kazemi attribution, because no note proposes it. The correction is stated in the Poisson's-ratio row and in NDR-0037, and whether it becomes an erratum is the authors' decision.
- NDR-0076 still cites Kazemi, Vaziri & Norato 2020, for multi-material design (the article's ref. [39]). That citation does not concern Eq. 54, so it was kept.
