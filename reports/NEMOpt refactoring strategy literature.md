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
| Finite-strain computational-homogenization mechanician | Preserve the mechanics exactly: Bloch-type total-displacement boundary conditions, the fixed center node (full domain) or six prescribed nodes (half domain), the symmetric total-field χ-form of 𝔸_M, and the parity rules of the point-symmetric half cell (odd fluctuations, χ fields and adjoints; even densities). A 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear load case, but only on the symmetric branch. Full-versus-half agreement ends at the first symmetry-breaking bifurcation. | Article §3, Eqs. 41–50 ([Güven & Özdemir 2026](https://doi.org/10.1016/j.finel.2026.104599)); Wallin & Tortorelli Eqs. 25, 32–35 and footnote 3 ([2020](https://doi.org/10.1016/j.mechmat.2020.103324)); [Ohno et al. 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) |
| Density-based topology-optimization methodologist | The result is a function of the optimization path. Rounding the objective, counting continuation per function evaluation, and switching between MMA87 and GCMMA each change the path, so none of them may change silently. Count continuation in design-loop iterations, as every source does, and compare designs after thresholding with a common evaluator. | [Stolpe & Svanberg 2001](https://doi.org/10.1007/s001580050177); [Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w); [Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2) |
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
| L4 full runs | J, Vf, ν, κ histories; iterates; MMA counters; final design | Bitwise trajectory replay | Inside an ensemble of round-off-perturbed legacy runs, on pre-registered metrics | Article values at printed precision; thresholded re-analysis, measure of non-discreteness (M_nd), intersection-over-union (IoU); Voigt-type energy bound at F_M; Hashin–Shtrikman bulk bound only on the linear-limit tangent 𝔸_M(F_M = I) of thresholded designs (NDR-0064) | [NASEM 2019](https://doi.org/10.17226/25303); [Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2); [Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x) |

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

These come from [article §5.2–5.3](https://doi.org/10.1016/j.finel.2026.104599). From Eq. 55, the two ν targets imply objective values of about 0.984 and 0.753 (inference). One plausibility check supports the code's bulk definition: κ = 28.4368 is 98.6% of the two-dimensional Hashin–Shtrikman upper bound of 28.85 for a porous solid at 60% solid (inference, from λ + μ = 96.15). This comparison is meaningful only because §5.3 case 1 (F_M = 1.0001·I) is essentially linear. Hashin–Shtrikman bounds are theorems for linear elastic moduli, so no such bound applies to the finite-strain tangent values of the other cases. Three further caveats apply:
- the ε = 5e−4 floor raises the bound to 28.905 (computed);
- p = 3 SIMP itself exceeds the bulk bound by up to 1.6% in gray elements under plane strain (computed);
- applying the isotropic bulk bound to a square-symmetric cell's hydrostatic-strain modulus is an assumption still to be confirmed.

NDR-0064 records how bounds are used from now on.

Starting bands of ±0.02 on ν and 1–2% on κ are experience values. They must be replaced by the measured spread of the legacy ensemble before any comparison is run, because the thesis itself notes that a different initial design can lead to a different optimum topology ([thesis §7.1.2, pp. 54–55](https://gcris.iyte.edu.tr/handle/11147/13904)). Following Sigmund's point that fair evaluation is the authors' responsibility, final designs should be re-analyzed with one frozen evaluator, thresholded and compared on M_nd and IoU ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2)).

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
| J and ν rounded to 3 decimals; the rounded J is returned to MMA and the rounded ν enters ∇J = 200(ν − ν₀)∇ν ([L79–80, L104, L117](../Full_Domain/NEMOpt_main.jl#L79)) | Verified likely defect; in no source | **A.** The returned value carries ±5e−4 noise, far above the 1e−7 conservativeness tolerance in Svanberg's GCMMA code; NonconvexMMA's check has no tolerance at all ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python)). The gradient factor error is up to 5e−4/\|ν+1\|, about 0.5% at ν = −0.90. With kkt = 1e−8, KKT termination becomes practically unattainable (inference). | Keep | Unrounded by default; round only in logs |
| p/β continuation advanced per objective evaluation, inside `Objective`, before `Objective_Grad` runs in the same rrule. The full-domain schedule (`β < 32`, then +2 from odd values) overshoots to β = 33 ([L45, L90–102, L213–216](../Full_Domain/NEMOpt_main.jl#L90)). | Verified likely defects; the article, thesis and Wallin & Tortorelli all count design iterations and state terminal values | **A.** At each step, J uses (p_old, β_old) but the gradient uses (p_new, β_new). This happens on about 50 of 300 half-domain evaluations (computed). The schedule also depends on how often the optimizer rejects iterations. | Keep, including the overshoot in the as-shipped preset | Per MMA outer iteration, applied between iterations; schedules stored as `{istart, max, every, delta}` data with explicit terminal values ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)) |
| Newton control: residual tested before the update (one extra correction per load step); `Rep_MU = \|max(F_M − I)\|·Lx`, which is 0 for any compressive F_M and then loops forever; no iteration cap or line search; 10 fixed load steps from U = 0 ([NEMOpt.jl L372–396](../Full_Domain/NEMOpt.jl#L372)) | Verified; Wallin & Tortorelli describe no solver; differs from thesis Algorithm 5.2 | **A** for moving the check, changing the scale or warm-starting. **B** for guards that fire only where legacy would hang or diverge. For F_M = 1.0001·I the absolute target is 2e−11, a plausible but unverified hang risk. | Keep, plus an explicit error instead of a hang | Check after the update; a scale that is never zero; tighter tolerance; iteration cap; bisection fallback on failure |
| Optimizer: MMA87 (full domain) versus MMA02 (half domain); `s_init` 0.25 and `maxinner` 40, against library defaults of 0.5 and 10 ([L235–237](../Full_Domain/NEMOpt_main.jl#L235); [DR L234](../Domain_Reduction/NEMOpt-DR_main.jl#L234)) | Verified. MMA87 matches Wallin & Tortorelli's plain MMA with "the numerical parameter values provided in Svanberg (1987)", wording the article echoes. MMA02 follows the thesis. No source documents `s_init` or `maxinner`. | **A.** The two shipped drivers cannot produce the "same optimum" of Fig. 8 under different optimizers. | Keep per preset | Explicit configuration; identical optimizer for any full-versus-half comparison |
| Stiffness interpolation `max(H^p, ε)` with zero derivative below ε = 5e−4 ([L59–60](../Full_Domain/NEMOpt.jl#L59)) | Verified reference formulation: Wallin & Tortorelli Eq. 71 (ε = 1e−5 there). The code's ε is stated in neither the article nor the thesis. | **A.** The clip creates a non-C¹ kink and a zero-sensitivity zone, while GCMMA's convergence theory assumes continuously differentiable functions ([Svanberg 2002](https://doi.org/10.1137/S1052623499362822); assumption not re-verified). Stiffness differs from modified SIMP by at most ε. | Keep, with ε recorded | **Keep the clip as the default.** Offer modified SIMP ε + (1 − ε)H^p ([Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)) only as an opt-in deviation behind the `stiffness_interp` switch, recorded as a deviation from both Wallin & Tortorelli and FEAD 2026, not as a bug fix. |
| Void switch Ĥ(ρ̃^p) with β̂ = 200, ω̂ = 0.05 ([L297–299](../Full_Domain/NEMOpt.jl#L297)) | Verified reference formulation. Matches the thesis and Wallin & Tortorelli Eq. 73 in argument and threshold; steepness is 200 against their 500, with no rationale anywhere. | **A** | Keep | Keep as default. Offer β̂ = 500 and a physical-density switch, Ĥ(H(ρ̃)^p) with (500, 0.01) ([Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4)), as exploratory options. Add distortion diagnostics. |
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
| Poisson's ratio | Eq. 54 (after Kazemi et al.), an objective toward −1 | Eq. 7.5 | Eq. 54 | Tangent uniaxial-stress −S₂₂₁₁E₁ (Eqs. 69, 96), as a constraint ν ≤ ν₀ | Eq. 54, a deliberate and documented deviation. On their uniaxial state, Eq. 54 gives 0.5596 against their 0.3561 (computed), so the two must never be compared numerically. |
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
| Assemble and factor the Helmholtz filter once per run; share it across objective, constraint and gradients | top88's PDE variant factors once outside the loop ([Andreassen et al. 2011](https://doi.org/10.1007/s00158-010-0594-7)); [Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072) | Reassembled on every call, and repeated inside `Volume_Constraint` ([L128–130](../Full_Domain/NEMOpt_main.jl#L128)) | B (same factorization algorithm) | Removes every per-evaluation filter assembly and factorization |
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
- the mandatory "Replication of results" section of *Structural and Multidisciplinary Optimization* ([Haftka et al. 2019](https://doi.org/10.1007/s00158-019-02298-4));
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
| NDR-0010 | Full-run acceptance: a legacy ensemble (size, perturbations, and metrics J, ν, κ, Vf, M_nd, thresholded ν, IoU) plus the article values at printed precision | P | — | Pre-registered bands | proposed |
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
| NDR-0033 | Stiffness interpolation. The clip `max(H^p, ε)` is Wallin & Tortorelli's reference law (Eq. 71, ε = 1e−5 there). | A (option only) | `max(H^p, ε)` with ε = 5e−4 → **unchanged default**; modified SIMP ε + (1 − ε)H^p as an opt-in deviation behind `stiffness_interp`, recorded as departing from both sources | Article presets re-run when the option is used; deviation row | proposed |
| NDR-0034 | Void (linear/nonlinear) switch, which matches the thesis and Wallin & Tortorelli Eq. 73 in argument and threshold. Its soft-element band is inherited from that scheme. | A (option only) | Ĥ(ρ̃^p) with (200, 0.05) → unchanged default; exploratory options β̂ = 500 and Ĥ(H(ρ̃)^p) with (500, 0.01); distortion diagnostics | Controlled study before any default change | proposed |
| NDR-0035 | Optimizer variant and options per preset. Full-domain MMA87 matches Wallin & Tortorelli's plain MMA; half-domain MMA02 follows the thesis. | A | MMA87 / MMA02, `s_init` 0.25, `maxinner` 40, kkt 1e−8 → explicit; identical optimizer in any full-versus-half comparison; new MMA implementations only by new record | Trajectory comparison | proposed |
| NDR-0036 | Initial-design ingestion and RNG | B/A | legacy crop, resize and rescale → hashed design-vector fixtures (B); corrected loader optional (A); StableRNGs for random starts | Fixture hashes | proposed |
| NDR-0037 | Effective-property definitions: keep ν (Eq. 54, a deliberate deviation from Wallin & Tortorelli's uniaxial-stress ν), κ and the shear sum. Rename the logged G to G₁₂ = 𝔸₁₂₁₂ (Wallin & Tortorelli's μ₁₂) and G_Σ = 𝔸₁₂₁₂ + 𝔸₂₁₂₁ (Eq. 53). Document them as tangent measures with ν in (−1, 1). Add an optional uniaxial-stress ν₁₂ diagnostic; never compare it numerically with Eq. 54. | P | — | Documentation and tests | proposed |
| NDR-0038 | Re-enable Eqs. 53 and 56 as objective types, with δ𝔸₂₂₁₁ = δ𝔸₁₁₂₂ by major symmetry | A (new capability) | commented out → selectable | Numerical symmetry check; gradient check | proposed |
| NDR-0039 | Filter length as a physical length with an explicit convention (Helmholtz r versus cone R ≈ 2√3·r); Wallin & Tortorelli's r = 2 elements equals the article's text value | P | r·Lx/Nelx with r = 1 element | Preset tests | proposed |
| NDR-0040 | Box constraint [0, 1] as coded; document the article's ρ_min | P | [0, 1] | Documentation | proposed |
| NDR-0041 | Provenance on every result file (git describe and patch, Manifest hash, `versioninfo()`, BLAS, solver path, profile flags, parameters, design hash); `CITATION.cff`; Zenodo DOIs for the legacy, bridge and refactored releases | P | — | Audit of the figure-to-file map | proposed |
| NDR-0042 | Keep in-house structured-grid assembly; defer Ferrite until a release untangles master–slave chains and a T-class record justifies the move | P | — | — | proposed |
| NDR-0043 | Add Wallin & Tortorelli §4.1 as an independent golden test. Requires a pluggable energy (their 2D form with the printed B = μ/8 typo kept, since the published numbers depend on it), a pure-SIMP mode S = ν³P with filter and projection bypassed, element-group densities assigned by centroid, and the index map (1, 2, 3, 4) → (1, 4, 3, 2). Port the evidence scripts as reference implementations. | B (new test modes; FEAD2026 unaffected) | — | Case A and case B S_M and 𝔸_M (including 𝔸₁₂₁₂ = −0.9861, F₂₂ = 0.872199, ν₁₂ = 0.356115) to printed digits; Table 1 sensitivities by symmetric adjoint and central differences (1e−5) to about 1e−6 relative; degree-3 Euler identity | proposed |
| NDR-0044 | Record in the profile and in the next paper's methods the FEAD2026 values that no source documents or justifies: ε = 5e−4 (against 1e−5), β̂ = 200 (against 500), `s_init` 0.25, `maxinner` 40, TOL 1e−7 scaled by Rep_MU, and 10 load steps. Optionally add a "Wallin & Tortorelli formulation" preset (β̂ = 500, ε = 1e−5, r = 2 elements, doubling β schedule) for formulation-sensitivity studies only, without claiming to reproduce their numbers, since energy and kinematics differ. | P (preset: A) | — | Documentation; preset tests | proposed |

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

The decision log continues there as **NDR-0045 to NDR-0063**, all proposed, and NDR-0064 is recorded in the next section. Two items in the addendum change records here:
- **The panel-table sentence** "Full-versus-half agreement ends at the first symmetry-breaking bifurcation" is amended. The *guarantee* of agreement ends there, but agreement itself may continue on an unstable symmetric branch. Agreement between the full and half cell is therefore never evidence of stability (addendum NDR-0051).
- **NDR-0022's "logged LDLᵀ/LU fallback"** is refined. UMFPACK LU serves the solve after a Cholesky failure, and LDLᵀ runs only as a logged diagnostic, because it is unpivoted and 10–21× slower than Cholesky (addendum NDR-0047).

## Decisions recorded after the addendum

### NDR-0064: elastic bounds are checks only where they hold

**Hashin–Shtrikman bounds are theorems for linear elastic moduli** ([Hashin & Shtrikman 1963, JMPS 11:127](https://garfield.library.upenn.edu/classics1988/A1988Q027300001.pdf)). They therefore cannot validate the finite-strain tangent values (𝔸_M, ν, κ at F_M ≠ I) that NEMOpt reports.

**Nonlinear extensions exist but do not supply a tangent check:**
- Talbot & Willis (1985) extended the variational principles to nonlinear media ([Milton, *Theory of Composites*](https://resolve.cambridge.org/core/books/theory-of-composites/variational-principles-and-inequalities/F55B22AA23653A981F56C1D0C9F22D78)).
- Ponte Castañeda (1991, JMPS 39:45) gave Hashin–Shtrikman-type bounds through a linear comparison composite, and Ponte Castañeda (1989, Proc. R. Soc. A 422:147) extended the approach to finite elasticity ([summary](https://apps.dtic.mil/sti/pdfs/ADA231637.pdf); verified only through secondary sources).
- These bound the effective *energy*, not tangent moduli.
- They require convexity that finite-strain hyperelastic energies lack. Once instabilities appear, the homogenized energy may need a multi-cell infimum ([Müller 1987, ARMA 99:189](https://link.springer.com/doi/10.1007/s00205-016-0992-0), cited from a reference list).
- The finite-strain results for porous rubbers are estimates for random microstructures, not bounds for periodic optimized cells ([Lopez-Pamies & Ponte Castañeda 2004, J. Elasticity 76:247](https://dialnet.unirioja.es/ejemplar/108150)).

**Computed against the FEAD2026 constants** (plane strain, E = 100, ν = 0.3, so κ_s = 96.15 and μ_s = 38.46; [script](../research_notes/NEMOpt%20elastic%20bounds%20check/hs_simp_check.py)):
- **2D Hashin–Shtrikman upper bulk bound for voids:** 15.38, 21.37 and 28.85 at solid fraction 0.4, 0.5 and 0.6. Modeling the void as the ε = 5e−4 material raises these to 15.44, 21.43 and 28.905.
- **p = 3 SIMP exceeds the bulk bound** by up to 1.6% at a density of about 0.93. The excess vanishes at p = 3.5, which matches the 2D compatibility condition usually attributed to Bendsøe & Sigmund (1999), p ≥ max{2/(1 − ν₀), 4/(1 + ν₀)}, evaluated with the in-plane Poisson's ratio ν₀ = 0.4286 ([Bendsøe & Sigmund 1999](https://link.springer.com/doi/10.1007/s004190050248); the formula's exact form in the paper was not confirmed).
- **Consequence:** a design that still contains gray elements can exceed the bound without any bug.

**The one rigorous and simple finite-strain bound is Voigt-type.** Periodic homogenization defines the effective energy as a minimum over periodic fluctuations. Taking zero fluctuation as the trial field gives W̄(F̄) ≤ (1/|Ω|)∫ W_ρ(x, F̄) dV (standard derivation). The bound holds only if the computed state is the energy minimizer, so a violation flags either a bug or a non-minimizing, possibly unstable, branch. In linear elasticity, exceeding the Voigt bound is likewise tied to loss of stability ([Kochmann & Milton 2014](https://arxiv.org/pdf/1401.4142)).

| ID | Proposed decision | Class | FEAD2026 value → proposed default | Acceptance evidence | Status |
|---|---|---|---|---|---|
| NDR-0064 | Use elastic bounds only where they hold. **(a) Linear limit:** apply the 2D Hashin–Shtrikman upper bound on the hydrostatic-strain bulk modulus κ = ¼(𝔸₁₁₁₁ + 𝔸₂₂₂₂ + 𝔸₁₁₂₂ + 𝔸₂₂₁₁) only to 𝔸_M(F_M = I). Apply it to thresholded designs, with the void modeled as the ε-material. Never apply it to raw gray designs (p = 3 SIMP exceeds the bound by up to 1.6%) or to finite-strain tangents. Do not apply the isotropic shear bound to anisotropic cells. **(b) Finite strain:** add a homogenized-energy output (new output; FEAD2026 arithmetic untouched) and check the Voigt-type bound W̄(F̄) ≤ (1/\|Ω\|)∫W_ρ(x, F̄) dV. Treat a violation as a bug or a non-minimizing branch, cross-referenced with the stability records NDR-0045 to NDR-0050. **(c)** Cite nonlinear Hashin–Shtrikman-type bounds and estimates for comparison only, never as acceptance tests. **(d)** This record amends the L4 oracle list and qualifies the §5.3 κ plausibility statement. | P (checks); B (energy output) | No bound checks; κ compared with the isotropic bound in prose → (a) and (b) as automated tests | Bound formulas reproduce the computed values (15.38 / 21.37 / 28.85; 28.905 with ε); F_M = I evaluations of thresholded regenerated designs stay below the bound within the measured noise band; the Voigt-type energy check passes on all golden-master states; before (a) becomes a hard gate, confirm in the composites literature (e.g., Milton 2002) that the bulk bound holds for the hydrostatic-strain modulus of anisotropic cells | proposed |

