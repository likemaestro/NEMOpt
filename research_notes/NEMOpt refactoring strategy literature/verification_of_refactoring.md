# Verifying that a refactored NEMOpt reproduces the published results, and recording refactoring decisions for a new TO research paper

Scope: how to define and show "same results" when refactoring NEMOpt (Julia; Güven & Özdemir, FEAD 260 (2026) 104599, doi:10.1016/j.finel.2026.104599) at four levels: (a) kernels, (b) single deterministic evaluations, (c) sensitivities, (d) full optimization runs. Also how to record every refactoring decision so a new, separate topology-optimization research paper that uses the refactored code can cite exactly what was run, how it was verified against FEAD 2026, and which formulation choices were deliberately changed. Researched October 2026.

Source labels: **[PR]** peer-reviewed journal or proceedings; **[Book]**; **[Report]** technical report, standard, or consensus report; **[Docs]** official software documentation; **[Grey]** blog, website, template, or secondary summary. **(mem)** means the bibliographic details come from prior knowledge and were not re-fetched in this session. Check the DOI of every (mem) item before final citation. **[repo]** means a direct observation of the NEMOpt source at the stated file and line (read-only inspection, Oct 2026).

---

## Q1. What do the V&V and scientific-software-testing literatures recommend for code verification and regression testing, and which metamorphic relations apply to NEMOpt?

### Takeaway
The V&V literature separates *code verification* (is the math implemented correctly? Use exact or manufactured solutions and observed order of accuracy) from *solution verification* (how large are the numerical errors in this run?). The software-testing literature adds that scientific codes have an oracle problem. The standard answers are pseudo-oracles (independent implementations, here the legacy code), analytical solutions, regression comparisons with tolerances, and metamorphic relations (MRs). For NEMOpt, several strong MRs and analytical oracles exist: the homogeneous-solid patch test, the small-strain limit, power-of-two scaling, periodic translation, 90° rotation and reflection, and full-domain vs. domain-reduction equivalence for symmetric designs. Together they cover what a single golden-master comparison cannot.

### Cited Findings
- **Oracle problem, systematically reviewed.** Kanewala & Bieman (2014) reviewed 62 studies. They group testing challenges into (1) those from the nature of scientific software, such as oracle problems, and (2) cultural differences between scientists and software engineers. — [PR] Kanewala, U., Bieman, J.M. (2014). Testing scientific software: A systematic literature review. *Information and Software Technology* 56(10):1219–1232. [doi:10.1016/j.infsof.2014.05.006](https://doi.org/10.1016/j.infsof.2014.05.006); full text [arXiv:1804.01954](https://ar5iv.labs.arxiv.org/html/1804.01954)
- **Oracle techniques the review lists:** pseudo-oracles ("an independently developed program that fulfills the same specification as the program under test"; caveat: "it is difficult to decide whether the oracle or the program is faulty when the answers do not agree"), analytical solutions (for example, comparing a theoretical convergence rate with the code's output), experimental results, simplified data, statistical oracles, reference data sets, metamorphic testing, and "golden run" traces. — [PR] Kanewala & Bieman 2014, [arXiv copy](https://ar5iv.labs.arxiv.org/html/1804.01954)
- **Round-off and tolerances are a named challenge:** "Faults can be masked by round-off errors, truncation errors and model simplifications"; "Choosing suitable tolerances for an oracle when testing numerical programs is difficult"; "Execution of some paths in scientific software are dependent on results of floating point calculations." Kelly et al. found that tightening the oracle tolerance raises the fault-detection rate. Clune et al. found that testing steps independently reduces the compounding of truncation and round-off errors. — [PR] Kanewala & Bieman 2014, [arXiv copy](https://ar5iv.labs.arxiv.org/html/1804.01954)
- **Regression testing in practice:** several reviewed studies compare "the current output to previous outputs", and one tool lets testers set "tolerances for comparisons". The review also notes "a lack of automated regression and acceptance testing in some scientific programs." — [PR] Kanewala & Bieman 2014, [arXiv copy](https://ar5iv.labs.arxiv.org/html/1804.01954)
- **Oracle-problem survey (general software engineering):** covers specified, derived (including metamorphic), and implicit oracles, and the human oracle as the fallback. — [PR] Barr, E.T., Harman, M., McMinn, P., Shahbaz, M., Yoo, S. (2015). The oracle problem in software testing: A survey. *IEEE TSE* 41(5):507–525. [doi:10.1109/TSE.2014.2372785](https://doi.org/10.1109/TSE.2014.2372785); OA at [UCL Discovery](https://discovery.ucl.ac.uk/1471263/)
- **Metamorphic testing for scientific software:** Kanewala & Chen argue that MT suits scientists because it tests programs "without oracles" and without extra tooling, and that it is "not very well utilized by scientific software developers." — [PR] Kanewala, U., Chen, T.Y. (2019). Metamorphic testing: A simple yet effective approach for testing scientific software. *Computing in Science & Engineering* 21(1):66–72. [doi:10.1109/MCSE.2018.2875368](https://doi.org/10.1109/MCSE.2018.2875368)
- A 2021 follow-up case study (Storm Water Management Model) found that MT adoption needs explicit MR formulation, scientific-domain examples, correlation with regression testing, and CI integration. — [PR] Peng, Kanewala, Niu (ESEM 2021, emerging results), [listing](https://conf.researchr.org/details/esem-2021/esem-2021-emerging-results-and-vision-papers/4/Contextual-Understanding-and-Improvement-of-Metamorphic-Testing-in-Scientific-Softwar)
- MT reviews (mem): Chen, T.Y. et al. (2018). Metamorphic testing: A review of challenges and opportunities. *ACM Computing Surveys* 51(1):4, [doi:10.1145/3143561](https://doi.org/10.1145/3143561); Segura, S., Fraser, G., Sanchez, A.B., Ruiz-Cortés, A. (2016). A survey on metamorphic testing. *IEEE TSE* 42(9):805–824, [doi:10.1109/TSE.2016.2532875](https://doi.org/10.1109/TSE.2016.2532875). — [PR] (mem)
- **V&V canon (mem):** Oberkampf & Roy distinguish code verification (exact/manufactured solutions, order-of-accuracy tests), solution verification (discretization and iterative error estimation), and validation. — [Book] Oberkampf, W.L., Roy, C.J. (2010). *Verification and Validation in Scientific Computing*. Cambridge University Press. [doi:10.1017/CBO9780511760396](https://doi.org/10.1017/CBO9780511760396) (mem)
- **Method of manufactured solutions (mem):** choose an analytical solution, derive the matching source term, and check the observed convergence order against the theoretical one. — [PR] Roache, P.J. (2002). Code verification by the method of manufactured solutions. *J. Fluids Eng.* 124(1):4–10, [doi:10.1115/1.1436090](https://doi.org/10.1115/1.1436090); [Report] Salari, K., Knupp, P. (2000). *Code Verification by the Method of Manufactured Solutions*. Sandia SAND2000-1444, [doi:10.2172/759450](https://doi.org/10.2172/759450) (mem)
- **Historical evidence that independent scientific codes disagree (mem, older work, 1994):** Hatton & Roberts found that nominally equivalent seismic-processing packages agreed to only one or two significant figures. — [PR] Hatton, L., Roberts, A. (1994). How accurate is scientific software? *IEEE TSE* 20(10):785–797, [doi:10.1109/32.328993](https://doi.org/10.1109/32.328993) (mem)
- **Analytical/pseudo-oracle for small-strain homogenization:** a widely used, independent periodic-homogenization implementation exists for linear elasticity. — [PR] Andreassen, E., Andreasen, C.S. (2014). How to determine composite material properties using numerical homogenization. *Computational Materials Science* 83:488–495, [doi:10.1016/j.commatsci.2013.09.006](https://doi.org/10.1016/j.commatsci.2013.09.006) (mem)
- **[repo] Facts that shape the MRs:** the Newton loop stops on `R_Norm / Rep_MU > TOL` with `Rep_MU = abs(maximum(F_M - I)) * Lx` (Full_Domain/NEMOpt.jl:375,388; Domain_Reduction/NEMOpt-DR.jl:592,606). The residual is normalized by a displacement-like scale, not by a stiffness, so the number of Newton iterations depends on the magnitude of E. The Helmholtz filter has its own periodic constraint matrix `D1_Filter` (Full_Domain/NEMOpt.jl:215). The full-domain design vector has `N_Elements/4` entries (Full_Domain/NEMOpt.jl:124), so symmetry is imposed through the design map.

### Inferences
- **Analytical patch test (strongest oracle).** For a uniform solid design (ρ = 1 everywhere), the exact periodic solution is the homogeneous deformation F_M with zero fluctuation. The homogenized P and A must then equal the compressible neo-Hookean first Piola–Kirchhoff stress and first elasticity tensor at F_M, which can be computed in closed form. Expect agreement to a few ulps × condition growth. The same holds for uniform ρ < 1, with moduli scaled by the SIMP law ε + (1-ε)ρ^p (or whatever interpolation the code uses; confirm in `constitutiveEqs`).
- **Small-strain limit MR.** As F_M → I (e.g., F_M = I + 10⁻⁶ H), the homogenized tangent of any design must approach the small-strain homogenized tensor of the same design with Lamé constants λ, μ derived from (E, ν). The Andreassen & Andreasen (2014) code is an independent pseudo-oracle; the expected agreement is O(‖F_M − I‖) plus the discretization differences between the two codes.
- **Power-of-two scaling MR (exact in binary floating point).** Multiplying E by 2 should multiply A and P by exactly 2 and leave the displacement field and ν₁₂ bit-identical. The one exception is the Newton stopping test, because the residual normalization (Rep_MU) is not stiffness-scaled ([repo] above). That makes this MR useful in two ways: it checks linearity in the moduli, and it exposes absolute-tolerance logic. Power-of-two scaling is exact barring overflow or underflow (standard IEEE property; see Goldberg 1991 / Higham 2002 under Q3).
- **Periodic translation MR.** Cyclically shifting the design by whole elements in x or y must leave A, P, ν₁₂, K, G unchanged up to round-off, provided both the FE problem and the Helmholtz filter are periodic (the filter's periodic map `D1_Filter` suggests it is). This checks the periodic boundary-condition and filter code paths.
- **Rotation/reflection MRs.** Rotating a design by 90° should swap the 1111 and 2222 components of A (and the matching P components) and leave ν₁₂ as defined in `Macroscopic_Quantities` unchanged when F_M is equibiaxial (F_M = 1.025 I in the driver). Reflections should leave the orthotropic components unchanged.
- **Full-domain vs. domain-reduction MR.** For designs that have the symmetry the domain-reduction formulation assumes, Full_Domain and Domain_Reduction should give the same A, P, and gradients (mapped) to within solver round-off. This is a cross-implementation pseudo-oracle that already exists inside the repository.
- **Order-of-accuracy check instead of a full MMS.** NEMOpt has no body-force term, so a classical MMS needs an extra source-term hook. A cheaper code-verification test is the observed convergence order of A for a fixed smooth design under mesh refinement (Roache-style grid convergence). The refactored and legacy codes should show the same observed order and the same extrapolated limit.
- **Use MRs as a complement to the golden master, not a replacement.** A golden-master test only shows "same as before", including any old bugs. The MRs and analytical oracles show "correct", and they stay valid after deliberate numerical changes, such as moving from LU to Cholesky, when bitwise golden masters stop matching.

### Gaps
- No published metamorphic-relation catalog was found specifically for nonlinear periodic homogenization or topology-optimization codes. The MRs above come from mechanics, not from a cited MR study.
- Oberkampf & Roy, Roache, and Salari & Knupp were not re-fetched in this session (mem). Chapter and page pointers for specific recommendations are not given.

---

## Q2. Characterization / golden-master / approval testing, Fowler-style refactoring, and incremental (strangler-fig) migration, adapted to floating-point outputs

### Takeaway
Feathers' "characterization tests" record what the code *actually* does, bugs included, before any change. Refactoring then proceeds in small, behavior-preserving steps with those tests as a safety net. For numerical codes, the golden master must be stored at full precision and compared with explicit, justified tolerances or, where possible, bitwise. Any behavior change must be a separate, deliberate, recorded step ("two hats": refactor, *or* change behavior, never both in one commit). NEMOpt contains several legacy behaviors that a naive refactor would silently "fix": a rounded objective, a rounded value inside the gradient, iteration-count-driven continuation, and dead option updates. Each of these needs a legacy switch.

### Cited Findings
- **Characterization test (definition and purpose):** Feathers coined the term for tests that characterize the *actual* behavior of existing code; they are rerun after changes to catch unintended behavior changes. He defines legacy code as code without automated regression tests. — [Book] Feathers, M. (2004). *Working Effectively with Legacy Code*. Prentice Hall, ISBN 978-0-13-117705-5 (mem for ISBN); definition as summarized by [InfoQ (2007)](https://www.infoq.com/news/2007/03/characterization-testing) and [DaedTech](https://daedtech.com/characterization-tests/) [Grey]
- The technique is also called **golden-master** or snapshot testing. Because it locks in current behavior, it can preserve bugs; the team then decides explicitly whether to keep a behavior or change it. — [Grey] [DaedTech](https://daedtech.com/characterization-tests/); [Learn Agile Practices](https://learnagilepractices.substack.com/p/characterization-tests-for-legacy)
- **Refactoring (mem):** Fowler defines refactoring as a change to internal structure that does not change observable behavior, done in small steps with tests run after each. He separates refactoring from adding functionality (Kent Beck's "two hats"). — [Book] Fowler, M. (2018). *Refactoring: Improving the Design of Existing Code*, 2nd ed. Addison-Wesley, ISBN 978-0-13-475759-9 (mem)
- **Strangler fig (mem):** incrementally replace parts of a legacy system behind a stable interface, so old and new coexist and can be compared until the old part is retired. — [Grey] Fowler, M. "StranglerFigApplication" (2004, updated 2024), [martinfowler.com](https://martinfowler.com/bliki/StranglerFigApplication.html) (mem)
- **Julia golden-master tooling:** ReferenceTests.jl's `@test_reference` compares a value with a stored reference file, creates missing references on first run, and offers an interactive update on mismatch. In non-interactive `pkg> test`, updates are forced only with `JULIA_REFERENCETESTS_UPDATE=true` or by deleting references. It supports `.sha256` hash references for low-storage exact checks. The README page did not document tolerance-based comparison. — [Docs] [ReferenceTests.jl README](https://github.com/JuliaTesting/ReferenceTests.jl)
- **Julia `isapprox` semantics (critical for tolerance tests):** `isapprox(x, y; atol=0, rtol=atol>0 ? 0 : √eps, nans=false[, norm])`, true if `norm(x-y) <= max(atol, rtol*max(norm(x), norm(y)))`. The default `rtol` is √eps (≈1.5e-8 for Float64), and arrays use `LinearAlgebra.norm` (2-norm) by default. — [Docs] [Julia Base.isapprox](https://docs.julialang.org/en/v1/base/math/)
- **[repo] Legacy behaviors that a refactor would otherwise change silently:**
  - The objective returned to MMA is rounded: `J, Vf = round(J, digits=3), …` then `return J` (Full_Domain/NEMOpt_main.jl:79,104; Domain_Reduction/NEMOpt-DR_main.jl:78,104).
  - `v_12` is rounded to 3 digits (NEMOpt_main.jl:80), and the gradient then uses this global `v_12` in `∇J = 200 * (v_12 - v_0) * ∇v_12` (NEMOpt_main.jl:117). The gradient is therefore not the exact derivative of the unrounded J.
  - Continuation is driven by the objective-call counter `iter += 1` (NEMOpt_main.jl:45), with `if iter % 5 == 0 && p < 3` (l. 90). Full_Domain raises β every 10 calls (l. 97) and Domain_Reduction every 5 (NEMOpt-DR_main.jl:97). With `maxinner=40` (l. 237), MMA inner iterations that call the objective can advance the schedule.
  - `options = MMAOptions(s_init=0.5/(β+1), …)` inside `Objective` (l. 100) assigns a function-local variable. The `options` passed to `optimize` at l. 236–241 is not changed, so this update appears to have no effect. This follows standard Julia scoping; confirm by instrumenting.
  - Reference outputs are written already rounded: `writedlm(…, round.(A, digits=4))` (l. 85–86). Existing output files are therefore too coarse to serve as golden masters for 1e-12-level comparisons.
  - The random initial design path `rand(Uniform(0.0, 1.0), …)` (Full_Domain/NEMOpt.jl:124) runs with `RANDOM_SEED` commented out (NEMOpt_main.jl:28). The published runs use PNG initial designs (`INITIAL_DESIGN_NAME = "ITER_1.png"`).

### Inferences
- **Step 0: capture before touching anything.** Pin the legacy environment (Manifest.toml, Julia version, BLAS vendor and thread count). Add non-invasive instrumentation that dumps *full-precision* (e.g., JLD2 Float64 or hex) values of ρ_design, ρ_filtered, U, χ₁..χ₄, P, A, unrounded J, ν₁₂, K, G, the gradient, the Newton iteration counts and residual histories, and the MMA iterates for every objective call of each published case. These are the characterization data. The existing 4-digit text outputs are not enough.
- **Two-mode refactor.** Keep a `legacy=true` mode in the refactored code that reproduces each legacy behavior listed above: the rounded J returned to MMA, the rounded ν₁₂ in the gradient, call-count continuation, and the same linear-solver path and operation order where feasible. Pure refactoring commits must match the legacy run *bitwise* in this mode (hash comparison, e.g., ReferenceTests `.sha256`). Each deliberate change (removing rounding, outer-iteration-based continuation, symmetric assembly with Cholesky, threading) then becomes its own commit with its own decision record (Q7) and its own tolerance-based evidence (Q3–Q5).
- **Strangler-fig ordering for NEMOpt.** Replace one kernel at a time behind the existing function signatures, such as `Helmholtz_Filter`, `Assembly`, `solve`, `χ_Solver`, `Homogenizer_VolumeIntegral`, and `Stiffness_Sensitivity_*`. Run old and new side by side on the same inputs (a "dual-run" pseudo-oracle) before deleting the old path. The DR/FD pair offers a second, built-in cross-check.
- **Do not use default `isapprox` blindly.** The default rtol √eps ≈ 1.5e-8 is far looser than what a pure refactor should achieve (≈1e-14 to 1e-12) and possibly tighter than what a solver change on an ill-conditioned K can achieve. Its norm-wise test can also hide large relative errors in small components, such as near-zero shear entries of A. State `rtol`, `atol`, and the norm explicitly per quantity, using the scales from Q3.
- **Treat the legacy rounding as a separate concern.** The 3-decimal quantization of J makes the objective MMA sees piecewise constant at the 1e-3 level. A refactor that returns unrounded J will change MMA's inner-loop acceptance decisions and so the trajectory, even with bit-identical physics. That is a results-changing decision, not a refactor.

### Gaps
- Feathers 2004 and Fowler 2018 were not accessed directly in this session. The definitions come via secondary summaries and prior knowledge.
- ApprovalTests (Falco) was not researched. No peer-reviewed study was found on golden-master testing specifically for floating-point scientific outputs.
- Whether Nonconvex.jl's MMA87 calls the objective inside inner iterations (and so advances `iter` there) was not verified against the Nonconvex.jl source. This should be confirmed by logging `iter` against MMA's outer iteration count.

---

## Q3. Floating-point reproducibility and nondeterminism: sources, literature, and how to choose tolerances scientifically

### Takeaway
Bitwise reproducibility is achievable only when operation order, libraries, thread counts, and hardware paths are all fixed. Any change of algorithm (LU → Cholesky), summation order, FMA contraction, or BLAS threading legitimately changes the last bits. The National Academies endorse relaxing to "reproducible results within an accepted range of variation" once the numerics are understood. Tolerances should be justified in one of two ways: (i) a priori, from conditioning (forward error ≲ κ·u·c), or (ii) empirically, from the original code's own sensitivity to round-off-sized perturbations. Climate modeling has done (ii) for decades: perturb the trusted code at round-off level, and accept a port/refactor if its divergence from the trusted run is no larger than that.

### Cited Findings
- **NASEM consensus (2019).** For computations "one may expect that the two results be identical (i.e., obtaining a bitwise identical numeric result)". The report notes that "floating-point addition is commutative but not associative", that in parallel runs "the final result will be different depending on the order in which the partial results are gathered", and that there is "a tension between computational performance and strict numerical reproducibility". It also notes that fused multiply-add cannot preserve operation order, and that it "may make sense to relax the requirement of bitwise reproducibility and settle on reproducible results within an accepted range of variation". FINDING 4-3: "Computational reproducibility, within the range of thoughtfully assessed uncertainties, can be expected…". RECOMMENDATION 4-1: convey "clear, specific, and complete information about any computational methods and data products", including the computational environment. — [Report] National Academies of Sciences, Engineering, and Medicine (2019). *Reproducibility and Replicability in Science*, ch. 4. [doi:10.17226/25303](https://doi.org/10.17226/25303); [chapter 4 online](https://www.nationalacademies.org/read/25303/chapter/7). Definition of reproducibility (Conclusion 3-1) in [chapter 3](https://www.nationalacademies.org/read/25303/chapter/6).
- **Port validation by perturbation growth (older work, 1997, still the canonical idea):** Rosinski & Williamson documented the growth of rounding differences in NCAR CCM2 and used it as the basis for validating ports across machines and compilers. — [PR] Rosinski, J.M., Williamson, D.L. (1997). The accumulation of rounding errors and port validation for global atmospheric models. *SIAM J. Sci. Comput.* 18(2):552–564. [doi:10.1137/S1064827594275534](https://doi.org/10.1137/S1064827594275534)
- The NCAR CCM3 FAQ turned this into an operational test (as summarized by search; direct fetch returned HTTP 401). Two simulations on different machines or compilers "should diverge no faster than" two simulations on a trusted machine where one has a minimal initial perturbation, with a suggested perturbation of 1e-14 (PERTLIM). — [Grey/Docs] [NCAR CCM3 FAQ: port verification](https://www2.cgd.ucar.edu/cms/ccm3/ccm-users-faq/port.html) (content via search summary only)
- **Statistical consistency testing when bit-for-bit is lost:** CESM's ensemble consistency test builds a statistical distribution from an ensemble of runs. It then tests whether a new run (new compiler, hardware, or code change) is statistically distinguishable from it, as an objective substitute for bit-for-bit comparison. — [PR] Baker, A.H. et al. (2015). A new ensemble-based consistency test for the Community Earth System Model (pyCECT v1.0). *Geosci. Model Dev.* 8:2829–2840. [doi:10.5194/gmd-8-2829-2015](https://doi.org/10.5194/gmd-8-2829-2015). Follow-up (mem): Milroy, D.J. et al. (2018). Nine time steps: ultra-fast statistical consistency testing of the Community Earth System Model (pyCECT v3.0). *GMD* 11:697–711, [doi:10.5194/gmd-11-697-2018](https://doi.org/10.5194/gmd-11-697-2018); Massonnet, F. et al. (2020). Replicability of the EC-Earth3 Earth system model under a change in computing environment. *GMD* 13:1165–1178, [doi:10.5194/gmd-13-1165-2020](https://doi.org/10.5194/gmd-13-1165-2020) (mem)
- **Random-rounding diagnosis in an industrial FE code:** EDF's Verrou (a Valgrind tool) randomizes floating-point rounding, using CESTAC / Monte Carlo Arithmetic ideas, to find numerically unstable code. Févotte & Lathuilière applied it to code_aster, a structural FE code, to understand "non-reproducibility between test computers". Verrou's localization found the source of the instability. — [PR] Févotte, F., Lathuilière, B. (2017). Studying the numerical quality of an industrial computing code: A case study on code_aster. NSV 2017, LNCS 10381 ([PDF](https://www.cs.ox.ac.uk/conferences/NSV17/NSV_fevotte.pdf); DOI 10.1007/978-3-319-63501-9_5 (mem)); follow-up on confidence intervals for stochastic arithmetic [arXiv:1807.09655](https://arxiv.org/pdf/1807.09655); [FOSDEM 2024 Verrou slides](https://archive.fosdem.org/2024/events/attachments/fosdem-2024-2020-verrou-a-valgrind-tool-dedicated-to-floating-point-error-diagnosis/slides/22471/presentation_ugL7yuB.pdf) [Grey]. Related (mem): Denis, C., de Oliveira Castro, P., Petit, E. (2016). Verificarlo: checking floating point accuracy through Monte Carlo Arithmetic. ARITH-23, [doi:10.1109/ARITH.2016.31](https://doi.org/10.1109/ARITH.2016.31)
- **Reproducible summation / BLAS:** Ahrens, Nguyen & Demmel give algorithms that compute a reproducible floating-point sum independent of summation order. The "binned" accumulator uses one read-only pass and one parallel reduction. ReproBLAS implements BLAS-1, GEMV, and GEMM reproducibly. — [Report] Ahrens, P., Nguyen, H.D., Demmel, J. (2015). *Efficient Reproducible Floating Point Summation and BLAS*. UC Berkeley EECS-2015-229, [PDF](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2015/EECS-2015-229.pdf); [ReproBLAS SCAN'14 slides](https://bebop.cs.berkeley.edu/reproblas/docs/talks/SCAN14.pdf). Journal version (mem; verify): Demmel, J., Nguyen, H.D. (2015). Parallel reproducible summation. *IEEE Trans. Computers* 64(7):2060–2070, doi:10.1109/TC.2014.2345391.
- **Vendor-level reproducibility controls:** Intel oneMKL's Conditional Numerical Reproducibility (CNR) gives "bitwise reproducible results from run to run … on a fixed number of threads for a specific Intel instruction set architecture". A "strict CNR" mode covers a subset of routines even when thread counts vary. Intel cites performance costs of 10–20%, and in some code-path-restricted cases more than half. — [Docs] [Intel oneMKL CNR training note](https://www.intel.com/content/www/us/en/developer/archive/training/conditional-numerical-reproducibility-cnr.html); [oneMKL developer reference: CNR control](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2026-0/conditional-numerical-reproducibility-control.html)
- **Threaded BLAS:** threaded GEMM combines partial sums in a thread-dependent order, so results can differ in the last bits. Setting `OPENBLAS_NUM_THREADS=1` gives bit-identical runs on one machine, but bitwise equality across BLAS vendors, CPU microarchitectures, or library versions is not guaranteed. — [Grey/Docs] [sparho reproducibility notes](https://sparho.readthedocs.io/en/latest/reproducibility.html); thread-dependent OpenBLAS discrepancy reports: [OpenBLAS #1071](https://github.com/OpenMathLib/OpenBLAS/issues/1071), [#2126](https://github.com/xianyi/OpenBLAS/issues/2126)
- **Julia FMA semantics:** `muladd` "may be implemented as an fma"; "The result can be different on different machines and can also be different on the same machine due to constant propagation or other optimizations". `fma` computes x*y+z "without rounding the intermediate result". — [Docs] [Julia Base math docs](https://docs.julialang.org/en/v1/base/math/)
- **Julia sparse `\` is a data-dependent polyalgorithm.** For a square `SparseMatrixCSC`, SparseArrays checks triangularity, then `ishermitian(A)` (an *exact* test), and dispatches Hermitian Float64 matrices to `cholesky(A; check=false)` (CHOLMOD), falling back to `ldlt!` if Cholesky fails. Otherwise it uses `lu(A)` (UMFPACK). — [Docs/source] [SparseArrays.jl src/linalg.jl](https://raw.githubusercontent.com/JuliaSparse/SparseArrays.jl/main/src/linalg.jl). The LinearAlgebra docs note that sparse `lu` uses UMFPACK with a control vector "modified from the UMFPACK defaults to disable iterative refinement", and that `LDLt` "does not use pivoting during the numerical factorization and therefore the procedure can fail even for invertible matrices". — [Docs] [Julia LinearAlgebra stdlib](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- **[repo] Which solver NEMOpt actually uses depends on exact symmetry.** Full_Domain builds `K_Reduced` from `Kmm + Kms*D1 + D1'*Ksm + D1'*Kss*D1` and the corresponding blocks, then solves `K_Reduced \ -R_Reduced` (Full_Domain/NEMOpt.jl:358–367). Element matrices are `B' * A_rho_2D * B` (l. 333). Neither is guaranteed bitwise symmetric in floating point, so the legacy run's choice of UMFPACK or CHOLMOD is data-dependent. Domain_Reduction builds `K_Upper = hcat(Kpm + Kps * D1, Kpi)` (NEMOpt-DR.jl:572, 664, 804). That appears to apply the constraint only on the column side, which suggests a non-symmetric reduced system by construction. This should be confirmed at runtime with `ishermitian(K_Reduced)`.
- **Julia RNG streams are not stable across versions:** "a minor release of Julia (e.g. 1.3 to 1.4) *may change* the sequence of pseudorandom numbers generated from a specific seed". For tests, the docs recommend storing data or using StableRNGs.jl. The default is Xoshiro256++ with per-Task state. — [Docs] [Julia Random stdlib](https://docs.julialang.org/en/v1/stdlib/Random/)
- **Foundations (mem):** Goldberg, D. (1991). What every computer scientist should know about floating-point arithmetic. *ACM Computing Surveys* 23(1):5–48, [doi:10.1145/103162.103163](https://doi.org/10.1145/103162.103163). Higham, N.J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, [doi:10.1137/1.9780898718027](https://doi.org/10.1137/1.9780898718027). Higham covers the rounding-error analysis of summation (error bounds growing with n·u and depending on ordering) and backward-stable linear solvers, where forward error ≲ condition number × backward error. — [PR]/[Book] (mem)
- **Reproducibility in scientific computing (overview):** Bailey, Borwein & Stodden discuss numerical variability among the causes of irreproducibility and summarize remedies, including high-precision arithmetic. — [Book chapter] Bailey, D.H., Borwein, J.M., Stodden, V. (2016). Facilitating reproducibility in scientific computing: principles and practice. In Atmanspacher & Maasen (eds.), *Reproducibility: Principles, Problems, Practices, and Prospects*, Wiley ([preprint](https://carmamaths.org/jon/reprod.pdf); DOI 10.1002/9781118865064.ch9 (mem)). Related workshop report: Stodden et al. (2013) "Setting the default to reproducible", [ICERM report](https://faculty.washington.edu/rjl/pubs/icerm2012/icerm_report.pdf) [Report]
- **Noise-aware tolerances (mem):** Moré & Wild give a method (ECnoise) to estimate the computational noise level of a simulation output from a few nearby evaluations, and use it to choose finite-difference steps. — [PR] Moré, J.J., Wild, S.M. (2011). Estimating computational noise. *SIAM J. Sci. Comput.* 33(3):1292–1314, [doi:10.1137/100786125](https://doi.org/10.1137/100786125); (2012) Estimating derivatives of noisy simulations. *ACM TOMS* 38(3):19, [doi:10.1145/2168773.2168777](https://doi.org/10.1145/2168773.2168777) (mem)

### Inferences
- **Two tolerance regimes for NEMOpt:**
  - **Regime B (bitwise).** Applies to pure refactors in legacy mode with the same Julia/SuiteSparse/BLAS versions, `BLAS.set_num_threads(1)`, no `@fastmath`, and no new `muladd`/`@simd` reassociation in summation loops. Accept only bit equality: hash the full-precision dumps.
  - **Regime T (tolerance).** Applies to deliberate numerical changes: solver, assembly order, threading, quadrature loop order, or library upgrades. The acceptance band must be justified and recorded per quantity.
- **How to set Regime-T bands (recommended order):**
  1. **Measure the legacy noise floor empirically.** Run the *legacy* code on the same input with round-off-level perturbations: ±1 ulp on ρ_design, a different BLAS thread count, forcing `lu` vs. `cholesky` on the same K, and/or Verrou random rounding if feasible under Julia's JIT (not verified). The spread of U, A, J, and ∇J across these runs is the band within which "the original code itself cannot distinguish results". This is the Rosinski–Williamson / CESM-ECT logic applied to a single evaluation.
  2. **Check against an a priori bound.** For linear solves, relative forward error ≲ c·κ(K_Reduced)·u with u ≈ 1.1e-16. Estimate κ with a 1-norm condition estimator. Void elements with ε = 5e-4 make κ large, so expect solver-swap differences well above 1e-14 in U, while integrated quantities such as A are usually less sensitive.
  3. **Account for the Newton stopping tolerance.** The state is converged only to `R_Norm/Rep_MU ≤ 1e-7` ([repo]). A refactor that changes round-off can change the number of Newton iterations near the threshold. Quantities can then jump by up to the size of the last Newton correction (small under quadratic convergence, but not zero). Record Newton iteration counts as a discrete, comparable characteristic, and treat an iteration-count change as an expected Regime-T event, not as a failure by itself.
- **Use relative, component-aware, scaled comparisons.** Compare A component-wise with `atol` scaled to ‖A‖ (some entries are near zero), U with a norm relative to Rep_MU, and ν₁₂ and J absolutely, because they are O(1) and J is itself quantized at 1e-3 in legacy mode. Report max-norm and 2-norm errors and ulp distances (e.g., via the bit patterns of `reinterpret(Int64, x)`) for kernel-level tests.
- **The LU-vs-Cholesky change is already latent in the legacy code.** Because `\` dispatches on an exact `ishermitian` check ([Docs] above), a refactor that assembles K symmetrically (upper triangle only, or wrapped in `Symmetric`) will silently move the solve from UMFPACK to CHOLMOD. If the tangent is ever indefinite, it would fall back to unpivoted LDLᵀ. The legacy solver path should be logged during the characterization run, and the refactored code should choose the factorization explicitly rather than relying on `\`.
- **Seed and RNG discipline.** The published cases use PNG initial designs, so the RNG matters only for random-start runs. Those should use `StableRNGs.jl` or a stored initial vector, because Julia does not guarantee stream stability across minor versions.

### Gaps
- No source was found that quantifies round-off-induced differences for UMFPACK vs. CHOLMOD specifically on FE stiffness matrices with SIMP void regions. The κ-based expectation above is inference.
- Verrou's compatibility with Julia's JIT-compiled code (and with libsuitesparse via Valgrind) was not verified.
- The NCAR CCM3 FAQ could not be fetched directly (HTTP 401). Its content is from a search summary.
- The IEEE TC version of Demmel & Nguyen was not confirmed in this session; the Berkeley technical report was.

---

## Q4. Sensitivity of nonconvex topology-optimization outcomes to small perturbations, and how to report and compare optimized designs

### Takeaway
Topology optimization is nonconvex with many local minima, and penalization/continuation does not guarantee a unique or even black-and-white limit. Round-off-level differences can therefore legitimately produce different final topologies after hundreds of iterations. "Same results" at the full-run level should be defined by (i) a bitwise legacy-mode replay for pure refactors, and (ii) for numerically changed code, statistical and property-based equivalence:
- the objective and homogenized properties of final designs re-analyzed with a common evaluator, after thresholding;
- volume fraction;
- the measure of non-discreteness (M_nd);
- geometric similarity (IoU/Hamming) of thresholded designs;
- the iteration and continuation history.

Divergence should be judged against the spread of the legacy code under round-off-level perturbations.

### Cited Findings
- **Fair comparison is the author's responsibility.** Sigmund (2022) states that "how to fairly evaluate and compare new concepts and ideas to existing ones is an open question". Without practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors". The note lists recommendations, tools, and concepts for good scientific practice. — [PR] Sigmund, O. (2022). On benchmarking and good scientific practise in topology optimization. *Struct. Multidisc. Optim.* 65(11):315. [doi:10.1007/s00158-022-03427-2](https://doi.org/10.1007/s00158-022-03427-2) ([DTU Orbit record](https://orbit.dtu.dk/en/publications/on-benchmarking-and-good-scientific-practise-in-topology-optimiza/))
- Secondary sources citing Sigmund (2022) report two of its practices: converting all solutions to binary before comparing designs, and judging success by whether the threshold designs are close. — [PR, secondary] 2025 CMES benchmarking study of the Zhou–Rozvany problem ([PDF](https://www.sciopen.com/local/article_pdf/10.32604/cmes.2025.065935.pdf); [HTML](https://www.techscience.com/CMES/v143n3/62838/html)); hyperparameter-tuning study noting that "slight differences in parameter settings can lead to the identification of a different local minimum" ([SMO 2024, doi:10.1007/s00158-024-03850-7](https://link.springer.com/article/10.1007/s00158-024-03850-7))
- **Penalization trajectories can be discontinuous.** Stolpe & Svanberg showed on chosen examples that the global trajectory of penalized optima can be discontinuous even when the 0–1 problem has a unique solution. They also show that penalization with continuation can fail to produce a black-and-white design however large the penalty. — [PR] Stolpe, M., Svanberg, K. (2001). On the trajectories of penalization methods for topology optimization. *Struct. Multidisc. Optim.* 21(2):128–139. [doi:10.1007/s001580050177](https://doi.org/10.1007/s001580050177)
- **Classic instabilities (older work, 1998, still canonical; mem):** checkerboards, mesh dependence, and local minima, with continuation as a common but non-guaranteeing remedy for local minima. — [PR] Sigmund, O., Petersson, J. (1998). Numerical instabilities in topology optimization: A survey on procedures dealing with checkerboards, mesh-dependencies and local minima. *Struct. Optim.* 16(1):68–75. [doi:10.1007/BF01214002](https://doi.org/10.1007/BF01214002) (mem)
- **Benchmarking solvers with performance profiles (mem):** Rojas-Labanda & Stolpe compared optimizers (including MMA) on a large set of TO problems using performance profiles, reporting objective quality, iterations/evaluations, and the final design's discreteness. — [PR] Rojas-Labanda, S., Stolpe, M. (2015). Benchmarking optimization solvers for structural topology optimization. *Struct. Multidisc. Optim.* 52(3):527–547. [doi:10.1007/s00158-015-1250-z](https://doi.org/10.1007/s00158-015-1250-z) (mem); performance profiles: Dolan, E.D., Moré, J.J. (2002). *Math. Program.* 91:201–213, [doi:10.1007/s101070100263](https://doi.org/10.1007/s101070100263) (mem)
- **Measure of non-discreteness.** Sigmund (2007), Eq. 41, defines a gray-level indicator: 1 when all densities are 1/2, 0 when all are 0 or 1. It is commonly written M_nd = (1/n)·Σ 4ρ_e(1−ρ_e) × 100%; the explicit formula is mem, consistent with the topoptlab description. — [PR] Sigmund, O. (2007). Morphology-based black and white filters for topology optimization. *Struct. Multidisc. Optim.* 33(4–5):401–424, [doi:10.1007/s00158-006-0087-x](https://doi.org/10.1007/s00158-006-0087-x); implementation note: [topoptlab.design_analysis docs](https://topoptlab.readthedocs.io/en/latest/topoptlab.design_analysis.html) [Docs] (topoptlab is described in [JOSS 10.21105/joss.09105](https://joss.theoj.org/papers/10.21105/joss.09105))
- **Projection/continuation (mem):** robust projection formulations and the effect of β continuation on convergence and final design. — [PR] Wang, F., Lazarov, B.S., Sigmund, O. (2011). On projection methods, convergence and robust formulations in topology optimization. *Struct. Multidisc. Optim.* 43(6):767–784, [doi:10.1007/s00158-010-0602-y](https://doi.org/10.1007/s00158-010-0602-y) (mem)
- **Geometric similarity of designs (IoU).** Recent ML-for-TO work compares thresholded designs with intersection-over-union (|S₁∩S₂|/|S₁∪S₂| of solid-element sets) and Dice, plus MAE on the continuous field before thresholding. It notes that thresholding can change performance by different amounts for different design families. — [Grey/preprint] Trajectory-aware flow matching for topology optimisation (arXiv 2607.14652, 2026), [PDF](https://arxiv.org/pdf/2607.14652); Meta-neural TO ([arXiv 2502.01830](https://arxiv.org/pdf/2502.01830))
- **[repo] Published-case configuration that defines the trajectory:**
  - F_M = 1.025 I, 80×80 elements, E = 100, ν = 0.3, r = 1.0 element, ε = 5e-4.
  - β₀ = 1, ω = 0.5, plus a second projection (β₂ = 200, ω₂ = 0.05).
  - N_ITER = 120, V_MAX = 0.5, MMA87 with kkt tol 1e-8, s_init = 0.25, maxinner = 40, objective J = 100(ν₁₂ − (−1))².
  - Initial design from PNG (Full_Domain/NEMOpt_main.jl:9–30, 104–117, 234–241).

### Inferences
- **Full-run acceptance protocol for NEMOpt (d):**
  1. **Legacy-mode replay.** The refactored code in `legacy=true` mode must reproduce the full trajectory bitwise on the reference machine and environment: every J, Vf, ν₁₂, K, G, design vector, and the iteration count. This is the only rigorous proof of "pure refactor" for a nonconvex continuation run.
  2. **Characterize the legacy code's own chaos.** Run the legacy code N times (e.g., 10–30) for each published case with round-off-level perturbations: ρ₀ ± 1 ulp per element, alternate solver path, or BLAS threads. Record the per-iteration divergence ‖ρᵏ − ρᵏ_ref‖ and the final-design statistics. This gives an empirical distribution of "legacy-equivalent" outcomes, as in Rosinski–Williamson and CESM-ECT (Q3).
  3. **Accept numerically changed code if:** its per-iteration divergence from the reference grows no faster than the legacy ensemble's, and its final-design metrics fall inside the legacy ensemble's range (or a pre-registered band).
  4. **Compare final designs with a common evaluator.** Re-analyze final designs from both codes with one frozen evaluator (the legacy analysis at full precision). Compare J, ν₁₂, K, G, Vf, M_nd, the thresholded design's (ρ > 0.5) ν₁₂, and IoU/Hamming distance between thresholded designs.
  5. **Compare to the article's reported numbers at the article's reported precision.** The 2026 article's values (to the digits printed) are the external acceptance target.
- **Topology differences are not automatically failures.** When perturbed legacy runs themselves reach different local optima (low IoU but similar ν₁₂), a refactored run reaching yet another optimum with ν₁₂ inside the legacy spread is "same result" in the scientific sense. This must be stated explicitly with the ensemble evidence, following Sigmund's (2022) burden-of-fairness principle.
- **Legacy design choices that amplify perturbations.** Several legacy behaviors make bitwise trajectory reproduction fragile once anything changes: the J quantization to 1e-3, the β step every 10 (FD) or 5 (DR) *objective calls*, and the β₂ = 200 second projection. Any change to how many objective calls MMA makes per outer iteration shifts the continuation schedule and so changes the trajectory. A refactor that moves continuation to MMA outer iterations is a results-changing decision (Q7).
- **Report the continuation details exactly** (p schedule, β schedule, trigger counter definition, MMA s_init/s_incr/s_decr, maxinner, stopping rule). Following Stolpe & Svanberg, Sigmund & Petersson, and Wang et al., the final design is a function of the path, not only of the problem.

### Gaps
- The full text of Sigmund (2022) was not accessible (Springer authentication redirect). Its specific recommendations (for example, on M_nd reporting, post-evaluation, and code sharing) are known here only through its abstract and secondary citations. Read it directly before quoting.
- No peer-reviewed TO paper was found that establishes IoU/Hamming as a standard design-similarity metric. Its use here is borrowed from ML-for-TO preprints.
- No study was found that quantifies trajectory divergence from round-off perturbations in density-based TO specifically. The climate-model analogy is an inference.

---

## Q5. Gradient/adjoint verification best practice: finite differences, complex step, Taylor remainder tests, and dot-product tests

### Takeaway
Verify sensitivities in two ways:
- **Independently of the legacy code:** Taylor remainder convergence (second order with a correct gradient) and, where feasible, the complex-step method (no subtractive cancellation, so tiny steps give near machine-precision derivatives).
- **Against the legacy code:** compare the gradient vector at frozen designs using the Q3 tolerance bands.

Also verify adjoint/linearization consistency with the dot-product identity. For NEMOpt, these tests must be run on the *unrounded* objective. The legacy gradient deliberately uses a rounded ν₁₂, so it is not the exact gradient of either the rounded or the unrounded J.

### Cited Findings
- **Taylor remainder test.** The first-order remainder |Ĵ(m + hδm) − Ĵ(m)| should converge at order 1. The second-order remainder |Ĵ(m + hδm) − Ĵ(m) − h∇Ĵ·δm| should converge at order 2 if the gradient is correct. In practice, "repeatedly halving h and checking that the result decreases by a factor of 4". If convergence orders are off, the docs suggest smaller perturbations. Non-differentiable operations or manual modification of assembled values break the test. — [Docs] [dolfin-adjoint verification docs](https://dolfin-adjoint-doc.readthedocs.io/en/latest/documentation/verification.html); primary method paper: [PR] Farrell, P.E., Ham, D.A., Funke, S.W., Rognes, M.E. (2013). Automated derivation of the adjoint of high-level transient finite element programs. *SIAM J. Sci. Comput.* 35(4):C369–C393, [doi:10.1137/120873558](https://doi.org/10.1137/120873558)
- **Complex-step derivative approximation.** Martins, Sturdza & Alonso present the method, improvements to its accuracy and robustness, and automatic implementations for Fortran and C/C++. (Mem: the method has no subtractive cancellation, so f'(x) ≈ Im f(x + ih)/h is accurate to machine precision for very small h, e.g., 1e-20 to 1e-30. Non-analytic functions such as `abs`, `min`/`max`, and comparisons need complex-safe redefinitions.) — [PR] Martins, J.R.R.A., Sturdza, P., Alonso, J.J. (2003). The complex-step derivative approximation. *ACM TOMS* 29(3):245–262, [doi:10.1145/838250.838251](https://doi.org/10.1145/838250.838251) ([MDO Lab record](https://mdolab.engin.umich.edu/bibliography/Martins2003a.html))
- **UMFPACK supports ComplexF64.** Julia's sparse `lu` "only supports sparse matrices with Float64 or ComplexF64 elements", and CHOLMOD supports "real or complex types in single or double precision". A complex-step pass through NEMOpt's sparse solves is therefore technically possible. — [Docs] [Julia LinearAlgebra stdlib](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- **Finite-difference step choice (mem):** forward differences have optimal h ≈ √u·|x| with error ≈ √u (~1e-8). Central differences have optimal h ≈ u^(1/3)·|x| with error ≈ u^(2/3) (~1e-11). When the function carries computational noise ε_f (for example, from iterative solver tolerances), replace u by ε_f. — [Book] Nocedal, J., Wright, S.J. (2006). *Numerical Optimization*, 2nd ed., Springer, §8.1, [doi:10.1007/978-0-387-40065-5](https://doi.org/10.1007/978-0-387-40065-5) (mem); noise-aware steps: Moré & Wild (2011, 2012) cited in Q3 (mem)
- **Dot-product (adjoint consistency) test (mem):** for a linearization F′ and its adjoint, ⟨ȳ, F′ẋ⟩ = ⟨F′ᵀȳ, ẋ⟩ must hold to round-off for random ẋ, ȳ. This is the standard consistency check between tangent-linear and adjoint code. — [Book] Griewank, A., Walther, A. (2008). *Evaluating Derivatives*, 2nd ed., SIAM, [doi:10.1137/1.9780898717761](https://doi.org/10.1137/1.9780898717761) (mem)
- **[repo] Sensitivity code paths:** `Stiffness_Sensitivity_Unfiltered` (Full_Domain/NEMOpt.jl:632) and `Stiffness_Sensitivity_Filtered` (l. 518) compute δA components. The driver combines them as ∇ν₁₂ = 2δA₁₁₂₂/(A₁₁+A₂₂) − 2A₁₂(δA₁₁₁₁+δA₂₂₂₂)/(A₁₁+A₂₂)² and ∇J = 200(ν₁₂ − ν₀)∇ν₁₂, with ν₁₂ already rounded to 3 decimals (NEMOpt_main.jl:80, 111–117).

### Inferences
- **Sensitivity verification protocol for NEMOpt (c):**
  1. **Test targets.** Expose unrounded pure functions: ρ_design ↦ A (all components), ρ_design ↦ ν₁₂, and ρ_design ↦ J_unrounded. Run the tests on these, never on the legacy `Objective`, whose output is quantized at 1e-3 and would make finite differences at h < ~1e-3 meaningless.
  2. **Taylor test.** Run it on a few random directions δm and on element-local directions at several designs: the initial PNG design, a mid-continuation iterate, and a near-binary final iterate. Expect slope ≈ 2. Note that the inner Newton tolerance (1e-7 relative residual) puts a noise floor on J that flattens the remainder curve at small h. Pick h ranges above that floor, or tighten TOL for the test only.
  3. **Complex step.** Feasible because UMFPACK/CHOLMOD accept ComplexF64 and the model functions (tanh projection, polynomial SIMP, neo-Hookean log/det terms) are analytic. Code must avoid `abs`, `real`-only `max`, `round`, and branch-on-value constructs, or redefine them for complex inputs. A complex-step gradient of A can serve as the near-exact reference for both the legacy and refactored adjoints.
  4. **Legacy vs. refactored gradient.** Compare at frozen designs with relative-norm error and cosine similarity, using bands from Q3. In legacy mode the gradient must be bitwise identical.
  5. **Dot-product tests** for each linear operator and its transpose: the Helmholtz filter (K_Filter and its adjoint use), the periodicity maps D1 and D2, and the reduced tangent K_Reduced and its use in the adjoint solves.
- **Rounded ν₁₂ in ∇J is a legacy bug to record, not to keep silently.** It is not the exact derivative of the unrounded objective. It also differs from the derivative of the rounded objective, which is zero almost everywhere. The relative gradient error is |ν₁₂_rounded − ν₁₂|/|ν₁₂ − ν₀|, which is small but nonzero, at most ~5e-4/|ν₁₂+1|. Correcting it is a results-changing decision to record (Q7).

### Gaps
- Martins et al. (2003) was not read in full this session. The recommended step range and the list of non-analytic functions to redefine are from prior knowledge.
- No literature was found on Taylor tests specifically for path-dependent (load-stepped) nonlinear homogenization sensitivities. The noise-floor guidance is inference.

---

## Q6. Reproducibility practices for computational research: rules, environment pinning, seeds/versions, FAIR4RS, and artifact badges

### Takeaway
The consensus practices are: track how every result was produced, version-control everything, archive exact versions of all dependencies, record seeds, keep raw data behind plots, and provide public access. For Julia, that means a committed Project.toml *and* Manifest.toml, a recorded Julia version and BLAS configuration, git-commit-tagged outputs (e.g., DrWatson `tagsave`), and archived releases with DOIs. Reproducibility terms are not standardized; define them explicitly in the article (e.g., ACM/NASEM sense).

### Cited Findings
- **Ten simple rules:**
  1. For every result, keep track of how it was produced.
  2. Avoid manual data manipulation.
  3. Archive exact versions of external programs.
  4. Version-control all custom scripts.
  5. Record intermediate results.
  6. Note underlying random seeds.
  7. Store raw data behind plots.
  8. Generate hierarchical analysis output.
  9. Connect textual statements to underlying results.
  10. Provide public access to scripts, runs, and results.

  — [PR] Sandve, G.K., Nekrutenko, A., Taylor, J., Hovig, E. (2013). Ten simple rules for reproducible computational research. *PLoS Comput. Biol.* 9(10):e1003285, [doi:10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285) ([PMC3812051](https://pmc.ncbi.nlm.nih.gov/articles/PMC3812051/))
- **Reproducibility spectrum (mem):** Peng places publications on a spectrum from "publication only" to "full replication", with code and data, linked and executable, as the reproducible end. — [PR] Peng, R.D. (2011). Reproducible research in computational science. *Science* 334(6060):1226–1227, [doi:10.1126/science.1213847](https://doi.org/10.1126/science.1213847) (mem)
- **Best / good-enough practices (mem):** use version control, automate workflows, write tests, and document environments. — [PR] Wilson, G. et al. (2014). Best practices for scientific computing. *PLoS Biol.* 12(1):e1001745, [doi:10.1371/journal.pbio.1001745](https://doi.org/10.1371/journal.pbio.1001745); Wilson, G. et al. (2017). Good enough practices in scientific computing. *PLoS Comput. Biol.* 13(6):e1005510, [doi:10.1371/journal.pcbi.1005510](https://doi.org/10.1371/journal.pcbi.1005510) (mem)
- **Computational methods policy (mem):** Stodden, V. et al. (2016). Enhancing reproducibility for computational methods. *Science* 354(6317):1240–1241, [doi:10.1126/science.aah6168](https://doi.org/10.1126/science.aah6168). — [PR] (mem)
- **FAIR4RS:** the FAIR Principles were adapted for research software by the FAIR4RS Working Group, convened by ReSA, FORCE11, and RDA. Adaptation was needed because of software's "executability, composite nature, and continuous evolution and versioning". — [PR] Barker, M., Chue Hong, N.P., Katz, D.S., et al. (2022). Introducing the FAIR Principles for research software. *Scientific Data* 9:622, [doi:10.1038/s41597-022-01710-x](https://doi.org/10.1038/s41597-022-01710-x) ([PMC9562067](https://pmc.ncbi.nlm.nih.gov/articles/PMC9562067/)); principles document: Chue Hong, N.P. et al. (2022). FAIR Principles for Research Software (FAIR4RS Principles) v1.0, RDA, [doi:10.15497/RDA00068](https://doi.org/10.15497/RDA00068) (mem for the DOI); [RDA v0.3 draft PDF](https://www.rd-alliance.org/system/files/FAIR4RS_Principles_v0.3_RDA-RFC.pdf)
- **Terminology (mem):** "reproducibility" and "replicability" are used with opposite meanings across fields. The SMO "Replication of results" editorial cites Plesser on this lack of consensus (see Q7). — [PR] Plesser, H.E. (2018). Reproducibility vs. replicability: a brief history of a confused terminology. *Front. Neuroinform.* 11:76, [doi:10.3389/fninf.2017.00076](https://doi.org/10.3389/fninf.2017.00076) (mem); NASEM (2019) definitions in Q3.
- **ACM artifact badges (mem; the page returned HTTP 403):** ACM v1.1 defines Repeatability (same team, same setup), Reproducibility (different team, same setup), and Replicability (different team, different setup). Badges are "Artifacts Available", "Artifacts Evaluated – Functional/Reusable", "Results Reproduced", and "Results Replicated". — [Docs] [ACM Artifact Review and Badging, current version](https://www.acm.org/publications/policies/artifact-review-and-badging-current) (not fetched; mem)
- **Run provenance in Julia:** DrWatson's `tagsave` adds a `:gitcommit` field "to establish reproducibility of results using Git". When the repository is dirty and `storepatch=true`, it adds `:gitpatch` with the diff. `@tagsave` adds a `:script` field with the calling script path and line. `gitdescribe` returns `git describe` output or the commit with a `-dirty` suffix. `produce_or_load` caches results keyed by parameters. — [Docs] [DrWatson.jl saving tools](https://juliadynamics.github.io/DrWatson.jl/stable/save/); paper (mem): Datseris, G. et al. (2020). DrWatson: the perfect sidekick for your scientific inquiries. *JOSS* 5(54):2673, [doi:10.21105/joss.02673](https://doi.org/10.21105/joss.02673)
- **Julia environment pinning (mem):** Pkg's `Manifest.toml` records the exact resolved versions (and git-tree hashes) of all direct and indirect dependencies. `Project.toml` records direct dependencies and compat bounds. — [Docs] [Pkg.jl: Project and Manifest](https://pkgdocs.julialang.org/v1/toml-files/) (mem)
- **[repo]** The README lists dependencies without versions and installs them with `Pkg.add([...])`. No Project.toml or Manifest.toml exists in the repository tree (README.md "Setup"; repository listing, Oct 2026). Outputs go to a timestamped folder (`createWorkspace`, Full_Domain/NEMOpt.jl:17) and record `Parameters.txt` and `Description.txt` (l. 38–39), but no git commit or package versions.

### Inferences
- **First reproducibility action, before any refactor:** reconstruct and commit the environment the published runs used (Julia version, package versions, and OpenBLAS/SuiteSparse versions from `versioninfo()` and `Pkg.status(mode=PKGMODE_MANIFEST)`). Tag the legacy commit (e.g., `v1.0.0-fead2026`) and archive it with a DOI (Q7). Without this, even the legacy code may not reproduce its own published numbers on a newer Julia (RNG streams, SuiteSparse/OpenBLAS updates).
- **Every output folder should carry its provenance:** git describe and dirty patch, Manifest hash, `versioninfo()`, BLAS vendor and thread count, solver path (UMFPACK/CHOLMOD), legacy-mode flags, and full-precision results (JLD2 via `tagsave`). This satisfies Sandve rules 1, 3, 5, and 6 and NASEM Rec. 4-1.
- **Define terms in the new article.** Use "reproduced" for same code and data giving consistent results within a stated tolerance. Use "replicated" for an independent re-implementation (the refactored code arguably sits between the two). Cite NASEM (2019) or ACM v1.1 to avoid the terminology confusion the SMO editorial notes.

### Gaps
- The ACM badging page could not be fetched (HTTP 403). Definitions are from prior knowledge.
- Peng 2011, Wilson 2014/2017, Stodden 2016, Plesser 2018, and the Pkg docs were not re-fetched this session (mem).

---

## Q7. (Scope addition) Recording every refactoring decision so the refactored code can support a NEW topology-optimization research paper

### Takeaway
The next publication is a TO *research* paper that uses the refactored code, not a paper about the code. The decision trail therefore has to serve that paper's methods and reproducibility statements. It needs six things:
- **(a)** Every reported result traces to an exact formulation, parameter set, and archived code release with a DOI.
- **(b)** A citable verification record shows that the refactored code reproduces the FEAD 2026 results, and to what tolerance.
- **(c)** An explicit register lists each deliberate deviation from the 2026 formulation, with rationale and citations.
- **(d)** The target journal's replication and data rules are met: SMO's mandatory "Replication of results" section, the fair-comparison expectations of Sigmund (2022), and Elsevier's data options for FEAD/CMAME.
- **(e)** Run-level provenance (commit, parameters, seeds, environment) is stored with every result.
- **(f)** If computational cost is reported, Hoefler–Belli and Chen–Revels practice is followed.

Lightweight, immutable, numbered decision records (Nygard / MADR), extended with numerics fields, are the vehicle that ties (a)–(f) together.

### Cited Findings
**Decision records: the vehicle**
- **Nygard ADRs (2011).** Five parts: Title, Context, Decision, Status ("proposed", "accepted"; later "deprecated" or "superseded" with a reference to the replacement), and Consequences. "ADRs will be numbered sequentially and monotonically. Numbers will not be reused." "If a decision is reversed, we will keep the old one around, but mark it as superseded." "The whole document should be one or two pages long." Store them "in the project repository under doc/arch/adr-NNN.md". — [Grey, foundational practitioner source] Nygard, M. (15 Nov 2011). Documenting Architecture Decisions. [Cognitect blog](https://www.cognitect.com/blog/2011/11/15/documenting-architecture-decisions)
- **Tyree & Akerman (2005).** Argue that explicitly documenting major architecture decisions makes the architecture transparent and clarifies its rationale. They provide a decision template; its full field list was not verified here; secondary sources describe it as having extensive traceability to requirements, principles, and related artifacts. — [PR] Tyree, J., Akerman, A. (2005). Architecture decisions: demystifying architecture. *IEEE Software* 22(2):19–27, [IEEE Xplore 1407822](https://ieeexplore.ieee.org/document/1407822) (DOI 10.1109/MS.2005.27 (mem))
- **Jansen & Bosch (2005).** Treat software architecture as a set of architectural design decisions, which helps trace decisions and keep conceptual integrity (content via secondary slides). — [PR] Jansen, A., Bosch, J. (2005). Software architecture as a set of architectural design decisions. *WICSA'05*, pp. 109–120, IEEE; citation as listed in [arXiv:1704.04798 references](https://arxiv.org/pdf/1704.04798) (DOI 10.1109/WICSA.2005.61 (mem))
- **MADR 4.0.0 (2024-09-17).** "Markdown Architectural Decision Records". Front matter: status, date, decision-makers, consulted, informed. Body: Title, Context and Problem Statement, Decision Drivers, Considered Options, Decision Outcome, Consequences, *Confirmation*, Pros and Cons of the Options, More Information. MADR says it can record any decision, not only architectural ones. Dual-licensed MIT OR CC0-1.0. — [Docs] [MADR](https://adr.github.io/madr/)
- **Empirical evidence on ADR use.** A GitHub mining study found ADR use growing each year but adoption still low. About 50% of repositories with ADRs contain only 1–5 of them. Sustained use is a team effort over a long period, and most projects use Nygard's template. — [PR] Buchgeher, G., Schöberl, S., Geist, V., Dorninger, B., Haindl, P., Weinreich, R. (2023). Using architecture decision records in open source projects—an MSR study on GitHub. *IEEE Access* 11:63725–63740 ([JKU record](https://research.jku.at/en/publications/using-architecture-decision-records-in-open-source-projects-an-ms/); [DOAJ](https://doaj.org/article/f838f2a538ea48e7a3b15ce99befe612))
- An ICSA 2026 study of 921 repositories and about 5,800 ADRs reports that about 63% were created directly as "accepted", skipping the deliberation ADRs are meant to capture. — [PR, conference listing] [ICSA 2026](https://conf.researchr.org/details/icsa-2026/icsa-2026-papers/34/Architecture-Decision-Records-Adoption-Impact-and-Developer-Engagement-in-Open-Sou) (abstract only)
- An industrial action-research study notes "little empirical evidence and research on how Architecture Decision Records are used in practice". — [PR] [Chalmers research record 538920](https://research.chalmers.se/en/publication/538920)
- **ISO/IEC/IEEE 42010 (mem).** The architecture-description standard (2011; revised 2022) includes recording decisions and their rationale. — [Std] [ISO/IEC/IEEE 42010:2022](https://www.iso.org/standard/74393.html) (mem; clause wording not verified)

**(a) A methods section that pins formulation, parameters, and code version**
- **NASEM Finding 4-2:** "the methods section of a traditional scientific paper is insufficient to convey the necessary information" for complex computational results. **Recommendation 4-1:** "researchers should convey clear, specific, and complete information about any computational methods and data products", covering input data, computational steps, and the computational environment. — [Report] [NASEM 2019, ch. 4](https://www.nationalacademies.org/read/25303/chapter/7), [doi:10.17226/25303](https://doi.org/10.17226/25303)
- **Software citation must be version-specific.** The *Specificity* principle says citations should identify "the exact version used" (version or revision numbers). The *Accessibility* principle says they should give access to the software and its metadata. — [PR] Smith, A.M., Katz, D.S., Niemeyer, K.E., FORCE11 SCWG (2016). Software citation principles. *PeerJ Comput. Sci.* 2:e86, [doi:10.7717/peerj-cs.86](https://doi.org/10.7717/peerj-cs.86)
- **Persistent identifiers and core metadata are essential.** The guide covers citing software one developed oneself. — [PR] Katz, D.S., Chue Hong, N.P., et al. (2021). Recognizing the value of software: a software citation guide. *F1000Research* 9:1257 (v2), [doi:10.12688/f1000research.26932.2](https://doi.org/10.12688/f1000research.26932.2)
- **FAIR4RS:** research software needs principles adapted to its "continuous evolution and versioning". — [PR] Barker, M., Chue Hong, N.P., Katz, D.S., et al. (2022). *Scientific Data* 9:622, [doi:10.1038/s41597-022-01710-x](https://doi.org/10.1038/s41597-022-01710-x)
- **Rules for traceability:** Sandve et al. Rule 1 (track how every result was produced), Rule 3 (archive exact versions of external programs), Rule 9 (connect textual statements to underlying results). — [PR] [Sandve et al. 2013, doi:10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285)
- **Release archiving (mem):**
  - GitHub's Zenodo integration mints a DOI per tagged release plus a concept DOI. — [Docs] [GitHub Docs: Referencing and citing content](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content) (mem)
  - CITATION.cff gives machine-readable citation metadata. — [Docs] [Citation File Format](https://citation-file-format.github.io/) (mem)
  - Keep a Changelog 1.1.0 groups changes under Added/Changed/Deprecated/Removed/Fixed. — [Docs] [keepachangelog.com](https://keepachangelog.com/en/1.1.0/) (mem)
  - Semantic Versioning 2.0.0. — [Docs] [semver.org](https://semver.org/) (mem)

**(b) Citable evidence that the refactored code reproduces FEAD 2026**
- **NASEM Finding 4-3:** "Computational reproducibility, within the range of thoughtfully assessed uncertainties, can be expected". The report accepts relaxing bitwise identity to "reproducible results within an accepted range of variation" once numerical issues are understood. — [Report] [NASEM 2019, ch. 4](https://www.nationalacademies.org/read/25303/chapter/7)
- **Elsevier's data options.** Its research-data guidelines define journal options from A (encourage deposit and citation) through B (encouraged deposit, citation, and linking, or a data availability statement), C (deposit, citation, and linking *required*, unless a statement explains why data cannot be shared), and D (required, without the statement alternative), to E (C plus peer review of the data). Each journal's Guide for Authors states which option it uses. — [Docs] [Elsevier research data guidelines](https://www.elsevier.com/authors/author-resources/research-data/data-guidelines); [Elsevier data statement](https://www.elsevier.com/authors/author-resources/research-data/data-statement)
- **Elsevier's research data policy** names "software, algorithms, protocols, methods" among research data and states that data "should be made available free of charge … wherever possible and with minimal reuse restrictions". — [Docs] [Elsevier research data policy](https://www.elsevier.com/about/policies-and-standards/research-data)
- The verification methods that back a citable tolerance claim (bitwise legacy mode, legacy noise-floor ensembles, metamorphic and analytical oracles, Taylor/complex-step tests) are in Q1–Q5 and Q8 of these notes.

**(c) Deliberate deviations from the 2026 formulation: what exists in the legacy code, and literature to justify changes**
- **[repo] Legacy behaviors that any "cleaned-up" formulation would change:**
  - The objective returned to MMA is rounded to 3 decimals (Full_Domain/NEMOpt_main.jl:79,104; Domain_Reduction/NEMOpt-DR_main.jl:78,104).
  - The gradient uses the rounded ν₁₂ (NEMOpt_main.jl:80,117).
  - p and β continuation is triggered by the objective-call counter `iter` (l. 45, 90–102). β steps every 10 calls in Full_Domain and every 5 in Domain_Reduction (NEMOpt-DR_main.jl:97).
  - An MMA `s_init` update inside `Objective` appears to have no effect, because it assigns a local variable (l. 100).
  - Exported A and P are rounded to 4 digits (l. 85–86).
  - A second steep projection uses β₂ = 200, ω₂ = 0.05 (l. 21).
  - The Newton stopping test is normalized by a displacement scale, `R_Norm/Rep_MU > TOL` (Full_Domain/NEMOpt.jl:375,388).
- **Continuation and penalization are path-defining choices.** Trajectories of penalized optima can be discontinuous, and continuation may not yield 0–1 designs. — [PR] [Stolpe & Svanberg 2001, doi:10.1007/s001580050177](https://doi.org/10.1007/s001580050177). Projection and β-continuation affect convergence and robustness. — [PR] Wang, Lazarov & Sigmund 2011, [doi:10.1007/s00158-010-0602-y](https://doi.org/10.1007/s00158-010-0602-y) (mem). Local minima and continuation. — [PR] Sigmund & Petersson 1998, [doi:10.1007/BF01214002](https://doi.org/10.1007/BF01214002) (mem)
- **MMA reference for describing optimizer settings (mem):** Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. *Int. J. Numer. Meth. Eng.* 24(2):359–373, [doi:10.1002/nme.1620240207](https://doi.org/10.1002/nme.1620240207). — [PR] (mem)
- **Describing the method relative to the prior paper (text recycling).**
  - TRRP Best Practices v1.1 (2021): recycle text "where consistency of language is needed for accurate communication" (mainly methods); check permissions for substantial reuse of published material; disclose recycling to editors (cover letter with citations) and readers (e.g., a footnote); do not misrepresent novelty through duplicate publication or "salami slicing". — [Grey/Report] [TRRP Best Practices for Researchers](https://textrecycling.org/resources/best-practices-for-researchers/)
  - COPE / BioMed Central guidelines: cited overlap in methods may be acceptable; undisclosed overlap, or overlap in results, discussion, or conclusions, is unlikely to be. — [Report] [COPE text recycling guidelines](https://nairobi2021.iamcr.org/publicationethics.org/files/Web_A29298_COPE_Text_Recycling.pdf)

**(d) Journal replication and reporting expectations for TO papers**
- **SMO.** The submission guidelines list a mandatory "Replication of results" section: "It is desirable that papers, especially on new methodology, should help readers reproduce results presented in the paper". The 2019 editorial says the section asks for a description "detailed enough" to replicate results consistent with the paper's claims, notes the terminology debate (citing Plesser 2018), and says almost all SMO papers concern computational algorithms. Published papers fill the section in different ways: code on request, a description of the implementation software, or a claim that the method description suffices. — [PR, editorial] Haftka, R.T., Zhou, M., Queipo, N.V. (2019). Replication of results. *Struct. Multidisc. Optim.* 60(2):405–409, [doi:10.1007/s00158-019-02298-4](https://doi.org/10.1007/s00158-019-02298-4); [SMO submission guidelines](https://link.springer.com/journal/158/submission-guidelines) (wording via search snippet; direct fetch blocked by Springer authentication)
- **Sigmund (2022).** Fair evaluation against existing methods is "an open question". Where common benchmarks are lacking or impractical, "the responsibility of fair evaluation of contributions is up to authors". The note lists recommendations, tools, and concepts for good practice. Secondary sources cite it for converting designs to binary before comparison and checking closeness of thresholded designs. — [PR] Sigmund, O. (2022). *Struct. Multidisc. Optim.* 65(11):315, [doi:10.1007/s00158-022-03427-2](https://doi.org/10.1007/s00158-022-03427-2); secondary: [2025 CMES benchmarking study](https://www.sciopen.com/local/article_pdf/10.32604/cmes.2025.065935.pdf)
- **FEAD and CMAME.** Both are Elsevier journals that follow one of the A–E research-data options above. Their specific option and any code-availability wording could not be retrieved: the ScienceDirect Guide for Authors returned HTTP 403. — [Docs] [Elsevier research data guidelines](https://www.elsevier.com/authors/author-resources/research-data/data-guidelines)
- *(One line, per scope.)* Software-paper venues (JOSS; SoftwareX) review the code itself. JOSS, for example, requires an OSI license, automated tests, and at least six months of public development history for recently public repositories. They are not the target here. — [Docs] [JOSS review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)

**(e) Run-level provenance for every reported result**
- DrWatson `tagsave` adds `:gitcommit` "to establish reproducibility of results using Git". It adds `:gitpatch` (the diff) when the repository is dirty and `storepatch=true`. `@tagsave` adds `:script` (path and line). `gitdescribe` returns tag/commit plus a `-dirty` suffix. `produce_or_load` keys results by parameters. — [Docs] [DrWatson.jl saving tools](https://juliadynamics.github.io/DrWatson.jl/stable/save/); (mem) Datseris et al. (2020) *JOSS* 5(54):2673, [doi:10.21105/joss.02673](https://doi.org/10.21105/joss.02673)
- Julia RNG streams "*may change*" between minor releases. Store data or use StableRNGs.jl for reproducible tests. — [Docs] [Julia Random stdlib](https://docs.julialang.org/en/v1/stdlib/Random/)
- Sandve Rule 6 (record random seeds) and Rule 5 (record intermediate results). — [PR] [Sandve et al. 2013](https://doi.org/10.1371/journal.pcbi.1003285)
- Optional machine-readable provenance: W3C PROV (entities, activities, agents). — [Std] [W3C PROV Overview](https://www.w3.org/TR/prov-overview/) (mem). Lab-notebook practice: Schnell (2015) *PLoS Comput. Biol.* 11(9):e1004385, [doi:10.1371/journal.pcbi.1004385](https://doi.org/10.1371/journal.pcbi.1004385) (mem)
- **[repo]** Current outputs record `Parameters.txt` and `Description.txt` in a timestamped folder (Full_Domain/NEMOpt.jl:17,38–39) but no git commit, package versions, BLAS configuration, or solver path. The repository has no Project.toml or Manifest.toml.

**(f) Performance reporting, only if computational cost is reported**
- **Hoefler & Belli (SC'15), twelve rules** (verified from the paper):
  1. Report whether the speedup base case is one parallel process or the best serial execution, plus its absolute performance; "never report ratios without absolute values".
  2. Justify reporting subsets; "report all results, not just the best".
  3. Arithmetic mean only for costs, harmonic mean for rates.
  4. Avoid summarizing ratios; if unavoidable, use the geometric mean.
  5. Report whether measurements are deterministic, with confidence intervals otherwise.
  6. Do not assume normality without diagnostics.
  7. Use statistically sound comparisons (e.g., non-overlapping CIs, ANOVA).
  8. Check whether mean or median is meaningful.
  9. Document all varying factors and the complete setup (software, hardware, techniques).
  10. Report measurement, synchronization, and summarization for parallel timings.
  11. Show upper performance bounds where possible.
  12. Plot enough to interpret results.

  — [PR] Hoefler, T., Belli, R. (2015). *SC'15*, Art. 73, [doi:10.1145/2807591.2807644](https://doi.org/10.1145/2807591.2807644) ([PDF](https://htor.inf.ethz.ch/publications/img/hoefler-scientific-benchmarking.pdf))
- **Chen & Revels (2016).** The minimum is the most robust run-time estimator under timer error and OS jitter, because "variations never reduce the run time". BenchmarkTools.jl automates the evaluations-per-sample choice. A 2026 critique notes the method assumes memoryless, sparse noise. — [Preprint] [arXiv:1608.04295](https://arxiv.org/pdf/1608.04295); critique [arXiv:2606.17261](https://arxiv.org/pdf/2606.17261)

### Inferences
- **What the new research paper should contain, item by item:**
  - **(a) Methods.** A short "Computational implementation" paragraph that:
    - cites the FEAD 2026 article for the base method;
    - cites the refactored code by release DOI (e.g., NEMOpt vX.Y.Z, Zenodo DOI) and the legacy release that produced FEAD 2026 (tagged and archived retroactively, e.g., `v1.0.0-fead2026`);
    - gives a parameter table per result: F_M, mesh, E, ν, r, ε, SIMP law, p/β schedules and their trigger (per MMA outer iteration or per objective call), β₂/ω₂, MMA variant and settings, stopping rule, Newton TOL and load step, and initial design (file hash or seed).
  - **(b) Verification statement.** One paragraph, plus a supplementary table or archived verification report with its own DOI, stating:
    - "In legacy-compatibility mode, NEMOpt vX reproduces the FEAD 2026 runs bit-for-bit on [environment]";
    - "In default mode, single-evaluation quantities agree to [band], gradients to [band], and full runs fall within the legacy round-off ensemble on [metrics]".

    Cite NASEM Finding 4-3 for the tolerance-based definition. This is citable evidence that later readers and reviewers can check.
  - **(c) Deviation register.** A table "Differences from Güven & Özdemir (2026)" with one row per results-changing decision record. Columns: what changed, why (with citations), legacy switch name, and effect on the FEAD 2026 benchmark (re-run value vs. published value). Candidate rows, from the [repo] facts:
    - remove the 3-decimal rounding of J returned to MMA;
    - use the unrounded ν₁₂ in ∇J;
    - trigger continuation per MMA outer iteration instead of per objective call;
    - make the β-dependent MMA `s_init` update actually take effect, or delete it;
    - any SIMP/projection modification;
    - explicit solver choice (LU vs. Cholesky);
    - Newton tolerance normalization.
  - **(d) Journal fit.** For SMO, the "Replication of results" section can point to the release DOI, the parameter table, the deviation register, and the verification report. For FEAD or CMAME, choose a data statement consistent with the journal's Elsevier option (C or stricter means deposit and link). Deposit final designs, full-precision histories, and run scripts with the code DOI. Following Sigmund (2022), compare against the FEAD 2026 designs fairly: same mesh, volume fraction, and evaluator; thresholded designs re-analyzed; M_nd reported.
  - **(e) Provenance.** Every figure and table value is generated by a script from `tagsave`d result files. Each file holds: git describe, dirty patch, Manifest.toml hash, `versioninfo()` (Julia, OS, CPU), BLAS vendor and thread count, SuiteSparse version, solver path taken, legacy-mode flags, parameter dict, seed or initial-design hash, wall time, and full-precision outputs. The paper's supplementary material maps each figure/table to result-file IDs (Sandve Rule 9).
  - **(f) Cost (only if reported).** Report absolute legacy and refactored times on the same machine and environment (Rule 1). Separate per-evaluation cost (implementation speedup at a fixed number of objective/gradient calls) from time to solution (which includes algorithmic changes to iteration counts). Report BenchmarkTools minima and distributions over repeated runs (Chen & Revels; Rules 5–7) and the full setup (Rule 9). Do not average speedup ratios arithmetically (Rule 4).
- **Recommended minimal decision-record template for NEMOpt** (`docs/decisions/NNNN-short-title.md`): Nygard numbering and immutability, MADR's "Considered Options" and "Confirmation", plus numerics and paper-traceability fields:

  ```markdown
  ---
  id: NDR-0007                       # monotonic, never reused (Nygard)
  title: Trigger p/β continuation per MMA outer iteration
  date: 2026-11-03
  status: proposed | accepted | superseded-by NDR-00xx | rejected
  decision-makers: [M. Güven, ...]
  supersedes: []
  related: [NDR-0003]
  ---
  ## Context
  Legacy behavior with file:line (e.g., Full_Domain/NEMOpt_main.jl:45,90–102: schedule
  advances per objective call, so MMA inner iterations shift it); decision drivers
  (correctness, fidelity to FEAD 2026, speed, maintainability).
  ## Considered options
  1. Keep per-call trigger (legacy)  2. Per outer iteration  3. Adaptive (e.g., on KKT/change)
  Pros/cons per option.
  ## Decision
  Chosen option; scope (Full_Domain / Domain_Reduction / both).
  ## Rationale and citations
  Primary sources with DOIs (e.g., Stolpe & Svanberg 2001; Wang et al. 2011).
  ## Effect on numerical results
  affects-results: yes | no | bitwise-only
  deviation-from-FEAD2026-formulation: yes | no   # if yes → row in paper's deviation table
  legacy-switch: `opts.continuation_trigger = :objective_call`
  expected/observed magnitude: e.g., trajectory diverges from iter k; final ν12 within legacy ensemble [a, b]
  ## Verification evidence (Confirmation)
  - Test IDs: K-012 (kernel), E-004 (single evaluation), S-003 (Taylor slope), R-002 (full-run ensemble)
  - Tolerance bands and their derivation (noise-floor run IDs, κ estimate)
  - FEAD 2026 Table/Figure X: reproduced to N digits | changed (value, reason)
  ## Performance evidence (only if cost will be reported)
  Hardware, OS, Julia/BLAS/SuiteSparse versions, threads; BenchmarkTools min & distribution;
  absolute legacy and new times (Hoefler–Belli rules 1, 5, 9)
  ## Paper traceability
  New-paper section/figure/table affected; result-file IDs; release tag + Zenodo DOI
  ## Consequences
  Positive/negative, follow-ups, risks.
  ## Links
  commit(s)/PR, CHANGELOG entry
  ```

- **Process.** Following Buchgeher et al. and the ICSA 2026 finding, open each record as "proposed", with the alternatives and the pre-registered acceptance criteria, *before* running verification. Then accept it with the evidence attached. Write one umbrella record per batch of pure refactors that were verified bitwise in legacy mode, listing the commits. Each results-changing record becomes one row of the paper's deviation table.

### Gaps
- FEAD's and CMAME's specific Elsevier data option and code-availability wording could not be retrieved (Guide for Authors returned HTTP 403). Check before submission.
- The full text of Sigmund (2022) was not accessible. Its exact recommendations (e.g., on reporting M_nd, post-evaluation, and code sharing) should be read directly before citing specifics.
- The current SMO guideline wording was seen only via a search snippet (mandatory section with "desirable" phrasing); the 2019 editorial's full text was not read.
- The copyright/license terms of the FEAD 2026 article, which govern figure and text reuse in the new paper, are unknown.
- No empirical study was found on decision documentation in computational-science or TO research codes specifically. Tyree & Akerman's template fields and ISO 42010 clause wording were not verified (paywalled; mem).
- Zenodo/GitHub, CITATION.cff, Keep a Changelog, SemVer, W3C PROV, and Svanberg 1987 details are from prior knowledge (mem).

---

## Q8. Synthesis: an operational definition of "same results" for NEMOpt at levels (a)–(d)

### Takeaway
Define "same results" as a hierarchy of increasingly statistical acceptance criteria, each justified by a cited principle:
- **(a) Kernels:** bitwise in legacy mode, otherwise ulp/condition-bounded.
- **(b) Single evaluations:** bitwise in legacy mode, otherwise within the legacy code's own measured round-off noise band and κ-based bounds, plus analytical and metamorphic oracles.
- **(c) Sensitivities:** bitwise vs. legacy in legacy mode, *and* independently verified by Taylor/complex-step/dot-product tests.
- **(d) Full runs:** bitwise trajectory replay in legacy mode; otherwise statistical equivalence to an ensemble of round-off-perturbed legacy runs, compared via re-evaluated properties, M_nd, and thresholded-design similarity, and checked against the FEAD 2026 published values at their printed precision.

### Cited Findings
- Bitwise equality is the natural expectation for identical computations; relaxation to "an accepted range of variation" is legitimate once numerical issues are understood. — [Report] [NASEM 2019, ch. 4](https://www.nationalacademies.org/read/25303/chapter/7)
- Port and refactor acceptance by comparing divergence against round-off-perturbed trusted runs. — [PR] [Rosinski & Williamson 1997](https://doi.org/10.1137/S1064827594275534); [PR] [Baker et al. 2015](https://doi.org/10.5194/gmd-8-2829-2015)
- Oracles for scientific software: pseudo-oracles, analytical solutions, metamorphic relations, regression with tolerances. — [PR] [Kanewala & Bieman 2014](https://doi.org/10.1016/j.infsof.2014.05.006); [PR] [Kanewala & Chen 2019](https://doi.org/10.1109/MCSE.2018.2875368)
- Characterization tests before change; behavior-preserving small steps. — [Book] Feathers 2004; Fowler 2018 (mem; see Q2)
- Gradient verification by Taylor remainder (order 2) and complex step. — [PR] [Farrell et al. 2013](https://doi.org/10.1137/120873558); [Docs] [dolfin-adjoint](https://dolfin-adjoint-doc.readthedocs.io/en/latest/documentation/verification.html); [PR] [Martins et al. 2003](https://doi.org/10.1145/838250.838251)
- TO outcomes are path-dependent and non-unique; fair comparison is the author's responsibility; compare thresholded designs. — [PR] [Stolpe & Svanberg 2001](https://doi.org/10.1007/s001580050177); [PR] [Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2); [PR] [Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x)
- Sparse `\` in Julia silently picks Cholesky vs. LU from an exact symmetry check. — [Docs] [SparseArrays.jl source](https://raw.githubusercontent.com/JuliaSparse/SparseArrays.jl/main/src/linalg.jl)

### Inferences
- **Level (a), kernels.** Element routines (`constitutiveEqs`, element K and F_int, filter element matrices), assembly, periodic reduction (D1, D2), the projection/SIMP maps, and `Macroscopic_Quantities`.
  - Legacy mode: bitwise equality on random and stored inputs.
  - New implementation (e.g., in-place, static arrays, different loop order): max ulp distance ≤ a small multiple of the summation length per entry (Higham-style n·u bound).
  - Plus property tests: symmetry of element K to ulp level, zero-energy rigid modes, objectivity of the neo-Hookean stress (rotation invariance).
- **Level (b), single evaluations.** Use a frozen design set: each published case's initial PNG design, every 10th iterate from the legacy trajectory, the final design, a uniform-solid design, a laminate, and a random design. For each design, compare U, χ₁..₄, P, A, unrounded J, ν₁₂, K, G, Vf, and Newton iteration counts. Pass if any of:
  - (i) bitwise equality in legacy mode;
  - (ii) within the legacy noise band measured by solver-swap, thread, or ulp perturbations;
  - (iii) within the κ-based a priori bound.

  In addition, the analytical patch test, power-of-two scaling, periodic-translation, rotation, and FD↔DR metamorphic relations must hold.
- **Level (c), sensitivities.**
  - Legacy-mode gradient is bitwise equal to the legacy gradient.
  - Default-mode gradient agrees within the band from (b), propagated.
  - Taylor second-order slope ≈ 2 over at least two decades of h above the Newton-noise floor.
  - Complex-step agreement ~1e-10 relative or better where feasible.
  - Dot-product identities hold to ~1e-13 relative.
- **Level (d), full runs.** Use the five-step protocol in Q4. Pre-register (in the decision record) the ensemble size, perturbation type, metrics (J, ν₁₂, K, G, Vf, M_nd, thresholded ν₁₂, IoU), and pass thresholds *before* running the comparison, to avoid post hoc tolerance tuning.
- **Separate the published-article target from the legacy-code target.** First prove that the legacy code in its pinned environment reproduces the FEAD 2026 numbers. If it does not, because the published runs used an unrecorded environment, that discrepancy must itself be documented as a decision record. Only then can the refactored code be held to "same as legacy".

### Gaps
- The FEAD 2026 article's own reported values (tables and figures, number of printed digits, run environment) were not examined in this note. Another workstream ("Establish reproduction baseline from article and code") is responsible for that.
- The numeric pass thresholds above (e.g., 1e-10 for complex step, 1e-13 for dot products) are typical values from experience, not cited standards. They should be replaced with measured values from the legacy noise-floor runs.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. The items below update it after open-access full texts were sought on 2026-10-12. Supporting note: [open_access_optimization_papers.md](../NEMOpt%20supplied%20papers%20verification/open_access_optimization_papers.md), called "optimization note" below. Line numbers refer to this file as of 2026-10-12.

- **L160 and L395 (Gaps: Sigmund 2022 full text not accessible).** Append: "No legal open copy exists as of 2026-10-12 (Springer paywall; DTU Orbit record without file; no arXiv)." The quoted abstract sentence is confirmed. It is conditional: "In the lack or impracticality of common benchmarks". The thresholding and binary-comparison practices (L135) remain secondary. Source: optimization note, §3 edit 4 and claims S-a and S-e.
- **L296 (SMO "Replication of results"; "(wording via search snippet; direct fetch blocked by Springer authentication)").** Corrected parenthetical: "(guideline wording via search snippet, unverified; the editorial was read in full from the publisher's free PDF on 2026-10-12 and does not contain this sentence)". The sentence meant is the quoted "It is desirable that papers …".
  - The editorial confirms the "detailed enough" description, the Plesser 2018 terminology point and the remark that almost all SMO papers concern computational algorithms (p. 405).
  - It prescribes no fixed content. It recommends source code with a tagged version and DOI, all tuning and convergence parameters, plot data, intermediate results and repeated runs to show scatter (pp. 406–409).
  - Source: optimization note, §4 edit 4 and claims H-c and H-e.
- **L396 (Gaps: "the 2019 editorial's full text was not read").** Resolved: the editorial (Haftka, Zhou & Queipo 2019, SMO 60:405–409) was read in full from the publisher's free-access PDF. Only the current guideline wording remains snippet-based. Source: optimization note, §4.
