# Julia architecture, tooling and ecosystem choices for refactoring NEMOpt (nonlinear FE topology optimization of periodic microstructures)

Scope: target architecture and toolchain for refactoring NEMOpt, the code behind Güven & Özdemir, "Topology optimization of nonlinear elastic micro-structures by a domain reduction technique," *Finite Elements in Analysis and Design* 260 (2026) 104599, doi:[10.1016/j.finel.2026.104599](https://doi.org/10.1016/j.finel.2026.104599). Metadata checked via Crossref (published Sep 2026).

Status date: 2026-10-12. Maintenance data comes from `git clone --bare --filter=tree:0` of each public repository on 2026-10-12 (last default-branch commit, newest tags, and commit count since 2025-10-12). The GitHub REST API was blocked in this session, so star and issue counts are not available. Bibliographic data was checked against Crossref and the arXiv API.

Evidence labels used below: **[PR]** peer-reviewed, **[Pre]** preprint, **[Doc]** official docs, README or changelog, **[Repo]** repository metadata or code observation, **[Comm]** community forum or blog (weak evidence), **[Ctx]** a measurement supplied in the task context that I did not reproduce here.

---

## Q1. Julia performance model and design rules (type stability, globals, closures, small fixed-size arrays) and the diagnostic tools

### Takeaway
Julia is fast only when the compiler can infer concrete types. The dominant NEMOpt costs match the pathologies that the papers and the official manual warn about: untyped globals, captured variables in closures, and heap-allocated small tensors. The fix is well established: pass all state through function arguments and immutable typed structs, keep kernels type-stable, use stack-allocated `Tensor`/`SVector` types for 2×2 algebra, and enforce these properties in CI with JET (inference), AllocCheck (static allocations), and BenchmarkTools (timing).

### Cited Findings
- **[PR]** Julia's design aims at "machine performance without sacrificing human convenience." It does this through "specialization and abstraction," where "multiple dispatch … picks the right algorithm for the right circumstance." Source: Bezanson, Edelman, Karpinski & Shah, "Julia: A Fresh Approach to Numerical Computing," *SIAM Review* 59(1):65–98, 2017, doi:[10.1137/141000671](https://doi.org/10.1137/141000671) (abstract via Crossref).
- **[PR]** Julia combines dynamic typing, automatic memory management, rich type annotations and multiple dispatch with "a specializing just-in-time compiler to eliminate much of the overhead of those features." The paper explains how these design choices affect performance and usability. Source: Bezanson, Chen, Chung, Karpinski, Shah, Vitek & Zoubritzky, "Julia: dynamism and performance reconciled by design," *Proc. ACM Program. Lang.* 2(OOPSLA), Article 120, 1–23, 2018, doi:[10.1145/3276490](https://doi.org/10.1145/3276490).
- **[PR]** Users "are expected to adhere to a coding discipline that enables so-called type stability." Informally, a function is type-stable if its output type depends only on its input types, not on their values. The paper formalizes type stability and a stronger property, *type groundedness*, which "enables compiler optimizations." It also presents a corpus analysis of real packages. Source: Pelenitsyn, Belyakova, Chung, Tate & Vitek, "Type stability in Julia: avoiding performance pathologies in JIT compilation," *Proc. ACM Program. Lang.* 5(OOPSLA), 1–26, 2021, doi:[10.1145/3485527](https://doi.org/10.1145/3485527); extended version [arXiv:2109.01950](https://arxiv.org/abs/2109.01950).
- **[Doc]** The *Performance Tips* section "Avoid untyped global variables" says: "The value of an untyped global variable might change at any point, possibly leading to a change of its type. This makes it difficult for the compiler to optimize code… Variables should be local, or passed as arguments to functions, whenever possible." It recommends `const` or typed globals and says "Passing arguments to functions is better style. It leads to more reusable code and clarifies what the inputs and outputs are." Source: [Julia manual, Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) (read from `doc/src/manual/performance-tips.md`, JuliaLang/julia master 30cd5fb, 2026-10-12).
- **[Doc]** The same manual has dedicated sections on several other pitfalls:
  - "Performance of captured variable" (closures)
  - "Write 'type-stable' functions"
  - "Avoid changing the type of a variable"
  - "Separate kernel functions (aka, function barriers)"
  - "Pre-allocate outputs"
  - "Be aware of when Julia avoids specializing"
  - "The dangers of abusing multiple dispatch (aka, more on types with values-as-parameters)"

  Source: [Julia manual, Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/).
- **[Doc]** "If your application involves many small (`< 100` element) arrays of fixed sizes … consider using StaticArrays.jl." These types avoid "unnecessary heap allocations" and let the compiler "specialize code for the *size* of the array, e.g. by completely unrolling vector operations … and storing elements in CPU registers." The manual's example is "computations with 2d geometries." Source: [Julia manual, Performance Tips, StaticArrays section](https://docs.julialang.org/en/v1/manual/performance-tips/#Consider-StaticArrays.jl-for-small-fixed-size-vector/matrix-operations).
- **[Doc]** The manual's "Tools" section lists Profiling/ProfileView, JET ("can help you find common performance problems"), `@time`/`@allocated`, `--track-allocation=user`, and `@code_warntype`. It advises: "If you don't see another reason for the allocations, suspect a type problem." Source: [Julia manual, Performance Tips → Tools](https://docs.julialang.org/en/v1/manual/performance-tips/#tools).
- **[Ctx]** Measured in this project: one current constitutive call takes about 40 µs and allocates 21.9 KB. A type-stable Tensors.jl version (stress and tangent from AD) takes about 0.08 µs and gives identical results, roughly a 500× speedup. *I did not reproduce this measurement.*
- **[Repo]** Code observation in the NEMOpt scripts:
  - Each `*_main.jl` declares `global 𝜌_Unfiltered, U_Nonlin, χ, 𝜌_Filtered, K_Filter, iter, Prev_J, p, β, P, A, v_12, v_0, Vf` inside the objective function (line 44).
  - The solver files use `@tullio` 11 times each.
  - There are no `const` declarations.

  Source: /home/user/NEMOpt/Full_Domain/NEMOpt_main.jl, /home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl, NEMOpt.jl, NEMOpt-DR.jl.
- **[Pre]** BenchmarkTools uses a benchmarking strategy that is "robust in the presence of timer error, OS jitter and other environmental fluctuations, and is insensitive to the highly nonideal statistics produced by timing measurements." It is "used in production continuous integration (CI) pipelines for developing the Julia language and its ecosystem." Source: Chen & Revels, "Robust benchmarking in noisy environments," [arXiv:1608.04295](https://arxiv.org/abs/1608.04295), 2016 (not peer-reviewed as far as I could verify). Maintenance: BenchmarkTools.jl v1.8.0 (2026-04-15), last commit 2026-10-08, 32 commits in the last 12 months ([repo](https://github.com/JuliaCI/BenchmarkTools.jl)).
- **[Doc]** JET has two analyses. `@report_opt` detects type instabilities and runtime dispatch. `@report_call` detects type errors. Maintenance: JET.jl v0.12.3 (2026-10-03), last commit 2026-10-12, 218 commits in 12 months. Source: [JET.jl README](https://github.com/aviatesk/JET.jl).
- **[Doc]** Cthulhu.jl is used for interactive descent into inferred code. Its README warns: "This package relies on internal APIs of the Julia base compiler … the latest version … does not guarantee compatibility with older Julia versions." Latest tag v2.16.7 (2026-03-03), last commit 2026-08-22. Source: [Cthulhu.jl](https://github.com/JuliaDebug/Cthulhu.jl).
- **[Doc]** AllocCheck "statically checks if a function call may allocate by analyzing the generated LLVM IR." Its entry point `@check_allocs f(x) = …` throws if the call may allocate or dynamically dispatch. This is "different from measuring allocations using … `@time` or `@allocated`." Maintenance: v0.2.7 (2026-10-12), 21 commits in 12 months. Source: [AllocCheck.jl](https://github.com/JuliaLang/AllocCheck.jl).
- **[Doc]** Aqua.jl automates checks for:
  - method ambiguities
  - undefined exports
  - unbound type parameters
  - stale dependencies
  - missing `compat` entries
  - "obvious" type piracy
  - persistent tasks

  Maintenance: v0.8.18 (2026-09-17), 27 commits in 12 months. Source: [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl).
- **[Doc]** Julia versions on 2026-10-12: current stable is **v1.13.1 (Sep 25, 2026)** and current LTS is **v1.10.12 (Aug 16, 2026)**. Source: [julialang.org manual downloads](https://julialang.org/downloads/manual-downloads/). Master reports `1.14.0-DEV` ([Repo]).

### Inferences
- The 40 µs vs 0.08 µs gap fits the documented pathologies: abstractly typed captured or global state, heap-allocated `Array` temporaries for 2×2 tensors, and Tullio-generated loops over small arrays. The fix is architectural, not micro-tuning. Every kernel (constitutive update, element residual and tangent, homogenization integrand, sensitivity integrand) should take explicit, concretely typed arguments, for example `Tensor{2,2,Float64}`/`SymmetricTensor{2,2}` from Tensors.jl or `SMatrix{2,2}`.
- Suggested CI "performance contract":
  - `JET.@test_opt` (or `report_opt`) on the top-level kernels.
  - `AllocCheck.@check_allocs` on the constitutive and element kernels, which should be zero-allocation by design.
  - BenchmarkTools `@benchmark` regression baselines for the constitutive call, the element routine, global assembly, and one full optimization iteration.
  - Aqua on the package.
- Pin Julia in CI to the current stable release (1.13.x) plus the LTS (1.10.x) only if the authors need LTS. Record the exact Julia version used in published runs, because compiler or LLVM changes can alter floating-point instruction selection (for example FMA contraction under `@fastmath`/`muladd`). This last point is general knowledge, not checked against a specific source here.
- Avoid encoding mesh size (`Nelx`, `Nely`) as type parameters. The manual warns against "abusing multiple dispatch" with values-as-parameters, because each new value triggers recompilation.

### Gaps
- I did not find peer-reviewed quantitative studies of the Tullio-with-closures pattern on very small (2×2) arrays. The claim that Tullio is suboptimal for 2×2 tensor algebra rests on the project's own measurement [Ctx] plus the manual's general StaticArrays guidance.
- I could not confirm that Chen & Revels (2016) was ever published in a peer-reviewed venue. Treat it as a preprint.

---

## Q2. Using multiple dispatch for extensibility (cell-symmetry types, material models, objectives)

### Takeaway
Dispatch on small singleton or parametric types is the idiomatic Julia replacement for NEMOpt's two near-duplicate scripts. Gridap's peer-reviewed software-design paper and the core Julia papers document this pattern for FE codes. Shared code stays generic, and only the symmetry-specific pieces become methods specialized on `FullCell`/`HalfCell` (constraint construction, design-to-element mapping, homogenization domain weights).

### Cited Findings
- **[PR]** Gridap "provides a low-level API that is modular and extensible via the multiple-dispatch paradigm of Julia." It "implements this API without introducing a DSL plus a compiler of variational forms. Instead, it leverages the Julia just-in-time compiler to build efficient code, specialized for the concrete problem at hand." The paper's benchmarks show assembly performance "comparable to FEniCS, demonstrating that the new software design does not compromise performance." Source: Verdugo & Badia, "The software design of Gridap: A Finite Element package based on the Julia JIT compiler," *Computer Physics Communications* 276 (2022) 108341, doi:[10.1016/j.cpc.2022.108341](https://doi.org/10.1016/j.cpc.2022.108341); abstract from [arXiv:2109.12818](https://arxiv.org/abs/2109.12818).
- **[PR]** Gridap is "an extensible Finite Element toolbox in Julia." Source: Badia & Verdugo, *Journal of Open Source Software* 5(52):2520, 2020, doi:[10.21105/joss.02520](https://doi.org/10.21105/joss.02520).
- **[PR]** Multiple dispatch "picks the right algorithm for the right circumstance," and generic programming captures "what remains the same after differences are stripped away." Source: Bezanson et al. 2017, *SIAM Review* 59:65, doi:[10.1137/141000671](https://doi.org/10.1137/141000671).
- **[Comm/talk]** Karpinski, "The Unreasonable Effectiveness of Multiple Dispatch," JuliaCon 2019. The talk argues that the unusually high degree of code sharing between Julia packages comes from multiple dispatch, which ties methods to generic functions rather than classes and makes double-dispatch patterns unnecessary. Source: [JuliaCon 2019 talk page](https://pretalx.com/juliacon2019/talk/BCYWZJ/). This is a talk, not peer-reviewed.
- **[PR]** The SciML ecosystem's founding paper describes a performant, feature-rich ecosystem built in Julia. Source: Rackauckas & Nie, "DifferentialEquations.jl – A Performant and Feature-Rich Ecosystem for Solving Differential Equations in Julia," *J. Open Research Software* 5(1):15, 2017, doi:[10.5334/jors.151](https://doi.org/10.5334/jors.151). *I verified only the bibliographic metadata; I did not review the content on dispatch-based architecture.*
- **[Doc]** The Julia manual warns against values-as-type-parameters abuse, function barriers that hide type instability, and cases where Julia "avoids specializing" on `Function`, `Type` and `Vararg` arguments unless a type parameter forces it. Source: [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing).
- **[Repo]** In NEMOpt, the full-domain and domain-reduction drivers differ mainly in symmetry handling:
  - `NEMOpt.jl` has 856 lines and `NEMOpt-DR.jl` has 1,000 lines.
  - The drivers are 249 and 255 lines and are nearly identical.
  - The full-domain driver uses `MMA87()`; the DR driver uses `MMA02()`.
  - Both drivers use 1/4 of the elements as design variables (`addvar!(model, zeros(Int(N_Elements / 4)), …)`).

  Source: /home/user/NEMOpt/*/NEMOpt*_main.jl.

### Inferences
- Proposed type hierarchy (design sketch, not from a source):
  - `abstract type CellSymmetry end`, with `struct FullCell <: CellSymmetry end` and `struct HalfCell <: CellSymmetry end`. Methods: `periodic_constraints(::CellSymmetry, grid)`, `design_map(::CellSymmetry, mesh)` (reduced design vector → element densities), `homogenization_weights(::CellSymmetry, …)`, and `macro_load_cases(::CellSymmetry)`.
  - `abstract type HyperelasticMaterial end`, with `struct NeoHookean{T} <: HyperelasticMaterial; μ::T; λ::T; end`. One method `ψ(m, C)` (strain energy). Generic `stress_tangent(m, F)` computes stress and tangent through `Tensors.hessian(…, :all)`. Hand-coded specialized methods can be added later as optional fast paths, with tests against the AD path.
  - `abstract type Interpolation end` (SIMP, RAMP), `abstract type Filter end` (density, Helmholtz), `abstract type Projection end` (Heaviside with β).
  - `abstract type Objective end`, for example `struct PoissonTarget{T}; ν0::T; end`, `struct BulkMax end`, `struct ShearMax end`. Methods `J(obj, A)` and `∂J∂A(obj, A)`. The commented-out alternative objectives in the drivers (`J = 1/(A[3,3]+A[4,4])`, bulk modulus) become selectable types instead of comment toggles.
  - `struct Continuation` (schedules for `p` and `β`) and an immutable `Problem`/`Params` struct loaded from TOML. These replace all globals.
- Dispatch on singleton types is resolved at compile time when the types are concrete, so it costs nothing at runtime. A `Union{FullCell,HalfCell}` field is also fine because the compiler splits small unions; this is a general Julia property and not checked against a specific source here.
- Keep the type hierarchy shallow. Use dispatch only at the true variation points (symmetry, material, objective, filter, optimizer backend). Deep hierarchies are a known maintainability risk. This is my opinion, not sourced.

### Gaps
- I found no peer-reviewed study that measures maintainability gains from dispatch-based design in scientific Julia codes. The evidence is design papers (Gridap CPC) and talks.
- I did not review SciML's design documentation in detail. I verified only the 2017 JORS paper's metadata.

---

## Q3. FE ecosystem fit: Ferrite.jl, Tensors.jl, Gridap.jl, GridapTopOpt.jl, TopOpt.jl

### Takeaway
**Ferrite.jl + Tensors.jl** is the best fit for 2D structured periodic unit cells with finite-strain hyperelasticity, and it is very actively maintained (Ferrite v1.7.0, 2026-08-31/09-04, 198 commits in 12 months). It provides `PeriodicDirichlet`, general `AffineConstraint`, hyperelasticity and computational-homogenization tutorials, and deterministic colored threaded assembly.

One limitation matters for the half-cell variant. Released Ferrite (≤ v1.7.0) cannot resolve constraints where a DOF is both master and slave. Untangling was merged to master on 2026-09-28 but was **unreleased** on 2026-10-12.

Gridap is mature and very active, but it is a different (weak-form, lazily evaluated) paradigm. GridapTopOpt is level-set based, and TopOpt.jl is linear-elastic only. Any library adoption changes DOF numbering, quadrature ordering and summation order, so bitwise reproduction of legacy results is not achievable. Verification must be tolerance-based.

### Cited Findings

**Ferrite.jl** ([repo](https://github.com/Ferrite-FEM/Ferrite.jl), [docs](https://ferrite-fem.github.io/Ferrite.jl/stable/))
- **[Repo]** Release history:
  - v1.0.0 2024-09-30
  - v1.1.0 2025-05-01
  - v1.2.0 2025-11-14
  - v1.3.0 2026-01-09
  - v1.4.0 2026-04-20
  - v1.5.0 2026-07-13
  - v1.6.0 2026-08-02
  - v1.7.0 2026-08-31 (git tag dated 2026-09-04)

  Last commit 2026-10-09; 198 commits since 2025-10-12. CI includes benchmark tracking and downstream tests. Formatting uses Runic via pre-commit (`Runic@1.11`). Source: [CHANGELOG.md](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md).
- **[Doc]** `PeriodicDirichlet(u, facet_mapping, components)` constrains DOFs on a "mirror" facet to an "image" facet. Variants accept a rotation matrix `R` or an inhomogeneity function. Rotations are "currently only supported when all components are periodic." Source: `src/Dofs/ConstraintHandler.jl` docstring ([repo](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/src/Dofs/ConstraintHandler.jl)).
- **[Doc]** The computational homogenization tutorial covers the 2D unit cell setup:
  - It uses `collect_periodic_facets(grid, "left", "right")` and `collect_periodic_facets!(…, "bottom", "top")`, then `PeriodicDirichlet(:u, periodic_facets, [1, 2])`.
  - It explains that `PeriodicDirichlet` is "simply just a more convenient way of" adding `AffineConstraint(u_m, [u_i => 1], 0)` constraints.
  - A newer "stress-driven computational homogenization" tutorial (v1.6.0) treats macroscopic strain as an unknown algebraic variable.

  Source: [computational homogenization tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/computational_homogenization/); CHANGELOG v1.6.0.
- **[Doc]** The hyperelasticity tutorial (keywords: "hyperelasticity, finite strain, large deformations, Newton's method") computes `∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)`. It also cross-checks the first Piola–Kirchhoff tangent against AD (`Tensors.hessian(x -> Ψ(tdot(x), mp), F, :all)`). Source: [hyperelasticity tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/hyperelasticity/).
- **[Doc] Limitation in released versions.** Ferrite v1.7.0 docs say: "Ferrite currently cannot untangle constraints where a DoF is both *master* and *slave*." The v1.7.0 source raises `error("nested affine constraints currently not supported")`. Source: v1.7.0 `docs/src/topics/constraints.md` and `src/Dofs/ConstraintHandler.jl` ([constraints docs](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/constraints/)).
- **[Repo] Change on master, unreleased.** Commit "Untangle affine constraints in close! (#1327)" landed on 2026-09-28. The CHANGELOG lists it under *Unreleased*: "Untangle the tangled `AffineConstraints` in `close!`". Dev docs say "In `close!` the tangled affine constraints are untangled" and point to the internal `Ferrite._istangled`. Also unreleased: "Constraint matrices correctly account for affine constraints with prescribed masters," and sparsity patterns with `keep_constrained=false` now include entries for affine constraints whose masters are in other cells, "e.g. periodic boundary conditions. Previously this resulted in a `SparsityError`" (#1470, #1546). Source: [CHANGELOG](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md).
- **[Doc] Correctness fix to note.** v1.6.0 fixed "bug applying the transpose operation in condensation of `AffineConstraints`. This bug gave silently wrong results when used on non-symmetric system matrices" (#1426). Source: CHANGELOG v1.6.0.
- **[Doc] Threaded assembly.** v1.6.0 added `start_assemble(K, f; atomic = true)` for lock-free concurrent assembly, with "non-deterministic (but correct up to summation order) results." v1.7.0 made `create_coloring` 3–10× faster serially plus 2–3× with 8 threads on 200k–400k-cell grids. The "resulting coloring is independent of the number of threads." Unreleased: `task_local_copy` replaces `copy` for per-task scratch data. Source: CHANGELOG; [threaded assembly howto](https://ferrite-fem.github.io/Ferrite.jl/stable/howto/threaded_assembly/).

**Tensors.jl** ([repo](https://github.com/Ferrite-FEM/Tensors.jl))
- **[PR]** Tensors.jl provides tensor computations in Julia. Source: Carlsson & Ekre, "Tensors.jl — Tensor Computations in Julia," *J. Open Research Software* 7(1):7, 2019, doi:[10.5334/jors.182](https://doi.org/10.5334/jors.182).
- **[Doc]** "`Tensors` supports forward mode automatic differentiation (AD) … by exploiting the `Dual` number defined in `ForwardDiff.jl`." With its own API (`gradient(f, A)`, `hessian(f, A)`, and the `:all` variant that also returns lower-order results), "everything will be stack allocated so there is no need to preallocate." Source: [Tensors.jl AD docs](https://ferrite-fem.github.io/Tensors.jl/stable/man/automatic_differentiation/).
- **[Repo]** Tensors.jl v1.17.1 (2026-06-26), last commit 2026-07-18, 17 commits in 12 months.

**Gridap.jl / GridapTopOpt.jl**
- **[PR]** Gridap papers: JOSS 2020 doi:[10.21105/joss.02520](https://doi.org/10.21105/joss.02520); CPC 2022 doi:[10.1016/j.cpc.2022.108341](https://doi.org/10.1016/j.cpc.2022.108341). The CPC paper lists "non-linear solid mechanics" among its applications.
- **[Repo]** Gridap.jl v0.20.10 (2026-09-21), last commit 2026-10-05, **680 commits** in 12 months. Version is still 0.x, so breaking changes are allowed under semver.
- **[Repo]** `CartesianDiscreteModel` supports periodicity through `isperiodic` ([source](https://github.com/gridap/Gridap.jl/blob/master/src/Geometry/CartesianDiscreteModels.jl)). The official [Tutorials](https://github.com/gridap/Tutorials) repo (v0.20.7, 2026-05-13) includes `hyperelasticity.jl` and `validation_DrWatson.jl`.
- **[PR]** GridapTopOpt is "an extendable framework for **level set-based** topology optimisation." It uses Gridap assembly and PETSc solvers and implements AD for sensitivities. Source: Wegert, Manyer, Mallon, Badia & Challis, "GridapTopOpt.jl: a scalable Julia toolbox for level set-based topology optimisation," *Structural and Multidisciplinary Optimization* 68:22, 2025, doi:[10.1007/s00158-024-03927-3](https://doi.org/10.1007/s00158-024-03927-3) ([arXiv:2405.10478](https://arxiv.org/abs/2405.10478)). Follow-up: Wegert et al., *CMAME* 445 (2025) 118203, doi:[10.1016/j.cma.2025.118203](https://doi.org/10.1016/j.cma.2025.118203) (unfitted FE plus automatic shape differentiation).
- **[Repo]** GridapTopOpt v0.5.1 (2026-09-16), 157 commits in 12 months. The README warns: "there have been several breaking releases since the above publication." The repo includes hyperelastic (including neo-Hookean) compliance and inverse-homogenization example scripts (`scripts/Examples/BoundarySmoothing/...`), and `NonlinearFEStateMap` / `ReverseNonlinearFEStateMap`. There is only one density-based example (`Density/ElectromagneticScattering.jl`). Source: [repo](https://github.com/zjwegert/GridapTopOpt.jl).

**TopOpt.jl**
- **[Doc]** TopOpt.jl "supports structural (**linear elasticity**) and heat-transfer problems on continuum and truss ground meshes in 2D and 3D, with automatic differentiation through every objective and constraint." Its dependencies include Ferrite, Tensors, Nonconvex, NonconvexMMA, NonconvexPercival, Zygote and ChainRulesCore. Source: [TopOpt.jl README / Project.toml](https://github.com/JuliaTopOpt/TopOpt.jl).
- **[Repo]** TopOpt.jl went from v0.12.0 to v0.14.0 between 2026-08-12 and 2026-08-18 and added a level-set (OpenLSTO port) module. It has 74 commits in 12 months, and its last commit is 2026-08-18.
- **[PR]/[Conf]** Citations requested by the TopOpt.jl README:
  - Tarek & Ray, "Adaptive continuation solid isotropic material with penalization for volume constrained compliance minimization," *CMAME* 363 (2020) 112880, doi:[10.1016/j.cma.2020.112880](https://doi.org/10.1016/j.cma.2020.112880).
  - Tarek, "TopOpt.jl: An efficient and high-performance package for topology optimization of continuum structures in the Julia programming language," Proc. 13th WCSMO, 2019 (no DOI found).
  - Huang & Tarek, "TopOpt.jl: Truss and Continuum Topology Optimization, Interactive Visualization, Automatic Differentiation and More," WCSMO 2021 ([PDF](https://web.mit.edu/yijiangh/www/papers/topopt_jl_WCSMO2021.pdf)).

### Inferences
- **Recommended FE layer:** Ferrite.jl (≥ 1.6, for the condensation transpose fix and atomic/colored assembly) plus Tensors.jl. Reasons:
  - It is imperative and low-level, so the authors keep full control of the element loop, the Newton solver and the adjoint. This matters for hand-coded sensitivities.
  - Its tutorials match the problem class (finite-strain hyperelasticity, periodic homogenization).
  - It is already the substrate of TopOpt.jl.
- **Half-cell (domain reduction) risk.** If the half-cell constraints chain mirror or symmetry relations onto periodic relations, so that one DOF is a slave in one relation and a master in another, the constraints are tangled. Released Ferrite rejects this. Three options:
  - (a) Pin Ferrite to a master commit after 2026-09-28 (or the next release) and test with `Ferrite._istangled`, which is internal.
  - (b) Pre-resolve chains by hand into untangled `AffineConstraint`s, each slave written in terms of free masters only.
  - (c) Keep NEMOpt's own master–slave reduction (a sparse transformation matrix `T` with u = T û) and use Ferrite only for assembly.

  Option (c) best preserves the legacy algebra. Whether chains actually occur depends on the DR paper's specific symmetry mapping. *Verify against the article.*
- **Gridap** is better suited to rapid weak-form prototyping and parallel (PETSc) scaling than to bit-level control. Its 0.x versioning and very high commit rate mean more API churn for a small research code. GridapTopOpt's level-set formulation does not match NEMOpt's density/SIMP + filter + Heaviside-projection formulation.
- **TopOpt.jl** is linear-elastic only, so it cannot host NEMOpt's finite-strain homogenization. It is still a useful reference for SIMP, filters, and the MMA/Nonconvex integration on Ferrite.
- **Bitwise reproduction.** Any of these libraries changes DOF numbering (Ferrite's `DofHandler` vs the hand-rolled numbering), element and quadrature traversal order, sparse pattern, and floating-point summation order. Ferrite's own howto stresses that summation order affects results (coloring is deterministic; locking and atomics are not). Topology optimization with continuation is nonconvex and path-dependent, so 1e-15 differences can grow into different final topologies. Recommended two-tier verification:
  1. **Pure restructuring** (no library change, same operation order) is checked against legacy outputs bitwise or to ≤1e-14 relative.
  2. **Library swaps** are checked with tolerance tests on single-call quantities (constitutive stress/tangent, element K and f, homogenized A and P for a fixed density field, sensitivities vs finite differences), and then on early optimization iterates. Final designs are compared on objective, volume fraction and ν₁₂ within agreed tolerances, not by bitwise equality.

### Gaps
- The Ferrite #1327 untangling has no release date yet. Re-check the [Ferrite releases](https://github.com/Ferrite-FEM/Ferrite.jl/releases) before adopting.
- I did not verify whether Gridap supports general affine (non-identity, sign-flipping) multipoint constraints for symmetry mappings beyond `isperiodic`.
- I did not find any published nonlinear (finite-strain) inverse-homogenization density-based topology optimization implemented on Ferrite. Feasibility is inferred from the tutorials.
- Star and issue counts and bus-factor data could not be gathered because the GitHub API was blocked.

---

## Q4. Optimization layer: MMA implementations, maintenance, hand-coded gradients, continuation

### Takeaway
NonconvexMMA (MMA87/MMA02) was revived in 2026 (v1.1.2, 2026-08-08) after about 27 months without commits. It still pulls in Zygote, ChainRulesCore, ForwardDiff and Optim. Its documented hand-gradient route is `CustomGradFunction(f, g)`, not hand-written `rrule`s.

NLopt (`:LD_MMA`/`:LD_CCSAQ`, an independent CCSA implementation by S. G. Johnson after Svanberg 2002) is the lowest-friction option for hand-coded gradients, because it uses an in-place `(x, grad)` callback with no AD plumbing. Its MMA is not the same algorithm as NonconvexMMA's MMA87, so switching changes optimization trajectories. A small in-house MMA/GCMMA, written from Svanberg's papers, is the only route to full control of continuation and determinism.

Continuation should be an outer loop with explicit restarts (or an Optimization.jl callback that halts), not a mutation of hidden state inside the objective.

### Cited Findings
- **[PR]** Primary algorithms:
  - Svanberg, "The method of moving asymptotes—a new method for structural optimization," *IJNME* 24(2):359–373, 1987, doi:[10.1002/nme.1620240207](https://doi.org/10.1002/nme.1620240207).
  - Svanberg, "A Class of Globally Convergent Optimization Methods Based on Conservative Convex Separable Approximations," *SIAM J. Optim.* 12(2):555–573, 2002, doi:[10.1137/S1052623499362822](https://doi.org/10.1137/S1052623499362822).
- **[Repo]** Maintenance on 2026-10-12:

| Package | Latest release | Last commit | Commits since 2025-10-12 | Notes |
|---|---|---|---|---|
| [NonconvexMMA.jl](https://github.com/JuliaNonconvex/NonconvexMMA.jl) | v1.1.2 (2026-08-08); v1.1.0 2026-02-09 | 2026-08-08 | 8 | Before Feb 2026, last commit 2023-11-15. Deps: ChainRulesCore, ForwardDiff, NonconvexCore, **Optim**, Parameters, Reexport, **Zygote** |
| [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) | **v3.0.0 (2026-08-11)**, v2.2.0 2026-03-09 | 2026-08-11 | 14 | v3 "require NonconvexUtils 1.0" (breaking) |
| [NonconvexCore.jl](https://github.com/JuliaNonconvex/NonconvexCore.jl) | v1.6.0 (2026-08-11) | 2026-08-11 | 12 | |
| [NonconvexNLopt.jl](https://github.com/JuliaNonconvex/NonconvexNLopt.jl) / [NonconvexIpopt.jl](https://github.com/JuliaNonconvex/NonconvexIpopt.jl) | v0.1.9 (2026-02-09) / v0.6.0 (2026-08-11) | 2026-02-10 / 2026-08-11 | 4 / 5 | |
| [NLopt.jl](https://github.com/jump-dev/NLopt.jl) (jump-dev) | v1.2.1 (2025-07-31) | 2026-08-27 | 4 | Wraps the C library |
| [nlopt](https://github.com/stevengj/nlopt) (C) | v2.11.0 (2026-06-05) | 2026-06-06 | 14 | |
| [Optimization.jl](https://github.com/SciML/Optimization.jl) (SciML monorepo) | many sub-package tags on 2026-10-12 | 2026-10-11 | 543 | Has `lib/OptimizationNLopt`, `OptimizationIpopt`, `OptimizationMOI`, `OptimizationMadNLP`, `OptimizationAuglag`, …; **no MMA-specific sub-package** besides via NLopt |
| [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) | v1.16.0 (2026-09-08) | 2026-09-08 | 32 | |
| [mohamed82008/MMA.jl](https://github.com/mohamed82008/MMA.jl) | none | 2018-12-12 | 0 | Predecessor of NonconvexMMA; abandoned |
| [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) | v1.26.1 (2026-04-02) | 2026-06-22 | 5 | |
| [Zygote.jl](https://github.com/FluxML/Zygote.jl) | v0.7.13 (2026-08-31) | 2026-08-31 | 60 | |

- **[Doc]** NonconvexMMA offers "the original MMA algorithm from the 1987 paper" (`MMA87()`/`MMA()`) and "the globally convergent MMA (GCMMA) algorithm from the 2002 paper" (`MMA02()`/`GCMMA()`). These support only inequality constraints. Source: [Nonconvex MMA docs](https://julianonconvex.github.io/Nonconvex.jl/stable/algorithms/mma/).
- **[Doc]** `MMAOptions` fields and defaults:
  - `maxiter=1000`, `outer_maxiter=10^8`, `maxinner=10` (MMA02 only)
  - `tol`, `s_init=0.5`, `s_incr=1.2`, `s_decr=0.7`
  - `store_trace`, `show_trace`, `auto_scale`, `keep_best`, `verbose`
  - `dual_options::Optim.Options`
  - `convcriteria=KKTCriteria()`

  The dual subproblem is solved with Optim.jl, and the constructor default is `dualoptimizer = Optim.GradientDescent()`. **No per-iteration user callback field** exists. Source: `src/mma_algorithm.jl` in [NonconvexMMA.jl](https://github.com/JuliaNonconvex/NonconvexMMA.jl/blob/master/src/mma_algorithm.jl).
- **[Doc]** The documented way to supply hand-derived gradients in Nonconvex is `F = CustomGradFunction(f, g)`. "Whenever `ForwardDiff` or any `ChainRules`-compatible AD package such as `Zygote` is used to differentiate `F`, the custom gradient/Jacobian function `g` will be used." There is also `CustomHessianFunction(f, g, h)` with an `hvp` option. Source: [Nonconvex docs, user-defined gradients](https://julianonconvex.github.io/Nonconvex.jl/stable/gradients/user_defined/).
- **[Doc]** NLopt's MMA (`NLOPT_LD_MMA`) is Johnson's implementation of "the globally-convergent method-of-moving-asymptotes (MMA) algorithm," as described in Svanberg 2002. It is "an improved CCSA … variant of the original MMA algorithm published by Svanberg in 1987." It is "completely independent of Svanberg's, whose code we have not examined." Johnson writes: "I generally recommend trying `NLOPT_LD_CCSAQ` first rather than `NLOPT_LD_MMA`." Tunable parameters: `inner_maxeval`, `dual_algorithm`, `dual_ftol_rel` (default 1e-14), `rho_init`, `inner_gradients`, `always_improve`. Source: [NLopt Algorithms, MMA and CCSA](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/) (read from `doc/docs/NLopt_Algorithms.md`, stevengj/nlopt master 2026-06-06).
- **[Doc]** The NLopt.jl API takes hand-coded gradients directly:

  ```julia
  function f(x::Vector, grad::Vector)
      if length(grad) > 0
          grad .= ...
      end
      return value
  end
  opt = NLopt.Opt(:LD_MMA, n)
  NLopt.min_objective!(opt, f)
  NLopt.inequality_constraint!(opt, c, tol)
  NLopt.optimize(opt, x0)
  ```

  Source: [NLopt.jl README](https://github.com/jump-dev/NLopt.jl).
- **[Doc]** In Optimization.jl, "When supported by `alg`, `callback` is called after an optimization step as `callback(state, objective)` … returning `true` stops the optimization." `maxiters` is a common keyword. Source: [Optimization.jl solve API](https://docs.sciml.ai/Optimization/stable/API/solve/).
- **[Pre]** Svanberg's MATLAB MMA code (`mmasub.m`, `subsolv.m`) "is freely available for download under the GNU GPLv3 license." Source: Alexandersen, "A detailed introduction to density-based topology optimisation of fluid flow problems with implementation in MATLAB," [arXiv:2207.13695](https://arxiv.org/abs/2207.13695). Julia ports of that code exist, for example [pollinico/topopt_jl](https://github.com/pollinico/topopt_jl), but I did not check their maintenance or licensing.
- **[Repo] How NEMOpt does continuation today.** The objective closure handles continuation itself:
  - It increments a global `iter` on every objective evaluation.
  - It raises `p` by 0.5 every 5 evaluations until p = 3, then raises `β` by 2 every 10 evaluations up to 32.
  - It assigns `options = MMAOptions(s_init=0.5/(β+1), …)` *inside* the objective function. That `options` is not declared `global` there; the global declaration on line 44 omits it.
  - It rounds `J` (`J, Vf = round(J, digits=3), …`) before `return J`.
  - It calls `Plot_Filtered_Topology` on every evaluation.

  Source: /home/user/NEMOpt/Full_Domain/NEMOpt_main.jl lines 44–106 and 225–240; the DR driver is identical at lines 44–100.

### Inferences
- **Probable latent bugs that must be preserved or consciously fixed** (code reading, *not* run-verified):
  1. The `options = MMAOptions(...)` assignment inside `Objective` creates a function-local binding. It cannot affect the `options` object already passed to `optimize`, so the intended `s_init = 0.5/(β+1)` schedule is probably a no-op.
  2. The value returned to MMA is rounded to 3 decimals. This matters for MMA02's conservativeness checks and for KKT/`f` tolerances, so a "faithful" refactor must reproduce the rounding to match legacy trajectories.
  3. `iter` counts objective *evaluations*, not outer iterations. Under MMA02 (inner iterations) the continuation schedule therefore depends on how many inner evaluations occur. This differs between MMA87 (full domain) and MMA02 (DR).

  These should be confirmed by the reproduction-baseline work and documented as "legacy semantics."
- **Hand-coded gradient plumbing.** If Nonconvex is kept, replace the custom `rrule`s with `CustomGradFunction(f, g)`, which is documented and simpler. The cleanest decoupling is to define the optimizer interface as plain `(x, grad) -> f`. NLopt uses this natively, Optimization.jl accepts user gradients, and an in-house MMA can use it directly. AD frameworks are then needed only for tests.
- **Optimizer choice for reproduction vs future work.** To reproduce the paper, keep NonconvexMMA (pin the exact version in `Manifest.toml`), because NLopt's CCSA-MMA is a different algorithm with different inner-loop logic. For future work, consider an in-house MMA87/GCMMA of about 200–400 lines, written clean-room from Svanberg 1987/2002 to avoid GPL entanglement. That gives deterministic control, explicit continuation hooks and no AD dependencies. Alternatively use NLopt `:LD_MMA`/`:LD_CCSAQ` through a thin wrapper.
- **Continuation design.** Use an explicit outer loop over continuation stages (`for (p, β) in schedule`), each running `optimize` for N iterations warm-started from the last design. Alternatively use an optimizer with an iteration callback (Optimization.jl `callback(state, obj)` returning `true` to stop, then re-solve with updated parameters). The trade-off: restarts reset MMA asymptote history, while in-objective mutation makes the objective non-stationary within a run, which conflicts with GCMMA's conservativeness logic and KKT convergence tests. Keep a "legacy mode" that reproduces the evaluation-count-driven schedule for paper reproduction.
- Move plotting and I/O out of the objective into a per-iteration callback or logger. Per-evaluation plotting also inflates wall time.

### Gaps
- I did not verify the exact keyword (`grad=`) for supplying hand-written gradients in Optimization.jl's `OptimizationFunction`. Check the [OptimizationFunction docs](https://docs.sciml.ai/Optimization/stable/).
- I did not verify whether NLopt.jl exposes `force_stop!` for clean mid-run termination.
- I found no maintained, registered, standalone native-Julia MMA package other than NonconvexMMA. A web search turned up nothing for 2025–2026, and the 2018 `MMA.jl` is abandoned.
- The reasons for NonconvexMMA's 2023-11 to 2026-02 dormancy and its bus factor are unknown.

---

## Q5. Automatic differentiation for constitutive laws and verification of hand-derived tangents

### Takeaway
For 2×2 finite-strain constitutive laws, use forward-mode AD on stack-allocated tensors (`Tensors.gradient`/`Tensors.hessian(..., :all)`, built on ForwardDiff duals). It is allocation-free, the Ferrite hyperelasticity tutorial already uses it, and the project measured it at about 0.08 µs per call. Keep hand-derived tangents and adjoint sensitivities only where they buy speed, and test them against AD and against finite differences. Enzyme is the option for reverse-mode or whole-program adjoint checks, but it adds heavier dependency and compile-latency costs.

### Cited Findings
- **[Pre]** ForwardDiff is "a Julia package for forward-mode automatic differentiation (AD) featuring performance competitive with low-level languages like C++." It "takes advantage of just-in-time (JIT) compilation to transparently recompile AD-unaware user code," and supports higher-order differentiation and custom number types. Source: Revels, Lubin & Papamarkou, "Forward-Mode Automatic Differentiation in Julia," [arXiv:1607.07892](https://arxiv.org/abs/1607.07892), 2016. Maintenance: ForwardDiff.jl v1.4.6 (2026-09-11), 28 commits in 12 months ([repo](https://github.com/JuliaDiff/ForwardDiff.jl)).
- **[PR]** Enzyme is "a high-performance automatic differentiation (AD) compiler plugin for the LLVM compiler framework capable of synthesizing gradients of statically analyzable programs expressed in the LLVM intermediate representation." Sources:
  - Moses & Churavy, "Instead of Rewriting Foreign Code for Machine Learning, Automatically Synthesize Fast Gradients," *Advances in Neural Information Processing Systems 33* (NeurIPS 2020), [arXiv:2010.01709](https://arxiv.org/abs/2010.01709).
  - Follow-up: Moses et al., "Reverse-mode automatic differentiation and optimization of GPU kernels via Enzyme," SC'21, doi:[10.1145/3458817.3476165](https://doi.org/10.1145/3458817.3476165).

  Maintenance: Enzyme.jl v0.13.216 (2026-10-12), **655 commits** in 12 months, still 0.x ([repo](https://github.com/EnzymeAD/Enzyme.jl)).
- **[Doc]** Tensors.jl AD: "it is simpler to use `Tensors` own AD API … This does not require any conversions and everything will be stack allocated." `hessian(f, A, :all)` returns the Hessian, gradient and value. Source: [Tensors.jl AD docs](https://ferrite-fem.github.io/Tensors.jl/stable/man/automatic_differentiation/); Carlsson & Ekre 2019, doi:[10.5334/jors.182](https://doi.org/10.5334/jors.182).
- **[Doc]** The Ferrite hyperelasticity tutorial computes `∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)`. It also checks the analytic `∂P∂F` against `Tensors.hessian(x -> Ψ(tdot(x), mp), F, :all)`. Source: [Ferrite hyperelasticity tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/hyperelasticity/).
- **[Ctx]** The project's Tensors.jl AD constitutive prototype gave identical results to the legacy implementation at about 0.08 µs per call.
- **[Doc]** Nonconvex's `CustomGradFunction` makes ForwardDiff or Zygote use the user gradient. Source: [Nonconvex docs](https://julianonconvex.github.io/Nonconvex.jl/stable/gradients/user_defined/).

### Inferences
- Use a layered verification strategy:
  - **(i) Constitutive.** Compare the hand tangent `∂P/∂F` (or `∂S/∂C`) to `Tensors.hessian` at randomized F with det F > 0, including near-reference and large-stretch states, to rtol ≈ 1e-12. Check the symmetries (major symmetry for hyperelastic laws).
  - **(ii) Element.** Compare the element tangent to a ForwardDiff Jacobian of the element residual with respect to element DOFs.
  - **(iii) Global sensitivities.** Run central finite-difference checks of dJ/dρ on a coarse mesh (for example 8×8) with random densities, at a few step sizes, and report the relative error. This is the standard "gradient check."
  - **(iv) Optional.** Use Enzyme reverse-mode on a small-mesh whole-pipeline function as an independent adjoint check.
- Keep AD packages in the test environment (`[extras]`/`test/Project.toml`) where possible. Runtime then depends only on Tensors.jl, which itself uses ForwardDiff duals.
- AD and closed-form expressions are not guaranteed to be bitwise identical in general, because operation order differs. The "identical results" observation [Ctx] should be re-checked across the full range of states, and tests should use tolerances.

### Gaps
- I found no peer-reviewed benchmark comparing Tensors.jl AD vs hand-coded tangents vs Enzyme for small constitutive kernels. Evidence is the project's own measurement plus documentation.
- I did not verify the NeurIPS proceedings page URL. The citation is given via arXiv.

---

## Q6. Package and workflow engineering (environments, extensions, configuration, testing, CI, docs, formatting, latency)

### Takeaway
Make NEMOpt a proper package. It needs:
- `Project.toml` with `[compat]`, and a committed `Manifest.toml` for the paper-reproduction environment
- `src/` modules and `scripts/` drivers
- TOML run configurations
- plotting moved to a package extension
- tests with `Test` + TestItemRunner, golden-file regression via ReferenceTests (low activity) or plain JLD2/DelimitedFiles comparisons
- Aqua + JET in CI
- GitHub Actions (`julia-actions/setup-julia@v3`, `cache@v3`, `julia-runtest`)
- Documenter docs
- an automatic formatter (Runic, as used by Ferrite, or JuliaFormatter with BlueStyle)
- PrecompileTools to reduce time-to-first-run

DrWatson is optional but well suited to parameter sweeps and result naming.

### Cited Findings
- **[Doc]** A `Manifest.toml` "includes exact information about (direct and indirect) dependencies of the project. Given a `Project.toml` + `Manifest.toml` pair, it is possible to instantiate the exact same package environment, which is very useful for reproducibility." Manifests also record registry metadata so "the exact source of each package version can be identified." Source: [Pkg docs, Project and Manifest](https://pkgdocs.julialang.org/v1/toml-files/) (read from JuliaLang/Pkg.jl master).
- **[Doc]** Package extensions use `[weakdeps]` and `[extensions]` sections in `Project.toml` for code that loads only when an optional dependency is present. Source: [Pkg docs, Creating packages → Conditional loading (Extensions)](https://pkgdocs.julialang.org/v1/creating-packages/).
- **[PR]** DrWatson is a "scientific project assistant" for managing simulations, simulation parameters, file locations, scripts and reproducibility. Source: Datseris, Isensee, Pech & Gál, "DrWatson: the perfect sidekick for your scientific inquiries," *JOSS* 5(54):2673, 2020, doi:[10.21105/joss.02673](https://doi.org/10.21105/joss.02673). Maintenance: v2.19.1 (2025-08-25), last commit 2026-09-28, 9 commits in 12 months ([repo](https://github.com/JuliaDynamics/DrWatson.jl)). Gridap's official tutorials include a `validation_DrWatson.jl` example ([gridap/Tutorials](https://github.com/gridap/Tutorials)).
- **[Doc]** TestItemRunner "runs `@testitem` tests." The framework is documented with the VS Code extension but "can be used entirely without VS Code." Maintenance: v1.3.2 (2026-08-31), 64 commits in 12 months. Source: [TestItemRunner.jl](https://github.com/julia-vscode/TestItemRunner.jl).
- **[Doc]** ReferenceTests.jl adds macros "for testing values against reference files, which in turn the package can help create and update." It builds on FileIO. Maintenance is **low**: last release v0.10.5 (2024-11-08), last commit 2025-11-30, 2 commits in 12 months. Source: [ReferenceTests.jl](https://github.com/JuliaTesting/ReferenceTests.jl).
- **[Doc]** PrecompileTools "allows you to reduce the latency of the first execution of Julia code" and allows "dramatic reduction in 'time to first execution' (TTFX)" on Julia ≥ 1.9 for package developers. Maintenance: v1.3.4 (2026-05-12). Source: [PrecompileTools.jl](https://github.com/JuliaLang/PrecompileTools.jl). The Julia manual has a section on "Execution latency, package loading and package precompiling time" ([Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)).
- **[Repo]** Maintenance of the remaining tooling on 2026-10-12:

| Tool | Latest release | Activity in last 12 months |
|---|---|---|
| [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) | v1.19.0 (2026-09-02) | 154 commits |
| [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) | v2.14.0 (2026-09-06) | 184 commits |
| [Runic.jl](https://github.com/fredrikekre/Runic.jl) | v1.13.0 (2026-10-06) | 43 commits; used by Ferrite via pre-commit |
| [invenia/BlueStyle](https://github.com/invenia/BlueStyle) (style guide) | none | last commit 2026-02-19 |
| [setup-julia](https://github.com/julia-actions/setup-julia) | v3 (2026-05-10) | last commit 2026-10-07 |
| [julia-actions/cache](https://github.com/julia-actions/cache) | v3.5.0 (2026-10-05) | |
| [julia-runtest](https://github.com/julia-actions/julia-runtest) | v1.11 / v1 tags (2026-07-02) | last commit 2026-10-01 |
| [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) | v0.6.7 (2026-09-16) | 32 commits |
| [Plots.jl](https://github.com/JuliaPlots/Plots.jl) | Plots v1.41.7 (2026-08-20) | 98 commits |
| [Makie.jl](https://github.com/MakieOrg/Makie.jl) | v0.24.15 (2026-09-19) | 208 commits |
| [Tullio.jl](https://github.com/mcabbott/Tullio.jl) | v0.3.9 (2026-02-12) | 4 commits |

- **[Repo]** NEMOpt currently has no `Project.toml`/`Manifest.toml` and no tests. Its `using` lines load Nonconvex, ChainRulesCore, Dates, Random, Distributions, SparseArrays, LinearAlgebra, Plots, Colors, Images, JLD2, DataFrames, Printf, DelimitedFiles, Tullio, FFMPEG and PlotlyJS. Source: /home/user/NEMOpt/*/NEMOpt*.jl.

### Inferences
- Suggested layout:

  ```
  NEMOpt/
    Project.toml (+ [compat], [weakdeps]/[extensions] for Plots/Makie/FFMPEG)
    Manifest.toml (committed for the paper environment; or under reproduce/FEAD2026/)
    src/NEMOpt.jl
    src/{mesh,constraints,material,assembly,newton,homogenization,sensitivity,filter,objective,optimizer,io}.jl
    ext/NEMOptPlotsExt.jl
    scripts/run.jl (reads config/*.toml)
    config/{full_domain,domain_reduction}.toml
    test/{runtests.jl, constitutive, element, homogenization, sensitivity_fd, regression}
    docs/
    reproduce/
  ```

- Move heavy optional dependencies (Plots, PlotlyJS, Images, FFMPEG, Colors, DataFrames) out of the core into extensions or scripts. The core then loads fast and precompiles cleanly. Distributions and Random should be needed only if the initial design is random. In that case, seed it explicitly through the config (Taschuk & Wilson Rule 10, see Q8).
- TOML (Julia's stdlib `TOML`) config files plus an immutable `Params` struct replace the hard-coded `Lx, Ly, Nelx, Nely, N_ITER, p/β schedules` and the comment-toggled objectives. This directly implements "Do not comment and uncomment sections of code to control a program's behavior" (Wilson et al. 2017, Q8).
- For golden-file regression, ReferenceTests' low activity argues for a thin in-house helper: compare against stored JLD2 or CSV arrays with explicit `rtol`/`atol`. ReferenceTests remains usable because its API is stable.
- CI matrix: Julia 1.13 (stable), plus 1.10 (LTS) only if supported. Jobs: Aqua, JET `report_package`/`test_opt` on kernels, a formatter check, fast unit tests on every push, and the long regression run on a schedule or by manual dispatch.

### Gaps
- I did not verify current julia-actions recommended workflow templates (for example, whether `julia-buildpkg` is still advised alongside `julia-runtest`).
- I did not check the JuliaFormatter vs Runic community-adoption trend beyond Ferrite's choice.

---

## Q7. Multithreading options for FE assembly and typical speedups

### Takeaway
Ferrite's officially documented options are deterministic grid coloring with task-local scratch buffers (via OhMyThreads + TaskLocalValues) and, since v1.6, atomic assembly. The atomic path is faster to set up but makes summation order non-deterministic. Polyester offers lower-overhead static scheduling with caveats. For NEMOpt's small 2D meshes (for example 80×80), threading assembly will likely matter less than removing allocations and reusing factorizations. Coloring is required to keep runs reproducible.

### Cited Findings
- **[Doc]** Ferrite's threaded-assembly howto lists the strategies and their trade-offs:
  - **Locking:** "can lead to lock contention … results will not be deterministic since floating point operations are neither associative nor commutative."
  - **Assembler task:** has task-switching overhead and is also non-deterministic.
  - **Grid coloring:** within each color no two elements share DOFs. It "gives predictable results because … we will always add the red first, then the blue, and finally the green."
  - **Atomic accumulation:** "results are not deterministic since the order of the additions depend on the task scheduling."

  The howto implements coloring with `OhMyThreads.@tasks` and `DynamicScheduler(; ntasks)` (or `GreedyScheduler` with chunking) plus TaskLocalValues for scratch data. Source: [Ferrite threaded assembly howto](https://ferrite-fem.github.io/Ferrite.jl/stable/howto/threaded_assembly/).
- **[Doc]** Ferrite v1.7.0 notes: `create_coloring` speedups of "3-10x serial, and another 2-3x with 8 threads" for 200k–400k-cell grids. The "resulting coloring is independent of the number of threads." Sparsity-pattern construction is 3–8× faster and multithreaded. Source: [Ferrite CHANGELOG](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md).
- **[Doc]** Polyester "provides low-overhead multithreading … `@batch` … implements static scheduling (c.f. `@threads :static`) and has minimal overhead because it manages and re-uses a dedicated set of Julia tasks." Caveats: `@batch` turns arrays into `StrideArraysCore.PtrArray`s, so "slices will create `view`s by default(!)", and threads are not pinned (use ThreadPinning.jl). Maintenance: v0.7.19 (2026-02-10), only 2 commits in 12 months. Source: [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl).
- **[Repo]** OhMyThreads.jl v0.8.7 (2026-10-02), 15 commits in 12 months ([repo](https://github.com/JuliaFolds2/OhMyThreads.jl)).
- **[Comm]** A community report: LowLevelFEM.jl threaded assembly on about 290,000 tetrahedra took 6.99 s on 1 thread and 4.61 s on 2 threads (about 1.5×) on a 2-core/4-thread laptop, with little gain from hyper-threading. Source: [Julia Discourse, LowLevelFEM v1.13.28 announcement](https://discourse.julialang.org/t/ann-lowlevelfem-jl-v1-13-28-help-test-multithreaded-fem-assembly/138710). Anecdotal and not peer-reviewed.
- **[Doc]** The Julia manual has a section on "Multithreading and linear algebra" (interaction between Julia threads and BLAS threads). Source: [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra).

### Inferences
- For NEMOpt, profile first. With 6,400 Q4 elements, after the constitutive fix (about 500× [Ctx]), assembly will likely be a small fraction of the time, and sparse factorizations will dominate:
  - the Newton tangent at each load step
  - the four χ (corrector) solves
  - the adjoint solves

  The biggest wins are probably reusing the symbolic factorization (same sparsity pattern), solving the four χ load cases as a multi-RHS solve against one factorization, and allocation-free element loops. Threading comes after that.
- If threading is used, prefer Ferrite-style coloring for deterministic, thread-count-independent results. Avoid atomic or lock-based accumulation in reproduction runs. Fix BLAS thread counts and record them.
- Parallelizing over independent work such as the four macroscopic load cases, or finite-difference gradient checks, is simpler and deterministic.

### Gaps
- I found no authoritative, peer-reviewed speedup figures for Ferrite or Polyester threaded assembly on 2D meshes of NEMOpt's size.
- I did not check whether SuiteSparse (UMFPACK/CHOLMOD) results are bitwise reproducible across thread counts or BLAS builds.

---

## Q8. Research software engineering guidance

### Takeaway
The RSE literature maps directly onto NEMOpt's gaps:
- version control and versioned releases
- explicit dependencies (Project/Manifest)
- configuration instead of comment toggles
- no duplication (merge FullCell and HalfCell)
- a small test set
- identical results for identical inputs (seeds, pinned environment, deterministic assembly)
- documentation with a README, examples and citation info
- a license and DOI archiving

### Cited Findings
- **[PR]** Wilson et al., "Best Practices for Scientific Computing," *PLoS Biology* 12(1):e1001745, 2014, doi:[10.1371/journal.pbio.1001745](https://doi.org/10.1371/journal.pbio.1001745). Practices:
  1. Write programs for people, not computers.
  2. Let the computer do the work.
  3. Make incremental changes.
  4. Don't repeat yourself (or others).
  5. Plan for mistakes.
  6. Optimize software only after it works correctly.
  7. Document design and purpose, not mechanics.
  8. Collaborate.
- **[PR]** Wilson, Bryan, Cranston, Kitzes, Nederbragt & Teal, "Good enough practices in scientific computing," *PLOS Comput. Biol.* 13(6):e1005510, 2017, doi:[10.1371/journal.pcbi.1005510](https://doi.org/10.1371/journal.pcbi.1005510). Software practices:
  - Place a brief explanatory comment at the start of every program.
  - Decompose programs into functions.
  - Be ruthless about eliminating duplication.
  - Always search for well-maintained software libraries that do what you need.
  - Test libraries before relying on them.
  - Give functions and variables meaningful names.
  - Make dependencies and requirements explicit.
  - **Do not comment and uncomment sections of code to control a program's behavior.**
  - Provide a simple example or test data set.
  - Submit code to a reputable DOI-issuing repository.

  The paper also recommends a project layout (`src/`, `data/`, `results/`, `doc/`, `bin/`) and version control.
- **[PR]** Taschuk & Wilson, "Ten simple rules for making research software more robust," *PLOS Comput. Biol.* 13(4):e1005412, 2017, doi:[10.1371/journal.pcbi.1005412](https://doi.org/10.1371/journal.pcbi.1005412). Rules:
  1. Use version control.
  2. Document your code and usage.
  3. Make common operations easy to control.
  4. Version your releases.
  5. Reuse software (within reason).
  6. Rely on build tools and package managers for installation.
  7. Do not require root or other special privileges.
  8. **Eliminate hard-coded paths.**
  9. **Include a small test set that can be run to ensure the software is actually working.**
  10. **Produce identical results when given identical inputs.**
- **[PR]** Lee, "Ten simple rules for documenting scientific software," *PLOS Comput. Biol.* 14(12):e1006561, 2018, doi:[10.1371/journal.pcbi.1006561](https://doi.org/10.1371/journal.pcbi.1006561). Rules:
  1. Write comments as you code.
  2. Include examples (and lots of them).
  3. Include a quickstart guide.
  4. Include a README file with basic information.
  5. Include a help command for command line interfaces.
  6. Version control your documentation.
  7. Fully document your application programming interface.
  8. Use automated documentation tools.
  9. Write error messages that provide solutions or point to your documentation.
  10. Tell people how to cite your software.
- **[PR]** Jiménez et al., "Four simple recommendations to encourage best practices in research software," *F1000Research* 6:876, 2017, doi:[10.12688/f1000research.11407.1](https://doi.org/10.12688/f1000research.11407.1). Recommendations:
  1. Make source code publicly accessible from day one.
  2. Make software easy to discover by providing software metadata via a popular community registry.
  3. Adopt a licence and comply with the licence of third-party dependencies.
  4. Define clear and transparent contribution, governance and communication processes.
- **[Repo]** NEMOpt's drivers hard-code output paths (for example `writedlm("Homogenized_Stiffnesses/A_ITER_$iter.txt", …)`) and switch objectives by commenting lines out (`#J = 1 / (A[3,3] + A[4,4])`). Source: /home/user/NEMOpt/Full_Domain/NEMOpt_main.jl.

### Inferences
- Map each rule to a concrete refactor item:
  - DRY / "be ruthless about duplication" → merge the two scripts through `CellSymmetry` dispatch (Q2).
  - "Do not comment and uncomment" and "Make common operations easy to control" → TOML config with an `Objective` type (Q6).
  - "Eliminate hard-coded paths" → an output directory taken from the config, for example via DrWatson `datadir()`/`savename`.
  - "Small test set" → a coarse-mesh (8×8 or 16×16) end-to-end regression that runs in seconds.
  - "Identical results for identical inputs" → seeded RNG, committed Manifest, deterministic assembly, and a recorded Julia version and thread counts.
  - "Tell people how to cite" → `CITATION.cff` with the FEAD DOI, plus Zenodo archiving of the tagged release that reproduces the paper.
  - License → choose one; also check GPL implications if any Svanberg-derived MMA code is ever incorporated.
- "Optimize only after it works correctly" argues for building the regression baseline (legacy outputs) *before* performance work.

### Gaps
- None material. These are well-established guidance papers, but they are not specific to Julia or FE codes.

---

## Q9. Unicode identifier pitfalls in Julia (𝜌 U+1D70C vs ρ U+03C1)

### Takeaway
Julia normalizes identifiers with **NFC** plus a small hard-coded table of look-alike mappings. It does **not** apply NFKC compatibility folding. As a result, `𝜌` (U+1D70C, MATHEMATICAL ITALIC SMALL RHO) and `ρ` (U+03C1, typed with `\rho<TAB>`) are **different identifiers**. A contributor who types `\rho` will get an `UndefVarError`, or worse, silently create a new variable. NEMOpt uses 𝜌 throughout (about 270 occurrences, zero ρ), so the refactor should either migrate to ASCII/standard names or enforce 𝜌 with a lint test.

### Cited Findings
- **[Doc]** "Julia identifiers are NFC-normalized." The non-standard equivalences are:
  - `ɛ` (U+025B) and `µ` (U+00B5) are treated as the corresponding Greek letters.
  - the middle dot `·` (U+00B7) and the Greek interpunct `·` (U+0387) are treated as the dot operator `⋅` (U+22C5).
  - the minus sign `−` (U+2212) is treated as hyphen-minus `-`.

  Source: [Julia manual, Variables → Allowed Variable Names](https://docs.julialang.org/en/v1/manual/variables/#man-allowed-variable-names) (read from `doc/src/manual/variables.md`, master 2026-10-12).
- **[Repo]** The parser's `julia_charmap.h` table contains exactly these mappings: U+025B→U+03B5, U+00B5→U+03BC, U+00B7→U+22C5, U+0387→U+22C5, U+2212→U+002D. It also contains **U+210F (ℏ) → U+0127 (ħ)** (#48870), which the manual paragraph does not mention. There is no mapping for mathematical alphanumeric symbols. Source: [src/flisp/julia_charmap.h](https://github.com/JuliaLang/julia/blob/master/src/flisp/julia_charmap.h).
- **[Repo]** In the REPL LaTeX completion table, `\rho` → "ρ" (U+03C1) and `\varrho` → "ϱ". The mathematical italic rho "𝜌" is produced by the `\it` prefix (`italic*"rho"`, i.e. `\itrho<TAB>`). Bold `\bfrho` → 𝛒 and bold-italic `\birho` → 𝝆 are also distinct characters. Source: [stdlib/REPL/src/latex_symbols.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/REPL/src/latex_symbols.jl) (lines with `italic*"rho"`, `const italic = "\\it"`).
- **[Standard]** U+1D70C has decomposition `<font> 03C1`, a compatibility (not canonical) decomposition. So NFC(𝜌) ≠ ρ while NFKC(𝜌) = ρ. Checked with Python `unicodedata` (Unicode 15.1) against the [Unicode Character Database](https://www.unicode.org/Public/UCD/latest/ucd/UnicodeData.txt).
- **[Repo]** Occurrences in NEMOpt: NEMOpt.jl has 98 `𝜌` and 0 `ρ`; NEMOpt-DR.jl has 104 and 0; each driver has 35 and 0. Source: grep over /home/user/NEMOpt.

### Inferences
- Choose one of three options:
  - **(a)** Rename to ASCII in public APIs (`rho`, `rho_filtered`, `rho_phys`), keeping Unicode only in local math-heavy kernels if desired. This is the most robust choice for collaborators and for grep.
  - **(b)** Standardize on plain `ρ` (U+03C1), which is what `\rho` produces.
  - **(c)** Keep `𝜌` but add a test that fails if any source file contains U+03C1, or if any identifier mixes look-alikes. Document `\itrho<TAB>` in CONTRIBUTING.

  A mechanical rename is behavior-preserving, but it must be done atomically across all files, since the two characters are distinct identifiers.
- Similar hazards apply to `β` vs `ϐ`, `ε` vs `ϵ` (U+03F5, `\epsilon`) vs `ɛ` (mapped to ε), and `χ` vs `𝜒`. A lint test can scan for characters in the Mathematical Alphanumeric Symbols block (U+1D400–U+1D7FF) and for Greek look-alikes.
- Formatters (Runic, JuliaFormatter) do not normalize identifiers. This is an assumption I did not verify, so the lint test is the actual safeguard.

### Gaps
- I did not verify whether JET or any linter warns about confusable identifiers.
- I did not check whether Julia ≥ 1.13 emits warnings for mathematical-alphanumeric identifiers (I found none in the manual).

---

## Q10. Synthesis: recommended target architecture and migration path (for the report writer)

### Takeaway
Refactor in two tracks.

**Track A, reproduction-faithful.** Package the current numerics unchanged:
- remove globals by threading state through typed structs
- merge FullCell and HalfCell via dispatch
- add regression tests against legacy outputs
- pin NonconvexMMA and every other dependency in a committed Manifest
- keep the legacy continuation and rounding semantics behind a `legacy=true` mode

**Track B, modernization.** Do this only after Track A's golden tests exist:
- Tensors.jl kernels (stress and tangent via AD; hand tangents tested against AD)
- optionally Ferrite for DOF handling and assembly
- a gradient-only optimizer interface (in-house MMA, NonconvexMMA via `CustomGradFunction`, or NLopt)
- outer-loop continuation
- tolerance-based verification

### Cited Findings
- **[PR]** Type stability and groundedness enable compiler optimization (Pelenitsyn et al. 2021, doi:[10.1145/3485527](https://doi.org/10.1145/3485527)). Untyped globals block optimization ([Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)).
- **[PR]/[Doc]** Dispatch-based extensibility is the documented design of Gridap (Verdugo & Badia 2022, doi:[10.1016/j.cpc.2022.108341](https://doi.org/10.1016/j.cpc.2022.108341)).
- **[Doc]** Ferrite provides periodic and affine constraints, a hyperelasticity tutorial, homogenization tutorials and deterministic colored threading. Released versions reject master-and-slave chains; untangling is unreleased as of 2026-10-12 ([CHANGELOG](https://github.com/Ferrite-FEM/Ferrite.jl/blob/master/CHANGELOG.md)).
- **[Doc]** NonconvexMMA (v1.1.2, 2026-08-08) supports `CustomGradFunction` for hand gradients and has no iteration callback. NLopt's MMA is a different, independent CCSA implementation ([NLopt docs](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/)).
- **[PR]** RSE rules call for identical results for identical inputs, small test sets and no hard-coded paths (Taschuk & Wilson 2017, doi:[10.1371/journal.pcbi.1005412](https://doi.org/10.1371/journal.pcbi.1005412)), and for not controlling behavior by commenting code (Wilson et al. 2017, doi:[10.1371/journal.pcbi.1005510](https://doi.org/10.1371/journal.pcbi.1005510)).

### Inferences
- **Recommended stack (2026-10):**

| Layer | Choice |
|---|---|
| Language | Julia 1.13.x |
| Tensors / constitutive | Tensors.jl (+ ForwardDiff underneath) |
| FE | Ferrite.jl (≥ 1.6; ≥ the release containing #1327 if tangled DR constraints arise), or the retained in-house structured-grid assembly with a master–slave transformation |
| Linear algebra | SparseArrays + SuiteSparse (factorization reuse) |
| Optimizer | NonconvexMMA pinned (Track A); in-house MMA/GCMMA or NLopt `:LD_MMA`/`:LD_CCSAQ` (Track B) |
| I/O | JLD2 + DelimitedFiles/CSV |
| Plotting | Plots or Makie in a package extension |
| Config | TOML (+ optional DrWatson) |
| Testing | Test + TestItemRunner, Aqua, JET, AllocCheck, BenchmarkTools |
| Docs | Documenter |
| Formatter | Runic or JuliaFormatter (BlueStyle) |
| Latency | PrecompileTools |
| CI | GitHub Actions (setup-julia v3, cache v3, julia-runtest) |

- **Migration order:**
  1. Freeze the legacy code and capture golden outputs at fixed iterations: densities, U, A, P, J, gradient, for both variants on a small mesh, plus the paper meshes.
  2. Create the package skeleton (Project/Manifest) and move code into functions with explicit arguments. Golden tests should pass bitwise.
  3. Merge the variants via `CellSymmetry` dispatch. Golden tests should still pass bitwise.
  4. Swap the constitutive kernel for Tensors.jl. Expect bitwise or near-bitwise equality ([Ctx]: identical); test to tolerance.
  5. Restructure the optimizer interface and continuation, keeping a legacy mode.
  6. Optionally adopt Ferrite assembly and constraints, with tolerance-level verification and re-run paper cases for qualitative and quantitative agreement.
  7. Add threading with coloring.
- **Decision trade-off:** adopting Ferrite gives maintainability and features but loses bitwise reproducibility. Keeping in-house assembly on a structured grid gives bit-level continuity, at the cost of maintaining constraint and assembly code. For a 2D structured unit-cell code, the in-house option is defensible if it is refactored into typed, tested functions.

### Gaps
- Whether DR constraints are tangled (master and slave chains) cannot be determined without analyzing the DR mapping in the article and code. This should be handed to the reproduction-baseline and numerical-method researchers.
- The exact number of objective evaluations per MMA02 outer iteration in the legacy runs (which drives the legacy continuation schedule) needs instrumentation of the legacy code.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. Source: [S&S note](../NEMOpt%20supplied%20papers%20verification/svanberg_2007_sigmund_2022.md), C.2.
- **L196.** Svanberg 2007 (§2, p. 3) advises solving again from the latest solution when problem data change. The notes also specify a complete subproblem solver (§5), useful for an in-house implementation.
- **L259, item 2.** For MMA87 the value enters only as a constant (Svanberg 2007, §5), and NonconvexMMA's `KKTCriteria` stopping test does not read it; the rounded ν in the gradient is what changes the path.
