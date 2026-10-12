# Buckling and Stability in Topology Optimization: Formulation Options for a Stability-Aware NEMOpt

Scope: density-based topology optimization (TO) with buckling or stability criteria. The emphasis is on periodic microstructures, finite-strain (hyperelastic) analysis, spurious modes, repeated eigenvalues and NEMOpt-specific interactions. Current as of 2026-10-12.

Evidence labels used below:
- **[PR]**: peer-reviewed, read in full text or an author preprint of the same paper.
- **[PR-abs]**: peer-reviewed, only the abstract or metadata was verified.
- **[2nd]**: claim about paper X taken from paper Y's description of X.
- **[press]**: non-scholarly secondary source.

Bibliographic metadata (authors, journal, volume, pages, DOI) was checked against Crossref for every DOI cited. Where a number comes from an arXiv version, this is stated, because it may differ from the journal version.

NEMOpt facts used here were read from the code in `/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR.jl`, lines 49–56 and 490–517:
- **Interpolation and clip:** `g(ρ) = max(H(ρ)^p, ϵ)` and `g_Deriv = 0` below the clip, where `H` is the tanh projection.
- **Stress blend:** `P_rho = g(ρ)·(H_Steep(ρ^p)·P + (1 − H_Steep(ρ^p))·S_Linear)`, with `S_Linear = ℂ:∇u` (small-strain linear elasticity).
- **Tangent blend:** `A_rho = g·(H_Steep·𝔸 + (1 − H_Steep)·ℂ)`.
- **Third-derivative tensor:** `Z_rho = g·H_Steep·Z`, where `Z = ∂𝔸/∂F` is hand-coded for the compressible neo-Hookean model `P = μF − (μ − λ ln J)F⁻ᵀ`.

---

## Q1. Linear (linearized) buckling TO and its pitfalls: spurious modes, interpolation of stress stiffness, aggregation, repeated eigenvalues

### Takeaway
Linearized buckling TO is mature: the eigenproblem is (K₀ + λK_σ)φ = 0, the stress stiffness K_σ is interpolated without the void floor, KS or p-norm aggregation handles many eigenvalues, and adjoint-consistent sensitivities are required. It has two known artifact classes:
- spurious modes in low-density regions, which separate K/K_σ interpolation cures;
- spurious localized modes in solid regions at stress concentrations, which Ferrari & Sigmund (2023) cure with stress filtering and erosion.

Repeated eigenvalues are non-differentiable. They are handled by Seyranian-type directional derivatives, symmetric functions, or aggregation.

### Cited Findings
**Origins (Neves et al. 1995, 2002).**
- Neves, Rodrigues & Guedes, "Generalized topology design of structures with a buckling load criterion," *Structural Optimization* 10(2):71–78 (1995), doi:10.1007/BF01743533, is the early density-based linear-buckling TO paper [PR-abs] — [DOI](https://doi.org/10.1007/BF01743533).
- Lindgaard & Dahl (2013) attribute to Neves et al. (1995) the cutoff remedy: the stress-stiffness contribution of elements below a density cutoff ρ_cutoff is ignored, using SIMP for K₀ and a cut-off scheme for K_σ. They note that the cutoff gives "abrupt changes" in objective and sensitivities and can cause oscillations [2nd] — [Lindgaard & Dahl postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf).
- Zhang, Khandelwal & Guo (2023) attribute to "Neves et al." two further points: spurious modes were suppressed by ignoring low-density contributions to the stress stiffness, and the non-differentiability of multiple eigenvalues was circumvented with a generalized-gradient concept [2nd] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911). (The page break hid whether their ref. [20] is the 1995 or the 2002 paper; it is not verified here.)
- Neves, Sigmund & Bendsøe, "Topology optimization of periodic microstructures with a penalization of highly localized buckling modes," *IJNME* 54(6):809–834 (2002), doi:10.1002/nme.449 [PR-abs; only title and metadata verified] — [DOI](https://doi.org/10.1002/nme.449). Thomsen et al. describe Neves et al. as the first to maximize the linearized buckling load of periodic microstructures by homogenization, for modes with the same periodicity as the cell. Later Bloch-wave work extended this to other wavelengths [2nd, via a thesis summary in search results] — [Figueiredo summary](https://fenix.tecnico.ulisboa.pt/downloadFile/1126295043838783/Resumo_MarianaFigueiredo.pdf).

**Bendsøe & Sigmund's smooth variant.**
- Lindgaard & Dahl reproduce the smooth version from Bendsøe & Sigmund (2003, *Topology Optimization: Theory, Methods and Applications*, Springer) [2nd] — [postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf):
  - E_eff(K₀) = [ρ_min + (1 − ρ_min)ρ_e^p]·E
  - E_eff(K_σ) = ρ_e^p·E

  The stress stiffness has no floor, so void elements contribute (almost) no geometric softening. A ρ_min that is too large makes the voids unrealistically stiff.

**Gao & Ma (2015).**
- Gao & Ma, "Topology optimization of continuum structures under buckling constraints," *Computers & Structures* 157:142–152 (2015), doi:10.1016/j.compstruc.2015.05.020 [PR-abs] — [DOI](https://doi.org/10.1016/j.compstruc.2015.05.020).
  - Zhang et al. (2023) say Gao & Ma "avoided spurious buckling modes by the eigenvalue shift and mode identification via computing modal strain energy ratio" [2nd] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).
  - Ferrari & Sigmund (2019) cite Gao & Ma for a strain-energy-density criterion that filters out modes deforming mainly low-density regions [2nd] — [arXiv 1809.05300](https://ar5iv.arxiv.org/html/1809.05300).

**Ferrari & Sigmund (2019), "Revisiting topology optimization with buckling constraints."** *SMO* 59(5):1401–1415, doi:10.1007/s00158-019-02253-3. Details below are from arXiv v2 (Sept 2018) [PR] — [arXiv](https://ar5iv.arxiv.org/html/1809.05300); [DOI](https://doi.org/10.1007/s00158-019-02253-3).
- **Interpolation:**
  - Elastic: h₁(x) = E₀ + x^p(E₁ − E₀) with E₀ = 10⁻⁶.
  - Stress stiffness: h₂(x) = x^p·E₁, with no E₀ term.
  - Penalty continuation: p from 1 to 6 in steps of 0.25 every 25 iterations.
  - The authors report this interpolation has "proven effective" against artificial modes (citing Gao & Ma 2015 and Thomsen et al. 2018) and that they never met such modes in their tests.
- **Aggregation:**
  - p-norm: M_ρ = (Σ(−μ_i)^ρ)^(1/ρ).
  - KS: M_KS = μ₁ + (1/ρ)ln Σ e^{ρ(μ_i − μ₁)}.
  - A scale factor s = μ₁/M[μ_i] is refreshed when p changes.
  - Recommended parameters: KS with ρ > 20, and p-norm with ρ ≥ 8. Doubling-ρ continuation works, but starting too low (ρ = 16) gave poor designs.
- **Eigenpairs:** 24 computed and the lowest 12 constrained on a 90×210 mesh; 32 computed and 24 constrained on a 180×420 mesh.
- **Repeated eigenvalues:**
  - Individual eigenvalue constraints use a separation factor α = 1.01.
  - Aggregated gradients are argued to be unique even for repeated eigenvalues, citing Gravesen et al. (2011) and Torii & de Faria (2017).
- **Sensitivity and the adjoint term:**
  - Formula: dλ_i/dx_e = φ_iᵀ(∂K/∂x_e + λ_i ∂K_σ/∂x_e)φ_i − λ_i vᵀ(∂K/∂x_e)u₀, with adjoint K v = φ_iᵀ[∇_u K_σ]φ_i.
  - Dropping the adjoint ("frequency-like" term only) gives "inconsistent sensitivities." Compliance worsened by about 0.4%, 1% and 6% at P̄c = 0.5, 0.75 and 1.0, and the optimization broke down by iteration 150 at P̄c = 2.0.
  - The dropped term is comparable in magnitude to the kept one. Alone, the frequency-like term would remove material from highly compressed regions.
- **Elements:**
  - On a 10×4 column mesh, Q4 overestimates the critical load by about 50%.
  - The incompatible-modes Q6 element is much more accurate on coarse meshes and converges from below. It is "generally recommended" but not essential.
- **Conclusions:** stiffness and stability compete strongly; many buckling modes coalesce at demanding constraints; the cost of eigenanalysis is the open issue.

**Ferrari, Sigmund & Guest (2021), 250-line Matlab code.** *SMO* 63(6):3045–3066, doi:10.1007/s00158-021-02854-x [PR] — [arXiv 2101.02973](https://arxiv.org/abs/2101.02973); [DOI](https://doi.org/10.1007/s00158-021-02854-x).
- Interpolation: E_K(x̂) = E_min + (E₀ − E_min)x̂^{p_K} and E_G(x̂) = E₀x̂^{p_G}.
- Formulation:
  - Eigenproblem: (G[x̂,u] + μK[x̂])φ = 0, with μ₁ = 1/λ₁.
  - KS aggregation of BLFs.
  - Updates: MMA-like approximations with an OC-like update.
  - Adjoint: w_i = K⁻¹[φ_iᵀ(∇_u G)φ_i].
  - The paper explicitly notes non-differentiability where μ_i are repeated (Seyranian et al. 1994).
- Efficiency: a vectorized K_σ (G) assembly and sensitivity removes the setup bottlenecks, so "problems of an interesting size" run on a laptop.

**Ferrari & Sigmund (2020), multilevel approach.** "Towards solving large-scale topology optimization problems with buckling constraints at the cost of linear analyses," *CMAME* 363:112911, doi:10.1016/j.cma.2020.112911 [PR-abs] — [arXiv 1909.13585](https://arxiv.org/abs/1909.13585).
- Preconditioned iterative solvers plus buckling modes taken from a coarse discretization.
- Computing modes on the coarse scale also alleviates artifacts from local stress concentrations and filters out less important local modes.

**Ferrari & Sigmund (2023), spurious modes in solid regions.** "A strategy for avoiding spurious localized buckling modes in topology optimization," *IJNME* 124(18):4118–4140, doi:10.1002/nme.7309 [PR-abs] — [DOI](https://doi.org/10.1002/nme.7309).
- Localized spurious modes also appear in solid regions at stress concentrations. The paper attributes these to inherent limitations of linearized buckling analysis for complex stress states, a case "hardly addressed" in the TO literature.
- Remedy: filtering and erosion applied to the stress field. Designs are verified by nonlinear analysis.

**Seyranian, Lund & Olhoff (1994).** "Multiple eigenvalues in structural optimization problems," *Structural Optimization* 8(4):207–227, doi:10.1007/BF01742705 [PR-abs] — [DOI](https://doi.org/10.1007/BF01742705).
- This is the standard reference for directional derivatives of N-fold eigenvalues (perturbation and sub-eigenproblem). Zhang et al. (2023) and Ferrari et al. (2021) use or cite it [2nd] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).

**Other cited contributions (via Zhang et al. 2023's reference list) [2nd]** — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911):
- Bruyneel, Colson & Remouchamps, *SMO* 35(2):181–186 (2008), doi:10.1007/s00158-007-0129-z. Enough representative buckling modes must be included to handle mode switching; otherwise convergence is poor.
- Du & Olhoff, *SMO* 34(2):91–110 (2007), doi:10.1007/s00158-007-0101-y. They used higher penalization of mass in low-density elements to suppress spurious vibration modes, and a two-loop procedure whose inner loop uses directional derivatives of multiple eigenvalues.

### Inferences
- In linear buckling, the "different interpolation for K₀ and K_σ" remedy works because the static problem and the eigenproblem are only weakly coupled (K_σ is linear in the prebuckling stress). Lindgaard & Dahl make this argument explicitly (see Q3). The remedy therefore does not transfer directly to a consistent finite-strain Newton solver.
- The consistent adjoint term (the ∇_u K_σ contraction) is not optional. Its finite-strain analogue is exactly the term that needs NEMOpt's third-derivative tensor Z (see Q3).
- Spurious solid-region localization at stress concentrations (Ferrari & Sigmund 2023) is likely relevant to periodic cells with sharp re-entrant corners, for example auxetic or re-entrant designs. This is untested for NEMOpt.

### Gaps
- Full text of Neves et al. (1995, 2002), Gao & Ma (2015) and Ferrari & Sigmund (2023) could not be retrieved (paywalls or 403). Their exact formulas are reported only second-hand. Details of the 2002 "penalization of highly localized buckling modes" (how localization is measured and penalized) remain unverified.
- The Bendsøe & Sigmund (2003) book was not consulted directly; its formula is quoted via Lindgaard & Dahl.

---

## Q2. Periodic materials: linearized bifurcation with Floquet–Bloch (Thomsen et al. 2018; Wang & Sigmund 2021, 2023; Christensen et al. 2023)

### Takeaway
State-of-the-art buckling-strength design of periodic materials (the DTU line) uses small-strain homogenization for the prebuckling stress, then a linearized eigenproblem [K₀(k) + λK_σ(k)]φ = 0 with Floquet–Bloch boundary conditions, sampled at a handful of wave vectors (IBZ vertices plus one small-|k| point for long-wave modes). KS aggregation over (k, mode) gives the objective. The results are hierarchical microstructures that trade stiffness for large strength gains.

This line is entirely linearized. Finite-strain, Bloch-aware TO of microstructures appears not to exist as a constraint formulation. The finite-strain Bloch machinery does exist, though: for analysis (Zhang, Feng & Khandelwal 2021) and for prestrained phononic TO (Dalklint et al. 2022).

### Cited Findings
**Thomsen, Wang & Sigmund (2018).** "Buckling strength topology optimization of 2D periodic materials based on linearized bifurcation analysis," *CMAME* 339:115–136, doi:10.1016/j.cma.2018.04.031 [PR-abs] — [DTU Orbit](https://orbit.dtu.dk/en/publications/buckling-strength-topology-optimization-of-2d-periodic-materials-/); [DOI](https://doi.org/10.1016/j.cma.2018.04.031).
- Method: homogenization plus a linearized stability criterion plus Floquet–Bloch theory, maximizing the buckling strength of the most critical failure mode.
- Cases: square, triangular and hexagonal cells under biaxial, uniaxial and shear loading, with failure surfaces computed in macroscopic stress space.
- Results: optimized designs are first-order hierarchical microstructures with substantially higher strength but lower stiffness. Geometric symmetry constraints can be used to control the shape of the failure surfaces.

**Wang & Sigmund (2021).** "3D architected isotropic materials with tunable stiffness and buckling strength," *JMPS* 152:104415, doi:10.1016/j.jmps.2021.104415 [PR] — [arXiv 2012.01359](https://arxiv.org/abs/2012.01359).
- Analysis:
  - LBA with Bloch–Floquet BCs captures both short- and long-wavelength modes.
  - It is solved as the largest τ_h = 1/λ_h.
  - Elements: 64³ incompatible H11 elements.
  - Interpolation: E₀ = 10⁻⁴E₁ for K₀, while K_σ uses the separate ρ^p interpolation "to suppress the spurious buckling modes."
- Wave-vector sampling:
  - Hydrostatic compression: five k-vectors, the four IBZ vertices plus k = [π/20, 0, 0] "for capturing global shear buckling modes."
  - Uniaxial compression: seven k-vectors, six IBZ vertices plus the same small k.
- Optimization settings:
  - Mode count: per k_l, the number of eigenvalues kept is those within 5% of max τ_h, updated every 20 iterations.
  - KS exponent: ζ = 100/max τ_h, updated every 100 iterations.
  - Projection continuation: β₁ doubled every 40 iterations up to 50.
  - Implementation: PETSc/SLEPc, MMA.
- Results versus the stiffness-optimal isotropic plate material: Young's modulus 79%–58%, buckling strength 180%–767%. Buckling strength "shows a weak dependence on the wave vector."

**Wang & Sigmund, "Architecting materials for extremal stiffness, yield and buckling strength."** *Programmable Materials* 1:e6 (2023), doi:10.1017/pma.2023.5. This is the journal version of arXiv 2210.00003 (v1 29 Sep 2022, v2 14 Apr 2023). Details are from arXiv v2 [PR] — [arXiv](https://arxiv.org/abs/2210.00003); [DOI](https://doi.org/10.1017/pma.2023.5).
- **Formulation:** the material buckling strength is σ_c = min over h and k of λ_h, with k ∈ [−π, π]². Sweeping k along the IBZ boundary is customary, and the IBZ is set by the shared symmetries of the geometry and the macroscopic stress state.
- **Interpolation:** E = ρ̄^p(E₁ − E₀) + E₀ for K₀; E = ρ̄^pE₁ for K_σ; ε-relaxed von Mises stress. E₀ = 10⁻⁵E₁ and p = 3.
- **Objective:** a single KS function aggregates elemental von Mises stresses and/or eigenvalues τ_h(k_l). The weighted objective is γ₁·KS + (1 − γ₁)·Ē⁻¹, with stiffness and yield constraints.
- **Sensitivity:** an adjoint K₀ψ_α = Σ_e φ_eᴴ(∂K_σᵉ/∂χ_α)φ_e. Aggregation gives "uniqueness of the gradient of eigenvalues, even when eigenvalues are repeated."
- **Robustness and meshes:**
  - One-case robust formulation with eroded/dilated designs, Δη = 0.05.
  - 512×512 Q6 elements at f* = 0.2, with resolution doubled for f* = 0.1 and 0.05.
  - 45° symmetry imposed.
  - Body-fitted COMSOL verification.
- **Findings:**
  - Designs shift from simple lattices (yield-dominated) to hierarchical lattices (buckling-dominated).
  - At low volume fraction buckling caps strength.
  - The strength–density exponent falls from about 2.3 (buckling-dominated TPU) to 1 (yield-dominated steel).
  - The paper cites prior 2D experimental verification of buckling superiority.
  - In the shown example (f = 20%), the critical mode is a global shear mode with σ_c/E₁ = 0.00060 and Ē/E₁ = 0.1074.

**Christensen, Wang & Sigmund (2023).** "Topology optimization of multiscale structures considering local and global buckling response," *CMAME* 408:115969, doi:10.1016/j.cma.2023.115969 [PR-abs] — [arXiv 2210.11477](https://arxiv.org/abs/2210.11477).
- Bloch–Floquet cell analyses give local buckling estimates. These are turned into density-dependent Willam–Warnke "yield" surfaces (local-buckling stress constraints) and combined with a global buckling criterion.
- De-homogenized structures confirm validity, with "huge structural gains as well as time savings" compared with single-scale approaches.

**Finite-strain Bloch machinery (analysis, not TO).** Zhang, Feng & Khandelwal, "A computational framework for homogenization and multiscale stability analyses of nonlinear periodic materials," *IJNME* 122(22):6527–6575 (2021), doi:10.1002/nme.6802 [PR-abs] — [arXiv 2010.02371](https://arxiv.org/abs/2010.02371).
- Homogenization results "lose their validity at the onset of 1st bifurcation."
- Microscale structural instability is computed by Bloch-wave analysis; macroscale material instability by rank-1 convexity checks on the homogenized tangent.
- Implementation details: selecting the wave-vector space, retrieving a real-valued mode from a complex Bloch wave, and three solution methods for the constrained eigenproblem (two condensation methods and a null-space projection).
- Long-wavelength microscopic buckling "can be equivalently detected by the loss of rank-1 convexity of the homogenized tangent moduli."

**Post-hoc Bloch check of finite-strain TO designs.** Zhang & Khandelwal, "Computational design of finite strain auxetic metamaterials via topology optimization and nonlinear homogenization," *CMAME* 356:490–527 (2019), doi:10.1016/j.cma.2019.07.027 [PR-abs] — [arXiv 1905.00482](https://arxiv.org/abs/1905.00482).
- Optimized finite-strain auxetic designs were checked afterward by Bloch-wave analysis and rank-one convexity.
- "Short and/or long wavelength instabilities can occur during the loading process, leading to a change of periodicity of the microstructure, which can affect the performance of an optimized design."

**Prestrained Bloch analysis inside TO.** Dalklint, Wallin, Bertoldi & Tortorelli, "Tunable phononic bandgap materials designed via topology optimization," *JMPS* 163:104849 (2022), doi:10.1016/j.jmps.2022.104849 [PR-abs] — [DOI](https://doi.org/10.1016/j.jmps.2022.104849).
- A finite macroscopic deformation is applied to the unit cell, and dispersion is computed with Floquet–Bloch BCs on the deformed state.
- Sensitivities use the adjoint method, and degenerate eigenvalues are handled with symmetric polynomials.
- Designs are verified on a conforming mesh.

### Inferences
- Dalklint et al. (2022) already assemble the prestrained tangent with Floquet–Bloch BCs, with sensitivities, in the same Wallin–Tortorelli finite-strain homogenization lineage that NEMOpt follows. A finite-strain Bloch buckling constraint mainly reuses that machinery with the inertia term dropped (or replaced by a pseudo-mass, see Q3). It is a natural novelty gap between the linearized DTU line and the finite-strain Lund/Illinois line.
- Cost scaling for Bloch-constrained designs: each sampled k needs one complex-Hermitian eigen-solve of cell size (several modes). DTU practice uses 5–7 k-points per iteration in 3D. Per-iteration cost is therefore roughly (n_k × n_modes) eigensolves plus one adjoint solve per KS aggregate. This is an inference from the reported sampling; the papers read here gave no wall-clock data.
- Zhang & Khandelwal (2019) is direct evidence that unconstrained finite-strain microstructure TO (the same problem class as NEMOpt's auxetic and Poisson's-ratio designs) can produce designs that bifurcate before the end of the prescribed macroscopic path. The homogenized 𝔸_M used in NEMOpt's objective is then no longer valid (Zhang, Feng & Khandelwal 2021).

### Gaps
- Thomsen et al. (2018) full text could not be retrieved: the DTU Orbit PDF returned 403 or a bot page. The number and placement of k-points during optimization, the failure-surface resolution (number of stress directions), mesh sizes, CPU cost and the quantitative strength gains are therefore not verified here.
- Neither Wang & Sigmund paper reported wall-clock costs in the sections read.
- No paper was found that imposes Floquet–Bloch stability constraints on finite-strain (hyperelastic) homogenized unit cells inside density-based TO. Searches of Semantic Scholar and the web in 2024–2026 returned none. This is an absence of evidence, not proof.

---

## Q3. Finite-strain stability in TO: formulations, eigenvalue computation and sensitivities, spurious-mode control

### Takeaway
The finite-strain literature evaluates the eigenvalues of the tangent stiffness K_T at the converged equilibrium. It does so in one of three ways:
- **Direct critical point:** Kemmler et al. 2005, using an extended system.
- **Extrapolation from a precritical step:** a linearized eigenproblem (Lindgaard & Dahl 2013; Dalklint et al. 2021).
- **Constraint at the target load:** the lowest eigenvalues of K_T at the target load are constrained to stay positive with a margin (Zhang, Khandelwal & Guo 2023).

Every variant needs an adjoint solve whose right-hand side contains ∂(φᵀK_Tφ)/∂u. That term requires ∂𝔸/∂F, the third derivative of the strain energy, which is exactly NEMOpt's Z tensor.

Spurious void modes are controlled by one of three devices:
- **Energy interpolation:** linear kinematics in void, so there is no geometric stiffness there (Wang et al. 2014; Dalklint et al. 2021; Zhang's adaptive version).
- **Pseudo-mass metric:** a pseudo-mass S in K_Tφ = λSφ that down-weights void DOFs (Zhang et al. 2023).
- **Element removal:** removing low-volume-fraction elements (Dalklint et al. 2020, frequency).

Separate interpolation of K_σ, the linear-buckling remedy, is inconsistent at finite strain (Lindgaard & Dahl).

### Cited Findings
**Kemmler, Lipka & Ramm (2005).** "Large deformations and stability in topology optimization," *SMO* 30(6):459–476, doi:10.1007/s00158-005-0534-0 [PR-abs] — [DOI](https://doi.org/10.1007/s00158-005-0534-0).
- Per Zhang et al.: they "exploited the so-called extended system of equations to directly calculate the critical point which was constrained from below to minimize compliance." The difficulty is that the current solution must be close enough to the critical point for convergence [2nd] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).

**Bruns, Sigmund & Tortorelli (2002).** "Numerical methods for the topology optimization of structures that exhibit snap-through," *IJNME* 55(10):1215–1237, doi:10.1002/nme.544 [PR-abs; title and metadata only] — [DOI](https://doi.org/10.1002/nme.544).
- Follow-up: Bruns & Sigmund, "Toward the topology design of mechanisms that exhibit snap-through behavior," *CMAME* 193(36–38):3973–4000 (2004), doi:10.1016/j.cma.2004.02.017 [PR-abs] — [DOI](https://doi.org/10.1016/j.cma.2004.02.017).

**Lindgaard & Dahl (2013).** "On compliance and buckling objective functions in topology optimization of snap-through problems," *SMO* 47(3):409–421, doi:10.1007/s00158-012-0832-2 [PR, postprint] — [postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf).
- **Critical-load approximation:** at a precritical step n, K_σ(D_c, γ_c) ≈ λK_σⁿ, giving (K₀ + K_Lⁿ)φ = −λK_σⁿφ and γ_c ≈ λγⁿ. The path is traced with an arc-length method.
- **Why the linear-buckling remedy fails:** applying separate interpolations to K₀, K_L and K_σ at finite strain makes the residual Q and the tangent K_T = ∂Q/∂D inconsistent. That causes slow or non-convergence, because in nonlinear buckling the static and eigen problems are strongly coupled.
- **Their remedy:** apply a modified SIMP with floor, [ρ_min + (1 − ρ_min)ρ^p]E, to all tangent terms and to the residual. This "seems to eliminate problems with artificial buckling modes" while keeping consistency, at the price of a slightly over-stiff structure.
- **Load level:** the load level at which the compliance objective is evaluated has "a tremendous effect" on the design. If the end load exceeds a limit point during optimization, the final design may be unstable.
- **Linear-buckling comparison:** linear-buckling optimization of the curved beam gave degenerate, gray, disconnected designs unless a strict compliance constraint was added.

**Dalklint, Wallin & Tortorelli (2020), eigenfrequencies.** "Eigenfrequency constrained topology optimization of finite strain hyperelastic structures," *SMO* 61(6):2577–2594, doi:10.1007/s00158-020-02557-9 [PR-abs] — [DOI](https://doi.org/10.1007/s00158-020-02557-9).
- Hyperelasticity with preload: a finite initial deformation followed by small oscillations.
- Adjoint gradients; simple and degenerate eigenfrequencies are both handled.
- "Excessive distortion and spurious eigenmodes in void regions are reduced by removing low volume fraction elements."
- Eigenfrequencies "drastically change with the load magnitude."

**Dalklint, Wallin & Tortorelli (2021), stability.** "Structural stability and artificial buckling modes in topology optimization," *SMO* 64(4):1751–1763, doi:10.1007/s00158-021-03012-z [PR-abs] — [DOI](https://doi.org/10.1007/s00158-021-03012-z).
- A "strain energy transition approach" removes artificial buckling modes.
- Analysis: nonlinear large-deformation hyperelastic simulation, Newton along the fundamental path, critical loads estimated by an eigenvalue analysis.
- Problem: minimize displacement subject to the lowest critical buckling loads and a volume limit, with a Helmholtz filter and MMA.
- Per Zhang et al., Dalklint et al. extend the Lindgaard–Dahl precritical approach by adopting the Wang et al. (2014) linear energy interpolation (to address mesh distortion) and an aggregation function for multiple-eigenvalue sensitivities [2nd] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).

**Wang, Lazarov, Sigmund & Jensen (2014).** "Interpolation scheme for fictitious domain techniques and topology optimization of finite strain elastic problems," *CMAME* 276:453–472, doi:10.1016/j.cma.2014.03.021. This is the original linear/nonlinear energy interpolation that removes void distortion [PR-abs] — [DOI](https://doi.org/10.1016/j.cma.2014.03.021).

**Zhang, Khandelwal & Guo (2023), "Finite strain topology optimization with nonlinear stability constraints."** *CMAME* 413:116119, doi:10.1016/j.cma.2023.116119. arXiv 2303.14911 (27 Mar 2023, 77 pp.), whose authors are Guodong Zhang, Kapil Khandelwal and Tong Guo [PR; details from arXiv v1] — [arXiv](https://arxiv.org/abs/2303.14911); [DOI](https://doi.org/10.1016/j.cma.2023.116119).
- **Stability criterion:**
  - Stable if K_T (with BCs) is positive definite.
  - Eigenproblem: (K_T − λ_q S_A)φ_q = 0 with S_A-orthonormal modes.
  - Classification at the critical point: φ_cᵀP̂ = 0 means a bifurcation point, ≠ 0 a limit point.
- **Formulation:**
  - Minimize end compliance subject to volume and to 1 − λ_q/λ̄ ≤ 0 for the first m distinct eigenvalue clusters of K_T at the final (target) load.
  - m must be large enough for mode switching; m = 6 in most examples.
  - One example had as many as 27 eigenvalues within the first six clusters in early iterations.
- **Pseudo-mass:**
  - S_A is a diagonal pseudo-mass matrix, equal to 1 at solid DOFs and small at void DOFs, which "degenerates into an identity matrix when all the elements are solid."
  - Nodal pseudo-density ϖ_i is a p-norm (q = 15) of the attached element densities.
  - Nodal mass: m̂_i = ε̂ + (1 − ε̂)ϖ_i^{p_s} for ϖ_i ≤ ϖ_l, a C¹ cubic for ϖ_l < ϖ_i < ϖ_u, and 1 for ϖ_i ≥ ϖ_u. The illustrated cut-offs are ϖ_l = 0.1 and ϖ_u = 0.2.
  - The idea follows mass interpolation for spurious vibration modes (Du & Olhoff 2007).
  - It removes spurious modes with "void and low-density elements" while leaving the physical eigenvalues almost unchanged.
- **Adaptive linear energy interpolation:**
  - Kinematics: F = I + η(ρ)∇u with η(ρ) = exp(βρ)/(exp(cβ) + exp(βρ)), β = 120 and c₀ = 0.08.
  - Internal force: ∫ηBᵀP + ∫(1 − η²)B_lᵀ[ℂ:ε]. Tangent: η²Bᵀ𝔸B + (1 − η²)B_lᵀℂB_l.
  - ℂ uses a higher SIMP exponent p_l > p "to discourage using low-density values to exploit small deformation kinematics."
  - c is increased by 0.05 whenever Newton–Raphson fails at minimum step size, then reset after convergence. The authors call this adaptivity "crucial."
- **Simple-eigenvalue sensitivity:**
  - dλ/dρ = ∂(φᵀK_Tφ)/∂ρ − λ∂(φᵀS_Aφ)/∂ρ + ηᵀ∂R/∂ρ.
  - The adjoint η solves ηᵀ(∂R/∂u) = −∂(φᵀK_Tφ)/∂u.
  - Appendix A gives ∂k_T/∂ρ and the needed ∂𝔸/∂F (eq. A3) for the neo-Hookean model.
- **Multiple eigenvalues:**
  - Seyranian et al. (1994) directional derivatives.
  - A two-loop procedure: the outer loop runs MMA on the design when eigenvalues are simple. The inner loop solves a linearized sub-problem for the design increment Δx when clusters are repeated, with extra equality constraints g_sk(Δx) = 0 on off-diagonal terms.
- **Conclusions:**
  - The stability constraints guarantee stability up to the target load. Raising λ̄ raises the critical load at the expense of stiffness.
  - Because of mode switching, the eigenmode of the smallest constrained eigenvalue may differ from the actual first buckling mode, which usually occurs above the target load.

**Han, Sun, Luo & Nishiwaki (2026).** "Material interpolation in finite strain topology optimization: a comprehensive assessment and a generalized energy-strain-combined interpolation scheme," *SMO* 69(5):122, doi:10.1007/s00158-026-04332-8 [PR-abs] — [DOI](https://doi.org/10.1007/s00158-026-04332-8).
- The paper compares parameter-based, energy-based, strain-based and combined interpolations.
- Schemes that look equivalent at small strain "can exhibit fundamentally different characteristics" at finite strain, in sensitivities, void distortion, stability and constitutive generality.
- It proposes an energy-strain-combined interpolation.

### Inferences
- **Third derivative.** For any eigenvalue of the finite-strain K_T, the adjoint right-hand side is ∂(φᵀK_Tφ)/∂u = Σ_e ∫ (∇φ ⊗ ∇φ) : (∂𝔸/∂F) : ∂F/∂u. This is exactly the hand-coded Z contracted twice with the mode gradient. NEMOpt already stores `Z_rho = g·H_Steep·Z`, which correctly gives zero third derivative for the linear branch. The new code is therefore an element routine computing the vector Z:(∇φ⊗∇φ), plus one extra solve with the Newton tangent (already factorized) per eigenvalue or per KS aggregate. This matches the structure of Zhang et al., eq. 24 and Appendix A.
- **Metric choice matters for NEMOpt.** With the plain eigenproblem K_Tφ = λIφ, NEMOpt's void elements have stiffness ε·ℂ. The lowest eigenvalues of K_T will be dominated by soft void-localized modes of order ε, even though they are positive and not instabilities. A constraint λ₁ ≥ λ̄ would then be governed by void modes. There are two ways out:
  - (a) Zhang's pseudo-mass S_A, which inflates the Rayleigh quotient of void-localized modes.
  - (b) A load-factor (buckling) form à la Lindgaard–Dahl/Dalklint, (K_mat + …)φ = −λK_geoφ. Here elements on the linear branch have K_geo = 0 and cannot enter the denominator.
- **Fully nonlinear soft band.** In NEMOpt, elements in the H_Steep ≈ 1 band are fully nonlinear but soft (g small). Their material and geometric stiffness both scale with g, so their local "buckling strain" is independent of g. If they carry large compressive strains, as they would when sandwiched between stiff struts, they can produce localized modes that neither the energy blend nor the clip removes. A pseudo-mass whose cut-off ϖ_u sits above this band would suppress them in metric (a). This is plausible from the cited mechanisms but not demonstrated in any source for this exact blend.
- **Possible material instability at extreme strains.** The NEMOpt neo-Hookean energy includes (λ/2)(ln J)². The second derivative of this term with respect to J is negative for J > e. Very large dilatations in soft nonlinear elements may therefore produce material-level (ellipticity) instabilities unrelated to the structure. This is analytical reasoning, not verified in a source; it is worth a check if extreme strains occur.

### Gaps
- The full text of Dalklint et al. (2021) could not be retrieved: Springer and Lund portal returned 403 or a bot page. Its exact energy-transition formula, the eigenproblem used, and whether it also uses KS or p-norm aggregation are reported only via the abstract and Zhang et al.
- Kemmler et al. (2005) and Bruns et al. (2002) were verified at metadata level only.
- No source found quantifies how many spurious modes arise from soft nonlinear (not void) elements specifically.

---

## Q4. Designs that exploit instabilities (snapping, multistable, auxetic, programmable) and whether to constrain or exploit

### Takeaway
TO that exploits instability exists, but evidence in the material- and unit-cell setting is thinner and newer than for constraining it:
- snap-through structures and mechanisms (Bruns et al. 2002, 2004);
- prestrain-tuned phononics (Dalklint et al. 2022);
- 2026 unit-cell TO with snap-through, bistability and contact (Aveline et al.).

Exploitation requires path-following past limit and bifurcation points, plus enlarged or Bloch-aware cells for period-changing modes. A stability-constrained finite-strain microstructure paper is the lower-risk, clearly novel next step for NEMOpt. Exploitation is a natural follow-up.

### Cited Findings
**Clausen et al. (2015).** Clausen, Wang, Jensen, Sigmund & Lewis, "Topology optimized architectures with programmable Poisson's ratio over large deformations," *Advanced Materials* 27(37):5523–5527, doi:10.1002/adma.201502485 [PR-abs; metadata only] — [DOI](https://doi.org/10.1002/adma.201502485).

**Wang (2018).** "Systematic design of 3D auxetic lattice materials with programmable Poisson's ratio for finite strains," *JMPS* 114:303–318, doi:10.1016/j.jmps.2018.01.013 [PR-abs] — [DTU Orbit](https://orbit.dtu.dk/en/publications/44a4c32d-bebb-4769-a26d-b929678c9c4c).
- Secant Poisson's ratio target, obtained by minimizing average lateral stress under prescribed deformation.
- TO followed by super-ellipsoid shape optimization, giving ν from about −0.78 to 0.00, nearly constant over 0–20% strain.
- The abstract does not mention buckling or stability analysis.

**Zhang & Khandelwal (2019).** Finite-strain auxetic TO, verified afterward by Bloch and rank-one convexity checks. It found short- and long-wavelength instabilities that change the periodicity and can affect performance [PR-abs] — [arXiv 1905.00482](https://arxiv.org/abs/1905.00482).

**Dalklint et al. (2022).** Deformation-tunable phononic bandgaps: the instability-adjacent regime is exploited for function (dispersion tuning), using a prestrained Bloch analysis [PR-abs] — [DOI](https://doi.org/10.1016/j.jmps.2022.104849).

**Aveline, Santer & Hewson (2026).** "A Topology Optimization Framework for the Inverse Design of Nonlinear Mechanical Metamaterials," *Advanced Engineering Materials* 28(9):e202503183, doi:10.1002/adem.202503183 [PR-abs via search snippet; Wiley returned 403] — [Wiley](https://advanced.onlinelibrary.wiley.com/doi/10.1002/adem.202503183).
- Unit cells are synthesized to hit prescribed nonlinear stress–strain points, including internal contact, snap-through buckling and bistability.
- Three designs were validated by mechanical testing.
- Press coverage reports about 8 h per optimization on an Apple M3 laptop and a negative-stiffness region in the bistable design [press] — [3D Printing Industry](https://3dprintingindustry.com/news/imperial-college-london-researchers-develop-topology-optimization-framework-for-nonlinear-mechanical-metamaterials-251514/).

**Lindgaard & Dahl (2013).** If the end load passes a limit point during optimization, compliance-based designs can end up unstable. Evaluating the objective at a fixed load level along a snapping path is therefore delicate [PR] — [postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf).

**Zhang, Feng & Khandelwal (2021).** First-order homogenization loses validity at the first bifurcation [PR-abs] — [arXiv 2010.02371](https://arxiv.org/abs/2010.02371).
- Inference: exploiting post-bifurcation microstructural pattern transformation therefore needs enlarged supercells or second-order or Bloch-informed homogenization, not NEMOpt's current single-cell first-order scheme.

### Inferences
- **Constrain-first is the stronger next paper for NEMOpt.**
  1. NEMOpt designs Poisson's-ratio and stiffness targets under prescribed F_M, exactly where Zhang & Khandelwal (2019) found post-hoc instabilities. The authors themselves list compressive-loading instability as open work.
  2. A stability-constrained formulation reuses Newton, the adjoint and Z.
  3. Finite-strain, Bloch-aware stability constraints for density-based homogenized unit cells appear unaddressed (Q2 gap).
- **Exploitation is a credible second paper.** Examples would be targeting a bifurcation at a prescribed F_M, or designing a desired post-buckling pattern. It would need:
  - arc-length or displacement-controlled path following through critical points;
  - branch switching using the critical mode;
  - supercell or Bloch sampling to capture period-doubling;
  - objectives defined on the post-bifurcated branch.

  The 2026 Aveline et al. work shows it is feasible at unit-cell level. Its specific method (density-based or not, path-following scheme) could not be verified here.

### Gaps
- Clausen et al. (2015): the abstract and full text were unavailable (Wiley 403). Whether its designs rely on, avoid or were checked for instability is not verified.
- Recent (2023–2026) density-based TO that exploits multistability or snapping in periodic cells: only Aveline et al. (2026) was identified, and only at abstract level. Other hits were curved-beam bistable cells or data-driven designs with unclear TO content, so they were not used. A targeted search of *SMO*, *CMAME*, *JMPS*, *Adv. Mater.* and *EML* for "pattern transformation" with "topology optimization" was rate-limited and is incomplete.

---

## Q5. Interaction with NEMOpt's specific choices (clip, linear/nonlinear blend, projection and continuation, half-cell symmetry)

### Takeaway
NEMOpt's linear branch in void already plays the role of the energy-interpolation remedy (Wang 2014; Dalklint 2021; Zhang 2023). Void elements on that branch have no geometric stiffness and cannot buckle by themselves. Three choices would still bias a stability analysis:
- **Identity-metric eigenvalues:** the plain K_Tφ = λφ spectrum is dominated by soft void modes. A pseudo-mass or load-factor form is needed.
- **Soft nonlinear band:** the ρ^p-switched band of soft but fully nonlinear elements can host localized modes.
- **Half-cell reduction:** point-symmetric half-cell reduction can only represent symmetry-preserving, cell-periodic modes. Stability must be evaluated on the full cell, or on both symmetry sectors, plus Bloch k.

Separately, the hard clip with zero gradient is non-smooth and blind to void reinforcement. The sources consistently use a smooth floor, E_min + (1 − E_min)ρ^p.

### Cited Findings
**Smooth floor rather than a clipped floor.**
- Every buckling-TO reference read uses a smooth floored SIMP for the elastic stiffness (E_min + (E₀ − E_min)x^p):
  - Ferrari & Sigmund 2019 [PR] — [arXiv](https://ar5iv.arxiv.org/html/1809.05300).
  - Ferrari, Sigmund & Guest 2021 [PR] — [arXiv](https://arxiv.org/abs/2101.02973).
  - Wang & Sigmund 2023 [PR] — [arXiv](https://arxiv.org/abs/2210.00003).
- Lindgaard & Dahl use [ρ_min + (1 − ρ_min)ρ^p] consistently in the residual and all tangent terms at finite strain [PR] — [postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf).

**The linear kinematics branch.**
- Zhang et al. (2023) use the linear branch with a higher penalization p_l > p, to stop the optimizer from exploiting cheap small-strain kinematics at low density.
- They make the blend parameter adaptive (c raised when Newton fails).
- They add a pseudo-mass with cut-offs (example values 0.1 and 0.2) so that both void and low-density nodes are down-weighted in the stability eigenproblem [PR] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).

**Low-volume-fraction elements.** Dalklint et al. (2020) remove low-volume-fraction elements to reduce distortion and spurious eigenmodes in finite-strain eigenvalue TO [PR-abs] — [DOI](https://doi.org/10.1007/s00158-020-02557-9).

**Continuation practice.**
- Wang & Sigmund double β every 40 iterations (to 50), refresh the eigenvalue count every 20 iterations and the KS exponent every 100. They use an eroded/dilated one-case robust formulation (Δη = 0.05) [PR] — [arXiv 2012.01359](https://arxiv.org/abs/2012.01359); [arXiv 2210.00003](https://arxiv.org/abs/2210.00003).
- Ferrari & Sigmund ramp p from 1 to 6 and rescale the aggregated constraint when p changes [PR] — [arXiv 1809.05300](https://ar5iv.arxiv.org/html/1809.05300).

**Symmetry.**
- The IBZ, and hence which k must be sampled, is fixed by the shared symmetries of the geometry and the macroscopic stress state [PR] — [arXiv 2210.00003](https://arxiv.org/abs/2210.00003).
- Thomsen et al. report that imposed geometric symmetry controls the shape of the failure surfaces [PR-abs] — [DOI](https://doi.org/10.1016/j.cma.2018.04.031).
- Instabilities can change the periodicity of finite-strain TO designs [PR-abs] — [arXiv 1905.00482](https://arxiv.org/abs/1905.00482).

**Repeated eigenvalues are expected in symmetric cells.** Zhang et al. saw up to 27 eigenvalues inside six clusters early in one optimization; Wang & Sigmund rely on KS for gradient uniqueness under repetition [PR] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911); [arXiv 2210.00003](https://arxiv.org/abs/2210.00003).

### Inferences
- **Clip `max(H^p, ϵ)` with zero gradient.**
  - It puts a kink in g at H^p = ϵ. Eigenvalue functions that are already non-smooth at crossings gain a second source of non-smoothness.
  - Void elements get exactly zero sensitivity, so a stability constraint cannot "ask" for material in void to brace a slender strut. Any re-entry must come through filter coupling with gray neighbors.
  - Recommendation: replace it with a smooth floor g = ϵ + (1 − ϵ)H^p, as in all cited sources. Keep ϵ around 10⁻⁶–10⁻⁴; the cited values are E₀/E₁ = 10⁻⁶, 10⁻⁵ and 10⁻⁴.
- **`H_Steep(ρ^p)` blend.**
  - The blend uses the unprojected filtered ρ, raised to p, while g uses the projected H(ρ).
  - Elements become fully nonlinear once ρ^p passes about ω₂, at which point g is already small. This creates the band of soft but fully nonlinear elements. After projection the band sits in the solid–void interface layer, roughly one filter radius wide.
  - Three remedies, none tested on NEMOpt:
    - (i) Cover the band with a pseudo-mass whose upper cut-off ϖ_u exceeds the band.
    - (ii) Shift the H_Steep switch so that nonlinearity is only switched on near the solid level, which means accepting more distortion risk.
    - (iii) Use a higher penalization on the linear branch (Zhang's p_l > p) to discourage designs that rely on the band.
- **Linear branch `S_Linear = ℂ:∇u`.** It contributes a constant tangent ε·ℂ and zero third derivative, which is correct for Z_rho. On its own it cannot create negative eigenvalues. Its non-objectivity (spurious stress under large rotations of void) adds positive stiffness and so cannot create instabilities, but it can mask or delay them slightly.
- **Projection and continuation.** Each β jump changes stiffness distribution, eigenvalues and mode ordering discontinuously. Recommendations:
  - Introduce the stability constraint only after the topology has stabilized (mid-continuation).
  - Re-scale the KS parameter and the constraint at each β step (Ferrari & Sigmund's s factor).
  - Refresh the number of tracked modes periodically (Wang & Sigmund).
  - Thin struts are where buckling governs, so an eroded-design evaluation (robust formulation) is advisable for the stability constraint specifically.
- **Half-cell (point-symmetry) reduction.**
  - For an operator invariant under 180° rotation, eigenmodes split into symmetric and antisymmetric families. A reduced half-cell with symmetric constraints represents only the symmetric family.
  - Bloch modes with k ≠ 0 (period-changing modes such as period doubling) are not representable in any single-cell periodic model.
  - Hence: solve equilibrium on the half-cell if desired, since the prebifurcation state is symmetric. Run the stability eigenproblem on the full cell (or on both symmetric and antisymmetric half-cell problems) and at a few k.
  - This is group-theoretic reasoning applied to NEMOpt; no source here tested domain reduction plus stability.
- **Prescribed F_M path.** NEMOpt is strain-driven, so instability may occur before the target F_M. Lindgaard & Dahl and Zhang et al. both show that the evaluation load level matters. The constraint should therefore cover the whole path, for example as a KS aggregate of λ_min over load increments, or at least be checked at intermediate increments. Otherwise a design can pass at the target F_M while having crossed a bifurcation earlier (a mode-switching risk noted by Zhang et al.).

### Gaps
- No source evaluated the stability of designs produced with a stress-blend (rather than energy-blend) linear/nonlinear interpolation exactly like NEMOpt's `H_Steep(ρ^p)` scheme.
- The effect of point-symmetry domain reduction on buckling-mode completeness in TO was not found in the literature.
- Element technology at finite strain was not researched here. Ferrari & Sigmund recommend Q6 for linear buckling, but the behavior of incompatible or enhanced elements in finite-strain compression (possible hourglassing) should be checked separately.

---

## Q6. Recommended formulation for a stability-aware NEMOpt paper (synthesis)

### Takeaway
Recommended: a stability-constrained finite-strain microstructure paper with a three-tier, KS-aggregated stability measure at (and along) the prescribed F_M path:
- **(T1) Macroscopic:** rank-one convexity of 𝔸_M. This is cheap and reuses existing derivatives.
- **(T2) Microscopic, cell-periodic (k = 0):** lowest eigenvalues of the full-cell K_T in a pseudo-mass metric.
- **(T3) Selected Bloch k-points:** IBZ vertices plus one small k.

Sensitivities are adjoint-consistent and use Z. Repeated eigenvalues are handled by KS (or symmetric-polynomial) aggregation. The clip is replaced by a smooth floor.

### Cited Findings (basis for each element)
- **T1 (macroscopic rank-one convexity).** Long-wavelength microscopic buckling is "equivalently detected by the loss of rank-1 convexity of the homogenized tangent moduli" [PR-abs] — [arXiv 2010.02371](https://arxiv.org/abs/2010.02371). NEMOpt already computes 𝔸_M and its design derivative via Z (code inspection).
- **T2 (pseudo-mass and constraint form).**
  - Pseudo-mass S_A and the constraint form 1 − λ_q/λ̄ ≤ 0 on the first m clusters at the target load, with the adjoint using ∂𝔸/∂F [PR] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911).
  - Alternative: a precritical extrapolated load factor (Lindgaard–Dahl; Dalklint et al. 2021) [PR] — [postprint](https://vbn.aau.dk/ws/files/311111755/LINDGAARD_DAHL_2012_On_compliance_and_buckling_objective_functions_in_topology_optimization_of_snap_through_problems_Postprint.pdf); [DOI](https://doi.org/10.1007/s00158-021-03012-z).
- **T3 (Bloch sampling).**
  - Sample the IBZ vertices plus k = [π/20, 0, 0] for global shear modes, keep the eigenvalues within 5% of the critical one, and aggregate with KS [PR] — [arXiv 2012.01359](https://arxiv.org/abs/2012.01359).
  - The finite-strain prestrained Bloch tangent with adjoint and symmetric polynomials is already demonstrated in TO [PR-abs] — [DOI](https://doi.org/10.1016/j.jmps.2022.104849).
- **Aggregation and repeated eigenvalues.**
  - KS with ρ > 20, or p-norm with ρ ≥ 8, performs as well as direct constraints when the parameter is high enough [PR] — [arXiv 1809.05300](https://ar5iv.arxiv.org/html/1809.05300).
  - Alternatives: Seyranian two-loop directional derivatives [PR] — [arXiv 2303.14911](https://arxiv.org/abs/2303.14911); symmetric polynomials [PR-abs] — [DOI](https://doi.org/10.1016/j.jmps.2022.104849).
- **Consistency.** Dropping the adjoint term yields wrong designs and breakdown at demanding constraints [PR] — [arXiv 1809.05300](https://ar5iv.arxiv.org/html/1809.05300).
- **Validation.** Verify designs by body-fitted (conforming) nonlinear analysis, as done in Wang & Sigmund (COMSOL), Dalklint et al. (2022, conforming mesh) and Ferrari & Sigmund (2023, nonlinear verification) [PR / PR-abs] — [arXiv 2210.00003](https://arxiv.org/abs/2210.00003); [DOI](https://doi.org/10.1016/j.jmps.2022.104849); [DOI](https://doi.org/10.1002/nme.7309).

### Inferences (concrete recommendation for NEMOpt)
1. **Problem statement.**
   - Keep NEMOpt's 𝔸_M-based objectives (Poisson's ratio, bulk and shear moduli).
   - Add g_stab = 1 − KS_ρ{λ_i(F_M^(n), k_l)}/λ̄ ≤ 0. The KS runs over modes i, the k-set {k_l}, and selected load increments n along F(γ) = I + γ(F_M − I).
   - Include T1 as min over unit vectors n of eig(Q(n)) with Q_ik(n) = 𝔸_M,ijkl n_j n_l, sampled over, say, 90–180 directions in 2D and KS-aggregated.
   - Start with T1 + T2 (k = 0, full cell). Add T3 at 3–6 k-points as the novelty.
2. **Eigenproblem.**
   - Use K_T(u; k)φ = λS_Aφ with Zhang's nodal pseudo-mass, built from the projected density with a smooth p-norm and cubic cut-offs. Choose ϖ_u above NEMOpt's `H_Steep` band.
   - Alternatively, use the load-factor form if a "critical F_M multiplier" γ_c is the desired design quantity (a more physical output for failure surfaces in F_M space, analogous to Thomsen's stress-space surfaces).
3. **Sensitivities.**
   - dλ/dρ = φᵀ(∂K_T/∂ρ)φ − λφᵀ(∂S_A/∂ρ)φ + adjoint term, where the adjoint solves K_T η = −Σ_e ∫ Bᵀ[Z_rho : (∇φ ⊗ ∇φ)].
   - This reuses NEMOpt's existing Z and Newton tangent, with one extra solve per aggregated constraint per load increment and per k.
   - For Bloch k the modes are complex; use φᴴ in the contractions.
4. **Interpolation changes.**
   - Make the floor smooth: g = ϵ + (1 − ϵ)H^p.
   - Consider a higher exponent on the linear branch, as in Zhang et al.
   - Optionally make the blend adaptive, as in Zhang et al., for robustness under compression.
5. **Symmetry.** Keep the half-cell only for equilibrium. Do eigenanalysis on the full cell, or on both symmetric and antisymmetric sectors. Do not impose symmetry on the stability modes.
6. **Continuation.**
   - Switch the stability constraint on after the first β stages.
   - Re-scale the KS and λ̄ normalization at each β change.
   - Refresh the tracked mode count periodically.
   - Consider an eroded-design evaluation of the stability constraint.
7. **Verification.**
   - Post-optimization: run full-path finite-strain analysis with a Bloch sweep along the IBZ boundary and a rank-one convexity scan (Zhang, Feng & Khandelwal 2021 protocol) on a body-fitted mesh.
   - Report the first critical F_M and whether it is microscopic (k ≠ 0), cell-periodic, or macroscopic.

These seven points are synthesis or inference from the sources above. No single paper implements this combination.

### Gaps
- Wall-clock cost per optimization iteration for finite-strain stability-constrained TO is not reported in the sections read: not for Zhang et al., Dalklint et al., or Wang & Sigmund. The relative cost of T2 and T3 versus the existing homogenization step is therefore unquantified.
- Whether the Seyranian two-loop scheme or KS aggregation works better with GCMMA in Nonconvex.jl has not been studied.

### Bibliography (verified against Crossref unless noted)
- Aveline C, Santer M, Hewson R (2026). *Adv. Eng. Mater.* 28(9):e202503183. doi:10.1002/adem.202503183
- Bruns TE, Sigmund O (2004). *CMAME* 193(36–38):3973–4000. doi:10.1016/j.cma.2004.02.017
- Bruns TE, Sigmund O, Tortorelli DA (2002). *IJNME* 55(10):1215–1237. doi:10.1002/nme.544
- Bruyneel M, Colson B, Remouchamps A (2008). *SMO* 35(2):181–186. doi:10.1007/s00158-007-0129-z (from Zhang et al. reference list, not Crossref-checked)
- Christensen CF, Wang F, Sigmund O (2023). *CMAME* 408:115969. doi:10.1016/j.cma.2023.115969
- Clausen A, Wang F, Jensen JS, Sigmund O, Lewis JA (2015). *Adv. Mater.* 27(37):5523–5527. doi:10.1002/adma.201502485
- Dalklint A, Wallin M, Tortorelli DA (2020). *SMO* 61(6):2577–2594. doi:10.1007/s00158-020-02557-9
- Dalklint A, Wallin M, Tortorelli DA (2021). *SMO* 64(4):1751–1763. doi:10.1007/s00158-021-03012-z
- Dalklint A, Wallin M, Bertoldi K, Tortorelli D (2022). *JMPS* 163:104849. doi:10.1016/j.jmps.2022.104849
- Du J, Olhoff N (2007). *SMO* 34(2):91–110. doi:10.1007/s00158-007-0101-y (from Zhang et al. reference list, not Crossref-checked)
- Ferrari F, Sigmund O (2019). *SMO* 59(5):1401–1415. doi:10.1007/s00158-019-02253-3
- Ferrari F, Sigmund O (2020). *CMAME* 363:112911. doi:10.1016/j.cma.2020.112911
- Ferrari F, Sigmund O, Guest JK (2021). *SMO* 63(6):3045–3066. doi:10.1007/s00158-021-02854-x
- Ferrari F, Sigmund O (2023). *IJNME* 124(18):4118–4140. doi:10.1002/nme.7309
- Gao X, Ma H (2015). *Comput. Struct.* 157:142–152. doi:10.1016/j.compstruc.2015.05.020
- Han J-K, Sun Z-Y, Luo Y-J, Nishiwaki S (2026). *SMO* 69(5):122. doi:10.1007/s00158-026-04332-8
- Kemmler R, Lipka A, Ramm E (2005). *SMO* 30(6):459–476. doi:10.1007/s00158-005-0534-0
- Lindgaard E, Dahl J (2013). *SMO* 47(3):409–421. doi:10.1007/s00158-012-0832-2
- Neves MM, Rodrigues H, Guedes JM (1995). *Struct. Optim.* 10(2):71–78. doi:10.1007/BF01743533
- Neves MM, Sigmund O, Bendsøe MP (2002). *IJNME* 54(6):809–834. doi:10.1002/nme.449
- Seyranian AP, Lund E, Olhoff N (1994). *Struct. Optim.* 8(4):207–227. doi:10.1007/BF01742705 (volume and issue from Zhang et al. reference list; DOI confirmed via Semantic Scholar)
- Thomsen CR, Wang F, Sigmund O (2018). *CMAME* 339:115–136. doi:10.1016/j.cma.2018.04.031
- Wallin M, Tortorelli DA (2020). *Mech. Mater.* 145:103324. doi:10.1016/j.mechmat.2020.103324
- Wang F (2018). *JMPS* 114:303–318. doi:10.1016/j.jmps.2018.01.013
- Wang F, Lazarov BS, Sigmund O, Jensen JS (2014). *CMAME* 276:453–472. doi:10.1016/j.cma.2014.03.021
- Wang F, Sigmund O (2021). *JMPS* 152:104415. doi:10.1016/j.jmps.2021.104415
- Wang F, Sigmund O (2023). *Programmable Materials* 1:e6. doi:10.1017/pma.2023.5 (arXiv 2210.00003)
- Zhang G, Khandelwal K (2019). *CMAME* 356:490–527. doi:10.1016/j.cma.2019.07.027
- Zhang G, Feng N, Khandelwal K (2021). *IJNME* 122(22):6527–6575. doi:10.1002/nme.6802
- Zhang G, Khandelwal K, Guo T (2023). *CMAME* 413:116119. doi:10.1016/j.cma.2023.116119 (arXiv 2303.14911)
