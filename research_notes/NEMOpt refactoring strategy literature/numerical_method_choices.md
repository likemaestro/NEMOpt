# Numerical-method choices in finite-strain density-based TO of periodic microstructures: what NEMOpt must preserve vs. what to change

Scope: the numerical choices in NEMOpt (Güven & Özdemir, *Finite Elem. Anal. Des.* 260 (2026) 104599, [doi:10.1016/j.finel.2026.104599](https://doi.org/10.1016/j.finel.2026.104599)) that materially affect results. Current as of 2026-10-12.

Source tags used below:
- **[PR]**: peer-reviewed. Its content was read in this session from the publisher abstract, the arXiv full text of the peer-reviewed version, or quoted equations.
- **[PR-bib]**: a peer-reviewed reference whose bibliographic data were verified through Crossref in this session. The content statement comes from prior knowledge and was **not re-verified** in this session, so treat it as unverified.
- **[Non-PR]**: software, documentation, preprint or secondary summary.
- **[Repo]**: NEMOpt source at commit `30d0842`, read in this session.
- **[Article]**: the authors' FEAD paper (the uploaded PDF), read in this session.

Recommendation codes:
- **KEEP**: preserve as the legacy default for reproduction.
- **SWITCH**: keep the current behavior behind a legacy flag and change the new default.
- **DOC**: document only.

Consolidated recommendation table. Details and evidence are in the sections below.

| Choice | Code | Why |
|---|---|---|
| SIMP with clipping `max(H^p, ε)`, ε=5e-4, zero derivative below ε | SWITCH, new default = modified SIMP ε+(1−ε)H^p | The clipping makes f non-C¹ and creates a zero-sensitivity dead zone. Modified SIMP is the community standard. |
| Linear/nonlinear stress blend, with Ĥ evaluated on ρ̃^p (β̂=200, ω̂=0.05) | SWITCH, new option = switch on physical (projected) density as in Wang et al. 2014 | The current switch leaves a band of soft but still hyperelastic elements at the end of continuation. |
| p/β continuation counted per function evaluation, applied inside the objective before the gradient is computed | SWITCH, new default = per outer iteration, applied between iterations | The gradient is inconsistent at continuation steps, and the schedule depends on the optimizer. |
| MMA87 (full domain) vs MMA02/GCMMA (half domain), s_init=0.25, maxinner=40 | KEEP per script, and make them explicit config. DOC the deviation from the "default parameters" | The scripts differ in the optimizer, so the full-vs-reduced comparisons are confounded. |
| Objective and ν rounded to 3 decimals before reaching MMA | SWITCH, new default = unrounded | The rounding error (≈5e-4) is far larger than GCMMA's conservativeness tolerance (1e-7). |
| ν = 2𝔸1122/(𝔸1111+𝔸2222) from the nominal tangent | KEEP + DOC | It is a tangent, in-plane, plane-strain measure, strictly inside (−1, 1) while the normal block is positive definite. |
| Bulk ¼(𝔸1111+𝔸2222+𝔸1122+𝔸2211), shear 𝔸1212+𝔸2121 | KEEP + DOC. Fix the inconsistent logged `G` in the DR script | These are tangent quantities, and the shear sum includes initial-stress terms. |
| Helmholtz filter, periodic BCs, `r` = 1 element (r²∇² form) | KEEP + DOC the convention (≈ cone radius 2√3·r) | The article text and its tables disagree on r. |
| Point-symmetry half-domain reduction plus quadrant-mirrored design | KEEP. Preserve the parity rules of each field when refactoring | It is valid for any F_M, including shear, as long as the solution stays on the symmetric branch. |

---

## 1. SIMP variants and the consequences of the `max(H(ρ̃)^p, ε)` clipping

### Takeaway
NEMOpt's Eq. 29 clipping, f = max(H^p, ε) with f′ = 0 below ε, is neither classic SIMP nor modified SIMP. Classic SIMP keeps a lower bound on the design variable. Modified SIMP, E_min + ρ^p(E0−E_min), is smooth. The clipping introduces a derivative jump and a zero-sensitivity zone. Results depend on it, so keep it under a legacy flag for reproduction and make smooth modified SIMP the new default.

### Cited Findings
- NEMOpt defines `g(ρ) = max(H(ρ)^p, ϵ)` and `g_Deriv(ρ) = H(ρ)^p > ϵ ? p*H^(p-1)*H′ : 0.0` with ϵ = 5e-4 and p starting at 1. [Repo] [Full_Domain/NEMOpt.jl L59–60](/home/user/NEMOpt/Full_Domain/NEMOpt.jl), [NEMOpt_main.jl L19](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl). The half-domain code is identical at [NEMOpt-DR.jl L55–56](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR.jl).
- The article states Eq. 29, f(ρ̃) = max((H_{α,β}(ρ̃))^p, ε), "to suppress 'zero stiffness' response that would lead to an ill-posed problem". It also writes the box constraint as ρ_min ≤ ρ ≤ 1 (Eqs. 53, 55). [Article] [Güven & Özdemir 2026, Sec. 2.5, 5.1–5.2](https://doi.org/10.1016/j.finel.2026.104599)
- The code's box constraint, however, is [0, 1]: `addvar!(model, zeros(...), ones(...))`. [Repo] [NEMOpt_main.jl L229](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR_main.jl L228](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl)
- Modified SIMP, E(x̂_e) = E_min + x̂_e^p (E0 − E_min), applied to the projected (physical) variable, is the form used in the DTU "top99neo" code (Eq. 6). [PR] [Ferrari & Sigmund 2020, SMO 62:2211–2228](https://doi.org/10.1007/s00158-020-02629-w), full text via [arXiv:2005.05436](https://arxiv.org/abs/2005.05436)
- In finite-strain compliant-mechanism TO, Kumar et al. use modified SIMP E_e(ρ̄_e) = E_v + (E_s − E_v)ρ̄_e^p with E_v = 10⁻⁶ E_s and p = 3, on the projected density. [PR] [Kumar, Schmidleithner, Larsen & Sigmund 2021, SMO 63:1351–1366](https://doi.org/10.1007/s00158-020-02764-4) ([arXiv:2008.08902](https://arxiv.org/abs/2008.08902))
- In finite-strain periodic microstructure TO, Behrou et al. use E_e = [ρ_min + ρ_e^p (1 − ρ_min)]E0 with ρ_min = 10⁻⁴, "to maintain positive definiteness". [Non-PR preprint] [Behrou et al., arXiv:2005.09111](https://arxiv.org/abs/2005.09111v1)
- Zhang, Khandelwal & Guo interpolate E(ρ) between E_min = εE and E (their Eq. 3), and use ε = 10⁻⁸ and p = 3 in the stability examples. [Non-PR preprint] [arXiv:2303.14911](https://arxiv.org/abs/2303.14911)
- Primary SIMP and interpolation references (bibliographic data verified):
  - Bendsøe 1989, "Optimal shape design as a material distribution problem", *Struct. Optim.* 1(4):193–202 [doi:10.1007/BF01650949](https://doi.org/10.1007/BF01650949).
  - Bendsøe & Sigmund 1999, "Material interpolation schemes in topology optimization", *Arch. Appl. Mech.* 69(9–10):635–654 [doi:10.1007/s004190050248](https://doi.org/10.1007/s004190050248).
  - Sigmund 2007, "Morphology-based black and white filters for topology optimization", *SMO* 33(4–5):401–424 [doi:10.1007/s00158-006-0087-x](https://doi.org/10.1007/s00158-006-0087-x).
  - Andreassen et al. 2011, "Efficient topology optimization in MATLAB using 88 lines of code", *SMO* 43(1):1–16 [doi:10.1007/s00158-010-0594-7](https://doi.org/10.1007/s00158-010-0594-7).
  - Stolpe & Svanberg 2001, "An alternative interpolation scheme for minimum compliance topology optimization" (RAMP), *SMO* 22(2):116–124 [doi:10.1007/s001580100129](https://doi.org/10.1007/s001580100129).
- [PR-bib] The following summaries are prior knowledge and were not re-verified:
  - Sigmund 2007 and Andreassen et al. 2011 introduced and popularized "modified SIMP" (E_min + x^p ΔE, E_min ≈ 10⁻⁹E0). It replaced the older approach of a lower bound x_min on the design variable with E = x^p E0, and it decouples void stiffness from the penalization.
  - Bendsøe & Sigmund 1999 discuss physical realizability of SIMP: p must exceed a Poisson-ratio-dependent bound (p ≥ 3 for ν0 = 1/3 in 2D).
  - RAMP, E = E0 ρ/(1+q(1−ρ)), has non-zero slope at ρ = 0.

### Inferences
- **Derivative jump and non-smoothness.** At H^p = ε the derivative of f drops from p·ε^{(p−1)/p}·H′ to 0.
  - For p = 3 this is ≈0.019·H′ → 0.
  - For p = 1, the first continuation stage, it is H′ → 0.
  - The objective is therefore only piecewise C¹. MMA/GCMMA theory assumes continuously differentiable functions (see Section 4), and the kink can cause oscillation of elements around the clip.
- **Dead zone and recovery.**
  - Below the clip, the stiffness-driven sensitivity is exactly zero. The only density dependence left is the blend term p ρ̃^{p−1} g H_Steep′(ρ̃^p) (see [NEMOpt.jl L618](/home/user/NEMOpt/Full_Domain/NEMOpt.jl)), which is scaled by g = ε.
  - The volume-constraint gradient (H′ > 0) still pushes these elements downward. Elements that fall below the clip are therefore hard to recover from objective information alone.
  - With modified SIMP the slope p H^{p−1}(1−ε)H′ is continuous. For p = 1 it is non-zero all the way to 0, so recovery is possible early in the continuation.
- **Magnitude.** The two laws differ by at most ε in stiffness, because ε+(1−ε)H^p − max(H^p, ε) ∈ [0, ε]. Final homogenized moduli therefore differ only slightly for the same design. The design trajectories can differ substantially, though, because of the dead-zone and kink effects.
  - ε = 5e-4 is large compared with 1e-6 to 1e-9 in linear and compliant-mechanism TO, and comparable to the 1e-4 used by Behrou et al.
  - The void stiffness contributes to 𝔸 and therefore to ν. This matters when the target is near the bound −1 (Section 5).
- **Box constraint.** Paper (ρ_min ≤ ρ) and code (0 ≤ ρ) differ. For reproduction, keep 0 (code). DOC the discrepancy.
- **Recommendation: SWITCH.**
  - `stiffness_interp = :clip_legacy` (Eq. 29, the default for paper-reproduction runs) vs `:modified_simp` (ε + (1−ε)H^p, the new default). RAMP is optional.
  - Keep ε as a separate parameter from the design lower bound.

### Gaps
- I could not confirm whether the authors' published results (Figs. 7–9, Tables 5–7) were produced with exactly ε = 5e-4. The article tables do not list ε. The code value is the only evidence.
- No paper found that studies `max(·, ε)` clipping specifically. The effects described above are derived, not cited.

---

## 2. Void-element distortion remedies and the definition of the switching variable (Eq. 30, Ĥ on ρ̃^p)

### Takeaway
Wang et al. (2014) switch low-density elements to small-strain energy using a smooth Heaviside γ. Per Kumar et al.'s restatement, γ is evaluated on the **physical (projected) density raised to p**, with β1 = 500 and η0 = 0.01. NEMOpt evaluates its switch Ĥ on the **filtered, unprojected** density raised to p, with β̂ = 200 and ω̂ = 0.05. As a result:
- the NEMOpt switch sits at ρ̃ = 0.05^{1/p} (0.05 → 0.37 during p-continuation), independent of β;
- at the end of continuation (p = 3, β = 32), a band ρ̃ ∈ (0.368, 0.462) remains that is still **hyperelastic but has only ε stiffness**, which is exactly the distortion-prone regime.

Keep the current switch for reproduction. Offer a Wang-style switch on the physical density as a new option.

### Cited Findings
- **NEMOpt implementation.**
  - `P_rho = g(ρ)*(H_Steep(ρ^p)*P + (1-H_Steep(ρ^p))*S_Linear)`
  - `A_rho = g(ρ)*(H_Steep(ρ^p)*A + (1-H_Steep(ρ^p))*C_Linear)`
  - `Z_rho = g(ρ)*H_Steep(ρ^p)*Z`
  - Here ρ is the filtered nodal field interpolated to Gauss points, with β₂ = 200 and ω₂ = 0.05.
  - [Repo] [Full_Domain/NEMOpt.jl L56–57, L297–299](/home/user/NEMOpt/Full_Domain/NEMOpt.jl), [NEMOpt_main.jl L21](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR.jl L507–509](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR.jl)
- **Article.** Eq. 30 is P_m = f(ρ̃)[Ĥ ∇ψ(F_m) + (1−Ĥ) C:∇u_m]. The article attributes the "modified stress expression" to its ref. [16] (Wallin & Tortorelli 2020) and states that the second term "prevents excessive deformation of elements with very low densities". The tables list (ω, β₂, ω₂) = (0.5, 200, 0.05). The article does not state that Ĥ's argument is ρ̃^p. [Article] [Sec. 2.5, Tables 3–4](https://doi.org/10.1016/j.finel.2026.104599)
- **Article conclusions.** "Excessive element distortion problem which typically emerges in very low density elements at large strains is a limiting barrier on the magnitude of macroscopic loading". [Article] [Sec. 6](https://doi.org/10.1016/j.finel.2026.104599)
- **Wang et al. 2014 energy interpolation**, as restated and used by Kumar et al.:
  - W̄_e(u_e) = [W_e(γ_e u_e) − W_e^L(γ_e u_e) + W_e^L(u_e)] E_e(ρ̄_e).
  - γ_e = [tanh(β1η0) + tanh(β1(ρ̄_e^p − η0))]/[tanh(β1η0) + tanh(β1(1−η0))].
  - ρ̄_e is the projected (physical) density. "We use β1 = 500 and η0 = 0.01 for all examples as suggested by (Wang et al., 2014)".
  - [PR] [Kumar et al. 2021, Eqs. 1, 4, 11–12](https://doi.org/10.1007/s00158-020-02764-4) ([arXiv full text](https://arxiv.org/abs/2008.08902)). Original: [PR-bib] [Wang, Lazarov, Sigmund & Jensen 2014, CMAME 276:453–472](https://doi.org/10.1016/j.cma.2014.03.021). The original full text was not accessed in this session.
- **Adaptive linear energy interpolation (Zhang, Alberdi & Khandelwal 2018).**
  - It interpolates the deformation gradient, F = I + η(ρ)∇u with η = exp(βρ)/(exp(cβ)+exp(βρ)), β = 120, c = c₀ = 0.08.
  - The internal force mixes a finite-strain term with a small-strain C:ε term. The linear part uses a **higher penalization p_L > p** "to discourage using low-density values to exploit small deformation kinematics".
  - c is increased by 0.05 whenever the FE analysis fails to converge, and reset afterwards. This adaptivity is described as "crucial for improving the robustness".
  - [Non-PR preprint describing PR method] [Zhang, Khandelwal & Guo, arXiv:2303.14911, Sec. 3.2](https://arxiv.org/abs/2303.14911). Original: [Zhang, Alberdi & Khandelwal 2018, IJNME 115(8):1015–1052, doi:10.1002/nme.5834](https://doi.org/10.1002/nme.5834). This is also the article's ref. [40].
- **Catalog of other remedies** (from the same source):
  - excluding void DOFs from the Newton convergence check (Buhl et al. 2000);
  - element removal and reintroduction (Bruns & Tortorelli 2003);
  - element connectivity parameterization with zero-length links (Yoon & Kim 2005);
  - a distortion-indicator-informed design update (Lahuerta et al. 2013);
  - linear energy interpolation (Wang et al. 2014; Zhang et al. 2018; Pajot 2006).
  - [Non-PR preprint] [arXiv:2303.14911](https://arxiv.org/abs/2303.14911)
  - Bibliographic data verified via Crossref: [Bruns & Tortorelli 2003, IJNME 57(10):1413–1430](https://doi.org/10.1002/nme.783); [Yoon & Kim 2005, IJSS 42(7):1983–2009](https://doi.org/10.1016/j.ijsolstr.2004.09.005).
  - Bibliographic data taken from the reference list of arXiv:2303.14911, not checked against Crossref: [Buhl, Pedersen & Sigmund 2000, SMO 19(2):93–104](https://doi.org/10.1007/s001580050089); [Lahuerta et al. 2013, Comput. Mech. 52(4):779–797](https://doi.org/10.1007/s00466-013-0843-x); Zhang, Alberdi & Khandelwal 2018 (doi:10.1002/nme.5834, which also matches the article's ref. [40]).
- **Other remedies named in the brief** (bibliographic data verified):
  - [van Dijk, Langelaar & van Keulen 2014, "Element deformation scaling for robust geometrically nonlinear analyses in topology optimization", SMO 50(4):537–560](https://doi.org/10.1007/s00158-014-1145-4).
  - [Luo, Wang & Kang 2015, "Topology optimization of geometrically nonlinear structures based on an additive hyperelasticity technique", **CMAME** 286:422–441](https://doi.org/10.1016/j.cma.2014.12.023). The brief lists this as IJNME, but Crossref gives CMAME.
  - [Klarbring & Strömberg 2013 (online 2012), "Topology optimization of hyperelastic bodies including non-zero prescribed displacements", SMO 47(1):37–48](https://doi.org/10.1007/s00158-012-0819-z).
  - [PR-bib] van Dijk et al. scale the deformation used in the nonlinear analysis of low-density elements. Luo et al. add a weak hyperelastic "additive" material everywhere to stabilize void regions. Neither summary was re-verified.
- **Wallin & Tortorelli 2020 abstract.**
  - Neo-Hookean constituents, Bloch-type BCs, MMA, "a periodic filter and thresholding".
  - Designs "for maximum load carrying capacity subject to negative Poisson's ratio constraints".
  - The adjoint analysis treats both displacement and incremental-displacement responses.
  - [PR] [Wallin & Tortorelli 2020, Mech. Mater. 145:103324](https://doi.org/10.1016/j.mechmat.2020.103324). Abstract via [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.mechmat.2020.103324?fields=title,abstract).
- A later work reports that linear void energy becomes problematic when voids are fully enclosed by solid, because the infinitesimal strain is not rotation-invariant. **Unverified**: this claim surfaced only in a search-engine summary of arXiv results (e.g., [arXiv:2607.13174](https://arxiv.org/abs/2607.13174)), and I could not confirm which paper makes it.

### Inferences
- **Computed switch locations.** These come from evaluating NEMOpt's H and Ĥ formulas (bisection, scratch computation). Here "Wang-type" means γ = H(H(ρ̃)^p; 500, 0.01).

  | p, β | NEMOpt switch (Ĥ=0.5) at ρ̃ | clip boundary H^p=ε at ρ̃ | Wang-type switch at ρ̃ | f at NEMOpt switch |
  |---|---|---|---|---|
  | 1, 1 | 0.050 | 0.001 | 0.012 | 4.4e-2 |
  | 2, 1 | 0.224 | 0.026 | 0.112 | 4.3e-2 |
  | 3, 1 | 0.368 | 0.090 | 0.231 | 4.6e-2 |
  | 3, 8 | 0.368 | 0.347 | 0.419 | 1.3e-3 |
  | 3, 16 | 0.368 | 0.423 | 0.460 | 5.0e-4 (clipped) |
  | 3, 32 | 0.368 | 0.462 | 0.480 | 5.0e-4 (clipped) |

- **What the table shows.**
  - Early in the run, NEMOpt linearizes elements that still carry about 4% stiffness. That is more aggressive than Wang et al.
  - From β ≈ 16 onward, elements with ρ̃ ∈ (0.37, ≈0.46) are fully hyperelastic (Ĥ ≈ 1) at the minimum stiffness ε. These are soft hyperelastic elements, which are prone to inversion and the negative detJ the code checks for at [NEMOpt.jl L265](/home/user/NEMOpt/Full_Domain/NEMOpt.jl).
  - A switch driven by the physical density moves toward 0.5 as β grows and never leaves such a band.
  - This is a plausible mechanism for the robustness limits the article reports. It is an inference and has not been tested.
- **Is ρ̃^p consistent with the literature?**
  - The *exponent p* matches Wang et al. (γ is driven by ρ̄^p).
  - Using the *filtered* rather than the *projected* density does not match.
  - Using β̂ = 200 and ω̂ = 0.05 in place of 500 and 0.01 is a further departure.
  - I could not verify which density Wallin & Tortorelli 2020 used (Gaps).
- **Blend form.** NEMOpt's blend is a convex combination of *stresses at the same F*, weighted by a density-only function. It derives from the potential f[Ĥψ(F) + (1−Ĥ)½ε:C:ε], so the tangent is symmetric and the code's `A_rho` and `Z_rho` are consistent. It is still not the Wang form, which evaluates W at the *scaled* displacement γu. Because Ĥ does not depend on u, the third-order term `Z_rho = gĤZ` is correct.
- **Expected effect of changing the switch.** It would alter which elements behave nonlinearly at intermediate β, and therefore both the optimization path and the final topology. It must therefore not be changed silently.
- **Rotating-hinge designs.** Rotating-hinge auxetic designs (the article's Table 5) rotate void regions. If the unverified enclosed-void caveat holds, the linear-void regions can produce spurious stresses under large rigid rotation. DOC this as a risk for large F_M.
- **Recommendation: SWITCH.**
  - `void_switch = :legacy_filtered_pow_p` (Ĥ(ρ̃^p), β̂ = 200, ω̂ = 0.05; the reproduction default).
  - `void_switch = :physical_pow_p` (Ĥ(H(ρ̃)^p) with β1 = 500 and ρ0 = 0.01 per Wang et al. 2014/Kumar et al. 2021; the recommended new default).
  - Optionally an adaptive scheme (Zhang et al. 2018), with a separate penalization for the linear term.
  - Add a per-solve diagnostic: minimum detJ and minimum principal stretch in elements with Ĥ > 0.5 and f ≤ 10ε.

### Gaps
- **Wallin & Tortorelli 2020 full text** (the source of Eq. 30) was not accessible: ScienceDirect returned 403 and the Elsevier API was rate-limited or invalid. Still unknown: their switch variable (filtered vs. thresholded density), parameter values, and whether the article's (200, 0.05) came from them.
- The original Wang et al. 2014 text was not read. The γ definition is taken from Kumar et al.'s restatement, which says the values are "as suggested by" Wang et al.
- I did not confirm the content of van Dijk et al. 2014, Luo et al. 2015 or Klarbring & Strömberg 2013 beyond titles.

---

## 3. Heaviside projection and p/β continuation (schedules, counting, MMA resets)

### Takeaway
NEMOpt's tanh projection matches the standard Wang–Lazarov–Sigmund (2011) form, and applying the volume constraint on the projected density is consistent with Ferrari & Sigmund (2020). The continuation implementation, however, has four reproduction-critical quirks:
- it counts **function evaluations**, so GCMMA inner iterations advance it;
- it updates p and β **inside the objective, before the gradient is computed in the same rrule**, so the gradient is inconsistent at those evaluations;
- the "reset" of MMA's s_init is dead code;
- the full-domain script's schedule differs from the article and does not reach β = 32 within N_ITER = 120.

### Cited Findings
- **Standard tanh projection.**
  - The form is H(x̃; β, η) = [tanh(βη) + tanh(β(x̃−η))]/[tanh(βη) + tanh(β(1−η))]. Kumar et al. attribute it to Wang, Lazarov & Sigmund 2011 (their Eq. 1) and use β from 1 to 32, "doubled each 60th MMA iteration".
  - They note that "as per β continuation, jumps in the convergence curves of the objectives and volume constraints can be noticed" and they use an external move limit of 0.1.
  - [PR] [Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4) ([arXiv](https://arxiv.org/abs/2008.08902)); [PR-bib] [Wang, Lazarov & Sigmund 2011, SMO 43(6):767–784](https://doi.org/10.1007/s00158-010-0602-y)
- **Ferrari & Sigmund 2020 (top99neo).**
  - Continuation is the structure `parCont = {istart, maxPar, isteps, deltaPar}`: it "starts when loop=istart and the parameter is increased by deltaPar each isteps, up to the value maxPar", counted in **design-loop iterations**.
  - Example: `penalCnt={50,3,25,0.25}`, `betaCnt={250,16,25,2}`. Acceleration is started only after "both continuations have finished".
  - With fixed η, "the volume constraint should be consistently applied on V(x̂), otherwise some violation or over-shooting of the constraint will happen". The paper also gives a Newton procedure for a volume-preserving η*.
  - [PR] [Ferrari & Sigmund 2020, SMO 62:2211–2228](https://doi.org/10.1007/s00158-020-02629-w) ([arXiv:2005.05436](https://arxiv.org/abs/2005.05436))
- **NEMOpt projection and volume constraint.**
  - Projection is applied at Gauss points to the bilinearly interpolated nodal filtered field.
  - The volume constraint integrates H(ρ̃) (Eq. 33), i.e., it is applied on the physical density.
  - [Repo] [NEMOpt.jl L53, L75–85](/home/user/NEMOpt/Full_Domain/NEMOpt.jl); [Article Eq. 33](https://doi.org/10.1016/j.finel.2026.104599)
- **Article schedules (stated per "optimization iteration").**
  - Sec. 5.1 (shear): p += 0.2 every 5 iterations to 3, then β += 1 every 5 iterations to 32.
  - Sec. 5.2 (Poisson): p += 0.1 every 5 iterations, then β += 1 every 5 iterations to 32. The article adds: "under relatively large strains, the optimization algorithm becomes sensitive to the parameters of the continuation procedure".
  - Sec. 5.3 (bulk): p += 0.5 every 5 iterations, then β += 4 every 5 iterations.
  - [Article] [Sec. 5](https://doi.org/10.1016/j.finel.2026.104599)
- **Code schedules.**
  - `Full_Domain/NEMOpt_main.jl`: p += 0.5 every 5 evaluations; then, once p == 3, β += 2 every 10 evaluations while β < 32. The comments say "0.2" and "by 1 every 10 until 40", which do not match the code.
  - `NEMOpt-DR_main.jl`: p += 0.1 every 5 evaluations; then β += 1 every 5 evaluations. This matches Sec. 5.2.
  - The counter `iter += 1` is incremented at the top of `Objective`.
  - [Repo] [NEMOpt_main.jl L45, L90–101](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR_main.jl L45, L91–101](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl)
- **What "iterations" mean in Nonconvex.jl.**
  - NonconvexMMA's GCMMA inner loop calls `value_jacobian(...)` on the objective and constraints at **every inner iteration**. The loop runs until the approximation is an upper bound or `maxinner` is reached. [Non-PR code] [NonconvexMMA.jl `mma_algorithm.jl`](https://raw.githubusercontent.com/JuliaNonconvex/NonconvexMMA.jl/master/src/mma_algorithm.jl)
  - The docs define `maxiter` as the "maximum number of inner iterations", with one inner iteration per outer iteration for MMA87. [Non-PR doc] [Nonconvex.jl MMA docs](https://julianonconvex.github.io/Nonconvex.jl/stable/algorithms/mma/)
- **Code path for the gradient.**
  - The objective's `rrule` calls `Objective(ρ)` and then `Objective_Grad(ρ)`. [Repo] [NEMOpt_main.jl L213–216](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl)
  - `Objective` updates the globals p and β at its end ([L90–101](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl)).
  - `Objective_Grad` then re-enters `constitutiveEqs`, `g`, `g_Deriv`, `H_Steep(ρ^p)` and `H_Deriv`, which read the global p and β ([NEMOpt.jl L612–625](/home/user/NEMOpt/Full_Domain/NEMOpt.jl)). It reuses `U_Nonlin`, `χ` and `A`, which were computed with the old values.
- The line `options = MMAOptions(s_init=0.5/(β+1), ...)` inside `Objective` assigns a *local* variable after `optimize(...)` has already received its options. [Repo] [NEMOpt_main.jl L100, L236–240](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl)
- **Further projection and continuation literature** (bibliographic data verified):
  - [Guest, Prévost & Belytschko 2004, IJNME 61(2):238–254](https://doi.org/10.1002/nme.1064)
  - [Xu, Cai & Cheng 2010, SMO 41(4):495–505](https://doi.org/10.1007/s00158-009-0452-7)
  - [Li & Khandelwal 2015, "Volume preserving projection filters and continuation methods in topology optimization", Eng. Struct. 85:144–161](https://doi.org/10.1016/j.engstruct.2014.10.052)
  - [Guest, Asadpoure & Ha 2011, "Eliminating beta-continuation from Heaviside projection and density filter algorithms", SMO 44(4):443–453](https://doi.org/10.1007/s00158-011-0676-1)
  - [Sigmund 2022, "On benchmarking and good scientific practise in topology optimization", SMO 65(11)](https://doi.org/10.1007/s00158-022-03427-2)
  - [PR-bib] Content summaries (not re-verified):
    - Guest et al. 2004 introduced regularized-Heaviside projection of nodal variables for minimum length scale, with β continuation.
    - Xu et al. 2010 shift the threshold to preserve volume.
    - Guest et al. 2011 propose replacing β-continuation by controlling MMA move limits and asymptotes as a function of β. NEMOpt's dead `s_init = 0.5/(β+1)` line resembles this idea.
    - Sigmund 2022 lists recommendations for reproducible TO. Per a search summary of its abstract, it puts "the responsibility of fair evaluation of contributions" on the authors.

### Inferences
- **Gradient inconsistency.** Every evaluation where `iter % 5 == 0`, or `iter % 10 == 0` for β in the full script, returns J(p_old, β_old) together with an approximate gradient computed with (p_new, β_new).
  - In the DR script that happens on about 50 of 300 evaluations.
  - For MMA87 this corrupts the linearization at that step.
  - For GCMMA, the next conservativeness checks compare against a different function.
  - This is a defect, not a modeling choice. It still shapes the published trajectories, so reproducing them bit for bit requires keeping it in legacy mode.
- **Schedule depends on the optimizer.** Continuation advances per evaluation, and MMA02 performs up to `maxinner = 40` evaluations per outer iteration. The p/β schedule *per outer iteration* therefore depends on how often the conservative approximation fails, which in turn is affected by rounding (Section 4).
  - The paper's "every 5 optimization iterations" holds only for MMA87.
  - In DR, with the code's schedule, p reaches 3 at evaluation 100 and β reaches 32 at evaluation 250 (computed), within `maxiter = 300`.
- **Full-domain script does not reproduce Sec. 5.1.**
  - p reaches 3 at evaluation 20, which also triggers β += 2 at that evaluation.
  - β then increases at 30, 40, …, 120, so β ≈ 23 after 120 evaluations (computed). If run long enough, β overshoots to 33 rather than stopping at 32.
  - Article Sec. 5.1 instead implies p = 3 at iteration 50 and β = 32 at about iteration 205.
  - For reproduction, the schedule must be data, not code.
- **No asymptote or GCMMA-ρ reset.** After a β jump, the MMA asymptotes and the oscillation-history factors carry over unchanged (the update rule is quoted in Section 4). The problem curvature changes abruptly with β, which tends to produce the "jumps" Kumar et al. report.
- **Recommendation: SWITCH.**
  - `continuation_counter = :function_evals_legacy` vs `:outer_iterations` (new default).
  - `continuation_timing = :inside_objective_legacy` vs `:between_iterations` (new default). Update p and β only after the gradient of the current iterate has been returned.
  - Express schedules as `{istart, max, every, delta}` tuples as in Ferrari & Sigmund 2020.
  - Provide per-case presets for paper Secs. 5.1, 5.2 and 5.3, plus a "code-as-shipped" preset.
  - Optional new behaviors, not in legacy: reset the MMA asymptotes after each β/p jump, and scale `s_init` with β as in Guest et al. 2011. These require Nonconvex.jl API support or an in-house MMA.

### Gaps
- **Sigmund 2022 full text** (Springer redirected to login) was not accessed. Its specific advice on continuation and reporting is unconfirmed.
- I found **no peer-reviewed source that verifies** the recommendation to reset MMA asymptotes or GCMMA ρ after continuation steps. Treat it as a heuristic.
- Li & Khandelwal 2015 and Guest et al. 2011 contents were not re-verified.

---

## 4. MMA87 vs MMA02 (GCMMA), inner iterations, and returning a rounded objective

### Takeaway
GCMMA (Svanberg 2002) accepts an iterate only if each conservative approximation is ≥ the true function value, within ~1e-7 in Svanberg's code. Each rejected inner iteration costs a full nonlinear solve plus gradients in Nonconvex.jl. NEMOpt returns J rounded to 1e-3 (error up to 5e-4), and its gradient uses ν rounded to 1e-3. In late iterations the rounding noise therefore exceeds both the predicted decrease and the acceptance tolerance. The full-domain and DR scripts also use *different* algorithms (MMA87 vs MMA02). Keep both per script for reproduction, and remove the rounding behind a legacy flag.

### Cited Findings
- **Algorithm references.**
  - MMA: [Svanberg 1987, IJNME 24(2):359–373](https://doi.org/10.1002/nme.1620240207).
  - GCMMA: [Svanberg 2002, "A class of globally convergent optimization methods based on conservative convex separable approximations", SIAM J. Optim. 12(2):555–573](https://doi.org/10.1137/S1052623499362822). These are the article's refs. [32]–[33].
  - [PR-bib] The global convergence results assume continuously differentiable problem functions, and convergence is to KKT points (not re-verified).
- **Svanberg's code conventions** (GPL Python port of his MATLAB code):
  - Asymptotes start at x ± asyinit·(xmax−xmin) for iterations ≤ 2, with asyinit = 0.5.
  - Afterwards the factor is asyincr = 1.2 if (x−x_old1)(x_old1−x_old2) > 0 and asydecr = 0.7 if < 0. Asymptote distances are bounded to [0.01, 10]·(xmax−xmin). Other defaults: move = 0.5, albefa = 0.1.
  - GCMMA `concheck` declares the approximations conservative only if f_app + epsimin ≥ f_new for the objective **and all constraints**, with epsimin = 1e-7.
  - Otherwise `raaupdate` raises ρ (raa) to 1.1·(raa + (f_new − f_app)/raacof), capped at 10× its previous value.
  - [Non-PR code] [GCMMA-MMA-Python `mma.py`](https://raw.githubusercontent.com/arjendeetman/GCMMA-MMA-Python/master/src/mmapy/mma.py)
  - The port's README states it derives from Svanberg's MATLAB code obtained from smoptit.se. [repo](https://github.com/arjendeetman/GCMMA-MMA-Python)
- **Nonconvex.jl implementation.**
  - `MMAOptions` defaults: `s_init = 0.5`, `s_incr = 1.2`, `s_decr = 0.7`, `maxinner = 10`, `maxiter = 1000`, `convcriteria = KKTCriteria()`.
  - Inner iterations re-evaluate **values and Jacobians** and accept when no exact value exceeds its approximation. Otherwise ρ is increased.
  - [Non-PR code] [NonconvexMMA.jl source](https://raw.githubusercontent.com/JuliaNonconvex/NonconvexMMA.jl/master/src/mma_algorithm.jl); [Nonconvex.jl docs](https://julianonconvex.github.io/Nonconvex.jl/stable/algorithms/mma/)
- **NEMOpt settings.**
  - The full-domain script uses `alg = MMA87()` and DR uses `alg = MMA02()`.
  - Both use `s_init = 0.25`, `s_incr = 1.2`, `s_decr = 0.7`, `maxinner = 40`, `tol = Tolerance(kkt = 1e-8, f = 0.0)`, `KKTCriteria()`.
  - [Repo] [NEMOpt_main.jl L235–240](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR_main.jl L234–239](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl)
  - The article says "The default parameters given by [32] are used for MMA". [Article] [Sec. 5](https://doi.org/10.1016/j.finel.2026.104599)
- **Rounding in the code.**
  - `J, Vf = round(J, digits=3), ...` and `v_12, K, G = round(v_12, digits=3), ...`, followed by `return J`.
  - The gradient is `∇J = 200*(v_12 − v_0)*∇v_12`, using the rounded global `v_12` and an unrounded ∇v_12.
  - The volume constraint is returned unrounded.
  - [Repo] [NEMOpt_main.jl L79–80, L104, L117](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR_main.jl L78–79, L104, L117](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl)
- **Benchmarking.** Rojas-Labanda & Stolpe benchmarked OC, MMA, GCMMA, IPOPT, FMINCON and SNOPT on compliance, volume and mechanism problems. They conclude that "general solvers are as efficient and reliable as classical structural topology optimization solvers". [PR] [Rojas-Labanda & Stolpe 2015, SMO 52(3):527–547](https://doi.org/10.1007/s00158-015-1250-z) (abstract via [DTU Orbit](http://orbit.dtu.dk/en/publications/benchmarking-optimization-solvers-for-structural-topology-optimization(472ad756-cb7d-49b4-bd24-80da05993ffc).html))
- **Link rot.** Svanberg's former KTH homepage (people.kth.se/~krille), which hosted the "MMA and GCMMA – two methods for nonlinear optimization" notes, now shows a notice that the page "is no longer in use" (fetched 2026-10-12). The notes could not be retrieved. [Non-PR] [people.kth.se/~krille](https://people.kth.se/~krille/mmagcmma.pdf)

### Inferences
- **Size of the quantization.** J = 100(ν+1)², so the rounding error |ΔJ| ≤ 5e-4. The gradient factor 200(ν_r+1) carries a relative error up to 5e-4/|ν+1|: about 0.5% at ν = −0.90 and about 5% at ν = −0.99.
- **Effect on GCMMA.** Late in the run, predicted decreases fall below ~1e-3. The ±5e-4 noise in the "true" objective is then ~5000× the 1e-7 tolerance in Svanberg's check, and Nonconvex.jl's check has no tolerance at all. Consequences:
  - spurious rejections, which grow ρ, shrink steps and add inner iterations, each a full nonlinear solve, four χ solves and the adjoints;
  - spurious acceptances, which break the conservativeness guarantee behind GCMMA's convergence;
  - because continuation counts evaluations, extra inner iterations also accelerate the p/β schedule (Section 3).
- **Effect on MMA87.** MMA87 has no acceptance test. Rounding mainly affects the KKT and convergence tests and the ν used in the gradient.
- **Termination.** With `kkt = 1e-8` and a quantized objective, KKT convergence is practically unattainable, so runs end at `maxiter`.
- **Confounded comparison.** The article reports that full and reduced domains give the same topology (Fig. 8). The shipped scripts, however, differ in optimizer, schedule and N_ITER, so a refactor's regression tests must pin each script's configuration rather than assume equivalence.
- **Recommendation.**
  - **SWITCH** for rounding: `round_objective_digits = 3` in legacy, `nothing` in the new default. Keep rounding for *logging* only.
  - **KEEP** the per-script algorithm (MMA87 full, MMA02 DR), `s_init = 0.25` and `maxinner = 40` as explicit, recorded config. DOC that these are not the library defaults (0.5 and 10).
  - Log the number of inner iterations per outer iteration.

### Gaps
- The Svanberg 2002 and 2007-notes statements on differentiability assumptions and the inner-loop termination tolerance could not be re-read, because the KTH host is gone and SIAM was not accessed. They are flagged [PR-bib].
- No literature found on returning quantized objectives to MMA. The effects above are derived from the acceptance logic in the code.

---

## 5. Poisson's ratio and effective properties under finite strain (Eq. 54, bulk and shear measures)

### Takeaway
ν = 2𝔸1122/(𝔸1111+𝔸2222) is an **incremental (tangent), in-plane** measure evaluated at the pre-stretched state F_M. For a plane-strain isotropic solid it equals λ/(λ+2μ) = ν3D/(1−ν3D), which is 0.4286 for the base material. For orthotropic cells it is the harmonic mean of the two compliance-based in-plane ratios, and it lies strictly inside (−1, 1) whenever the normal-normal 2×2 block of the tangent is positive definite. ν0 = −1 is therefore an unattainable bound. This differs from the **secant, uniaxial, engineering-strain** Poisson's ratios used by Clausen et al. 2015 and Wang 2018. KEEP the definition for reproduction, and DOC its meaning.

### Cited Findings
- **Article definitions.**
  - Eq. 54, ν = 2(A_M)1122/((A_M)1111+(A_M)2222), is attributed to Kazemi et al. [39]. "As opposed to linear case, Poisson's ratio definition is not unique for non-linear elasticity, [16]."
  - The linear comparison uses 2(C_M)1122/((C_M)1111+(C_M)2222).
  - Reported values: −0.9008 (load case 1, nonlinear), −0.9132 (case 2, with shear), −0.9218 (linear). The article says case 2 is "closer to the theoretical lower limit, i.e., −1".
  - Bulk moduli are computed from the tangent, and "care has to be taken since these values are based on macroscopic tangent stiffness tensor which deviates from small strain elasticity tensor as deformation gets larger".
  - [Article] [Sec. 5.2–5.3](https://doi.org/10.1016/j.finel.2026.104599); [PR-bib] [Kazemi, Vaziri & Norato 2020, CMAME 363:112895](https://doi.org/10.1016/j.cma.2020.112895)
- **Code definitions.**
  - `v_12 = 2A[1,2]/(A[1,1]+A[2,2])`; `K = ¼(A11+A22+A12+A21)`.
  - `G = A[3,3]+A[4,4]` in the full-domain script but `G = A[3,3]` in DR. This affects logging only; the shear objective `1/(A[3,3]+A[4,4])` is identical in both mains.
  - [Repo] [NEMOpt.jl L65–71](/home/user/NEMOpt/Full_Domain/NEMOpt.jl), [NEMOpt-DR.jl L61–67](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR.jl)
- **Clausen et al. 2015.** Designs with Poisson's ratios from −0.8 to 0.8 maintained over strains ≥ 20%, using a geometrically nonlinear model and "minimizing the error between the actual and the pre-defined value of Poisson's ratio over a range of discrete, nominal strain values up to 20%". [PR] [Clausen, Wang, Jensen, Sigmund & Lewis 2015, Adv. Mater. 27(37):5523–5527](https://doi.org/10.1002/adma.201502485) (abstract via [DTU Orbit](https://orbit.dtu.dk/en/publications/topology-optimized-architectures-with-programmable-poissons-ratio/))
- **Wang 2018.** "The secant Poisson's ratio is defined as the negative ratio between the lateral and the longitudinal engineering strains". [PR] [Wang 2018, JMPS 114:303–318](https://doi.org/10.1016/j.jmps.2018.01.013). Quoted from an abstract snippet in search results; full text not read.
- **Wang, Sigmund & Jensen 2014.** [PR-bib] [JMPS 69:156–174](https://doi.org/10.1016/j.jmps.2014.05.003). They design materials to follow prescribed nonlinear stress–strain (load-path) responses via finite-strain homogenization (not re-verified).
- **Wallin & Tortorelli 2020.** They discuss "domain symmetry of the unit cell and its effect on uniaxial loading conditions" and impose negative Poisson's ratio *constraints*. Their response functions include both the displacement and the incremental displacement. [PR] [Mech. Mater. 145:103324](https://doi.org/10.1016/j.mechmat.2020.103324) (abstract)
- **Finite-strain Poisson function and 2D bounds** (bibliographic data verified):
  - [Beatty & Stalnaker 1986, "The Poisson function of finite elasticity", J. Appl. Mech. 53(4):807–813](https://doi.org/10.1115/1.3171862)
  - [Thorpe & Jasiuk 1992, "New results in the theory of elasticity for two-dimensional composites", Proc. R. Soc. A 438:531–544](https://doi.org/10.1098/rspa.1992.0124)
  - [PR-bib] Beatty & Stalnaker define a stretch-dependent Poisson function for finite elasticity. Thorpe & Jasiuk treat 2D elasticity, in which isotropic Poisson's ratio is bounded by −1 ≤ ν ≤ 1 (not re-verified).
  - Ting & Chen 2005 (QJMAM 58:73, "Poisson's ratio for anisotropic elastic materials can have no bounds"): **DOI 10.1093/qjmam/hbh021 could not be verified** (Crossref failed).
- **Validity limit of homogenized tangents.** Homogenized results lose validity after microscopic bifurcation, and Bloch-wave analysis is the standard way to detect such instabilities. [Non-PR preprint] [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371)

### Inferences
These are derivations. I did not find a source that states them in this form.
- **Plane-strain isotropic solid.** C1111 = λ+2μ and C1122 = λ, so Eq. 54 gives λ/(λ+2μ) = ν/(1−ν). For E = 100 and ν = 0.3 this is 0.4286 (computed). The measure is the in-plane (2D) Poisson's ratio, not the 3D ν.
- **Orthotropic 2D, linear.** Under uniaxial in-plane stress, ν_xy = C12/C22 and ν_yx = C12/C11. Then 2C12/(C11+C22) = 2/(1/ν_yx + 1/ν_xy), the harmonic mean, which is exact for square symmetry.
- **Bound.**
  - Positive definiteness of the normal block gives C11C22 > C12², and C11+C22 ≥ 2√(C11C22). Hence |2C12/(C11+C22)| < 1.
  - J = 100(ν+1)² therefore never reaches 0, and its gradient never vanishes from the objective side. The optimum is set by the volume constraint and design freedom.
  - Under finite strain the bound holds only while the normal block of 𝔸 at F_M remains positive definite. Near a hinge-buckling instability it may not, and then ν_tangent < −1 is possible but physically meaningless.
- **Tangent vs secant.**
  - Eq. 54 is incremental at F_M. For diagonal F_M on an orthotropic cell (load case 1), normal–shear coupling terms vanish. The incremental uniaxial ratios are then 𝔸2211/𝔸2222 and 𝔸1122/𝔸1111, and Eq. 54 is their harmonic mean.
  - For F_M with shear (load case 2), coupling terms such as 𝔸1112 are generally non-zero, and Eq. 54 ignores them.
  - Secant, engineering-strain ratios over a loading path (Clausen; Wang 2018) are different quantities, and numbers are not directly comparable.
- **Bulk measure.** ¼(𝔸1111+𝔸2222+𝔸1122+𝔸2211) is the response of mean nominal stress to an equibiaxial increment. For linear square symmetry it reduces to the 2D bulk modulus (C11+C12)/2.
- **Shear measure.** 𝔸1212+𝔸2121 equals 2C1212 in linear elasticity (minor symmetry). The nominal tangent has only major symmetry and contains initial-stress terms (∂P/∂F ⊃ δ_ik S_JL), so under pre-stretch the sum also includes S11+S22-type contributions.
- **Recommendation.**
  - **KEEP** all three definitions (they define the published objectives).
  - **DOC** them as tangent, plane-strain, in-plane measures, and state the (−1, 1) bound.
  - Optionally add diagnostics, without changing objectives: secant ν from F_M, and incremental uniaxial ν via a static condensation of 𝔸 that includes the shear components.
  - Fix the inconsistent logged `G` in DR to match the full-domain definition, or name the two quantities differently.

### Gaps
- Kazemi et al. 2020's exact wording and justification of the formula were not read.
- The Clausen et al. 2015 and Wang 2018 full texts were not accessed: the Wiley PDF and the Harvard mirror returned HTML.
- Primary verification of the 2D Poisson bounds (Thorpe & Jasiuk; Ting & Chen) is still outstanding.

---

## 6. Periodic homogenization and point-symmetry domain reduction

### Takeaway
Point symmetry (180° rotation, Q = −I in 2D) leaves *every* macroscopic F_M invariant (QF_MQᵀ = F_M). The odd-fluctuation half-cell reduction is therefore valid for arbitrary F_M, including the shear load case, provided:
- the density is point-symmetric;
- the material law is identical at Y and −Y;
- the solution remains on the point-symmetric branch (no symmetry-breaking bifurcation).

Mirror-symmetry (quarter-cell) analysis reductions would *not* be valid under shear. A refactor must preserve the field-parity rules: the fluctuations ũ, χ and the adjoints are odd, and the densities ρ and ρ̃ are even. It must also keep the quadrant-mirrored design parameterization, which both scripts use.

### Cited Findings
- **Article: symmetry and reduction.**
  - Design symmetry is enforced by "taking density field associated with one quadrant as independent design variable and mirroring it to other three quadrants". "Not all point symmetric microstructures result in orthotropic response but all orthotropic microstructures possess point symmetry."
  - Fluctuations satisfy ũ_m = −ũ*_m (Eq. 42), with ũ = 0 at the corner and mid-edge points, and the top and middle edges are coupled antisymmetrically (Eqs. 43–44).
  - "Any periodically fluctuating field over an orthotropic unit cell can be solved for by using this domain reduction technique."
  - Load case 2 (biaxial stretch plus shear) is solved in the reduced domain.
  - The authors state that extension to staggered and point-and-staggered symmetries "is realizible, [20]".
  - [Article] [Secs. 2.3, 3, 5.2, 6](https://doi.org/10.1016/j.finel.2026.104599)
- **Code.** Both mains use quadrant design variables: `N_Elements/4` in the full domain, and `N_Elements/2` on the half mesh in DR, i.e., the same quarter. The DR code uses a separate filter mesh and constraint matrix (`generateMesh_Filter`, `D1_Filter`). [Repo] [NEMOpt_main.jl L229](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR_main.jl L37–38, L228](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR_main.jl)
- **Saavedra Flores & de Souza Neto 2010.**
  - They treat staggered-translational and point symmetries under three kinematic constraints: periodic fluctuations, linear boundary displacements, and minimum constraint (uniform traction).
  - They extend Ohno et al. to the latter two, and report speed-ups above 7× when both symmetries are present in coupled two-scale analyses.
  - [PR] [Eng. Comput. 27(4):551–575](https://doi.org/10.1108/02644401011044612). Abstract content via a search summary of the Emerald and Swansea Cronfa records ([Cronfa](https://cronfa.swan.ac.uk/Record/cronfa14268)). Whether finite strain is covered is unconfirmed.
- **Point-symmetry homogenization and reduced domains** (bibliographic data verified):
  - [Ohno, Matsuda & Wu 2001, "A homogenization theory for elastic–viscoplastic composites with point symmetry of internal distributions", IJSS 38(16):2867–2878](https://doi.org/10.1016/S0020-7683(00)00188-8)
  - [Matsuda, Nimiya, Ohno & Tokuda 2007, Compos. Struct. 79(4):493–500](https://doi.org/10.1016/j.compstruct.2006.02.008)
  - [PR-bib] Ohno et al. show that point symmetry allows half-cell analysis under arbitrary macroscopic loading (not re-verified).
- **Unit-cell and homogenization foundations** (bibliographic data verified; content not re-verified):
  - [Li & Wongsto 2004, Mech. Mater. 36(7):543–572](https://doi.org/10.1016/S0167-6636(03)00062-0)
  - [Li 2008, "Boundary conditions for unit cells from periodic microstructures and their implications", Compos. Sci. Technol. 68(9):1962–1974](https://doi.org/10.1016/j.compscitech.2007.03.035)
  - [Miehe & Koch 2002, Arch. Appl. Mech. 72:300–317](https://doi.org/10.1007/s00419-002-0212-2)
  - [Kouznetsova, Brekelmans & Baaijens 2001, Comput. Mech. 27:37–48](https://doi.org/10.1007/s004660000212)
  - [Blanco, Sánchez, de Souza Neto & Feijóo 2016, Arch. Comput. Methods Eng. 23:191–253](https://doi.org/10.1007/s11831-014-9137-5) (article ref. [23])
- **Wallin & Tortorelli 2020.** Bloch-type displacement BCs, equivalently a prescribed F_M, plus a "detailed discussion regarding domain symmetry of the unit cell and its effect on uniaxial loading conditions". [PR] [abstract](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Bifurcation limit.** Homogenization results "are only valid before the bifurcation occurs". Bloch-wave analysis or enlarged-RVE analysis detects micro-instabilities. [Non-PR preprint] [Zhang, Feng & Khandelwal, arXiv:2010.02371](https://arxiv.org/abs/2010.02371)

### Inferences
- **Why the reduction holds for shear.**
  - With Q = −I, QF_MQᵀ = F_M for any F_M, and the hyperelastic law depends only on C = FᵀF.
  - If ρ(−Y) = ρ(Y), then ũ(−Y) = −ũ(Y) solves the cell problem.
  - Periodicity plus oddness forces ũ = 0 at Y = ±L/2 (corners and mid-edges).
  - ∇ũ is even, so F_m, the tangent 𝔸_m, the χ-problem operators and the adjoint operators are all point-symmetric. This justifies reducing every fluctuation-type field.
  - Reflection R ≠ −I does not satisfy RF_MRᵀ = F_M when F12, F21 ≠ 0. Mirror-based quarter-cell analysis is therefore invalid for load case 2.
- **Orthotropy is stronger than needed.** The article's "orthotropic" condition exceeds what the reduction requires (point symmetry of ρ suffices), and the design is restricted further by double mirroring in *both* scripts.
  - A refactor can separate `design_symmetry` (`:quadrant_mirror` in legacy, or `:point`, `:none`) from `analysis_reduction` (`:full` or `:point_half`).
  - Legacy must keep quadrant mirroring and the specific gradient-summation mapping over the four mirrored elements ([NEMOpt.jl L688–695](/home/user/NEMOpt/Full_Domain/NEMOpt.jl)).
- **Parity rules a refactor must preserve.**
  - Even fields (ρ, ρ̃, the filter adjoint w*): mapping with the + sign (`D1_Filter`).
  - Odd fields (ũ, χ^(ij), η̃*): mapping with the − sign, plus zero constraints at the corner and mid-edge nodes.
  - Half-domain volume integrals: correct weighting when assembling 𝔸, P, the volume fraction and the sensitivities.
  - Mixing these up passes many tests silently, because results are still periodic. Regression tests should compare full vs DR fields for a fixed design under a *shear* F_M.
- **Branch-following risk.** Under large F_M the full cell may bifurcate to a non-point-symmetric mode, while DR enforces the symmetric (possibly unstable) branch. Full-vs-DR agreement is therefore expected only before bifurcation. Single-cell periodicity already excludes period-multiplying modes in both scripts. **DOC** this, and optionally add a smallest-eigenvalue check of the reduced and full tangent stiffness.

### Gaps
- The full text of Saavedra Flores & de Souza Neto 2010 was not accessible. Its explicit finite-strain statements, if any, are unconfirmed.
- Wallin & Tortorelli's specific remarks on symmetry under uniaxial loading are known only from the abstract.
- Wongsto & Li, Miehe & Koch, and Kouznetsova et al. were bibliographically verified only.

---

## 7. Helmholtz filter: parameter meaning and boundary conditions

### Takeaway
NEMOpt solves (r²∇ρ̃·∇w + ρ̃w) = ρw with periodic BCs and r = 1 element. The DTU convention links the PDE length r to the classic cone or density filter radius R by r = R/(2√3). NEMOpt's r = 1 element therefore corresponds to R ≈ 3.5 elements. The article is internally inconsistent: its tables give r = 1, while its text says "2.0 × Lx/100", which would correspond to R ≈ 6.9 elements. Periodic BCs are the right choice for unit cells. KEEP the code convention and DOC the mapping.

### Cited Findings
- **Article Eq. 27 and BCs.**
  - Eq. 27: ∫(r²∇w·∇ρ̃ + wρ̃) = ∫wρ, with r "the 'length scale' controlling the non-locality".
  - "Instead of using Neumann type boundary conditions … volume preserving … periodic boundary conditions are imposed on ρ̃."
  - Tables 3 and 4 list r = 1. Secs. 5.2–5.3 say "the filter radius is 2.0 × Lx/100" on 100 × 100 meshes.
  - Sec. 5.1 states "r = (1*Lx)/Nelx".
  - [Article] [Secs. 2.4, 5](https://doi.org/10.1016/j.finel.2026.104599)
- **Code.** Element matrix `r^2*B'B + N'N`, with `r*Lx/Nelx` passed in and `r = 1.0`. [Repo] [NEMOpt.jl L711, L788](/home/user/NEMOpt/Full_Domain/NEMOpt.jl), [NEMOpt_main.jl L18, L52](/home/user/NEMOpt/Full_Domain/NEMOpt_main.jl), [NEMOpt-DR.jl L416](/home/user/NEMOpt/Domain_Reduction/NEMOpt-DR.jl)
- **DTU 82-line PDE-filter code** (Nov 2010, corrected 2011): `Rmin = rmin/2/sqrt(3)`. It refuses rmin < 2.5 because "for PDE filter the mesh has to be able to resolve the PDE solution". [Non-PR code, DTU TopOpt group] [top82.m](https://www.topopt.mek.dtu.dk/-/media/6a047031ecc24368877a50cf6f8d2805.m)
- **Lazarov & Sigmund 2011.** The source of the approximate R_min = r_min/(2√3) relation, per secondary sources. One notes that the convolution filter regularizes more strongly for large radii, and another that the PDE filter's singular Green's function makes it less robust. [PR-bib] [Lazarov & Sigmund 2011, IJNME 86(6):765–781](https://doi.org/10.1002/nme.3072); secondary: [Non-PR] [arXiv:2501.12756](https://arxiv.org/abs/2501.12756), [arXiv:2507.16108](https://arxiv.org/abs/2507.16108)
- **Ferrari & Sigmund 2020.** Filter BCs are a user choice (`ftBC`: 'N' zero-Neumann or 'D' zero-Dirichlet), so the BC affects results. [PR] [Ferrari & Sigmund 2020](https://doi.org/10.1007/s00158-020-02629-w)
- **Wallin et al. 2020.** PDE-filter boundary effects motivate "padding" or an augmented PDE filter with a tunable boundary penalty. [PR] [Wallin, Ivarsson, Amir & Tortorelli 2020, SMO 62(3):1299–1311](https://doi.org/10.1007/s00158-020-02556-w)
- **Wallin & Tortorelli 2020** use "a periodic filter and thresholding" for unit cells. [PR] [abstract](https://doi.org/10.1016/j.mechmat.2020.103324)

### Inferences
- **Volume preservation is not the distinguishing property.** Testing the weak form with w = 1 gives ∫ρ̃ = ∫ρ for both Neumann *and* periodic BCs. The real difference is that Neumann BCs impose ∂ρ̃/∂n = 0 on the cell edges. That breaks translation invariance and distorts features crossing the cell boundary, whereas periodic BCs filter the cell as a torus. The article's wording ("volume preserving … periodic") should be documented accordingly.
- **Reproduction consequence.** If the published figures used r = 2 elements at 100 × 100 (article text) but the code ships r = 1 at 80 × 80, the minimum feature size differs by about 2×, and the topologies (hinge widths in particular) will differ. The baseline must fix (mesh, r) per figure.
- **Mesh resolution.** r = 1 element (R ≈ 3.46) satisfies the DTU resolution guidance (R > 2.5 elements).
- **Recommendation.**
  - **KEEP** periodic BCs and the r²∇² form.
  - Expose `filter_length` with explicit units (physical length) and an explicit convention (`:helmholtz_r` vs `:cone_R`, with r = R/(2√3)).
  - **DOC** the article's table/text inconsistency, and pin per-case values in the reproduction presets.

### Gaps
- The Lazarov & Sigmund 2011 derivation was not read directly. The 2√3 factor is verified via DTU code and secondary papers only.
- Which r and mesh the authors actually used for each published figure cannot be determined from the article alone.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. The items below correct it after the supplied papers and open-access sources were read in full. Each item gives the corrected statement and the verification note that supports it. All notes are in [NEMOpt supplied papers verification](../NEMOpt%20supplied%20papers%20verification/). Short names: "B&S note" ([bendsoe_sigmund_1999.md](../NEMOpt%20supplied%20papers%20verification/bendsoe_sigmund_1999.md)), "mechanics note" ([open_access_mechanics_sources.md](../NEMOpt%20supplied%20papers%20verification/open_access_mechanics_sources.md)) and "optimization note" ([open_access_optimization_papers.md](../NEMOpt%20supplied%20papers%20verification/open_access_optimization_papers.md)). Line numbers refer to this file as of 2026-10-12.

**Interpolation and Bendsøe & Sigmund (1999)**
- **L36 ("Classic SIMP keeps a lower bound on the design variable").** Corrected: "Classic SIMP keeps a lower bound on the design variable ("0 < ρmin ≤ ρ, is usually imposed, in order to avoid a singular FEM problem", Bendsøe & Sigmund 1999, p. 638). The clip equals classic SIMP evaluated at max(H, ε^{1/p}), with the floor applied to the stiffness only, not to the variable or the volume." The same line's "neither classic nor modified SIMP" and its "SWITCH … new default = modified SIMP" are superseded by the main report: the clip is Wallin & Tortorelli's Eq. 71, and NDR-0033 keeps it as the default. Source: B&S note §4.5(a).
- **L48 (Bendsøe & Sigmund 1999 bibliographic entry).** Volume, year and pages are confirmed from the PDF; the issue number "9–10" is not printed in the PDF. Source: B&S note §4.5(h), claim C4.
- **L53 ("Sigmund 2007 and Andreassen et al. 2011 introduced and popularized "modified SIMP" …").** Corrected: "Sigmund 2007 and Andreassen et al. 2011 popularized "modified SIMP" (E_min + x^p ΔE, E_min ≈ 10⁻⁹E0) as a void regularization (not re-verified). The same algebraic form appears earlier as the two-material SIMP C(ρ) = ρ^p C¹ + (1 − ρ^p)C² ([PR] Bendsøe & Sigmund 1999, Eq. 23, citing Sigmund & Torquato 1997). It replaced the older approach of a lower bound x_min on the design variable with E = x^p E0 ([PR] Bendsøe & Sigmund 1999, p. 638), and it decouples void stiffness from the penalization." Source: B&S note §4.5(b).
- **L54 ([PR-bib] "Bendsøe & Sigmund 1999 discuss physical realizability of SIMP: p must exceed a Poisson-ratio-dependent bound …").** Corrected and upgraded to [PR] (supplied PDF, read in full).
  - SIMP satisfies the Hashin–Shtrikman bounds for a void/solid composite if and only if p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} in 2D (Eq. 11, p. 639; the bulk and shear terms respectively).
  - In 3D the condition is p ≥ max{15(1 − ν⁰)/(7 − 5ν⁰), (3/2)(1 − ν⁰)/(1 − 2ν⁰)} (Eq. 13, p. 640).
  - In 2D, p\* ≥ 3, with equality "only for ν⁰ = 1/3". ν⁰ is the 2D Poisson's ratio (Eq. 9), so NEMOpt's plane-strain solid (ν⁰ = 0.4286) needs p ≥ 3.5; under plane stress (ν⁰ = 0.3) p\* = 3.08, set by shear.
  - Compatibility is necessary for realization. Whether every compatible tensor is realizable "is still an open problem" (p. 639), but microstructures were designed for ν⁰ = 1/3 with p = 3 (Fig. 3) and for ν⁰ = 0 and 0.5 with p = 4 (Fig. 4).
  - Source: B&S note §4.5(c).
- **L55 (RAMP).** Corrected: "RAMP, E = E0 ρ/(1+q(1−ρ)), has non-zero slope at ρ = 0 (Stolpe & Svanberg 2001; Bendsøe & Sigmund 1999 predates and does not discuss RAMP). With q = 2 it coincides with their 2D Hashin–Shtrikman Young's-modulus bound ρE⁰/(3 − 2ρ) (Eqs. 7, 17; algebra)." Source: B&S note §4.5(d).
- **After L68 (§1 Inferences), new inference: bound compatibility of the floors** (computed, plane strain).
  - Against void, εC⁰ exceeds the 2D HS bulk bound for H < 1.75e−3.
  - Against an ε-phase, the clip falls to 0.89 × the HS lower bulk bound at H = 0.0794, and modified SIMP to 0.99 × near H = 0.015.
  - Near H = 0.93 both exceed the upper bulk bound by about 1.6% at p = 3.
  - This is consistent with Bendsøe & Sigmund 1999, p. 647, and is irrelevant for black-and-white designs (p. 651).
  - Source: B&S note §4.5(e).
- **§3 Cited Findings (continuation; the note proposes it after L207), new finding.** [PR] Bendsøe & Sigmund 1999, p. 643: p-continuation "is often advisable", raising p "until the final design is arrived at for a power satisfying (11) or (13)". It is "a compromise, since initial designs will be analysed using an interpolation which is not realizable as a composite structure." NEMOpt's terminal p = 3 is below the plane-strain p\* = 3.5 (computed). Source: B&S note §4.5(f).
- **L332 (plane-strain isotropic solid: Eq. 54 gives the in-plane ratio 0.4286).** Append: "This in-plane ratio is the ν⁰ that enters Bendsøe & Sigmund's 2D SIMP condition (1999, Eqs. 9 and 11)." Source: B&S note §4.5(g).
- **L353 (Gaps: 2D Poisson bounds).** Append: "Bendsøe & Sigmund 1999 evaluate their 2D condition at ν⁰ = ±1 "(in 2-D)" (Eq. 12), consistent with the 2D range −1 ≤ ν ≤ 1. This is indirect support only." Source: B&S note §4.5(g).

**Eq. 54 and Kazemi et al. (2020)**
- **L310 and L314 (Eq. 54 "is attributed to Kazemi et al. [39]"; [PR-bib] Kazemi, Vaziri & Norato 2020).** Correction: the article's attribution, and thesis Eq. 7.5's "(Kazemi et al., 2020)", is incorrect.
  - Kazemi et al.'s only Poisson's-ratio formula is Eq. 24 (p. 8 of the version of record; arXiv 1910.13313 has the same numbering): ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)], for 3D linear homogenization of cubic lattices.
  - The paper contains no 2D formula and no (C₁₁₁₁ + C₂₂₂₂) form.
  - Even for the isotropic FEAD2026 solid in plane strain the two differ: Eq. 54 gives 0.4286, the Kazemi form restricted to in-plane indices 0.300.
  - NDR-0037's decision to keep Eq. 54 is unaffected; only the credit is wrong. Whether to issue an erratum is the authors' decision.
  - Source: mechanics note, Item 2.
- **L351 (Gaps: "Kazemi et al. 2020's exact wording and justification of the formula were not read").** Resolved: the paper was read (arXiv preprint and NSF PAR version of record). Eq. 24 is stated without justification and is not Eq. 54 (see the previous item). Source: mechanics note, Item 2.

**Sigmund (2022)**
- **L218 ("Per a search summary of its abstract …").** Corrected: "Its abstract (Springer page, read 2026-10-12) states that, in the absence of practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors"." Source: optimization note, §3 edit 3.
- **L243 (Gaps: Sigmund 2022 full text not accessed).** Append: "No legal open copy exists as of 2026-10-12 (Springer paywall; DTU Orbit record without file; no arXiv)." Source: optimization note, §3 edit 4.

**Svanberg's GCMMA tolerance**
- **L25 ("far larger than GCMMA's conservativeness tolerance (1e-7)").** Corrected: "far larger than the 1e-7 slack (epsimin) in the conservativeness check of the Python port of Svanberg's GCMMA." Source: optimization note, §5 edit 2.
- **L252 ("within ~1e-7 in Svanberg's code").** Corrected: "within epsimin = 1e-7 in the Python port of Svanberg's code (not checked against Svanberg 2007)". In the port, 1e−7 is primarily the subproblem tolerance, reused as the concheck slack and set by the GCMMA driver. Source: optimization note, §5 edit 3 and claim V-a.
- **L281 (link rot; Svanberg's notes not retrieved).** Append: "Wayback snapshots exist (20210204193223, 20220120004729, 20250316062505; gcmma07.pdf 20230603101120) but `web.archive.org` was unreachable from the research container on 2026-10-12." Do not cite Svanberg 2007 for any parameter value until a snapshot has been read. Source: optimization note, §5 edits 4 and 5.

**Filter length convention (Lazarov & Sigmund 2011)**
- **L420 (DTU convention r = R/(2√3); "R ≈ 3.5 elements").** Corrected: "The DTU convention links the PDE length r to the classic cone or density filter radius R by r = R/(2√3), from matching second central moments on an unbounded domain (Lazarov, Wang & Sigmund 2016, Eq. 10, citing Lazarov & Sigmund 2011). NEMOpt's r = 1 element therefore corresponds to R ≈ 3.46 elements." Source: optimization note, §6 edit 3.
- **L431 ("The source of the approximate R_min = r_min/(2√3) relation").** The letters are reversed against L420. Corrected: "The source of the approximate relation r = R/(2√3) (PDE length r, cone radius R; top82 writes `Rmin = rmin/2/sqrt(3)` with the letters reversed), per Lazarov, Wang & Sigmund 2016 (Eq. 10, ref. [93]); the 2011 full text has no legal open copy." Source: optimization note, §6 edit 4.

**Second pass: Lazarov & Sigmund (2011), read in full** (user-supplied version of record; [LKO note](../NEMOpt%20supplied%20papers%20verification/lazarov_sigmund_2011_kubo_matsuda_2016.md))
- **L420, and the log item for L420 above.** The source is now primary: "Lazarov & Sigmund 2011, Eqs. 30–34, p. 775: r_h = 2√3 r, from 1D kernels on an infinite domain; restated in Lazarov, Wang & Sigmund 2016, Eq. 10". The 2D second-moment match gives R ≈ 3.65 r, so R ≈ 3.46 elements is the convention, not an exact 2D equivalence. Source: LKO note §1.2 and §1.4 edit 6.
- **L431, and the log item for L431 above.** Lazarov & Sigmund 2011 write r_h = 2√3 r, with r the PDE length and r_h the hat radius (Eq. 34); "the 2011 full text has no legal open copy" is superseded, because the user supplied the version of record. The two robustness statements in L431's second sentence come from the arXiv preprints, not from Lazarov & Sigmund 2011, which report more contrast for the Helmholtz filter at equal length scale (p. 775) and do not discuss a singular Green's function. Source: LKO note §1.4 edits 6–7.
- **L437.** Lazarov & Sigmund 2011 prove volume preservation for Neumann BCs only (Eqs. 12–14, p. 769) and say that other BCs require normalization. They do not consider periodic BCs, for which the w = 1 argument holds. Source: LKO note §1.4 edit 8.
- **L446 (Gaps).** Resolved: the derivation was read (Eqs. 30–34). It is one-dimensional. Source: LKO note §1.4 edit 9.

**Second pass: Ohno, Matsuda & Wu (2001), read in full**
- **L383 ([PR-bib] "Ohno et al. show that point symmetry allows half-cell analysis under arbitrary macroscopic loading (not re-verified)").** Corrected and upgraded to [PR]: "Ohno et al. assume point-symmetric constituents, stress and strain at infinitesimal strain, and show that the half cell then holds under any uniform macroscopic loading ("the transformation does not affect macroscopic loading", p. 2870). They verify it against full-cell results under transverse and off-axial loading (p. 2876)." Source: LKO note §3.6 edit 10.
