# Wallin & Tortorelli 2020 ("Nonlinear homogenization for topology optimization") vs. NEMOpt code vs. Güven & Özdemir 2026: equation-level comparison

Scope: full text of M. Wallin, D.A. Tortorelli, *Mechanics of Materials* 145 (2020) 103324 ("W&T"), read from the user-supplied PDF (17 pp.; received 28 Mar 2019, accepted 14 Jan 2020; the PDF states it is "an open access article under the CC BY-NC-ND license", so the earlier notes' statement that no full text exists is superseded). It is compared with the NEMOpt repository (commit `30d0842` code; links relative to this file, `../../` = repo root), with Güven & Özdemir, *Finite Elem. Anal. Des.* 260 (2026) 104599 ("ART"), and, for targeted lookups only, with Güven's 2023 M.Sc. thesis ("THESIS", user-supplied PDF, 98 pp.; printed page = PDF page − 14). Current as of 2026-10-12.

Source key:
- WT = [Wallin & Tortorelli 2020, doi:10.1016/j.mechmat.2020.103324](https://doi.org/10.1016/j.mechmat.2020.103324). Page numbers are the journal page footers (1–17). Equations were read from 220-dpi renders, because the PDF text layer drops most math symbols.
- ART = [Güven & Özdemir 2026, doi:10.1016/j.finel.2026.104599](https://doi.org/10.1016/j.finel.2026.104599)
- THESIS = [Güven 2023, İYTE handle 11147/13904](https://gcris.iyte.edu.tr/handle/11147/13904)
- CALC = my own computations in this session. Scripts: `wt_bench.py`, `wt_fd.py` and `wt_adj.py` preserved in `evidence/wallin_tortorelli_2020/` (run with `python3 -I <script>` from that folder; Python 3.13, NumPy 2.5.3). `wt_bench.py` is also reproduced in full under Q7. The coordinator re-ran all three from the evidence folder and checked the output against the printed W&T values (p. 11 and Table 1, p. 12).

Classification labels used below:
- **MATCH**: NEMOpt matches W&T.
- **DEV-INT**: NEMOpt deviates, intentional per the article or thesis.
- **DEV-DEFECT**: NEMOpt deviates, likely a defect.
- **DEV-UNDOC**: NEMOpt deviates; the value is stated in neither the article nor the thesis, and it is not evidently a defect. This fourth label was added because the three requested labels do not fit these items.

**Consolidated item-by-item classification.** Evidence and quotes are in the sections below.

| # | Item | W&T 2020 | NEMOpt code | ART 2026 | THESIS | Class |
|---|---|---|---|---|---|---|
| 1 | Variable driving the linear/nonlinear blend | H_{β̂,ω̂}(ν̃^p): filtered volume fraction raised to p (Eq. 73, p. 9) | `H_Steep(𝜌^p)`, 𝜌 = filtered ρ̃ at the Gauss point ([NEMOpt.jl:297–299](../../Full_Domain/NEMOpt.jl#L297)) | Eq. 30 omits the argument | H_{β̂,ω̂}(ρ̃^p), Eqs. 5.23–5.24 | **MATCH** (the article omits the argument) |
| 2 | Blend steepness / threshold | β̂ = 500, ω̂ = 0.05 (p. 9) | β₂ = 200, ω₂ = 0.05 ([NEMOpt_main.jl:21](../../Full_Domain/NEMOpt_main.jl#L21)) | (α̂, β̂) = (200, 0.05), Table 1; (β₂, ω₂) = (200, 0.05), Tables 3–4 | β̂ = 200, ω̂ = 0.05, §5.6 | **DEV-INT** (threshold matches; steepness 200 vs 500, with no rationale given anywhere) |
| 3 | Linear stress in voids | D²Ψ̃(I)[F − I] = ℂ[∇u] (p. 9) | `S_Linear = C:∇u` ([NEMOpt.jl:287–291](../../Full_Domain/NEMOpt.jl#L287)) | ℂ:∇u_m, Eq. 30 | D²Ψ̂(I)·(F − I), Eq. 5.23 | **MATCH** |
| 4 | SIMP form | g = max((H_{β,ω}(ν̃))^p, ϵ), Eq. 71 | `g = max(H^p, ϵ)` ([NEMOpt.jl:59](../../Full_Domain/NEMOpt.jl#L59)) | Eq. 29, same form | Eq. 5.20, same form | **MATCH** (max form, not ε + (1−ε)H^p) |
| 5 | SIMP floor ε | ϵ = 10⁻⁵ (p. 13) | ϵ = 5e-4 ([NEMOpt_main.jl:19](../../Full_Domain/NEMOpt_main.jl#L19); DR main :17) | not stated | not stated (only the FD step 10⁻⁶) | **DEV-UNDOC** |
| 6 | Design derivative of the blended stress | δŜ = (g′/g)S_m + pν̃^{p−1} g H′(ν̃^p)(∇Ψ̃ − D²Ψ̃(I)[F − I]), Eq. 75; "δ𝔸_m … is similarly computed" | `Coeff1`, `Coeff2` ([NEMOpt.jl:618–622](../../Full_Domain/NEMOpt.jl#L618)), algebraically identical | text only ("additional terms … would appear") | Eq. 5.26, identical to W&T Eq. 75 | **MATCH** |
| 7 | Projection function and names | H_{β,ω}(x) = [tanh(βω) + tanh(β(x − ω))]/[tanh(βω) + tanh(β(1 − ω))]; β = steepness; ω = 0.5 (Eq. 72, p. 13) | `H(𝜌)` with β, ω = 0.5 ([NEMOpt.jl:53](../../Full_Domain/NEMOpt.jl#L53)) | Eq. 28 renames to α/β and has a typo | Eq. 5.21, cites W&T | **MATCH** |
| 8 | p continuation | 1 → 3, +0.05 "after every 5 design iterations" | FD: +0.5 every 5 evaluations; DR: +0.1 every 5 evaluations | §5.1 +0.2, §5.2 +0.1, §5.3 +0.5, each every 5 iterations | +0.5 every 5 iterations (Ch. 7) | **DEV-INT** |
| 9 | β continuation | β = 1 until p = 3, then "doubled after every 20 design iterations" to 40 | FD: +2 every 10 evaluations while β < 32 (overshoots to 33); DR: +1 every 5 to 32 | additive, +1 or +4 every 5, to 32 | +2 every 10 to 32 | **DEV-INT** (schedule); FD overshoot to 33 is **DEV-DEFECT** |
| 10 | Continuation counter | "design iterations" | objective evaluations, `iter += 1` in `Objective` ([NEMOpt_main.jl:45, 90–99](../../Full_Domain/NEMOpt_main.jl#L45)) | "optimization iteration" | "iterations" | **DEV-DEFECT** |
| 11 | Filter PDE and boundary conditions | ∫(∇w·r²∇ν̃ + wν̃) = ∫wν, periodic, Eq. 74 | r²B'B + N'N, periodic ([NEMOpt.jl:700–752](../../Full_Domain/NEMOpt.jl#L700)) | Eq. 27, periodic | Eq. 5.22 | **MATCH** |
| 12 | Filter radius | r = l/40 = 2 element widths (80 elements) | r = 1 element ([NEMOpt_main.jl:18, 52](../../Full_Domain/NEMOpt_main.jl#L18)) | Tables 3–4: r = 1; §5.2–5.3 text: "2.0 × Lx/100" | r = 1.5 elements (Ch. 7 tables) | **DEV-INT** (code = article tables; the article's text value equals W&T's 2-element choice) |
| 13 | Filter adjoint | one adjoint w per response θ (Eqs. 76–77) | one multi-RHS solve, one column per requested 𝔸 component ([NEMOpt.jl:632–698](../../Full_Domain/NEMOpt.jl#L632)) | "an adjoint field w* for each component" (Eqs. 31–32) | App. A | **MATCH** (equivalent) |
| 14 | Volume constraint | θ₁ = (1/\|Ω_m\|)∫H_{β,ω}(ν̃)dv − 0.4 ≤ 0 (Eq. 96) | (1/(Lx·Ly))∫H(ρ̃) − V_MAX ([NEMOpt.jl:75–89](../../Full_Domain/NEMOpt.jl#L75)) | Eq. 33 lacks 1/Ω_m; Eqs. 53, 55, 57 have it | Eqs. 7.1–7.6 | **MATCH** (the article's Eq. 33 is a typo) |
| 15 | Homogenized tangent | Eq. 24 (one-sided ≡ symmetric) and Eq. 35 (symmetric, total χ) | symmetric, total χ ([NEMOpt.jl:484](../../Full_Domain/NEMOpt.jl#L484)) | Eq. 16, one-sided | §5.4 | **MATCH** (code = W&T Eq. 35) |
| 16 | Unit-strain problems | total χ^(kl) with Bloch constraints (Eqs. 32–34); 9 problems in 3D | 4 total χ with Bloch jump ([NEMOpt.jl:401–462](../../Full_Domain/NEMOpt.jl#L401)) | Eq. 14, fluctuation form, 4 problems | §5.4 | **MATCH** |
| 17 | Rigid-translation fix | "prescribing zero displacement to the center of the unit cell", fn. 3 | center node fixed (`zeros(2, 1)` appended, [NEMOpt.jl:384, 454, 588](../../Full_Domain/NEMOpt.jl#L384)) | §4 | — | **MATCH** |
| 18 | 𝔸 sensitivity form | symmetric, total fields (Eq. 51) | symmetric ([NEMOpt.jl:624](../../Full_Domain/NEMOpt.jl#L624)) | **one-sided Eq. 20: an erratum** (CALC shows it is wrong by O(1)) | symmetric, Eq. 4.25 | **MATCH** (article erratum) |
| 19 | Displacement adjoint | one û per component (ij, rs), Eq. 50 | one `FindTrialU` per component, 3 for ν ([NEMOpt_main.jl:109](../../Full_Domain/NEMOpt_main.jl#L109)) | Eq. 19 | Eq. 4.21 | **MATCH** |
| 20 | Aggregation with ∂J/∂𝔸 weights | none (per component in §2.4; per response function in §4.2) | none | none | none | **MATCH** |
| 21 | Poisson's ratio | tangent, uniaxial-stress, compliance-based: ν₁₂ = −S₂₂₁₁E₁ = −du₂₂ᵖ/du₁₁ᵖ (Eqs. 69, 94), used as constraint ν ≤ ν₀ | `2A[1,2]/(A[1,1]+A[2,2])` at prescribed F_M ([NEMOpt.jl:67](../../Full_Domain/NEMOpt.jl#L67)), objective 100(ν − ν₀)² | Eq. 54, citing [39] | Eq. 7.5, citing Kazemi | **DEV-INT** |
| 22 | Shear measure | μ₁₂ = μ₂₁ = A₁₂₁₂ (Eq. 69) | FD logs A₁₂₁₂ + A₂₁₂₁ ([NEMOpt.jl:71](../../Full_Domain/NEMOpt.jl#L71)); DR logs A₁₂₁₂ | Eq. 53: 1/(A₁₂₁₂ + A₂₁₂₁) | Eq. 7.3: −½(A₁₂₁₂ + A₂₁₂₁) | **DEV-INT** (the DR log = W&T μ₁₂) |
| 23 | Bulk measure | none | ¼(A₁₁₁₁ + A₂₂₂₂ + A₁₁₂₂ + A₂₂₁₁) ([NEMOpt.jl:69](../../Full_Domain/NEMOpt.jl#L69)) | Eq. 56 | Eq. 7.4 | n/a (no W&T counterpart) |
| 24 | Optimizer | "MMA algorithm with the numerical parameter values provided in Svanberg (1987)"; nothing more | FD `MMA87`, DR `MMA02`; s_init = 0.25, maxinner = 40, KKT 1e-8 ([NEMOpt_main.jl:235–237](../../Full_Domain/NEMOpt_main.jl#L235)) | "The default parameters given by [32] are used for MMA" | GCMMA, KKT stop, "default MMA parameters presented by (Svanberg, 1987)", 250 iterations | FD MMA87: **MATCH**; DR MMA02: **DEV-INT** (thesis); s_init and maxinner: **DEV-UNDOC** |
| 25 | Newton solver | not described | 10 load steps, TOL 1e-7·Rep_MU, restart from 0, extra correction, no cap ([NEMOpt.jl:372–396](../../Full_Domain/NEMOpt.jl#L372)) | "incremental-iterative … Newton–Raphson", Eqs. 51–52 | Algorithm 5.2, no TOL or Δt values | no W&T counterpart (see Q6) |
| 26 | FD verification | central FD, 10⁻⁵, no filter/projection, explicit densities, 6 (5) significant figures | none in repo | forward FD, 10⁻⁶, full regularization, densities unpublished, \|Δ\| ≤ 5.6e−8 | Tables 5.1–5.2 (same as article) | **DEV-INT** (different test; see Q7 benchmark) |
| 27 | Energy | Ψ̂(C) = A[(J² − 1) − ln J²] + B[tr C̄ − 3], A = κ/4, B = μ/2 (Eq. 97); 2D test uses B = μ/8 | ψ = (μ/2)(I₁ − 2) − μ ln J + (λ/2)(ln J)² ([NEMOpt.jl:281](../../Full_Domain/NEMOpt.jl#L281)) | Eqs. 22–24 (with the 1/μ typo) | Eq. 5.4, "Following (Bonet & Wood, 2008)" | **DEV-INT** |
| 28 | Kinematic setting | 3D 80 × 80 × 1 eight-node bricks; out-of-plane macro stress free | 2D plane strain, Q4, 2 × 2 Gauss | "plane strain" | "plane strain conditions" | **DEV-INT** |
| 29 | Loading / objective | Bloch displacement u₁₁ᵖ = 0.2l (uniaxial stress); max (S₁₁⁽¹¹⁾S₂₂⁽²²⁾)^½; ν ≤ ν₀; V ≤ 0.4 | prescribed full F_M; min 100(ν + 1)²; V = 0.5 (FD) / 0.4 (DR) | Eq. 55, Table 4 | Eq. 7.6, V 0.5 | **DEV-INT** |
| 30 | Design symmetry and initial design | lower-left quadrant reflected; random initial values | quadrant mirrored ([NEMOpt.jl:822–834](../../Full_Domain/NEMOpt.jl#L822)); random or image | §2.3, §5.1 "same random initial densities" | — | **MATCH** |
| 31 | Density discretization | ν piecewise constant per element; ν̃ ∈ H̃ (nodal) | element ρ, nodal ρ̃ interpolated to Gauss points | — | — | **MATCH** |

---

## Q1. Stress interpolation for void distortion (source of article Eq. 30) and the SIMP factor

### Takeaway
Article Eq. 30 is W&T Eq. 73. W&T drive the linear/nonlinear switch with the **filtered volume fraction raised to p**, H_{β̂,ω̂}(ν̃^p), with **β̂ = 500, ω̂ = 0.05**. NEMOpt's `H_Steep(ρ̃^p)` therefore follows W&T in its argument and threshold and deviates only in steepness (200, stated in the article and thesis without justification). The SIMP factor in W&T is the **clipped** form g = max((H_{β,ω}(ν̃))^p, ϵ) with ϵ = 10⁻⁵. NEMOpt copies the form but uses ϵ = 5e-4, a value that appears in neither the article nor the thesis.

### Cited Findings
- **W&T ersatz material (p. 9, §3.3).**
  - Rationale: "To remedy this situation we take an ersatz material approach wherein a lower bound is placed on the volume fraction field, i.e. 'void' regions have volume fraction ϵ where ϵ > 0 is small."
  - Eq. 70: S_m(F, ν) = g(ν̃)∇Ψ̃_m(F).
  - Eq. 71: g(ν̃) = max((H_{β,ω}(ν̃))^p, ϵ).
  - Eq. 72: H_{β,ω}(x) = [tanh(βω) + tanh(β(x − ω))]/[tanh(βω) + tanh(β(1 − ω))], "is such that lim_{β→∞} H_{β,ω}(x) → H(x − ω). In this way, g(1) → 1 and g(0) → ϵ as β → ∞."
  - [WT p. 9, Eqs. 70–72](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T augmented stress (p. 9).**
  - Motivation: "Unfortunately the above ersatz material is not enough to prevent convergence problems as regions with small volume fractions exhibit large deformations which potentially result in tangled elements with nonphysical deformation gradients characterized by negative determinants. To circumvent these issues we augment the nonlinear strain energy function, Ψ̃_m, such that (70) now reads"
  - **Eq. 73: S_m(F, ν) = g(ν̃)( H_{β̂,ω̂}(ν̃^p) ∇Ψ̃_m(F) + (1 − H_{β̂,ω̂}(ν̃^p)) D²Ψ̃_m(I)[F − I] )**
  - "where β̂ = 500 and ω̂ = 0.05. Effectively we have replaced the nonlinear stress, ∇Ψ̃_m(F), with the 'linear stress,' D²Ψ̃_m(I)[F − I] = ℂ[∇u], in the regions with small penalized volume fraction."
  - [WT p. 9, Eq. 73](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T's stated provenance (p. 2).** "To circumvent the numerical degeneracies that appear when low volume fraction elements suffer extreme distortion we modify the strain energy similar to Wang et al. (2014a)." Wang et al. 2014a is "Interpolation scheme for fictitious domain techniques and topology optimization of finite strain elastic problems", CMAME 276, 453–472 (W&T reference list, p. 17). [WT p. 2, p. 17](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T sensitivity of the blended stress (Eq. 75, p. 9).**
  - Eq. 75: δŜ_m(F_m, ν̃; δν̃) = (g′(ν̃)/g(ν̃)) S_m(F, ν) + p ν̃^{p−1} g(ν̃) H′_{β̂,ω̂}(ν̃^p)(∇Ψ̃_m(F) − D²Ψ̃_m(I)[F − I]).
  - "The variation δ𝔸_m(F_m, ν̃; δν̃) is similarly computed."
  - [WT p. 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T values (p. 13).** "For our computations we use ϵ = 10⁻⁵ and ω = 0.5." [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **NEMOpt code.**
  - `P_rho = g(𝜌)*(H_Steep(𝜌^p)*P + (1 - H_Steep(𝜌^p))*S_Linear)`; `A_rho` uses the same blend with `C_Linear`; `Z_rho = g(𝜌)*H_Steep(𝜌^p)*Z` ([NEMOpt.jl:297–299](../../Full_Domain/NEMOpt.jl#L297)). The DR code is identical ([NEMOpt-DR.jl:507–509](../../Domain_Reduction/NEMOpt-DR.jl#L507)).
  - 𝜌 is the nodal filtered field interpolated to the Gauss point (`𝜌_x = N_Filter*𝜌_Filtered[Ce_Filter]`, e.g. [NEMOpt.jl:412–413](../../Full_Domain/NEMOpt.jl#L412)).
  - `H_Steep` uses (β₂, ω₂) = (200.0, 0.05) ([NEMOpt.jl:56](../../Full_Domain/NEMOpt.jl#L56), [NEMOpt_main.jl:21](../../Full_Domain/NEMOpt_main.jl#L21)).
  - `g(𝜌) = max(H(𝜌)^p, ϵ)`, `g_Deriv = H^p > ϵ ? p*H^(p-1)*H′ : 0.0` ([NEMOpt.jl:59–60](../../Full_Domain/NEMOpt.jl#L59)), with `p, ϵ = 1, 5e-4` ([NEMOpt_main.jl:19](../../Full_Domain/NEMOpt_main.jl#L19)).
  - Sensitivity coefficients: `Coeff1 = g′·H_s(ρ^p) + p·ρ^(p−1)·g·H_s′(ρ^p)`, `Coeff2 = g′ − Coeff1`, `δP = Coeff1*P + Coeff2*S_Linear`, `δA = Coeff1*A + Coeff2*C_Linear` ([NEMOpt.jl:618–622](../../Full_Domain/NEMOpt.jl#L618)).
- **Article.**
  - Eq. 29 is f(ρ̃) = max((H_{α,β}(ρ̃))^p, ϵ).
  - Eq. 30 is P_m(F_m, ρ̃) = f(ρ̃)[H_{α̂,β̂} ∇_F ψ(F_m) + (1 − H_{α̂,β̂}) ℂ : ∇u_m]. H_{α̂,β̂} has no argument.
  - The text: "As a remedy, [16] uses a modified stress expression which is a linear combination of chosen elasticity model and linear elastic response".
  - Table 1 gives (α̂, β̂) = (200.0, 0.05). The ε value is not given anywhere in the article.
  - [ART §2.5, Eqs. 29–30, Table 1](https://doi.org/10.1016/j.finel.2026.104599)
- **Thesis.**
  - Eq. 5.20: g(ρ̃) = max(H_{β,ω}(ρ̃)^p, ϵ). Eq. 5.21: H_{β,ω}, "defined as the smoothed Heaviside step function (Wallin & Tortorelli, 2020)" (p. 38, PDF p. 52).
  - §5.6 Eqs. 5.23–5.24 blend P_m and 𝔸_m with H_{β̂,ω̂}(ρ̃^p), "with β̂ = 200 and ω̂ = 0.05".
  - Eq. 5.26 is W&T Eq. 75 term by term (p. 40, PDF p. 54).
  - No numeric value of the SIMP ϵ appears in the thesis text: a search of the full text found only the FD step "𝜖 = 10⁻⁶" (p. 42).
  - [THESIS §§5.5–5.6](https://gcris.iyte.edu.tr/handle/11147/13904)

### Inferences
- **Algebra check.** NEMOpt's δP = Coeff1·P + Coeff2·S_lin expands to g′[Ĥ P + (1 − Ĥ)S_lin] + p ρ̃^{p−1} g Ĥ′ (P − S_lin). This equals W&T Eq. 75, because (g′/g)·S_m = g′[ĤP + (1 − Ĥ)S_lin]. The same holds for δ𝔸 with (𝔸, ℂ). The `Z_rho = gĤZ` third-order term is consistent, because Ĥ depends on ρ̃ only, not on F.
- **The earlier criticism needs re-attribution.** [numerical_method_choices.md §2](numerical_method_choices.md) recommended switching the void-switch variable to the physical (projected) density because NEMOpt's ρ̃^p departs from Wang et al. 2014 (as restated by Kumar et al. 2021: projected density^p, β₁ = 500, η₀ = 0.01). That departure is W&T's, not NEMOpt's: W&T Eq. 73 uses the filtered ν̃^p. W&T's β̂ = 500 coincides with Wang's β₁ = 500 per that restatement, which suggests W&T took β̂ from Wang 2014a and replaced η₀ = 0.01 by 0.05 and the projected density by ν̃. This is an inference, because Wang 2014a itself was not read.
- **Switch location and width.**
  - The switch center is Ĥ = 0.5 at ν̃^p ≈ 0.05, i.e. ν̃ = 0.05^{1/p}, for both β̂ = 200 and 500, since tanh(β̂ω̂) ≈ 1 in both cases.
  - The 10–90% transition width in ν̃^p is ≈ 2·atanh(0.8)/β̂: 0.011 for β̂ = 200 and 0.0044 for β̂ = 500 (my arithmetic). The 200-vs-500 change therefore mainly smooths the switch. It does not move it.
- **The "soft hyperelastic band" is inherited from W&T.**
  - At the end of W&T's continuation (p = 3, β = 40, ϵ = 10⁻⁵), the clip H = ϵ^{1/3} = 0.0215 lies at ν̃ = 0.5 − atanh(0.957)/40 ≈ 0.452.
  - Elements with ν̃ ∈ (0.368, 0.452) are therefore hyperelastic (Ĥ ≈ 1) at near-minimum stiffness.
  - For NEMOpt (β = 32, ϵ = 5e-4) the band is (0.368, 0.462). The mechanism flagged in numerical_method_choices.md §2 is thus a property of W&T's scheme, not something NEMOpt introduced (my arithmetic).
- **Effect of the larger floor.** NEMOpt's ϵ = 5e-4 is 50× W&T's floor. The clip boundary in H moves from ϵ^{1/p} = 0.0215 to 0.0794 at p = 3. The zero-gradient dead zone (`g_Deriv = 0` below the clip) is larger, and the void phase is stiffer, which raises 𝔸_M slightly and affects ν near −1.
- **Classification.**
  - Blend argument ρ̃^p, ω̂ = 0.05 and the linear stress ℂ:∇u: **MATCH**.
  - β̂ = 200: **DEV-INT** (article Table 1; thesis §5.6), with no stated reason.
  - SIMP max form: **MATCH**.
  - ϵ = 5e-4: **DEV-UNDOC**.
  - A reproduction preset "W&T 2020" should use β̂ = 500 and ϵ = 10⁻⁵. The legacy preset should keep 200 and 5e-4.

### Gaps
- W&T do not state g′ at the clip (Eq. 75 uses g′/g without defining the kink). NEMOpt's choice g′ = 0 below ϵ cannot be compared.
- W&T do not show the blended 𝔸_m or ∂𝔸_m/∂F explicitly ("similarly computed"). That NEMOpt's `A_rho` and `Z_rho` are what W&T used is inferred from Eq. 73, not stated.
- Why the article and thesis use β̂ = 200 instead of W&T's 500 is not explained in either source.

---

## Q2. Homogenized tangent and sensitivity formulation (article Eqs. 11–20): symmetric or one-sided, number of adjoints, total or fluctuation fields

### Takeaway
W&T give 𝔸_M in both one-sided and symmetric forms. These are equal at χ-equilibrium (Eq. 24), and the symmetric total-field form (Eq. 35) is exactly what NEMOpt computes. W&T's final sensitivity (Eq. 51) is **symmetric** in total fields, (E_ij + ∇χ̃^(ij))·δ𝔸[E_rs + ∇χ̃^(rs)] + δS·∇û, with **one displacement adjoint per 𝔸 component** (Eq. 50) and no aggregation by ∂J/∂𝔸. NEMOpt matches W&T Eq. 51 exactly. **Article Eq. 20 is an erratum**: it keeps only E^(mn):δ𝔸:(E^(rs) + ∇χ̃^(rs)). On W&T's own test case that form is wrong by O(1), while the symmetric form reproduces W&T Table 1 to every printed digit (CALC).

### Cited Findings
- **Unit-strain problems in fluctuation form (W&T p. 4).**
  - Eq. 20: ∫𝔸_m(F_m)[∇χ̃^(ij)]·∇û dv = −∫𝔸_m(F_m)[E_ij]·∇û dv for all û ∈ H̃ᵘ_m.
  - "We can satisfy the above for any dF_Mij by solving nine 'unit strain' problems". The "periodic cell … is subjected to the initial strains, E_ij, (or equivalently the incremental initial stresses, 𝔸_m(F_m)[E_ij],)".
  - Eq. 23: 𝔸_M(F_M) = (1/|Ω_m|)∫𝔸_m(F_m)(𝕀 + ∇χ̃^(ij) ⊗ E_ij)dv.
  - [WT p. 4](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Eq. 24 (p. 5), both forms.**
  - 𝔸_Mmnrs = (1/|Ω_m|)∫E_mn·𝔸_m(F_m)[E_rs + ∇χ̃^(rs)]dv = (1/|Ω_m|)∫(E_mn + ∇χ̃^(mn))·𝔸_m(F_m)[E_rs + ∇χ̃^(rs)]dv.
  - "The last line in (24) commonly appears in the literature, cf. e.g. Neves et al. (2000); it follows from the equality ∫∇χ̃^(mn)·𝔸_m(F_m)[E_rs + ∇χ̃^(rs)]dv = 0, cf. (20)."
  - "This result is consistent with the linear result wherein 𝔸_m is independent of F_M."
  - [WT p. 5, Eq. 24](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Total-field (Bloch) form (W&T §2.3, p. 5).**
  - Eq. 32: χ^(kl)(y) = E_kl y + χ̃^(kl)(y) = E_kl y + dũ_m(F_M, y; E_kl), "so that ∇χ^(kl) = E_kl + ∇χ̃^(kl)".
  - Eq. 33: χ^(kl)(y) − χ^(kl)(y − l e_i) = E_kl l e_i on Y₊ᵢ.
  - "here we directly solve (34) for χ^(kl) by imposing the Bloch boundary condition of (33); notably neither initial strain nor initial stresses are prescribed."
  - **Eq. 35: A_Mijkl(F_M) = (1/|Ω_m|)∫E_ij·𝔸_m(F_m)[∇χ^(kl)]dv = (1/|Ω_m|)∫∇χ^(ij)·𝔸_m(F_m)[∇χ^(kl)]dv.**
  - [WT p. 5, Eqs. 32–35](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Total displacement with a Bloch jump.** Eq. 25: u_m(y) − u_m(y − l e_i) = (F_M − I) l e_i; "we directly solve for u_m by imposing the Bloch boundary condition of (26)". [WT p. 5](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Rigid-body fix.** Footnote 3: "In practice this lack of uniqueness is resolved by prescribing zero displacement to the center of the unit cell, i.e. ũ(0) = 0 (De Souza Neto et al., 2015)." [WT p. 4](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Sensitivity derivation (W&T §2.4, pp. 6–7).**
  - Eq. 42 starts from the one-sided component functional θ_A(ν) = A_Mijrs(ν) = (1/|Ω_m|)∫E_ij·𝔸_m(F_m, ν)[E_rs + ∇χ̃^(rs)(ν)]dv, "for the sensitivity of the elasticity tensor we refer to Watts and Tortorelli (2016) where the variations of the individual components (24) are obtained".
  - Two zero integrals are added (Eqs. 39 and 44).
  - Eq. 47, the explicit part, contains **three** terms: E_ij·δ𝔸_m[E_rs + ∇χ̃^(rs)] + δS_m·∇û + δ𝔸_m[∇χ̃^(rs) + E_rs]·∇z.
  - Eq. 48: "this term is annihilated if we equate z = χ̃^(ij)".
  - [WT pp. 6–7, Eqs. 42–49](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Adjoint (W&T Eq. 50, p. 7).**
  - ∫∇δu_m·𝔸_m(F_m, ν)[∇û]dv = −∫∇δu_m·(∂𝔸_m(F_m, ν)/∂F_m)[∇χ̃^(ij) + E_ij, E_rs + ∇χ̃^(rs)]dv for all δu_m ∈ H̃ᵘ_m.
  - "we liken ∂𝔸_m/∂F_m[…] to an initial stress, cf. (20) and (21)."
  - [WT p. 7](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Final sensitivity (W&T Eq. 51, p. 7).**
  - **δθ_A = (1/|Ω_m|)∫[(E_ij + ∇χ̃^(ij))·δ𝔸_m(F_m, ν; δν)[E_rs + ∇χ̃^(rs)] + δS_m(F_m, ν; δν)·∇û]dv.**
  - "Note that if the microscopic stress, S_m, is a linear function of F_m, then ∂𝔸_m/∂F_m = 0 so û = 0 and we recover the result in Watts and Tortorelli (2016)."
  - [WT p. 7](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Per-component adjoint, worked example (W&T p. 7).**
  - For A_M1122: "we solve the linear equation div(𝔸_m[∇û] + ∂𝔸_m/∂F_m[∇χ̃^(11) + E₁₁, E₂₂ + ∇χ̃^(22)]) = 0 … for the periodic adjoint response û".
  - δA_M1122 = (1/|Ω_m|)∫[(E₁₁ + ∇χ̃^(11))·δ𝔸_m[E₂₂ + ∇χ̃^(22)] + δS_m·∇û]dv.
  - Footnote 7: "We do not need χ̃^(11) to evaluate A_M1122, but we do need it to evaluate its sensitivity δA_M1122."
  - [WT p. 7](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Per response function in the displacement formulation.**
  - For uniaxial-stress Poisson's-ratio responses, W&T solve the adjoint for (dû, β̂, dû₂₂ᵖ, dû₃₃ᵖ) (Eq. 92 analog) and for (û, λ̂, û₂₂ᵖ, û₃₃ᵖ) (Eq. 93), then use Eq. 95.
  - Conclusions: "it was shown that the use of the tangent response requires the solution to two adjoint problems. As per usual the adjoint problems are linear and add negligible overall computational cost to the optimization."
  - [WT pp. 12–13, Eqs. 90–95; p. 15](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Filter adjoint (W&T Eqs. 76–77, p. 9).**
  - "We now require the heretofore arbitrary w ∈ H̃ᵘ_m to solve the adjoint problem A(w, δν̃(ν; δν)) = δθ(ν̃(ν); δν̃(ν; δν)) … whereupon we obtain δθ(ν; δν) = B(δν, w)."
  - "As seen here the Helmholtz filter requires two solves."
  - [WT p. 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **NEMOpt code.**
  - χ is a total field: `χ_Slaves = D1*χ_Masters + H` with `H = D2*vcat(unitStrain(k,l)*e[1,:]*Lx, unitStrain(k,l)*e[2,:]*Ly)` ([NEMOpt.jl:451, 458](../../Full_Domain/NEMOpt.jl#L451)), solved for 4 (k, l) ([NEMOpt_main.jl:56](../../Full_Domain/NEMOpt_main.jl#L56)).
  - Tangent: `Aijkl += w*(∇χ⁽ⁱʲ⁾'*A_rho_2D*∇χ⁽ᵏˡ⁾)*detJ`, scaled by 1/(Lx·Ly) ([NEMOpt.jl:484, 491](../../Full_Domain/NEMOpt.jl#L484)).
  - Adjoint RHS: `Z_Cont = Z_rho·∇χ⁽ⁱʲ⁾·∇χ⁽ʳˢ⁾` (total fields); solve with `Bm = -(Fm + D1'*Fs)` and homogeneous periodic BCs ([NEMOpt.jl:542–545, 578–588](../../Full_Domain/NEMOpt.jl#L542)).
  - Sensitivity: `elA_δθ += w*N_Filter'*(∇χ⁽ⁱʲ⁾'*δA*∇χ⁽ʳˢ⁾ + δP'*∇TrialU)*detJ` ([NEMOpt.jl:624](../../Full_Domain/NEMOpt.jl#L624)).
  - One `FindTrialU` per (i, j, r, s) ([NEMOpt.jl:598–599](../../Full_Domain/NEMOpt.jl#L598)). The ν objective requests (1111), (2222), (1122): three displacement adjoints plus one three-column filter solve ([NEMOpt_main.jl:109](../../Full_Domain/NEMOpt_main.jl#L109), [NEMOpt.jl:664](../../Full_Domain/NEMOpt.jl#L664)).
- **Article.**
  - Eq. 16 is (𝔸_M)_mnrs = (1/Ω_m)∫E^(mn):𝔸_m:(E^(rs) + ∇χ̃^(rs))dΩ_m (one-sided). Eq. 19 is the adjoint, with Q = ∂𝔸_m/∂F_m:(∇χ̃^(mn) + E^(mn)):(∇χ̃^(rs) + E^(rs)).
  - **Eq. 20 is δ(𝔸_M)_mnrs = (1/Ω_m)∫[E^(mn):δ𝔸_m:(E^(rs) + ∇χ̃^(rs)) + (δP_m:∇η̃*)]dΩ_m.**
  - The text says δ_χ "can be annihilated, please see [16,28]" and that the article follows "primarily the formalism of [16]".
  - [ART §2.1–2.2, Eqs. 14–20](https://doi.org/10.1016/j.finel.2026.104599)
- **Thesis.** Eq. 4.21 is W&T Eq. 50. Eq. 4.25 is δθ_A = (1/|Ω|)∫∇χ^(ij)·δ𝔸_m·∇χ^(rs) + δP_m·∇û dv (symmetric, total fields). Eq. 4.26 is the linear reduction "which is the same result obtained by (Watts & Tortorelli, 2016)". [THESIS pp. 28–29 (PDF 42–43)](https://gcris.iyte.edu.tr/handle/11147/13904)
- **CALC: one-sided vs. symmetric sensitivity on W&T's §4.1 case.** This is the configuration of Q7 (F_M = diag(1.3, 0.8), W&T 2D energy, SIMP ν³, the adjoint of Eq. 50 with ∂𝔸/∂F by central differences). Group derivatives d𝔸/dν_i, i = 1…4:

  | Component | Symmetric (W&T Eq. 51 = NEMOpt) | One-sided (article Eq. 20) | W&T Table 1 |
  |---|---|---|---|
  | A₁₁₁₁ | 26.548821, 16.46246, 1.850108, 0.596271 | 10.146857, 12.574467, 5.745697, 5.193686 | 26.548821, 16.46246, 1.850108, 0.596271 |
  | A₂₂₂₂ | 67.338406, 16.639579, 46.683352, 25.566853 | 27.879547, 20.915624, 44.901299, 37.730266 | 67.338406, 16.639579, 46.683352, 25.566853 |
  | A₁₁₂₂ | 42.464146, 16.125462, 4.69197, −0.893395 | 16.410031, 13.813575, 9.008626, 5.718293 | 42.464146, 16.125462, 4.69197, −0.893395 |
  | A₁₂₁₂ | 3.886347, 8.312414, −5.343278, −5.044428 | 2.993385, 6.336012, −2.645841, −5.74678 | 3.886347, 8.312414, −5.343278, −5.044428 |

  Source: CALC `wt_adj.py` vs. [WT Table 1, p. 12](https://doi.org/10.1016/j.mechmat.2020.103324).

### Inferences
- **Article Eq. 20 is an erratum.** The missing term ∫∇χ̃^(mn):δ𝔸:(E^(rs) + ∇χ̃^(rs)) is the δ𝔸[∇χ̃^(rs) + E_rs]·∇z contribution of W&T Eq. 47, which survives when z = χ̃^(ij) is used to annihilate δ_χ. The article's Eq. 19 adjoint is the same û as W&T Eq. 50, so it does not absorb that term. The code, the thesis and W&T all use the symmetric form, and the article's Table 2 FD agreement (≤ 5.6e−8) was evidently produced by the symmetric code. The table should therefore be read as validating the code, not the printed Eq. 20.
- **Eq. 16 is harmless.** The one-sided 𝔸 (article Eq. 16) equals the symmetric one only to solver precision (W&T Eq. 24 identity). Golden masters should use the code's symmetric form, as the reproduction-baseline notes already advise.
- **Fields.** W&T derive in fluctuation notation but show that total-field Bloch solves (Eqs. 32–35) are equivalent. NEMOpt implements the total-field variant for u, χ and û (û with homogeneous periodic jump), which is W&T §2.3 verbatim.
- **No aggregation exists in either code.** Neither W&T nor NEMOpt aggregates the displacement adjoint with ∂J/∂𝔸 weights. Because Eq. 50 is linear in its RHS and the operator 𝔸_m is shared, a refactor could solve one adjoint û_J with RHS −Σ_c (∂J/∂A_c) Q_c and accumulate Σ_c (∂J/∂A_c) ∇χ^(c₁)·δ𝔸·∇χ^(c₂). That reduces the three ν adjoint solves (and three filter RHS) to one each. It changes round-off only.
  - W&T's §4.2 approach does the analogous thing at the modeling level: it treats ν₁₂ itself as the response function (2 adjoints per response).
  - This is my inference, not a W&T statement. It is a performance change (class B), not a behavior change.
- **Index order of ∂𝔸/∂F.** W&T write ∂𝔸_m/∂F_m[∇χ̃^(ij) + E_ij, E_rs + ∇χ̃^(rs)], and NEMOpt contracts Z_ijklmn with ∇χ^(ij)_kl ∇χ^(rs)_mn. The two agree because the third derivative of a hyperelastic energy is fully symmetric in its three index pairs. This also holds for the blend, since Z_rho = gĤZ.

### Gaps
- W&T give no cost figures for the adjoint solves beyond "negligible overall computational cost".
- W&T do not discuss the derivative of the max(·) clip, so the dead-zone behavior of NEMOpt's δ𝔸 below ϵ is not addressed by the source.

---

## Q3. Definitions of effective properties (Poisson's ratio, bulk, shear) and bounds

### Takeaway
W&T use a **tangent (incremental), uniaxial-stress, compliance-based** Poisson's ratio at a 20% uniaxially stretched state: ν₁₂ = −dF₂₂/dF₁₁ = −S₂₂₁₁E₁, with the lateral macro stresses free. It is imposed as an inequality **constraint** ν ≤ ν₀ (ν₀ = −0.8 … −0.2). NEMOpt and the article instead use 2𝔸₁₁₂₂/(𝔸₁₁₁₁ + 𝔸₂₂₂₂) at a fully prescribed F_M (Kazemi et al.) as a least-squares **objective** toward ν₀ = −1. This is a deliberate deviation (article Eq. 54; thesis Eq. 7.5) and gives different numbers for the same state. W&T define shear moduli as μ₁₂ = 𝔸₁₂₁₂, define no bulk measure, and discuss no bounds.

### Cited Findings
- **W&T framing (p. 2).**
  - "What is the Poisson's ratio for nonlinear case? One possibility is to define the Poisson's ratio via the secant of the contraction, i.e. the ratio between the total contraction and total elongation as done in Wang et al. (2014b). Another possibility is to consider the incremental rather than total response, i.e. define the Poisson's ratio as in the linear case, but replace the linearized elasticity tensor, ℂ_lin, with the elasticity tensor, 𝔸."
  - "For situations where material instability is important, constraints built upon 𝔸 are the only viable option."
  - "a design tailored for the total response may collapse if an incremental load is applied whereas a design tailored for the incremental response may give unexpected results."
  - Fig. 1 illustrates the tangent and secant definitions.
  - [WT p. 2](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T §4 (p. 9).** "The notion of Poisson's ratio is inherently linked to a specimen in the state of uniaxial stress." [WT p. 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Eq. 69 (p. 8), orthotropic constants from 𝔸 and its inverse.**
  - E₁ = 1/S₁₁₁₁ = A/(A₂₂₂₂A₃₃₃₃ − A₂₂₃₃²), and similarly E₂, E₃.
  - ν₁₂ = −S₂₂₁₁E₁ = (A₁₁₃₃A₂₂₃₃ − A₁₁₂₂A₃₃₃₃)/(A₂₂₃₃² − A₂₂₂₂A₃₃₃₃).
  - ν₁₃ = −S₃₃₁₁E₁ and ν₃₂ = −S₂₂₃₃E₃.
  - μ₂₃ = μ₃₂ = A₂₃₂₃, μ₁₃ = μ₃₁ = A₁₃₁₃, **μ₁₂ = μ₂₁ = A₁₂₁₂**, "where A = det𝔸/(A₁₂₁₂A₁₃₁₃A₂₃₂₃)".
  - "These parameters have the usual meanings, e.g. under a uniaxial stress load increment dS = dS₁₁e₁⊗e₁ we have E₁ = dS₁₁/dF₁₁, ν₁₂ = −dF₂₂/dF₁₁ and ν₁₃ = −dF₃₃/dF₁₁."
  - "the linearized elasticity tensors about such uniaxially deformed configurations are characterized by 3 Young's moduli, 6 Poisson's ratios and 3 shear moduli, albeit they are not all independent."
  - Notation: Eq. 69 writes S for the compliance that the text calls ℂ = 𝔸⁻¹.
  - [WT p. 8, Eq. 69](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T optimization responses (p. 13, Eq. 96).**
  - "To compute, e.g. ν₁₂ = −du₂₂ᵖ/du₁₁ᵖ we first apply the axial displacement, u₁₁ᵖ = u, and then the differential displacement, du₁₁ᵖ = du = 1."
  - Problem: min θ₀ = −(S_M11⁽¹¹⁾S_M22⁽²²⁾)^½ subject to θ₁ = (1/|Ω_m|)∫H_{β,ω}(ν̃)dv − 0.4 ≤ 0, θ₂ = ν₁₂⁽¹¹⁾ − ν₀ ≤ 0, θ₃ = ν₂₁⁽²²⁾ − ν₀ ≤ 0.
  - "where ν₀ is upper limit for the Poisson's ratio … In all cases the inequality constraints are active, hence they are equivalent to equality constraints."
  - [WT pp. 13–14](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T uniaxial-stress enforcement (§4.1–4.2, pp. 10–12).** F_M or the Bloch displacements u_iᵖ are partly unknown and are solved so that S_Mij = 0 for ij ≠ 11 (Eqs. 78–79, 86–87). For orthotropic cells, "we need only consider u₂₂ᵖ and u₃₃ᵖ as additional variables". [WT pp. 10–12](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T conclusion.** "Since the developed procedure accounts for finite deformation and non-linear material response the definition of the Possions' ratio is not unique. In our work we showed that it is feasible to design for constraints on the tangent definition of the Poisson's ratio." [WT p. 15](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Article.**
  - "As opposed to linear case, Poisson's ratio definition is not unique for non-linear elasticity, [16]. In this paper, the following definition of Poisson's ratio [39], ν = 2(𝔸_M)₁₁₂₂/((𝔸_M)₁₁₁₁ + (𝔸_M)₂₂₂₂) (54) is used which is forced to be negative." Objective min 100 × (ν − ν₀)², volume ≤ 0.4.
  - Bulk moduli "are based on macroscopic tangent stiffness tensor which deviates from small strain elasticity tensor as deformation gets larger".
  - [ART §5.2–5.3](https://doi.org/10.1016/j.finel.2026.104599)
- **Thesis.** Eq. 7.5 is ν₁₂ = 2A₁₁₂₂/(A₁₁₁₁ + A₂₂₂₂), "(Kazemi et al., 2020)" (p. 61, PDF 75). Shear objective Eq. 7.3: min −½(A₁₂₁₂ + A₂₁₂₁) (p. 56). Bulk Eq. 7.4: min −¼(A₁₁₁₁ + A₂₂₂₂ + A₁₁₂₂ + A₂₂₁₁), "Although there exists a different definition of bulk modulus in the context of topology optimization, the particular one used here is given below" (p. 59). [THESIS §§7.2–7.4](https://gcris.iyte.edu.tr/handle/11147/13904)
- **Code.** `v_12 = 2*A[1,2]/(A[1,1]+A[2,2])`, `K = 1/4*(A[1,1]+A[2,2]+A[1,2]+A[2,1])`, `G = A[3,3]+A[4,4]` (FD) ([NEMOpt.jl:65–71](../../Full_Domain/NEMOpt.jl#L65)). DR logs `G = A[3,3]` ([NEMOpt-DR.jl:67](../../Domain_Reduction/NEMOpt-DR.jl#L67)). Voigt order is (11, 22, 12, 21) ([NEMOpt.jl:44](../../Full_Domain/NEMOpt.jl#L44)), so A[3,3] = 𝔸₁₂₁₂.
- **CALC: the two definitions on W&T's uniaxial state.**
  - W&T report, at F_M = [1.3, 0; 0, 0.872199] (uniaxial): A₁₁₁₁ = 6.13408, A₂₂₂₂ = 22.4903, A₁₁₂₂ = 8.00914, and ν₁₂ = −E_M1 C_M2211 = 0.356115 ([WT p. 11](https://doi.org/10.1016/j.mechmat.2020.103324)).
  - Eq. 54 applied to the same 𝔸_M gives 2·8.00914/(6.13408 + 22.4903) = 0.5596 (my arithmetic).
  - W&T's value equals A₁₁₂₂/A₂₂₂₂ = 0.35612, because shear–normal coupling vanishes here.

### Inferences
- **What Eq. 54 is relative to W&T.** For a 2D orthotropic tangent without normal–shear coupling, W&T's ν₁₂ = A₁₁₂₂/A₂₂₂₂ and ν₂₁ = A₁₁₂₂/A₁₁₁₁. Eq. 54 is their harmonic mean 2/(1/ν₁₂ + 1/ν₂₁) (consistent with [numerical_method_choices.md §5](numerical_method_choices.md)). It is exact only for square symmetry (A₁₁₁₁ = A₂₂₂₂).
- **Differences in state.**
  - NEMOpt evaluates the tangent at a prescribed biaxial F_M, with lateral stresses not zero.
  - W&T evaluate it at a uniaxial-stress state reached by an unknown lateral stretch.
  - W&T's model is 3D one-layer with a free out-of-plane stretch; NEMOpt's is plane strain (Q8).
  - The ν values are therefore not comparable between the two papers even for identical topologies.
- **Shear.** The DR script's logged G = 𝔸₁₂₁₂ is exactly W&T's μ₁₂. The FD log and the article's Eq. 53 objective use 𝔸₁₂₁₂ + 𝔸₂₁₂₁. W&T Eq. 58 states that 𝔸 "has major symmetry, but no minor symmetries", so the two are different quantities under pre-stress.
- **Bounds.** W&T do not discuss Poisson bounds. Their constraints (ν₀ ≥ −0.8) never approach −1. NEMOpt's target ν₀ = −1 is the 2D isotropic bound, and the open-interval argument in numerical_method_choices.md §5 stands unchallenged by W&T.
- **Classification.**
  - Poisson definition and objective form: **DEV-INT** (article Eq. 54, thesis Eq. 7.5).
  - Shear measure: **DEV-INT** (article Eq. 53, thesis Eq. 7.3). The DR-vs-FD log inconsistency remains a separate logging defect (class B, as in the baseline notes).
  - Bulk: no W&T counterpart.

### Gaps
- W&T give no numeric ν, E or μ values for the optimized designs, only figures. Fig. 10b plots ν trajectories versus u/l, and the text concludes that "the constraints ν < ν₀ = −0.6 are fulfilled at the end state, but not in general at low load levels".
- W&T provide no bulk-modulus definition to compare with ¼(A₁₁₁₁ + A₂₂₂₂ + A₁₁₂₂ + A₂₂₁₁).

---

## Q4. Filter, projection and continuation schedule

### Takeaway
W&T use the same Helmholtz PDE (r² convention), **periodic** BCs, element-wise ν with nodal ν̃, and the same tanh projection (β steepness, ω = 0.5 threshold), applied in g and in the volume constraint. Their radius is **r = 2 element widths** (r = l/40 on 80 elements). Their continuation is **p: 1 → 3 by +0.05 every 5 design iterations; then β = 1 doubled every 20 design iterations to 40**. NEMOpt matches the filter and projection forms. It uses r = 1 element, and its continuation is additive, larger-stepped, ends at 32, and is counted per objective evaluation.

### Cited Findings
- **W&T filter rationale (p. 9).** "To limit the chatter, we filter, i.e. mollify or locally average, the design using the PDE filter proposed by Lazarov and Sigmund (2011) such that ν̃ ← ν". "The size of these intermediate regions is eliminated by thresholding the filtered volume fraction such that intermediate volume fraction values are pushed to either zero or one (Guest et al., 2004)." [WT p. 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Eq. 74 and BCs (p. 9).**
  - ∫_{Ω_m}(∇w·r²∇ν̃ + w·ν̃)dv = ∫_{Ω_m} w·ν dv for all w ∈ H̃ᵘ_m, "where r controls the chatter such that the chatter increases as r → 0".
  - "Rather than the usual homogeneous Neumann condition, ∇ν̃·n = 0, on ∂Ω_m we enforce periodicity on ν̃ due to the periodicity of Ω_m. The former boundary condition conserves volume, i.e. ∫ν̃ dv = ∫ν dv. But, the latter also conserves volume because of the antiperiodic flux ∇ν̃·n that results from the periodicity boundary condition."
  - [WT p. 9, Eq. 74](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T discretization (p. 13).** "The volume fraction field ν is discretized in the usual way as being piecewise uniform over the elements." ν̃ lives in H̃ᵘ_m, the periodic displacement-type space, i.e. it is nodal. [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T parameters (p. 13).**
  - "For our computations we use ϵ = 10⁻⁵ and ω = 0.5."
  - "The penalty exponent is initially equated to p = 1 and incremented by 0.05 after every 5 design iterations until it reaches its terminal value p = 3."
  - "The steepness parameter is initially equated to β = 1 until p = 3 after which time it is doubled after every 20 design iterations until it reaches the terminal value of β = 40."
  - "We apply all the 'tricks' discussed in Section 3.3. The Helmholtz filter is applied with r equal to twice the element width, i.e. r = l/40."
  - [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Code.**
  - Filter element matrix `r^2*B'*B + N'*N`, load `𝜌_Unfiltered[el]*N'`, periodic `D1_Filter` reduction ([NEMOpt.jl:700–752](../../Full_Domain/NEMOpt.jl#L700)); called with `r*Lx/Nelx`, `r = 1.0` ([NEMOpt_main.jl:18, 52](../../Full_Domain/NEMOpt_main.jl#L18)).
  - Projection `H` with `β, ω = 1.0, 0.5` ([NEMOpt.jl:53](../../Full_Domain/NEMOpt.jl#L53), [NEMOpt_main.jl:20](../../Full_Domain/NEMOpt_main.jl#L20)).
  - FD schedule: `iter % 5 == 0 && p < 3 → p += 0.5`; `p == 3 && iter % 10 == 0 && β < 32 → β += 2` ([NEMOpt_main.jl:90–99](../../Full_Domain/NEMOpt_main.jl#L90)).
  - DR schedule: `p += 0.1` every 5, `β += 1` every 5 while β < 32 ([NEMOpt-DR_main.jl:91–99](../../Domain_Reduction/NEMOpt-DR_main.jl#L91)).
  - `iter` is incremented on each `Objective` call ([NEMOpt_main.jl:45](../../Full_Domain/NEMOpt_main.jl#L45)).
- **Article.**
  - Eq. 27 is the same PDE. "Instead of using Neumann type boundary conditions … volume preserving … periodic boundary conditions are imposed on ρ̃" ([ART §2.4](https://doi.org/10.1016/j.finel.2026.104599)).
  - Tables 3–4: r = 1. §5.2–5.3 text: "filter radius is 2.0 × Lx/100".
  - Schedules: §5.1 p +0.2, §5.2 p +0.1, §5.3 p +0.5, all every 5 iterations, then β +1 (§5.1–5.2) or +4 (§5.3) every 5 iterations to 32 ([ART §5](https://doi.org/10.1016/j.finel.2026.104599)).
- **Thesis.**
  - Eq. 5.22 is the same PDE (p. 38).
  - Ch. 7: "Helmholtz radius values are presented in terms of element size, e.g., r = L/80", tables r = 1.5.
  - "For p, it is increased by 0.5 for every 5 iterations until p = 3; after that, β is increased by 2 for every 10 iterations until a terminal value of β = 32" (p. 49, PDF 63).
  - [THESIS Ch. 7](https://gcris.iyte.edu.tr/handle/11147/13904)

### Inferences
- **Filter radius.** Both W&T and NEMOpt multiply the Laplacian by r², so their r values are directly comparable. W&T's r = 2h equals the article's §5.2–5.3 text value "2.0 × Lx/100" (2h on 100 × 100). The article's prose thus matches W&T's convention, while its tables and the code use r = 1h. This strengthens the case that both r = 1 and r = 2 elements must be tested for the article's §5.2–5.3 reproductions.
  - Classification: form and BCs **MATCH**; radius **DEV-INT** (code = article tables; the thesis uses 1.5).
- **W&T's schedule as a timeline** (my arithmetic):
  - p reaches 3 after 40 increments × 5 = 200 design iterations.
  - β then takes 2 (≈220), 4 (≈240), 8 (≈260), 16 (≈280), 32 (≈300), and presumably caps at 40 (≈320). The cap is inferred, because doubling 32 gives 64.
  - W&T's Fig. 9 note, "The large oscillations around design iteration 225 are due to large design changes" ([WT p. 14](https://doi.org/10.1016/j.mechmat.2020.103324)), coincides with the first β doublings.
- **Comparison with NEMOpt.**
  - NEMOpt's p steps (0.5 in FD, 0.1 in DR) are 2–10× larger than W&T's 0.05.
  - NEMOpt ramps β additively to 32, while W&T double β to 40.
  - Classification: **DEV-INT**, because the article and thesis document their own schedules.
- **Counting unit.** All three sources speak of "design/optimization iterations". NEMOpt counts objective evaluations (inner iterations in MMA02, plus the setup evaluation). Classification: **DEV-DEFECT**, consistent with [numerical_method_choices.md §3](numerical_method_choices.md).
- **FD overshoot.** FD's β overshoot to 33 is also **DEV-DEFECT**: W&T, the article and the thesis all state terminal values.
- **Volume constraint.** W&T apply the projected, normalized volume (Eq. 96), which matches [NEMOpt.jl:75–89](../../Full_Domain/NEMOpt.jl#L75). The article's Eq. 33 (no 1/Ω_m) is the outlier.

### Gaps
- W&T do not say whether the doubling stops exactly at 40 (64 > 40). They also do not say whether "design iterations" are MMA outer iterations; plain MMA has no inner loop, so this is likely.
- W&T state no move limits or box bounds for ν.

---

## Q5. Optimizer: MMA or GCMMA, settings, resets, convergence criteria, iteration counts

### Takeaway
W&T say only "MMA … with the numerical parameter values provided in Svanberg (1987)". There is no GCMMA, no move-limit or asymptote settings, no reset at continuation steps and no stopping criterion, and the convergence plot runs 400 design updates. The article's sentence "The default parameters given by [32] are used for MMA" is essentially W&T's sentence. NEMOpt's FD script (MMA87) is consistent with W&T. The DR script's MMA02 (GCMMA) follows the thesis. Neither source documents `s_init = 0.25` or `maxinner = 40`.

### Cited Findings
- **W&T (p. 2).** "The optimization problem is formulated similar to that in Wang et al. (2014b) and solved via the Method of Moving Asymptotes (MMA) (Svanberg, 1987)." [WT p. 2](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T (p. 14).** "The optimization problems are solved using the MMA algorithm with the numerical parameter values provided in Svanberg (1987). In all cases the inequality constraints are active, hence they are equivalent to equality constraints." [WT p. 14](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T initial design (p. 13).** "To avoid uniform sensitivities which confuse the optimizer, the initial volume fraction values in this quadrant are randomly assigned." [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Fig. 9 (p. 14).**
  - The x-axis "Design update" runs 0–400.
  - "Initially, the optimization strives to fulfill the Poisson's ratio constraints which worsens the objective. However, upon satisfying the constraints, the objective reduces. The large oscillations around design iteration 225 are due to large design changes and after that only minor changes in the design takes place."
  - Read by eye, θ₀ settles near −0.08 after ≈250.
  - [WT p. 14, Fig. 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Article.** "The default parameters given by [32] are used for MMA" ([32] = Svanberg 1987) ([ART §5](https://doi.org/10.1016/j.finel.2026.104599)).
- **Thesis.** "the Globally Convergent Method of Moving Asymptotes (GCMMA) (Svanberg, 2002) is used … For the stopping criteria of the optimization, Karush-Kuhn-Tucker (KKT) optimality conditions are used. The default MMA parameters presented by (Svanberg, 1987) are used in MMA. … As an upper limit for the number of iterations, 250 is selected." [THESIS p. 49 (PDF 63)](https://gcris.iyte.edu.tr/handle/11147/13904)
- **Code.**
  - FD: `alg = MMA87()`. DR: `alg = MMA02()`.
  - Both: `MMAOptions(tol=Tolerance(kkt=1e-8, f=0.0), s_init=0.25, s_incr=1.2, s_decr=0.7, maxiter=N_ITER, maxinner=40)`, `convcriteria=KKTCriteria()`.
  - `N_ITER` = 120 (FD) and 300 (DR).
  - [NEMOpt_main.jl:25, 235–240](../../Full_Domain/NEMOpt_main.jl#L235); [NEMOpt-DR_main.jl:23, 234–239](../../Domain_Reduction/NEMOpt-DR_main.jl#L234)
  - The dead `options = MMAOptions(s_init=0.5/(β+1), …)` inside `Objective` is never used ([NEMOpt_main.jl:100](../../Full_Domain/NEMOpt_main.jl#L100)).

### Inferences
- The article's "default parameters" sentence is almost certainly inherited from W&T's wording. It cannot justify the code's non-default `s_init = 0.25` and `maxinner = 40` (NonconvexMMA defaults are 0.5 and 10, per [reproduction_baseline.md Q2](reproduction_baseline.md)). W&T give no information that would.
- Classification:
  - FD MMA87: **MATCH** (W&T plain MMA).
  - DR MMA02: **DEV-INT** (thesis GCMMA).
  - `s_init`, `maxinner`: **DEV-UNDOC**.
  - No MMA reset at continuation steps: consistent with W&T, which mention none.
- W&T's 400 design updates vs. NEMOpt's 120 (FD) and 300 (DR) evaluations: with W&T's slower schedule, continuation alone takes ≈320 iterations. NEMOpt's faster p steps are what let it finish within 120 or 300 evaluations.

### Gaps
- W&T report no run times, no per-iteration costs, no stopping criterion, and no move limits or asymptote parameters beyond the Svanberg (1987) reference.
- Svanberg's 1987 default values themselves were not re-read in this session.

---

## Q6. Nonlinear solver: Newton–Raphson details, load stepping, tolerance, warm start, non-convergence

### Takeaway
W&T do not describe the nonlinear solver at all: no Newton details, load stepping, tolerance, warm start or failure handling. Their only related statement is that small-volume-fraction regions cause "convergence problems", which motivates Eq. 73. NEMOpt's Newton scheme therefore has no W&T reference. The thesis' Algorithm 5.2 is the documented baseline for the code's control flow, and it differs from the code by one extra correction per load step (see [reproduction_baseline.md Q3 item 3, Q5](reproduction_baseline.md)).

### Cited Findings
- **W&T.** Full-text search of the PDF text layer for "Newton", "load step", "tolerance" and "increment" finds only the Poisson-ratio "incremental response" usage and the "convergence problems" sentence on p. 9. For Eq. 11 they say only "To compute S_M for a given F_M(x) we solve the microscopic equilibrium Eq. (11) for ũ_m." [WT pp. 4, 9](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Article.** "Resulting equilibrium equations are solved in an incremental-iterative way using the Newton–Raphson method, see [37]. … R = F^(ext) − F^(int) (51) and eventually the norm of R is expected to drop below a threshold value, i.e., convergence tolerance." [ART §4, Eqs. 51–52](https://doi.org/10.1016/j.finel.2026.104599)
- **Thesis.**
  - Algorithm 5.2: "1: Initialize U, Δt, TOL, |R|; 2: while t < 1 do; 3: t ← t + Δt; 4: H ← H(F_M, t); 5: U ← Reduction(U, H); 6: Assemble K, F; 7: while |R| > TOL do; 8: Solve K·dU = R; 9: U ← U + dU; 10: Assemble K, F; 11: end while; 12: end while; 13: return U".
  - Eq. 5.14: "R = F^Ext − F^Int … where F^Ext = 0 for the problem in hand". No numeric TOL or Δt is given.
  - [THESIS §5.3, p. 36 (PDF 50)](https://gcris.iyte.edu.tr/handle/11147/13904)
- **Code.**
  - `Displacement_Solver(𝜌_Filtered, 1e-7, 0.1)` ([NEMOpt_main.jl:54](../../Full_Domain/NEMOpt_main.jl#L54)).
  - `Rep_MU = abs(maximum(F_M - I))*Lx` ([NEMOpt.jl:375](../../Full_Domain/NEMOpt.jl#L375)).
  - The load loop is `while round(t, digits=2) < tFinal` with `H = D2*vcat(A(t)*L1, A(t)*L2)` and `A(t) = (t/tFinal)*(F_M - I)` ([NEMOpt.jl:373–379](../../Full_Domain/NEMOpt.jl#L373)).
  - The inner loop is `R_Norm = Inf; while R_Norm/Rep_MU > TOL` with the correction applied before the test ([NEMOpt.jl:386–394](../../Full_Domain/NEMOpt.jl#L386)).
  - It starts from U = 0 on each call, and there is no iteration cap.

### Inferences
- With no W&T counterpart, the classification follows the earlier notes:
  - **DEV-INT** for the load-stepped Newton itself (thesis Algorithm 5.2).
  - **DEV-DEFECT** for `Rep_MU = 0` when no entry of F_M − I is positive (an infinite loop) and for the missing iteration cap.
  - Legacy behavior to keep for the extra post-convergence correction (it differs from Algorithm 5.2; see the baseline notes).
- For the Q7 benchmark, any Newton variant converged to round-off gives the same equilibrium, because the hyperelastic cell problem is path-independent before bifurcation. In my reproduction, 10 load steps and ‖R‖ < 1e−12 reproduced W&T's printed digits.

### Gaps
- No source (W&T, article or thesis) gives a solver tolerance definition, a warm-start policy or failure handling. The code is the only specification.

---

## Q7. Finite-difference sensitivity verification, and an independent benchmark for the refactor

### Takeaway
W&T verify on a fully specified 4 × 4 periodic cell: four orthotropic density groups with **given** values, SIMP ν³, no filter or projection, F_M = diag(1.3, 0.8), and central FD with ϵ = 10⁻⁵. They report "agree to 6 significant figures" (non-uniaxial, Table 1) and "agree to 5 significant figures" (uniaxial, Table 2). I reproduced **all** published S_M, 𝔸_M, uniaxial F₂₂, S₁₁, ℂ_M, E₁, ν₁₂ and every Table 1 sensitivity to the printed digits. That required W&T's 2D energy **with B = μ/8 exactly as printed**, Q4 elements, 2 × 2 Gauss and a fixed center node. This makes W&T §4.1 a ready-made **independent benchmark** for a refactored NEMOpt with a pluggable energy and a "pure SIMP" mode. The article's own check (Table 1/2) is different: forward FD 10⁻⁶ with full regularization, and unpublished densities.

### Cited Findings
- **W&T setup (pp. 10–11).**
  - "To validate the sensitivity analysis we performed computations on the two-dimensional 2 × 2 (meshed) unit cell domain illustrated in Fig. 3. The unit cell material is the two-dimensional version of the Neo-Hookean material described in Appendix A, i.e. Ψ̂(C) = A[(J² − 1) − ln J²] + B[(J²)^{−1/2}trC − 2] where A = (λ + μ)/4 and B = μ/8, λ = νE/((1 + ν)(1 − 2ν)) and μ = E/(2(1 + ν)); E = 100 and ν = 0.3 are Young's modulus and Poisson's ratio."
  - "To simplify the sensitivity study the stress is represented as in the SIMP, i.e. S_m(F, ν) = ν³∇Ψ̃_m(F), (85) i.e. we do no thresholding, filtering, etc."
  - Fig. 3 shows the domain [−1, 1]², 4 × 4 elements and 25 nodes, numbered column-wise from the bottom-left: element 1 at x ∈ [−1, −0.5], y ∈ [−1, −0.5]; elements 1–4 form the left column.
  - "To enforce orthotropic symmetry, elements (1, 4, 13, 16) are assigned the uniform volume fraction ν = 0.3, elements (2, 3, 14, 15) are assigned ν = 0.5, elements (5, 8, 9, 12) are assigned ν = 0.7 and elements (6,7,10,11) are assigned ν = 0.9."
  - [WT pp. 10–11, Eq. 85, Fig. 3](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T results, case A (p. 11), F_M = [1.3 0; 0 0.8], not uniaxial.**
  - Formats: S as [S₁₁ S₂₁ S₁₂ S₂₂]ᵀ; 𝔸 rows and columns ordered (11, 21, 12, 22) (footnotes 13–14).
  - S_M = [0.8694 0.0000 0.0000 −1.8120]ᵀ.
  - 𝔸_M = [6.0092 0 0 7.7608; 0 1.1531 0.7871 0; 0 0.7871 −0.9861 0; 7.7608 0 0 28.0699].
  - "The analytical derivative values of 𝔸_M = DS_M, ∂S_M/∂ν_i and ∂𝔸_M/∂ν_i … agree to 6 significant figures with those computed via central finite difference computations using a perturbation of ϵ = 0.00001."
  - [WT p. 11](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T results, case B (p. 11), uniaxial with F_M11 = 1.3.**
  - F_M = [1.3 0 0 0.872199]ᵀ, S_M = [1.43773 0 0 0].
  - 𝔸_M = [6.13408 0 0 8.00914; 0 1.07234 −0.050092 0; 0 −0.050092 −0.074661 0; 8.00914 0 0 22.4903].
  - ℂ_M = 𝔸_M⁻¹ = [0.3047 0 0 −0.108508; 0 0.904204 −0.60665 0; 0 −0.60665 −12.9868 0; −0.108508 0 0 0.083105].
  - E_M1 = 1/0.3047 = 3.28191; ν₁₂ = −E_M1 C_M2211 = 0.356115.
  - Perturbing F_M11 by dϵ₁₁ = 10⁻⁶ gives −dF_M22/dF_M11 = 0.356114, "which agrees with the Poisson's ratio ν₁₂ to within 5 significant figures".
  - "These analytical derivative values agree to 5 significant figures with those computed via finite differences."
  - [WT pp. 11–12, Table 2](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Table 1 (p. 12), ∂/∂ν_i of S_M and 𝔸_M** (group derivatives; full 4 × 4 matrices printed):
  - ν₁: ∂S/∂ν₁ = [1.861827, 0, 0, 1.994601]; ∂𝔸 diagonal (26.548821, 1.959737, 3.886347, 67.338406), (11,22) = 42.464146, (21,12) = 0.857289.
  - ν₂: [3.471355, 0, 0, 2.075059]; diagonal (16.46246, 4.835893, 8.312414, 16.639579), (11,22) = 16.125462, (21,12) = 3.519132.
  - ν₃: [−0.339064, 0, 0, −5.389651]; diagonal (1.850108, 0.267006, −5.343278, 46.683352), (11,22) = 4.69197, (21,12) = 0.857714.
  - ν₄: [0.612492, 0, 0, −3.665842]; diagonal (0.596271, 0.296141, −5.044428, 25.566853), (11,22) = −0.893395, (21,12) = −0.284385.
  - [WT Table 1](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Table 2 (p. 12), uniaxial.**
  - ∂S/∂ν_i = [2.72865], [3.46339], [0.850051], [1.29763] in the S₁₁ entry, with zeros elsewhere.
  - ∂A₁₁₁₁/∂ν_i = 29.3634, 14.884, 1.8568, 0.946049.
  - ∂A₂₂₂₂/∂ν_i = 85.5031, 20.03, 30.0856, 11.9391.
  - ∂A₁₁₂₂/∂ν_i = 47.4552, 13.4836, 5.09922, −0.578267.
  - The ν₁ block has typesetting glitches ("047.4552", "085.5031", a missing zero column).
  - [WT Table 2](https://doi.org/10.1016/j.mechmat.2020.103324)
- **CALC reproduction** (`wt_bench.py`, `wt_fd.py`; B = μ/8, 2 × 2 Gauss, center node (0, 0) fixed, periodic fluctuation, 10 load steps, ‖R‖ < 1e−12):
  - Case A: S_M = [0.86937, 0, 0, −1.81203]. 𝔸_M = [6.0092 0 0 7.76077; 0 1.1531 0.78707 0; 0 0.78707 −0.98606 0; 7.76077 0 0 28.06994].
  - Case B: F₂₂ = 0.8721985, S₁₁ = 1.4377283. 𝔸_M = [6.13408 0 0 8.00913; 0 1.07234 −0.05009 0; 0 −0.05009 −0.07466 0; 8.00913 0 0 22.49033]. ℂ_M diagonal (0.3047, 0.9042, −12.98684, 0.0831), ℂ₁₁₂₂ = −0.10851. E₁ = 3.281912, ν₁₂ = 0.3561146.
  - Central FD (ϵ = 10⁻⁵) on the case-A configuration reproduces **every Table 1 entry to all printed digits** (e.g. ν₁: [1.861827, 0, 0, 1.994601]; ∂A₁₁₁₁ = 26.548821, ∂A₂₂₂₂ = 67.338406, ∂A₁₁₂₂ = 42.464146).
  - The adjoint of Eqs. 50–51 (`wt_adj.py`) gives the same values (Q2 table).
- **CALC controls.**
  - With B = μ/2 (the value consistent with A = (λ + μ)/4; see Q8), case A gives S_M = [2.37173, 0, 0, −6.84968], ν₁₂ = 0.17560. This does not match W&T.
  - With NEMOpt's energy (article Eq. 24, same E, ν), case A gives S_M = [2.86344, 0, 0, −5.93282]. 𝔸_M = [10.30737 0 0 4.82886; 0 4.75225 4.14311 0; 0 4.14311 −0.68347 0; 4.82886 0 0 46.11114]. Uniaxial: F₂₂ = 0.9570835, S₁₁ = 3.514007, E₁ = 8.996065, ν₁₂ = 0.1146009.
  - These NEMOpt-energy numbers are my computation, not published, so they are only a cross-model reference.
- **Euler check (my arithmetic on W&T's printed values).** Because S_m ∝ ν³, S_M and 𝔸_M at fixed F_M are homogeneous of degree 3 in (ν₁…ν₄), so Σᵢ νᵢ ∂(·)/∂νᵢ = 3(·).
  - Table 1: A₁₁₁₁ gives 0.3·26.548821 + 0.5·16.46246 + 0.7·1.850108 + 0.9·0.596271 = 18.0276 = 3 × 6.0092. A₂₂₂₂ gives 84.2098 ≈ 3 × 28.0699. S₂₂ gives −5.4361 ≈ 3 × (−1.8120).
  - Table 2 (F₂₂ is scale-invariant): S₁₁ gives 4.3132 = 3 × 1.43773; A₂₂₂₂ gives 67.471 ≈ 3 × 22.4903.
- **Article.** Table 1: F_M = [1.05 0.05; 0.05 1.05], (L₁, L₂) = (3, 3), (E, ν) = (100, 0.3), r = 2, p = 3, (α, β) = (32, 0.5), (α̂, β̂) = (200, 0.05), 4 × 4 mesh, four design variables (Fig. 6). "one of the design variables is perturbed by ϵ = 10⁻⁶ … first order finite difference". Table 2 lists adjoint − FD differences, maximum 5.59858e−8. [ART §4.1, Tables 1–2](https://doi.org/10.1016/j.finel.2026.104599)
- **Thesis.** Table 5.1 has the same parameters as the article. FD step "𝜖 = 10⁻⁶". "Negligible differences between adjoint sensitivities and finite difference sensitivities" (pp. 41–43, PDF 55–57). [THESIS §5.7](https://gcris.iyte.edu.tr/handle/11147/13904)
- **Repo.** No test or FD script exists in the repository (only `Full_Domain/`, `Domain_Reduction/`, `README.md`, `reports/`, `research_notes/`; local `ls` and `find -iname '*test*'`).

<details><summary>CALC script <code>wt_bench.py</code> (reproduces W&T §4.1; run <code>python3 wt_bench.py wt_mu8</code>; modes <code>wt_mu2</code>, <code>nemopt</code>)</summary>

```python
import numpy as np, sys
E, nu_m = 100.0, 0.3
lam = nu_m*E/((1+nu_m)*(1-2*nu_m)); mu = E/(2*(1+nu_m))
MODEL = sys.argv[1] if len(sys.argv) > 1 else "wt_mu8"
def P_of_F(F):                      # first Piola-Kirchhoff stress
    J = F[0,0]*F[1,1]-F[0,1]*F[1,0]
    FinvT = np.array([[F[1,1], -F[1,0]], [-F[0,1], F[0,0]]])/J
    if MODEL.startswith("wt"):      # W&T 2D: A[(J^2-1)-ln J^2] + B[(J^2)^(-1/2) trC - 2]
        A = (lam+mu)/4; B = mu/8 if MODEL == "wt_mu8" else mu/2
        trC = np.sum(F*F)
        return A*(2*J*J-2)*FinvT + B*(2*F/J - trC/J*FinvT)
    return mu*F - (mu - lam*np.log(J))*FinvT   # NEMOpt / article Eq. 24
def tangent(F):                     # complex-step dP/dF
    h = 1e-30; At = np.zeros((2,2,2,2))
    for k in range(2):
        for l in range(2):
            Fc = F.astype(complex); Fc[k,l] += 1j*h
            At[:,:,k,l] = np.imag(P_of_F(Fc))/h
    return At
n = 4; L = 2.0; h = L/n
xy = np.array([[-1+i*h, -1+j*h] for i in range(n+1) for j in range(n+1)])  # W&T Fig. 3 numbering
nid = lambda i, j: i*(n+1)+j
els = [[nid(i,j), nid(i+1,j), nid(i+1,j+1), nid(i,j+1)] for i in range(n) for j in range(n)]
dens = np.zeros(n*n)
for v, lst in {0.3:[1,4,13,16], 0.5:[2,3,14,15], 0.7:[5,8,9,12], 0.9:[6,7,10,11]}.items():
    for e in lst: dens[e-1] = v
rep = {nid(i,j): nid(i % n, j % n) for i in range(n+1) for j in range(n+1)}   # periodic fluctuation
indep = sorted(set(rep.values())); idx = {nd:k for k,nd in enumerate(indep)}
fixed = idx[nid(2,2)]; ndof = 2*len(indep)                                    # center node fixed
gp = np.array([-1,1])/np.sqrt(3)
def edofs(e): return [d for nd in els[e] for d in (2*idx[rep[nd]], 2*idx[rep[nd]]+1)]
def grads(e, xi, eta):
    dN = 0.25*np.array([[-(1-eta), -(1-xi)], [(1-eta), -(1+xi)], [(1+eta), (1+xi)], [-(1+eta), (1-xi)]])
    Jm = xy[els[e]].T @ dN; return dN @ np.linalg.inv(Jm), np.linalg.det(Jm)
def assemble(FM, u):
    R = np.zeros(ndof); K = np.zeros((ndof, ndof)); Sbar = np.zeros((2,2)); vol = 0
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta); F = FM + ue.T @ dNdX; s = dens[e]**3
                P = s*P_of_F(F); At = s*tangent(F); Sbar += P*dJ; vol += dJ
                for a in range(4):
                    for k in range(2):
                        R[d[2*a+k]] += (P[k,:] @ dNdX[a])*dJ
                        for b in range(4):
                            for m in range(2):
                                K[d[2*a+k], d[2*b+m]] += np.einsum('j,l,jl->', dNdX[a], dNdX[b], At[k,:,m,:])*dJ
    return R, K, Sbar/vol
def solve_red(K, rhs):
    keep = [i for i in range(ndof) if i not in (2*fixed, 2*fixed+1)]
    x = np.zeros(ndof); x[keep] = np.linalg.solve(K[np.ix_(keep,keep)], rhs[keep]); return x
def equilibrium(FM, nsteps=10):
    u = np.zeros(ndof)
    for t in np.linspace(0,1,nsteps+1)[1:]:
        Ft = np.eye(2) + t*(FM-np.eye(2))
        for it in range(50):
            R, K, _ = assemble(Ft, u)
            if np.linalg.norm(R) < 1e-12: break
            u -= solve_red(K, R)
    R, K, S = assemble(FM, u); return u, K, S
Es = [(0,0),(1,0),(0,1),(1,1)]      # W&T order 11, 21, 12, 22
def homog_A(FM, u, K):              # A_M = <(E_ij + grad chi_ij) : A : (E_kl + grad chi_kl)>
    chis = {}
    for (i,j) in Es:
        f = np.zeros(ndof)
        for e in range(n*n):
            d = edofs(e); ue = u[d].reshape(4,2)
            for xi in gp:
                for eta in gp:
                    dNdX, dJ = grads(e, xi, eta); At = dens[e]**3*tangent(FM + ue.T @ dNdX)
                    for a in range(4):
                        for k in range(2): f[d[2*a+k]] += (At[k,:,i,j] @ dNdX[a])*dJ
        chis[(i,j)] = solve_red(K, -f)
    AM = np.zeros((4,4)); vol = 0
    for e in range(n*n):
        d = edofs(e); ue = u[d].reshape(4,2)
        for xi in gp:
            for eta in gp:
                dNdX, dJ = grads(e, xi, eta); At = dens[e]**3*tangent(FM + ue.T @ dNdX); vol += dJ
                G = []
                for (i,j) in Es:
                    Eij = np.zeros((2,2)); Eij[i,j] = 1; G.append(Eij + chis[(i,j)][d].reshape(4,2).T @ dNdX)
                for p_ in range(4):
                    for q in range(4): AM[p_,q] += np.einsum('ij,ijkl,kl->', G[p_], At, G[q])*dJ
    return AM/vol
FM = np.array([[1.3,0],[0,0.8]]); u, K, S = equilibrium(FM)
print("S_M [S11 S21 S12 S22] =", [S[0,0],S[1,0],S[0,1],S[1,1]]); print(homog_A(FM, u, K))
f22 = 0.85                          # uniaxial: Newton on F22 so that S22 = 0
for it in range(30):
    FM = np.array([[1.3,0],[0,f22]]); u, K, S = equilibrium(FM); AM = homog_A(FM, u, K)
    if abs(S[1,1]) < 1e-11: break
    f22 -= S[1,1]/AM[3,3]
CM = np.linalg.inv(AM); print(f22, S[0,0], AM, CM, 1/CM[0,0], -CM[3,0]/CM[0,0])
```
</details>

### Inferences
- **Recommended benchmark (golden master independent of the authors' code).** Implement in the refactor:
  1. a pluggable energy, including W&T's 2D energy with B as a parameter;
  2. a "pure SIMP" interpolation S = ν³P (g = ν^p with no projection, no blend and no filter, i.e. a filter bypass);
  3. W&T's element-group densities.

  Then assert case A and case B to 4–6 significant figures, and Table 1 sensitivities against both the adjoint and a central FD (ϵ = 10⁻⁵) to ~1e−6 relative. This tests all of the following without any unpublished input:
  - the periodic Bloch kinematics and center-node fix;
  - χ and the symmetric 𝔸_M;
  - the Eq. 50 adjoint, including ∂𝔸/∂F;
  - the Eq. 51 sensitivity.

  Table 2 additionally exercises a uniaxial-stress mode: an outer solve for F₂₂ with S₂₂ = 0 and the derivative of that implicit F₂₂. NEMOpt does not have that mode today.
- **Index mapping is required.** W&T's (11, 21, 12, 22) corresponds to NEMOpt's Voigt (11, 22, 12, 21) as W&T index (1, 2, 3, 4) → NEMOpt (1, 4, 3, 2). For example, W&T 𝔸[1,4] = A₁₁₂₂ = NEMOpt `A[1,2]`. Also, W&T's 2 × 2 domain is centered at the origin, while NEMOpt's Generate_Mesh uses its own node numbering. Element-group assignment must be done by element centroid, not by index.
- **Negative tangent component.** W&T's case A has 𝔸₁₂₁₂ = −0.9861 < 0 (compressive F₂₂ = 0.8 gives an initial-stress term). This makes it a good regression test of sign handling in Z and in the tangent assembly.
- **Euler test.** The degree-3 homogeneity identity is a cheap invariant test for any pure-power-law SIMP mode of the refactor. It does not hold with g = max(H^p, ϵ) and the blend.
- **Comparison with the article's check.** W&T verify with central differences (O(ϵ²) truncation) and relative agreement. The article uses forward differences (O(ϵ) truncation) and absolute differences, with densities that cannot be reproduced. The W&T case is the stronger, fully reproducible test. The article's Table 1/2 case should still be kept as a regression test with documented chosen densities, as reproduction_baseline.md proposes.

### Gaps
- W&T do not state their element type, quadrature or Newton tolerance for §4.1. Q4 with 2 × 2 Gauss and a fixed center node reproduces the digits, which strongly suggests that setup but does not prove it.
- Table 2 (uniaxial) sensitivities were checked only through the Euler identity, not by FD. The uniaxial outer solve is in `wt_bench.py`, but `wt_fd.py` was run for Table 1 only.

---

## Q8. Constitutive model and kinematic setting (compressible neo-Hookean form; plane strain vs. plane stress)

### Takeaway
W&T use an isochoric/volumetric-split neo-Hookean energy, Ψ̂(C) = A[(J² − 1) − ln J²] + B[tr C̄ − 3], with A = κ/4 and B = μ/2 in 3D (E = 3, ν = 0.3 in §5). They solve a **3D one-element-thick** cell (80 × 80 × 1 eight-node bricks) with the out-of-plane macro stretch free. Their 2D verification energy prints, and evidently uses, B = μ/8. NEMOpt and the article/thesis use a different compressible neo-Hookean, ψ = (μ/2)(I₁ − 2) − μ ln J + (λ/2)(ln J)² ("Following Bonet & Wood"), in **plane strain**. This is an intentional deviation, and it means no NEMOpt result can match W&T numerically without swapping the energy.

### Cited Findings
- **W&T Appendix A (p. 15).**
  - Eq. 97: "Ψ̂(C) = A[(J² − 1) − ln J²] + B[tr C̄ − 3] = A[(J² − 1) − ln J²] + B[(J²)^{−1/3}trC − 3], where J² = det C so that C̄ = J^{−2/3}C is isochoric". "In this way B(trC̄ − 3) models the isochoric, i.e. volume preserving, deformation energy while A[(J² − 1) − ln J²] models the dilatational energy."
  - Eq. 103 and p. 16: "so together with ∇Ψ̂_{C=I} = 0 we have 𝔸_{F=I} = 4(A − B/3)I ⊗ I + 4B𝕀_S. For isotropy 𝔸 = 2μ𝕀_S + λI ⊗ I, thus we equate A = (3λ + 2μ)/12 = κ/4 and B = μ/2."
  - [WT pp. 15–16, Eqs. 97–106](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T §5 material (p. 13).** "The energy function of (97) is used with A = E/(12(1 − 2ν)) and B = E/(4(1 + ν)) where E = 3 and ν = 0.3." [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T §5 kinematics (p. 13).** "Although the formulation and computations are applicable to three-dimensional structures, for simplicity we consider the two-dimensional case by only considering a single layer of elements in the e₃ direction. Our l × l computational domain is discretized by a 80 × 80 × 1 mesh of Lagrangian eight node displacement based brick elements." The load cases use the §4.2 displacement approach, with u₂₂ᵖ and u₃₃ᵖ unknown so that the macro stress is uniaxial (Eqs. 86–87). [WT pp. 11–13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T 2D test energy (p. 10).** A = (λ + μ)/4, **B = μ/8**, E = 100, ν = 0.3 (quoted in Q7). [WT p. 10](https://doi.org/10.1016/j.mechmat.2020.103324)
- **W&T Eq. 58 (p. 7).** 𝔸(F) = DS(F) = 2𝕀⊠∇_CΨ̂(C) + 4(F⊠I)D²Ψ̂(C)(F⊠I)ᵀ. "Note that 𝔸 has major symmetry, but no minor symmetries." [WT p. 7](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Code / article / thesis.**
  - Code: `P = μ*F - (μ - λ*log(detF))*Finv'`, tangent A(i,j,k,l), λ and μ from E = 100, ν = 0.3 on a 2 × 2 F ([NEMOpt.jl:250, 281–284](../../Full_Domain/NEMOpt.jl#L250)).
  - Article: Eqs. 22–24, with ψ printed as "(1/μ)(I₁ − 2)", a typo for μ/2 per Eq. 24 and the thesis. "a number of plane strain microstructural optimization problems are solved" ([ART abstract, §2.3](https://doi.org/10.1016/j.finel.2026.104599)).
  - Thesis: "Following (Bonet & Wood, 2008), the strain-energy density function" (Eq. 5.4); "The analyzed unit cell is assumed to be a square made of a two-phase neo-Hookean material under plane strain conditions" (pp. 32–33, PDF 46–47) ([THESIS §5.1](https://gcris.iyte.edu.tr/handle/11147/13904)).

### Inferences
- **W&T's printed B = μ/8 is internally inconsistent but is what they computed** (my derivation plus CALC).
  - For the 2D energy, D²Ψ̂ at C = I is (A − B/2)I⊗I + B𝕀_S, so 𝔸_{F=I} = 4(A − B/2)I⊗I + 4B𝕀_S.
  - Matching 2μ𝕀_S + λI⊗I requires B = μ/2 and A = (λ + μ)/4. W&T's A is consistent with B = μ/2.
  - With B = μ/8, the small-strain shear modulus is μ/4 and the Lamé λ is λ + 3μ/4, so the test material is not "E = 100, ν = 0.3" at small strain.
  - Nevertheless the published numbers are reproduced only with μ/8 (Q7). A benchmark implementation must use the printed μ/8 and must not "correct" it.
- **W&T's §5 state is not plane strain.** With one periodic brick layer and u₃₃ᵖ free, F₃₃ is uniform and set by zero mean out-of-plane stress. This is a generalized-plane-strain state with zero net out-of-plane force, not plane strain (F₃₃ = 1) and not pointwise plane stress (my inference from Eqs. 86–87 and the mesh description). NEMOpt's plane strain therefore differs in kinematics as well as energy.
- **The two energies differ even at moderate stretch.** On W&T's 4 × 4 case, NEMOpt's energy gives S_M₁₁ = 2.863 vs. 0.869 (W&T, μ/8) and 2.372 (W&T form with μ/2) (CALC). Even W&T's form with the consistent μ/2 departs from NEMOpt's at 30% stretch: 2.372 vs. 2.863.
- **Classification.** Energy form and plane strain: **DEV-INT** (article Eqs. 22–24; thesis Eq. 5.4, §5.1).

### Gaps
- W&T do not discuss the μ/8 vs. μ/2 discrepancy. Whether it is a typo that propagated into their code cannot be determined from the paper.
- The out-of-plane thickness of W&T's one-layer cell is not stated.

---

## Q9. Example problems and reported results; what can serve as benchmark or verification

### Takeaway
W&T's design examples are qualitatively different from NEMOpt's:
- objective: maximize the geometric mean of two uniaxial-stress load capacities at 20% strain;
- constraints: two tangent Poisson's ratios ν ≤ ν₀ ∈ {−0.8, −0.6, −0.4, −0.2} and V ≤ 0.4;
- discretization: 3D one-layer bricks, E = 3.

They report figures only, with no numeric optimum values or timings. They are therefore not usable as quantitative NEMOpt benchmarks. The only quantitative, fully specified W&T case is the §4.1 verification (Q7), which I reproduced exactly.

### Cited Findings
- **Problem (Eq. 96, p. 13).**
  - "We now optimize the microstructure … to maximize its axial load carrying capacity in the e₁ and e₂ directions. To do this we apply the displacement, u₁₁ᵖ = u = 0.2l … To ensure that the structure can support loads in both directions we maximize the geometric mean, (S_M11⁽¹¹⁾S_M22⁽²²⁾)^½. Material failure due to excessive stress and microscopic- and macroscopic buckling is not considered at this time. Our design is subject to a resource constraint that limits the total volume of the material to 40% of the unit cell volume."
  - Design: lower-left quadrant reflected (orthotropic); 80 × 80 × 1 bricks; r = l/40; ϵ = 10⁻⁵; p and β schedules (Q4).
  - [WT p. 13](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Results (pp. 13–14).**
  - Figs. 4–7 show optimized designs for ν₀ = −0.8, −0.6, −0.4, −0.2: a 3 × 3 array, the undeformed cell, and the cell deformed in the vertical and horizontal directions.
  - Fig. 8 shows ν₀ = −0.6 optimized "under the small strain, i.e. 0.2%".
  - "The microstructural beams optimized for large deformations are initially curved and then straighten during the deformation, whereas the microstructural beams optimized for small strains are initially straight as their deformation is small. This observation is in accordance with the microstructures experimentally tested in Clausen et al. (2015)."
  - Fig. 9 shows objective and constraint histories over 400 design updates.
  - Fig. 10 shows expansion and ν trajectories to u/l = 0.2: "the lateral displacement is almost linear for vertical load case (red curve) but not so for the horizontal load case (blue curve)".
  - Conclusion: "The numerical examples show that it is possible to design microstructures that exhibit ν = −0.8 behaviour at a 20% macroscopic elongation."
  - [WT pp. 13–15](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Symmetry.**
  - "macroscopic shear stress generally arise when the cell is subject to a triaxial extension. To ensure a triaxial stress state, i.e. to eliminate the shear stress, we imposed domain symmetry when designing the unit cell microstructure such that the macroscopic material properties are orthotropic" (p. 14).
  - §3.2: reflections R_i = I − 2e_i⊗e_i make the macroscopic response orthotropic, and "points exhibiting orthotropic symmetry in the initial configuration exhibit linearized orthotropic symmetry with respect to deformed configurations that results from triaxial deformations" (p. 8).
  - [WT pp. 8, 14](https://doi.org/10.1016/j.mechmat.2020.103324)
- **Article cases.**
  - Shear (Eq. 53).
  - ν target −1 (Eq. 55), with F_M = 1.075·I and [1.025 0.025; 0.025 1.025]; reported ν = −0.9008 and −0.9132, and −0.9218 for the linear model.
  - Bulk (Eq. 56), reported κ = 28.4368, 27.9679, 24.1839.
  - Fig. 7 timings.
  - [ART §5](https://doi.org/10.1016/j.finel.2026.104599); details in [reproduction_baseline.md Q2](reproduction_baseline.md).

### Inferences
- **Usable benchmarks from W&T:**
  - §4.1 case A and case B and Table 1, which are exact and reproduced (Q7);
  - Table 2, if a uniaxial-stress outer solve is added;
  - qualitatively, the observation that large-strain designs have initially curved ligaments (Figs. 5 vs. 8). It could be a qualitative check of NEMOpt's large-F_M vs. small-F_M designs. It is not a numeric test.
- **W&T's design examples are not reproducible targets for NEMOpt.** Their objective, constraint type, Poisson definition, kinematics (3D layer), energy (E = 3, split form), radius and schedules all differ (rows 2, 5, 8–9, 12, 21, 27–29 of the table above), and no numeric optima are published.
- **NEMOpt's load case 2 is outside W&T's discussion.** NEMOpt's shear F_M is outside W&T's uniaxial-stress analysis. W&T's symmetry argument concerns triaxial (diagonal) deformations only. The article's domain reduction for shear relies on point symmetry instead (see numerical_method_choices.md §6). W&T is silent on this, so it provides no check.

### Gaps
- W&T publish no numeric optimized objective, ν values, iteration counts beyond Fig. 9, or run times, so no quantitative comparison with the article's −0.9008 / −0.9132 or κ values is possible.
- W&T do not analyze the full (non-reduced) vs. reduced domain. The article's domain-reduction claims have no counterpart in W&T.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. The item below corrects it after Kazemi, Vaziri & Norato (2020) was read in full (arXiv 1910.13313 and the NSF PAR version of record). Supporting note: [open_access_mechanics_sources.md](../NEMOpt%20supplied%20papers%20verification/open_access_mechanics_sources.md), Item 2. Line numbers refer to this file as of 2026-10-12.

- **L41 (row 21: "Eq. 54, citing [39]"; "Eq. 7.5, citing Kazemi"), L214 ("2𝔸₁₁₂₂/(𝔸₁₁₁₁ + 𝔸₂₂₂₂) at a fully prescribed F_M (Kazemi et al.)") and L244 (thesis Eq. 7.5, "(Kazemi et al., 2020)").** These record the article's and the thesis's citations correctly, but the attribution itself is incorrect.
  - Kazemi et al.'s Eq. 24 (p. 8) is ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)], a 3D linear formula for cubic lattices. The paper contains no 2D formula and no (C₁₁₁₁ + C₂₂₂₂) form.
  - Corrected statement: NEMOpt and the article use 2𝔸₁₁₂₂/(𝔸₁₁₁₁ + 𝔸₂₂₂₂) at a fully prescribed F_M (article Eq. 54, thesis Eq. 7.5). Both sources credit Kazemi et al. (2020), whose paper does not contain it. For a square-symmetric linear tensor, Eq. 54 equals the 2D uniaxial-stress ratio C₁₁₂₂/C₁₁₁₁.
  - The deliberate deviation from Wallin & Tortorelli's ν₁₂ stands, and NDR-0037 keeps Eq. 54. Whether to correct the attribution in an erratum is the authors' decision.
