# Michel, Lopez-Pamies, Ponte Castañeda & Triantafyllidis (2007): claim check and benchmark data

Scope: verify the NEMOpt reports' claims that rest on, or bear on, this paper. Extract benchmark data for the stability addendum (k = 0 check, ellipticity of 𝔸_M, Bloch–Floquet module, B-S2). Check NDR-0064's wording on second-order estimates. Recommend edits. Date: 2026-10-12.

## 0. Sources: read this first

**The supplied PDF is not Michel et al. (2007).** It is the companion paper by two of the same authors:

| | Supplied PDF (read in full, 25 pp.) | Michel et al. (2007), the paper the task names |
|---|---|---|
| Authors | O. Lopez-Pamies, P. Ponte Castañeda | J.-C. Michel, O. Lopez-Pamies, P. Ponte Castañeda, N. Triantafyllidis |
| Title | "Homogenization-based constitutive models for porous elastomers and implications for macroscopic instabilities: I—Analysis" | "Microscopic and macroscopic instabilities in finitely strained porous elastomers" |
| Journal | J. Mech. Phys. Solids **55** (2007) **1677–1701** (as printed on the PDF); issue 8 per Crossref | J. Mech. Phys. Solids **55**(5) (2007) **900–938** (Crossref; HAL cover page) |
| DOI | 10.1016/j.jmps.2007.01.007 (printed on the PDF) | 10.1016/j.jmps.2006.11.006 (Crossref; reference list of the supplied PDF) |
| Dates | Received 9 Oct 2006, revised 13 Jan 2007, accepted 18 Jan 2007 | Issued May 2007 (Crossref) |
| Content | 3D random isotropic porous elastomers (spherical voids); second-order (SOH) theory; no FE, no Bloch, no 2D periodic data | 2D plane strain; random, square and hexagonal arrays of circular holes; SOH vs FEM; Bloch analysis; onset-of-failure curves |

- **How Michel et al. was read.** The full text came from the open HAL deposit [hal-00111501v1](https://hal.science/hal-00111501), file `Michel2006.pdf`: a HAL cover page plus 39 manuscript pages numbered 1–39. The journal version runs 39 pages, 900–938. OpenAlex labels the HAL file "submittedVersion".
- **Page citations.** "HAL p. N" refers to the printed manuscript page number. Equation and figure numbers are the paper's own. The text may differ slightly from the journal version; recheck the figure values against the journal PDF before freezing any tolerance.
- **Digitization.** The figures are vector graphics. Curves were extracted from the PDF drawing paths (pdfplumber), with Bézier segments sampled. Axes were calibrated from the frame lines and tick labels. Each curve was intersected with radial load paths. A 600-dpi raster read agreed within 0.002 in log strain.
  - Uncertainty of a digitized point: ±0.002 in log strain (frame calibration plus line width).
  - Source resolution: load step Δλ = 10⁻³ on biaxial-compression paths (HAL p. 17).
  - The plotted first-bifurcation curve in Fig. 4b is a polyline with vertices only at φ = 180°, 195°, …, 270°, so values between those angles are interpolated.
  - The digitization scripts were temporary and are not stored in the repository.
- **Supplied PDF.** It is cited below as "LP&PC 2007-I, p. N" with journal page numbers.

## 1. Claim verification (Task 1)

Grep results:
- **"Michel":** no hits in either report. Five hits in the research notes: `instability_theory_detection.md` L106 and L118; `stability_implementation_verification.md` L225, L352 and L377.
- **"Lopez-Pamies":** one hit in the reports, the main report's NDR-0064 support bullet (L504). No hits in the notes beyond L352.
- **Related claims:** claims about Geymonat et al., Bloch sampling, B-S2 and NDR-0064 that this paper bears on are checked too.

| # | Where | Claim (condensed) | Verdict | Evidence |
|---|---|---|---|---|
| 1 | Main report L504 (NDR-0064 support) | Finite-strain porous-rubber results are "estimates for random microstructures, not bounds for periodic optimized cells" (LP&PC 2004, J. Elasticity 76:247) | **Supported in substance, incomplete** | Michel: the LCC uses a "Hashin–Shtrikman (H–S) estimate (not to be confused with the Hashin–Shtrikman bounds)" (HAL p. 15, Eq. 3.12). LP&PC 2007-I uses a "stationary variational statement" (p. 1681, Eq. 6). But the SOH is *not* limited to random media: Michel apply it with periodic square and hexagonal P tensors (App. B, Eqs. B.1–B.2), "for the first time here … in periodic porous elastomers" (HAL p. 1). Bibliographic data "J. Elasticity 76, 247–287" agree in both reference lists. LP&PC 2004 itself was not read. |
| 2 | Main report L503 | Multi-cell infimum after instability (Müller 1987, ARMA 99:189, "cited from a reference list") | **Supported; link wrong** | Michel: "a fully consistent definition (see Müller, 1987) of the homogenized energy … requires … fluctuations u′ that are kD-periodic" (HAL p. 8, Eq. 2.22). **The report's link resolves to DOI 10.1007/s00205-016-0992-0, which is Duerinckx & Gloria 2016, ARMA 221:1511–1584** (Crossref). Müller 1987 is 10.1007/BF00284506, ARMA 99(3):189–212 (Crossref). |
| 3 | Main report L501 | Ponte Castañeda (1989, Proc. R. Soc. A 422:147) extended the approach to finite elasticity | **Bibliography verified; content partly** | Michel ref. list: "The overall constitutive behavior of nonlinear elastic composites. Proc. R. Soc. A 422, 147–171". Michel: finite-strain homogenization was addressed in papers "starting with Ponte Castañeda (1989)" (HAL p. 2). The bound content itself is not checked. |
| 4 | NDR-0064(b), L509 and L515 | Voigt-type energy bound holds only if the computed state is the minimizer | **Consistent** | One-cell energy is a min over periodic fluctuations (Eq. 2.13); the solution must be a "stable equilibrium" (Eqs. 2.14–2.15); W̃ ≤ Ŵ (Eq. 2.22, HAL p. 9). LP&PC 2007-I: "Beyond the first microscopic instability W̃(F̄) ≤ Ŵ(F̄)" (p. 1681). |
| 5 | Addendum L13 and L56 | β ≤ β(k→0) ≤ B; micro precedes or coincides with macro; loss of rank-one convexity equals long-wave buckling. Currently "verified as relayed" by Zhang et al. only | **Supported; stronger second source** | Eq. 2.30: B̂ ≥ β₀ ≥ β; "microscopic stability … implies macroscopic stability" (HAL p. 10). "B̂(F) and β₀(F) always vanish simultaneously" (HAL p. 10). Michel is co-authored by Triantafyllidis, a co-author of Geymonat et al. 1993. |
| 6 | Addendum L58 | A long-wave mode is only approached as k → 0; the surface is discontinuous at the origin; it belongs to the 𝔸_M test | **Supported** | "the limit ω → 0 is a singular one depending on the ratio of the ω components" (HAL p. 10). A long-wave onset "is always detectable as a loss of strong ellipticity of the one-cell homogenized moduli" (HAL p. 10). |
| 7 | Addendum L3 and L13 | Homogenized moduli are meaningful only before the first bifurcation | **Supported** | After a finite-wavelength onset, "the one-cell homogenization is no longer physically meaningful" (HAL p. 11). |
| 8 | Addendum L46 (B-S1) | Assert only t_micro ≤ t_macro | **Consistent** | "microscopic bifurcation precedes the macroscopic loss of ellipticity for periodic microstructures" (HAL p. 32). Equality also occurs: "For parts of the mixed loading region … the first instability encountered is a long wavelength instability" (HAL p. 22). So "≤" is correct and a strict "<" would be wrong. |
| 9 | Addendum L66 and L3 | A strongly elliptic constituent cannot lose ellipticity under compression, so any loss of ellipticity of 𝔸_M is structural | **Consistent (general)** | Matrix laws are "poly-convex" (HAL p. 19). The homogenized energy shows "consistent loss (at adequately large strains) of the rank-one convexity … independently of the microgeometric pore arrangement" (HAL p. 31). |
| 10 | Addendum L129 | TRIM points cover the 2×1, 1×2 and 2×2 period-doubling modes | **Consistent; Michel's critical mode is a TRIM point** | Square, f₀ = π/16, θ = 0: "the critical wavenumbers are found to be either ω₁L₁ = ω₂L₂ = π or ω₁L₁, ω₂L₂ → 0" (HAL p. 17). |
| 11 | Addendum L131 | Grid scan needed; contour not guaranteed | **Consistent** | Microscopic stability "requires tedious and time consuming calculations since one has to scan using a fine grid" the full ω-cell (HAL p. 10). Michel use no contour. Their full Bloch scan was done "for only one periodic elastomer", on "a relatively coarse mesh" (HAL p. 17). |
| 12 | Addendum L127 and NDR-0050 | No center pin for θ ≠ 0 | **Not addressed; one notational trap** | Michel define D# as "zero-average displacement functions" (HAL p. 7) and reuse D# for the Bloch amplitude u′ (Eqs. 2.25–2.26). Taken literally, a zero-mean u′ at ω ≠ 0 is a constraint like the pin. It would also exclude the long-wave modes whose vanishing must coincide with B̂ = 0 (HAL p. 10), so it cannot be intended (inference). Do not transcribe it. |
| 13 | Addendum L179 (B-S2) (a) | First instability at the M point θ = (π, π), a 2×2 pattern | **Supported for porosity π/16 (≈ 0.196), not shown for 0.5–0.6** | M-point (2×2) onset with "an antisymmetric mode" (HAL p. 22). Compressive mode: "An alternating pattern … based on void ovalization" (HAL p. 22). No micro onset is computed for f₀ ≠ π/16. |
| 14 | Addendum L179 (B-S2) (b) | Ellipticity is lost later | **Supported** | Fig. 4b; "The loss of ellipticity … occurs at higher strains along the same path" (HAL p. 22). Values in §2. |
| 15 | Addendum L179 (B-S2) (c), and impl. note L368 (d) | "below porosity ≈ 0.34, localization first", read as macroscopic first | **Contradicted for Michel's material** | At porosity π/16 < 0.34 with κ/μ = 10, the 2×2 mode comes first under uniaxial-strain and equibiaxial compression: ε = −0.154 vs −0.177, and −0.0785 vs −0.0945 (Fig. 4b). "A bifurcation with a finite wavelength eigenmode is the first instability encountered under biaxial compression" (HAL p. 22). Bertoldi's "localization" (nearly incompressible Rivlin, finite specimens) is not tested here. |
| 16 | Addendum L179 (B-S2) | Quantitative acceptance rests on non-peer-reviewed slides | **Superseded** | Michel give peer-reviewed onset curves with a fully specified material and geometry (§2 below), in figures only. |
| 17 | `instability_theory_detection.md` L106 | Bibliographic data; "compares periodic FE Bloch analyses with second-order homogenization estimates of macroscopic loss of ellipticity" | **Bibliography verified; content imprecise** | Bibliography: Crossref and HAL. The SOH-vs-FEM comparisons are of *macroscopic loss of ellipticity*, with FEM moduli from Eq. 2.18 (Figs. 7–12). The first-bifurcation curves come from "a 2 × 2 cell aggregate subjected to periodic boundary conditions" (HAL p. 22). The Bloch scan was a one-off (HAL p. 17). |
| 18 | `instability_theory_detection.md` L118; `stability_implementation_verification.md` L225 and L377 | Full text not retrievable; critical values not extractable | **Superseded** | Read from HAL. Critical values exist in figures only (no tables) and are digitized in §2. |
| 19 | `stability_implementation_verification.md` L352 | JMPS 55(5):900–938 | **Verified** | Crossref; HAL cover page. |
| 20 | `stability_implementation_verification.md` L353 | Triantafyllidis, Nestorović & Schraad 2006: "failure surfaces of porous neo-Hookean solids with circular holes in strain space" | **Partly supported** | Michel: TNS 2006 applied onset-of-failure surfaces "to periodic porous and particulate elastomers" (HAL p. 3), and Michel's FEM algorithm "details … are given in Triantafyllidis et al. (2006)" (HAL pp. 16–17). Michel do not state TNS 2006's matrix law, so "neo-Hookean" remains unverified. |
| 21 | Supplied paper LP&PC 2007-I, p. 1679 (not cited in the reports) | Michel et al. "have shown that the second-order theory can deliver accurate estimates … of macroscopic instabilities in porous elastomers with periodic microstructures" | **Over-generalized relative to Michel's own results** | Michel find good agreement for square arrays at f₀ ≈ 0.2 but "an unsatisfactory overestimation for the hexagonal microstructure" at f₀ = 0.5 (HAL p. 27). With θ = π/10 loading, "the SOH method considerably overpredicts" (HAL p. 31). Do not cite LP&PC 2007-I as evidence of SOH accuracy for arbitrary periodic cells. |

## 2. Benchmark data (Task 2)

### 2.1 Common setup (Michel et al., §§3–4)

- **Geometry.** 2D plane strain with cylindrical holes of circular cross-section (Fig. 1b).
  - Square lattice: L₂ = L₁, one hole per cell. f₀ = π/16 means R = L/4.
  - Hexagonal lattice: rectangular period L₁ × √3·L₁ with two holes, close-packed rows along X₁ (Fig. 1b). FE cells are a body-fitted square and a hexagon.
  - Porosities: π/16 ≈ 0.196 (most results), 0.1, 0.3 and 0.5 (Figs. 8b–10).
- **Matrix material (Eq. 4.2, HAL p. 18).** Compressible neo-Hookean: W = (μ/2)[(I − 2) − 2 ln J] + ((κ − μ)/2)(J − 1)², with I = F_ij F_ij and J = det F, where "μ and κ are … the shear and bulk moduli of the solid at zero strain".
  - All onset curves use **κ/μ = 10** ("neo-Hookean matrix material with κ/μ = 10", HAL p. 21).
  - Linearizing Eq. 4.2 gives W ≈ μ|ε|² + ((κ − μ)/2)(tr ε)². So λ_Lamé = κ − μ and κ = λ + μ is the plane-strain bulk modulus. κ/μ = 10 therefore means **ν = 0.45** (inference, by linearization).
  - Gent variant: J_m = 50 (Fig. 11). Its compressive onset curves show "no appreciable difference" from neo-Hookean (HAL p. 30).
- **FE discretization (HAL p. 17).**
  - Elements: Q4 bilinear with 2 × 2 Gauss points.
  - Square cell: 1800 elements, 3840 DOF. Hexagonal cell: 2700 elements. Imperfect 4 × 4 aggregate: 59,970 DOF.
  - "further mesh refinement does not result in appreciably different onset-of-failure curves" (HAL p. 17).
  - Newton stopping criterion ‖S‖ ≤ 10⁻⁴‖S̄‖. Step Δλ = 10⁻³ (biaxial compression) and 10⁻² otherwise. 360 path angles.
- **Loading (Eqs. 3.1–3.2).** Strain-controlled, proportional in principal log-strain space: ln λ₁ = λ cos φ and ln λ₂ = λ sin φ, with principal axes at angle θ to the lattice (θ = 0 unless stated).
  - "Uniaxial" means **uniaxial strain** (ε₂ = 0, φ = π), not uniaxial stress ("uniaxial (ε₂ = 0) loading", HAL p. 19).
  - Equibiaxial compression is φ = 5π/4.
  - NEMOpt's ramp F(t) = I + t(F_M − I) passes through the same states for F_M = diag(s, 1) and F_M = s·I, since both are radial in log strain. In general, the principal-branch state at a given F does not depend on the path, as long as no bifurcation is crossed (inference).
- **Definitions.**
  - Macro onset: first zero of B̂ = min over unit a and n of a_i n_j L̂_ijkl a_k n_l (Eq. 2.29), with n in the reference configuration and L̂ = L^H from the χ-problems (Eqs. 2.18–2.20). These are NEMOpt's 𝔸_M and Q(n).
  - Micro onset: first zero of β (Eq. 2.26). In practice this was computed on a 2 × 2 aggregate with periodic boundary conditions (HAL p. 17).
  - Void-surface instability: an analytic half-space (Biot) check on the surface stretch (App. A.1). It is not captured by FE.

### 2.2 Square lattice, f₀ = π/16, κ/μ = 10 (FEM, Fig. 4b; SOH, Fig. 8a)

Values are the load parameter λ = radius in log-strain space, digitized to ±0.002.

| Path | First bifurcation (micro) | Mode | Loss of ellipticity of 𝔸_M (macro) | Void-surface instability | SOH estimate of macro onset (periodic P) | Order |
|---|---|---|---|---|---|---|
| Uniaxial strain, φ = π (ε₂ = 0) | **ε₁ = −0.154** (stretch 0.857; engineering −0.143) | M point (π, π), 2 × 2 antisymmetric; voids ovalize with alternating orientation | **ε₁ = −0.177** (stretch 0.838) | ε₁ = −0.170 | ε₁ = −0.171 | micro < surface < macro |
| Equibiaxial, φ = 5π/4 | **ε₁ = ε₂ = −0.0785** (λ = 0.111; stretch 0.9245) | Same M-point mode (Fig. 4b inset) | **ε₁ = ε₂ = −0.0945** (λ = 0.134; stretch 0.910) | ε = −0.123 | ε = −0.093 | micro < macro < surface |
| Intermediate φ (micro curve vertices) | λ = 0.127 (φ = 195°), 0.115 (210°), 0.111 (225°) | M point | λ = 0.149 (195°), 0.136 (210°), 0.134 (225°) | λ = 0.159 (195°), 0.160 (210°) | n/a | micro first throughout the compressive quadrant |

- **Direction at macro onset.** The critical reference normal satisfies "ζ = 0, π for π ≤ φ ≤ 5π/4 but ζ = π/2 for 5π/4 ≤ φ ≤ 3π/2" (HAL p. 25, Fig. 6b). So n ∥ X₁ for compression along X₁, which tests the ellipticity module's argmin.
- **One-cell and k ≠ M modes.** For this solid the only critical wavenumbers were (π, π) or ω → 0 (HAL p. 17). So Γ (one-cell), X and Y stay stable through the micro onset (inference from that statement).
- **Imperfect 4 × 4 aggregates.** Imperfection amplitude ξ = 0.1, three shapes (Fig. 5b). Their loss of ellipticity falls at ε₁ = −0.151 to −0.155 (uniaxial strain) and ε = −0.079 to −0.080 (equibiaxial). This "practically coincides" with the perfect-lattice micro onset (HAL p. 24).
- **Mixed loading (Fig. 4a).**
  - Micro and macro coincide (long-wave first) near φ ≈ 125°–135° and the mirror range 315°–325°. At φ = 135° (ε₂ = −ε₁): onset λ ≈ 0.566 ± 0.005, i.e. ε₁ ≈ −0.400, ε₂ ≈ +0.400.
  - A void-surface instability comes earlier there (λ ≈ 0.445; "the surface instability appears prior to any other failure", HAL p. 22).
  - By φ = 170°–180° the micro onset is clearly first (0.185 vs 0.205; 0.154 vs 0.177).
- **Equibiaxial tension (φ = 45°).** Micro λ = 1.383, macro 1.63. These strains are far outside NEMOpt's range.

### 2.3 Other cases (macro onset only; no micro data in the paper)

Values are log strain at loss of ellipticity, digitized to ±0.002.

| Case | Uniaxial strain along X₁ (φ = π) | Along X₂ (φ = 3π/2) | Equibiaxial (per axis) | Source |
|---|---|---|---|---|
| Square, f₀ = 0.5, FEM | **−0.090** (stretch 0.914) | −0.090 | **−0.060** (λ = 0.085) | Fig. 8b |
| Square, f₀ = 0.5, SOH | −0.108 | −0.108 | −0.0735 | Figs. 8b, 9a |
| Square, f₀ = 0.5, surface instability | −0.299 | −0.299 | −0.203 | Fig. 8b |
| Square, f₀ = 0.1 / 0.3, SOH only | −0.181 / −0.153 | same | −0.0905 / −0.089 | Fig. 9a |
| Hexagonal, f₀ = π/16, FEM | −0.168 | −0.198 | −0.126 | Figs. 5b, 8a |
| Hexagonal, f₀ = π/16, SOH | −0.241 | off chart | −0.146 | Fig. 8a |
| Hexagonal, f₀ = π/16, surface instability | **−0.153 (precedes macro)** | −0.186 | −0.140 | Fig. 8a |
| Hexagonal, f₀ = 0.5, FEM | −0.218 | −0.221 | −0.126 | Fig. 8b |
| Hexagonal, f₀ = 0.5, SOH | n/a | n/a | −0.254 | Fig. 8b |
| Random, f₀ = π/16, SOH | no loss of ellipticity on the axis within the plot; the curve spans only φ ≈ 185°–265° and ends near λ ≈ 0.36, where pore closure begins | same | −0.145 | Fig. 5 |

- **Porosity trend.** In strain space the square-lattice SOH curves are "nested one inside the other according to increasing initial porosity" (HAL p. 28). In stress space, critical stresses fall with porosity for every lattice (Fig. 10).
- **Stress data.** Normalized nominal stresses at onset (S₁/μ, S₂/μ) are given in Figs. 6a and 10. They were not digitized here.

### 2.4 What the paper does not provide

- **No tables.** All values are in figures.
- **Micro onset for one case only.** The first-bifurcation (micro) curve exists only for the square lattice at f₀ = π/16, θ = 0. **There is no micro onset at porosity 0.5–0.6**, the B-S2 range.
- **No uniaxial-stress paths.** Bertoldi-type uniaxial-stress compression with lateral traction-free boundaries does not appear.
- **No sensitivity studies.** No ν (κ/μ) sensitivity and no mesh-convergence numbers are given.

### 2.5 Usability for NEMOpt

| Issue | NEMOpt | Michel | Consequence |
|---|---|---|---|
| Volumetric term | λ/2 (ln J)² (code: P = μF − (μ − λ ln J)F⁻ᵀ, NEMOpt.jl L281) | (κ − μ)/2 (J − 1)² | Same linearization. Finite-strain pressure differs, e.g. ln J/J = −0.0314 vs J − 1 = −0.030 at J = 0.97 (inference). The shift in the critical strain is unquantified. Option: a **test-only** variant P = μF − (μ − λJ(J − 1))F⁻ᵀ with tangent μδ_ikδ_jl + λ(2J − 1)J F⁻ᵀ_ij F⁻ᵀ_kl + (μ − λJ(J − 1)) F⁻ᵀ_il F⁻ᵀ_kj (derivation, not checked in code). It must stay outside FEAD2026 (R-S0). |
| Compressibility | FEAD2026 ν = 0.3, so κ/μ = 2.5 | κ/μ = 10, i.e. ν = 0.45 | FEAD constants cannot reproduce Michel. Use ν = 0.45 in the Default profile. |
| Hole geometry | Staircase pixels; void = ρ = 0 on the linear branch with ε = 5e−4·E | Body-fitted circle, true void | At f₀ = π/16 the hole is well resolved: R = 20 px at 80×80, pixel porosity 0.1975 vs 0.1964, ligament 40 px. At f₀ = 0.5 the ligament is 0.20·L (16 px at 80×80). The surface-instability criterion is meaningless on staircase edges, and corners may introduce mesh-dependent local modes. The residual void stiffness delays onset slightly. All of these shifts are unquantified. |
| Element | Q4, full integration | Q4, 2 × 2 Gauss | Same element family. Mild volumetric locking at ν = 0.45 applies to both. |
| Cell shape | Square periodic cell | Square, and L × √3L for hexagonal | **Square lattice: directly usable.** Hexagonal: not representable on a square pixel cell (√3 aspect), so low usability. Random: not applicable. |
| Loading | F_M prescribed (strain control) | Strain control, radial log-strain paths | Matches for uniaxial strain and equibiaxial. Needs the compressive Newton options (NDR-0031 and NDR-0032). |
| Quantities | k = 0 inertia, ellipticity of 𝔸_M, TRIM, grid | Micro (2 × 2 = TRIM), macro (L^H), mode, argmin n | All B-S2a targets are computable on NEMOpt's one-cell principal branch: macro onset after micro onset is still defined on that branch, as Michel compute it. |

**Assessment.**
- **B-S2a (square, f₀ = π/16, κ/μ = 10, uniaxial strain and equibiaxial): high value.** It is a peer-reviewed, fully specified case that checks the M-point onset, the Γ/X/Y non-criticality, the macro onset, the micro-before-macro ordering and the argmin n.
  - Acceptance should combine Michel's digitized values (±0.002, plus the unknown source discretization error) with a NEMOpt pixel-refinement trend (80, 100, 160 and 200 px) and an independent body-fitted computation.
  - Expect deviations of a few percent from the volumetric term and the staircase boundary. That magnitude is an estimate, not a measurement.
- **B-S2b (square, f₀ = 0.5): medium value.** It has macro targets only (ε = −0.090 and −0.060), so it can assert t_micro ≤ t_macro but not the micro value.
- **Bertoldi porosity 0.5–0.6: qualitative only.**
- **Mixed path φ = 135°: low to medium value.** It is a long-wave-first check (β and B vanish together), but at large strains (ε ≈ ±0.40) with pixel holes.

## 3. SOH estimates vs FE analyses of periodic cells (Task 3)

1. **The comparison is not random-medium SOH against FE Bloch analyses.** Michel compare SOH estimates *specialized to the same periodic arrays* (square and hexagonal P tensors, App. B) with FEM loss of ellipticity of the one-cell homogenized moduli. Random-medium SOH appears only alongside, as a microstructure comparison. Micro (Bloch-type) onsets come only from FEM: "the FEM is the only method at this point that can be effectively used for the calculation of the microscopic instabilities" (HAL p. 32).
2. **Agreement varies widely (figures digitized).**

   | Case | SOH vs FEM | Paper's verdict |
   |---|---|---|
   | Square, f₀ = π/16 | −3% (uniaxial strain), −2% (equibiaxial) | "remarkably good in the biaxial compression quadrant" (HAL p. 27) |
   | Square, f₀ = 0.5 | +20% (uniaxial strain), +22% (equibiaxial) | "a reasonable estimate" (HAL p. 27) |
   | Hexagonal, f₀ = π/16 | +44% (φ = π), +15% (equibiaxial) | — |
   | Hexagonal, f₀ = 0.5 | ≈ +100% (equibiaxial) | "an unsatisfactory overestimation" (HAL p. 27) |
   | Larger porosities generally | — | SOH "overestimates the critical stresses" (square; HAL p. 28) and "substantially overestimates the critical strains" (hexagonal; HAL p. 29) |
   | Loading at θ = π/10 | — | "the SOH method considerably overpredicts" (HAL p. 31) |

   Overall verdict: "sometimes with great accuracy" (HAL p. 32).
3. **Random vs periodic.** In biaxial compression "always the failure along a given path happens at larger load parameters in the random medium" (HAL p. 32). Random curves lie "outside (and often at considerable distance)" of the periodic ones (HAL p. 23).
4. **Micro before macro.** For periodic cells the true onset is the micro bifurcation, which precedes macro onset (HAL p. 32). Even a perfect macro estimate therefore overstates the stable range of a periodic cell.
5. **LP&PC 2007-I adds two points.**
   - At small strain the second-order estimates "recover the classical HS upper bounds" (3D voids; p. 1689, §4.3).
   - The HS-type LCC estimate "may become inaccurate for large values of f₀" (pp. 1684–1685).
   - Its summary of Michel (p. 1679) is more favorable than Michel's own results (claim 21).

**Implication for NDR-0064.**
- The core wording ("estimates …, not rigorous bounds"; comparison only) is **supported and should be strengthened**, not reversed.
- "For random microstructures" is too narrow: SOH estimates exist for simple periodic arrays too. The decisive facts are these:
  - they are estimates;
  - they predict only macroscopic loss of ellipticity;
  - they can be non-conservative for periodic cells by large factors;
  - random-microstructure estimates fail later than periodic arrays.
- No closed-form estimate exists for NEMOpt's optimized topologies, because the LCC needs a microstructural tensor P for the actual geometry; Michel give P only for circular holes (App. B; inference about applicability).
- Edits are in §5.

## 4. Bloch method details and conflict check (Task 4)

**Method in Michel (§2.2, §3.3).**
- Coercivity constant (Eq. 2.26): β = inf over ω of the min over Bloch modes v = u′e^{iω·X} of the Rayleigh quotient ∫L v_{i,j}v̄_{k,l} / ∫ v_{i,j}v̄_{i,j}. The denominator is the **gradient norm** of v, not an L² or "mass" norm.
- Wave vectors: 0 ≤ ω_i L_i < 2π.
- k = 0: β_D (Eq. 2.15) is the one-cell (D-periodic) stability, i.e. ω = 0.
- Long-wave limit (Eq. 2.28): β₀ = lim inf as ω → 0, a singular limit (HAL p. 10). β₀ = 0 if and only if B̂ = 0, so the long-wave onset is detected by ellipticity of L̂ (HAL p. 10).
- Ordering (Eq. 2.30): B̂ ≥ β₀ ≥ β.
- **Sampling in practice.**
  - The method calls for a fine grid over the whole ω-cell (HAL p. 10).
  - Actual practice: one full Bloch scan on a coarse mesh, for one solid (square, π/16, θ = 0) (HAL p. 17). Then a "2 × 2 cell aggregate with the final refined mesh and periodic boundary conditions" was used for all first-bifurcation curves (HAL p. 17).
  - No contour, symmetry reduction or TRIM terminology appears. The 2 × 2 aggregate is equivalent to the four TRIM points (Γ, X, Y, M).
- **Micro vs macro comparison.**
  - Both are computed along the same radial paths.
  - Micro comes first in biaxial compression and in equibiaxial tension; the two coincide in parts of mixed loading (HAL p. 22).
  - A 4 × 4 aggregate with ξ = 0.1 imperfection turns the perfect-lattice bifurcation into a macroscopic loss of ellipticity of the aggregate at nearly the same strain (HAL p. 24).
- **Other limits on the principal solution** (App. A): void-surface instability (analytic), percolation (FEM threshold d_min/l_max ≤ 0.005), lock-up (moduli ratio ≥ 10) and pore self-contact.

**Conflict check against the planned NEMOpt implementation.**

| Plan item | Michel | Verdict |
|---|---|---|
| TRIM points (Level 1) | Critical k for the canonical case is M or k → 0; the 2 × 2 aggregate is the production method | **No conflict; supports it.** Michel justified this only after a full scan, for one geometry and aligned loading. NEMOpt's grid on final designs plays the same role. |
| Contour labeled "not guaranteed" (Level 2) | Not used | No conflict. |
| Grid scan on final designs (Level 3) | Required in principle ("fine grid") | No conflict. Michel's coarse-mesh scan followed by a refined supercell is a cheap precedent for confirming a Level-3 result at the critical k. |
| No center pin for k ≠ 0 | Not discussed. The literal D# ("zero-average") in Eqs. 2.25–2.26 would act like a pin. | **Flag:** do not copy Michel's D# into the Bloch problem; keep NDR-0050 as written. |
| k = 0 check (Tier 0a) | β_D, Eq. 2.15 | Consistent. |
| Long-wave limit via ellipticity of 𝔸_M | β₀ = 0 ⟺ B̂ = 0 (HAL p. 10); FEM L̂ = L^H | **Consistent.** |
| E-S4: "λ_min(K(θ))·N_el²/\|θ\|² → min eig Q(n) of 𝔸_M" | With gradient normalization only B̂ ≥ β₀ holds (inequality, simultaneous zeros) | **Not a conflict but a citation caveat.** The E-S4 equality is a Bloch-wave-homogenization statement under NEMOpt's uniform-node (Euclidean) metric, checked numerically in the addendum. It is not a Geymonat or Michel result. Cite it as such. |
| Argmin direction | n in the reference configuration; ζ = 0 at φ = π (square) | Add as a test of the ellipticity module's argmin (B-S2a). |
| Void-surface instability | Separate analytic check; FE cannot capture it | Not in the NEMOpt plan. It is not meaningful for pixel holes. Note it as a limitation of B-S2 on staircase boundaries. |

## 5. Recommended edits (Task 5; not applied)

**E1. Main report L504 (NDR-0064 support bullet).**
- Old: `- The finite-strain results for porous rubbers are estimates for random microstructures, not bounds for periodic optimized cells ([Lopez-Pamies & Ponte Castañeda 2004, J. Elasticity 76:247](https://dialnet.unirioja.es/ejemplar/108150)).`
- New: `- The second-order finite-strain results for porous rubbers use a Hashin–Shtrikman-type *estimate* for the linear comparison composite, "not to be confused with the Hashin–Shtrikman bounds", so they are estimates, not bounds ([Michel et al. 2007, §3.2.1](https://doi.org/10.1016/j.jmps.2006.11.006); verified from the HAL author version). Derived for random microstructures ([Lopez-Pamies & Ponte Castañeda 2004, J. Elasticity 76:247–287](https://dialnet.unirioja.es/ejemplar/108150), bibliographic data verified; [2007, JMPS 55:1677–1701](https://doi.org/10.1016/j.jmps.2007.01.007), read), they have also been specialized to square and hexagonal arrays of circular holes. There they estimate only macroscopic loss of ellipticity, and against finite elements they ranged from "remarkably good" (square, porosity ≈ 0.2) to "an unsatisfactory overestimation" (hexagonal, porosity 0.5). Random-microstructure estimates lose ellipticity later than periodic arrays of equal porosity ([Michel et al. 2007, §§4.3, 5](https://doi.org/10.1016/j.jmps.2006.11.006)). They are therefore not conservative for periodic cells, and no such estimate exists for NEMOpt's optimized topologies (inference).`

**E2. Main report L503 (Müller link).**
- Old: `([Müller 1987, ARMA 99:189](https://link.springer.com/doi/10.1007/s00205-016-0992-0), cited from a reference list)`
- New: `([Müller 1987, ARMA 99:189–212](https://doi.org/10.1007/BF00284506); restated in [Michel et al. 2007, Eq. 2.22](https://doi.org/10.1016/j.jmps.2006.11.006))`
- Reason: the old link is Duerinckx & Gloria 2016.

**E3. Main report L515 (NDR-0064 row, part c).**
- Old: `**(c)** Cite nonlinear Hashin–Shtrikman-type bounds and estimates for comparison only, never as acceptance tests.`
- New: `**(c)** Cite nonlinear Hashin–Shtrikman-type bounds and second-order estimates for comparison only, never as acceptance tests. The second-order results are estimates, not bounds; they predict only macroscopic loss of ellipticity, which a periodic cell's microscopic bifurcation precedes; and for periodic arrays they overestimated finite-element onset strains by up to about a factor of two (Michel et al. 2007, hexagonal array, porosity 0.5; digitized), so they never certify stability.`

**E4. Addendum L56.**
- Old: `Their 1993 full text could not be read, so these statements are **verified as relayed** by Zhang et al.'s restatement ([Geymonat et al. 1993](https://doi.org/10.1007/BF00380256); [Zhang et al. 2021, Eqs. 55–56](https://arxiv.org/abs/2010.02371)).`
- New: `Their 1993 full text could not be read, so these statements are **verified as relayed** by two restatements. Michel et al., a paper co-authored by Triantafyllidis, give B̂ ≥ β₀ ≥ β and state that B̂ and β₀ "always vanish simultaneously" ([Michel et al. 2007, Eq. 2.30 and §2.2](https://doi.org/10.1016/j.jmps.2006.11.006); HAL author version), and Zhang et al. give the same chain ([Zhang et al. 2021, Eqs. 55–56](https://arxiv.org/abs/2010.02371)); see [Geymonat et al. 1993](https://doi.org/10.1007/BF00380256).`
- Also, after `Microscopic instability therefore precedes macroscopic instability or coincides with it.` insert: `Porous elastomers show both cases: in a square array of circular holes the 2×2 mode comes first under biaxial compression, while in parts of mixed tension–compression loading the first instability is long-wave and the two onsets coincide ([Michel et al. 2007, §4.3, Fig. 4](https://doi.org/10.1016/j.jmps.2006.11.006)).`

**E5. Addendum L58.**
- Old: `so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); verified).`
- New: `so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); verified; [Michel et al. 2007, §2.2](https://doi.org/10.1016/j.jmps.2006.11.006): the limit is "a singular one depending on the ratio of the ω components").`

**E6. Addendum L131.** After the sentence ending `with three refined strips ([Zhang et al. 2021](https://arxiv.org/abs/2010.02371); verified).` insert:
- `Michel et al. used a similar economy for a square array of circular holes: one Bloch scan on a coarse mesh found critical wave numbers only at (π, π) or k → 0, after which a refined 2×2 aggregate, equivalent to the four TRIM points, produced every first-bifurcation curve. The precedent covers that geometry and aligned loading only ([Michel et al. 2007, §3.3](https://doi.org/10.1016/j.jmps.2006.11.006)). Their Bloch space is written with zero-average amplitudes; that constraint must not be transcribed, since at θ ≠ 0 it acts like the center pin (inference).`

**E7. Addendum L179 (B-S2 row; replace the whole row).**
- New: `| B-S2 | Square array of circular holes, strain-controlled. **B-S2a:** porosity π/16 (R = L/4); compressible neo-Hookean W = (μ/2)[(I − 2) − 2 ln J] + ((κ − μ)/2)(J − 1)², κ/μ = 10 (NEMOpt: ν = 0.45 and a test-only (J − 1)² volumetric variant); F = diag(s, 1) and F = s·I. **B-S2b:** porosity 0.5, same material and paths. **B-S2c:** porosity 0.5–0.6, qualitative | B-S2a: Γ, X, Y stay stable; the first instability is at the M point, a 2×2 antisymmetric mode with alternating void ovalization, at log strain ε₁ = −0.154 (uniaxial strain; stretch 0.857) and ε₁ = ε₂ = −0.0785 (equibiaxial); 𝔸_M loses ellipticity later, at ε₁ = −0.177 and ε₁ = ε₂ = −0.0945, with critical normal n ∥ X₁ for compression along X₁. B-S2b: 𝔸_M loses ellipticity at ε₁ = −0.090 and ε₁ = ε₂ = −0.060; assert t_micro ≤ t_macro only. Targets digitized from figures (±0.002 plus the source discretization). Tolerance from a NEMOpt pixel-refinement trend (80–200 px) and an independent body-fitted computation. B-S2c: M-point mode first, ellipticity later | [Michel et al. 2007, Figs. 4b, 6b, 8b; §§3.3, 4.1, 4.3](https://doi.org/10.1016/j.jmps.2006.11.006) (verified from the HAL author version; values digitized from vector figures); [Bertoldi et al. 2010](https://doi.org/10.1002/adma.200901956) (qualitative; its "localization first below 0.34" is not an expectation, since Michel et al. find the 2×2 mode first at porosity π/16); [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) (non-peer-reviewed, context only) |`

**E8. Addendum L182, new row after X-S1 (optional).**
- New: `| X-S2 | Cross-code targets, porous elastomer | Michel et al.'s square array (porosity π/16, κ/μ = 10): on the path ε₂ = −ε₁ the first instability is long-wave, so β and B vanish together at \|ε\| ≈ 0.400 ± 0.004 (a void-surface instability precedes it physically); hexagonal arrays (porosity π/16 / 0.5) lose ellipticity at ε₁ = −0.168 / −0.218, ε₂ = −0.198 / −0.221 and ε = −0.126 / −0.126 (equibiaxial); valid only with their material and, for hexagonal, an L × √3L cell | [Michel et al. 2007, Figs. 4a, 5, 8](https://doi.org/10.1016/j.jmps.2006.11.006) (digitized) |`

**E9. Addendum L306 (NDR-0050, literature column).**
- Old: `Literature: Åberg & Gudmundson 1997; Triantafyllidis et al. 2006; Zhang et al. 2021; Wang & Sigmund 2023`
- New: `Literature: Åberg & Gudmundson 1997; Triantafyllidis et al. 2006; Michel et al. 2007; Zhang et al. 2021; Wang & Sigmund 2023`

**E10. Addendum L310 (NDR-0054).**
- Old: `B-S2 against an independent body-fitted computation (NDR-0042 unaffected)`
- New: `B-S2a against the digitized Michel et al. 2007 targets (porosity π/16, κ/μ = 10) and an independent body-fitted computation; B-S2b asserting t_micro ≤ t_macro against the porosity-0.5 ellipticity targets; B-S2c qualitative (NDR-0042 unaffected; the (J − 1)² volumetric variant is test-only and outside FEAD2026)`
- If E8 is adopted, also change `K-S1–3, E-S1–6, B-S1–4, X-S1 and R-S0` to `K-S1–3, E-S1–6, B-S1–4, X-S1–2 and R-S0`.

**E11. Addendum L175 (E-S4, source column).**
- Old: `Re-computed, six digits`
- New: `Re-computed, six digits; a Bloch-homogenization statement under NEMOpt's uniform-node metric, not part of the Geymonat et al. chain (with gradient normalization only B̂ ≥ β₀ holds; Michel et al. 2007, Eq. 2.30)`

**E12. Research notes.**
- `instability_theory_detection.md` L106.
  - Old: `The abstract was not retrievable **[content not re-verified]**. As I recall it, the paper compares periodic FE Bloch analyses with second-order homogenization estimates of macroscopic loss of ellipticity.`
  - New: `Full text read from the HAL author version (hal-00111501). The paper compares second-order homogenization estimates of macroscopic loss of ellipticity (random, square and hexagonal arrays of circular holes) with finite-element one-cell homogenized moduli; first-bifurcation curves come from 2×2 aggregates after a single coarse Bloch scan of one case (square, porosity π/16). See ../NEMOpt supplied papers verification/michel_et_al_2007.md.`
- `instability_theory_detection.md` L118.
  - Old: `The abstracts and full texts of Miehe et al. (2002), Michel et al. (2007), Bertoldi et al. (2008, JMPS) and Overvelde & Bertoldi (2014) could not be retrieved (Elsevier paywall and 403s).`
  - New: `The abstracts and full texts of Miehe et al. (2002), Bertoldi et al. (2008, JMPS) and Overvelde & Bertoldi (2014) could not be retrieved (Elsevier paywall and 403s); Michel et al. (2007) was later read from HAL.`
- `stability_implementation_verification.md` L225.
  - Old: `and of Michel et al. 2007 and Triantafyllidis et al. 2006 (numerical algorithms and tabulated critical values) could not be retrieved`
  - New: `and of Triantafyllidis et al. 2006 (numerical algorithms and tabulated critical values) could not be retrieved (Michel et al. 2007 was later read from HAL; it reports critical values in figures only)`
- `stability_implementation_verification.md` L353.
  - Old: `(failure surfaces of porous neo-Hookean solids with circular holes in strain space)`
  - New: `(onset-of-failure surfaces for periodic porous and particulate elastomers, per Michel et al. 2007, §1; matrix law unverified)`
- `stability_implementation_verification.md` L368, expectation (d).
  - Old: `(d) for porosity below about 0.34, macroscopic first.`
  - New: `(d) no macroscopic-first expectation: at porosity π/16 with κ/μ = 10, Michel et al. 2007 find the 2×2 mode first under uniaxial-strain and equibiaxial compression.`
- `stability_implementation_verification.md` L377.
  - Old: `Exact published critical strains or stretches with full parameter sets for neo-Hookean hole arrays could not be extracted. Michel et al. 2007, Triantafyllidis et al. 2006,`
  - New: `Michel et al. 2007 give full parameter sets for compressible neo-Hookean hole arrays (κ/μ = 10), with critical strains in figures only (digitized in the supplied-papers note). Triantafyllidis et al. 2006,`

## 6. Gaps and open items

- **Unread primary sources.** LP&PC 2004 (J. Elasticity), the 2D second-order paper behind NDR-0064, was not read; only its bibliographic data were confirmed. Geymonat et al. 1993 remains unread; this note confirms it only through a co-author's restatement.
- **Version risk.** Michel values come from the HAL file, which OpenAlex labels as a submitted version. Verify Figs. 4b and 8 against the journal version before freezing tolerances.
- **Unquantified model-form shift.** The effect of NEMOpt's (ln J)² term, staircase holes and residual void stiffness on the B-S2a onsets is unknown. A body-fitted run with both volumetric laws would quantify it.
- **Bertoldi loading type.** Whether Bertoldi et al.'s "uniaxial compression" is uniaxial stress (lateral traction-free) was not checked. NEMOpt can impose only F_M, so Michel's uniaxial-strain path is the directly reproducible one.
- **Surface instability.** Michel's void-surface criterion (App. A.1) has no pixel-mesh analogue. On staircase holes, a mesh-dependent local mode is a possible artifact to watch for in B-S2.

---
*Coordinator note.* The PDF the user supplied as "Michel et al. 2007" is Lopez-Pamies & Ponte Castañeda (2007), JMPS 55:1677–1701, doi:10.1016/j.jmps.2007.01.007; Michel et al. (2007), JMPS 55(5):900–938, doi:10.1016/j.jmps.2006.11.006, was read from its open HAL deposit (hal-00111501, submitted version). Benchmark values were digitized from vector figures (about ±0.002 in log strain) and must be checked against the journal version before tolerances are fixed.
