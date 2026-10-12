# Triantafyllidis, Nestorović & Schraad (2006): full-text verification and benchmark extraction

**Paper.** N. Triantafyllidis, M. D. Nestorović & M. W. Schraad, "Failure surfaces for finitely strained two-phase periodic solids under general in-plane loading," *J. Appl. Mech.* 73(3):505–515 (May 2006), [doi:10.1115/1.2126695](https://doi.org/10.1115/1.2126695). Called "TNS" below.

**How it was read (2026-10-12).**
- Source: the user-supplied PDF, 11 pages, journal pp. 505–515. PDF page n is journal page 504 + n.
- All 11 pages were read as page images. The text layer (pdftotext) was used only to locate and check short quotes.
- Figures 3–9 are raster JPEGs at 300 ppi. Critical strains were digitized from them by pixel sampling along radial rays (Section 6). No value in the paper is tabulated.
- Manuscript received 9 Feb 2005, final version 1 Sep 2005 (p. 505 footnote). The online date is not printed in the PDF.

**Labels.** *verified* = read in the full text (page, equation or figure given); *digitized* = read from a figure with the stated uncertainty; *inference* = derived here; *unverified* = recalled or not checked against a primary source. Report line numbers refer to `reports/NEMOpt buckling and stability addendum.md` as of this session.

**Bottom line.**
- Every claim the reports and notes attribute to TNS is correct in substance. Three need correction in detail:
  - the "established workflow" in the theory note (TNS use a sign test, not an eigen-solve, and do not retrieve modes);
  - the note's description "porous neo-Hookean solids with circular holes" (they are two-phase solids, and the "voids" are soft inclusions);
  - NDR-0050's attribution (TNS do not use TRIM points, contours or time reversal, and their k = 0 constraint rule differs from ours).
- TNS give peer-reviewed, strain-controlled hole-array data that fits NEMOpt's full-F_M driver. There are two catches:
  - the porosity is π/16 ≈ 0.196 with soft inclusions, not 0.5–0.6 with empty holes;
  - every value must be digitized from figures.
- The strongest case for NEMOpt is equibiaxial compression of the square array with κ/μ = 9.8. There the first instability is at the M point θ = (π, π), a 2×2 antisymmetric mode, at principal log strain −0.089 ± 0.003, and ellipticity is lost later at −0.103 ± 0.003.
- TNS's energy has exactly NEMOpt's shear part, μ/2(I₁ − 2) − μ ln J. Only the volumetric term differs: κ/2(J − 1)² against λ/2(ln J)². Two changes make a true reproduction possible, and both leave FEAD2026 untouched:
  - a test-only volumetric option;
  - a two-phase material map that keeps both phases on the nonlinear branch.
- TNS's Bloch method agrees with the planned module: phases on boundary DOFs, a full-zone grid, and the long-wave limit left to the ellipticity of the homogenized tangent. One rule conflicts and should not be copied: TNS constrain a corner DOF whenever a single ω_i = 0, that is, along the zone axes, not only at k = 0.

---

## 1. Notation traps when reading TNS

| Item | What TNS write | Consequence for NEMOpt |
|---|---|---|
| Cell size L_i | Sec. 3 puts the cell at X_i ∈ [−L_i, L_i] (p. 509, Fig. 3), so L_i is the **half-width**. The square cell is 2L × 2L, the diagonal cell 2L × 4L, and "the volume fraction of the inclusions is (π/4)(R/L)²" (p. 508) | NEMOpt's cell, [−Lx/2, Lx/2]² with Lx = Ly = 2 ([NEMOpt_main.jl L10](../../Full_Domain/NEMOpt_main.jl#L10); [NEMOpt.jl L173](../../Full_Domain/NEMOpt.jl#L173)), is TNS's cell with L = 1, so R = 0.5 |
| Bloch phase | "V_r = exp(2iω₁L₁)V_l" (Eq. 3.9). Scan range "(ω₁L₁, ω₂L₂) ∈ {[0, π) × [0, π)}" (p. 509) | NEMOpt's phase per lattice vector is θ_i = 2ω_iL_i, and TNS's scan is the full zone θ ∈ [0, 2π)² |
| Inconsistent Sec. 2 notation | Eq. 2.6 writes (ω₁L₁, ω₂L₂) ∈ (0, 2π)², and Eq. 2.3 uses L₁ as the period (with the typo "X₁ + n₁L₁" in both slots) | Use the Sec. 3 convention. TNS's "(ω₁L₁)_c = (ω₂L₂)_c = π/2" (p. 511) is **θ = (π, π), the M point**, confirmed by "the repeating unit pattern of the eigenmode is a 2 × 2 cell region" (p. 512, Fig. 5). Reading it as θ = π/2 would be wrong |
| κ | W = μ/2(I₁ − ln I₂ − 2) + κ/2(√I₂ − 1)² (Eq. 3.3), with I₂ = det C. Hence W = μ/2(I₁ − 2) − μ ln J + κ/2(J − 1)²; "κ determines its compressibility" (p. 509) | κ is the small-strain Lamé constant, not a bulk modulus. ν₀ = κ/(2(κ + μ)) gives 0.49495 at κ/μ = 98 and 0.45370 at κ/μ = 9.8. NEMOpt's E = 100 and ν = 0.3 give λ/μ = 1.5 (inference, derived by linearizing at F = I) |
| Stress and moduli | "Π̇_ji = L_ijkl Ḟ_kl" with L_ijkl = ∂²W/∂F_ij∂F_kl (Eqs. 2.1, 3.5) | The same first elasticity tensor as NEMOpt's 𝔸, because Π is the transposed nominal stress. The rank-one test a_i n_j L^H_ijkl a_k n_l (Eq. 2.11) is (m⊗N):𝔸_M:(m⊗N) |
| Eigenvalue metric | ‖v‖² = ∫ v̄_m,n v_m,n dV (Eq. 2.4): an **H¹-seminorm** (gradient) metric, also used in Eq. 2.7 | The numerics test only signs, so the metric is irrelevant to them. It matters only for margins (Section 4.3, C5) |
| "Voids" | "The case of voids is approximated by soft inclusions (μ_f/μ_m = 0.02)" (p. 510), and "the compressibility ratio has the same value in each phase" (p. 511, footnote 4) | The "void" has κ_f = 98 × 0.02 μ_m = 1.96 μ_m (0.196 μ_m at κ/μ = 9.8). It is a soft, nearly incompressible filler, not an empty hole (inference from the stated parameters) |
| Line styles | The text calls the microscopic curves "dotted" (p. 511); the legends draw them dashed | Cosmetic |

---

## 2. Verification of claims (Task 1)

**Search.** "Triantafyllidis" occurs on 5 lines of the addendum (L13, L46, L56, L178, L306), on no line of `reports/NEMOpt refactoring strategy literature.md`, on 14 lines of `instability_theory_detection.md` and on 17 lines of `stability_implementation_verification.md`. No other file under `research_notes/` contained a match when searched, before this note was written.

### 2.1 Claims attributed to TNS 2006

| # | Location | Claim (short) | Verdict | Evidence in TNS |
|---|---|---|---|---|
| A1 | Addendum L13 (panel row "Stability theorist"; TNS is one of three anchors) | (a) homogenized moduli are meaningful only before the first bifurcation; (b) micro (any k) precedes or coincides with macro; (c) one-cell definiteness covers only k = 0; (d) Newton converges just as readily to unstable equilibria | (a) **needs narrowing** with respect to TNS; (b) **confirmed**; (c) **confirmed in substance**; (d) **not in the paper** (inference, fine as such) | (a) TNS tie the existence of L^H to uniqueness of the cell problem: L^H "requires uniqueness in the solution of the unit cell boundary value problem" (p. 508). The any-k statement is Geymonat's or Zhang's, so see the GMT note for the wording. (b) The microscopic surface "always lies inside its macroscopic counterpart" (p. 506), and λ_h = inf over ω → 0 of λ_m(ω) (Eq. 2.13) with λ_c = inf over ω of λ_m(ω) (Eq. 2.8). (c) ω = 0 is one point of the scan (Eqs. 2.6–2.8) |
| A2 | Addendum L306 (NDR-0050, "Literature: … Triantafyllidis et al. 2006 …") | Source for the Bloch module: phases on the master–slave map, no pin for θ ≠ 0, K(θ) as a sum, time-reversal halving, TRIM points, contour, half-zone grid, ball near θ = 0, condensation | **Partly confirmed; needs correction in detail** | TNS support the phase relations (Eq. 3.9), a uniform full-zone grid (p. 510), boundary condensation (Eqs. 3.7–3.12), an LDU sign test (p. 510) and the long-wave limit through L^H (Eq. 2.13). They use **no** TRIM, contour or time-reversal reduction ("not necessarily unique due to symmetry", p. 510). Their k = 0 rule, "V_Ai = 0" whenever "ω_i = 0" (p. 510), differs from ours, and they are **not** the source for "no pin at θ ≠ 0" (Zhang et al. are). See Section 4 |
| N1 | Theory note L34 (Zhang's sentence "micro … before or simultaneously with the macro", citing Geymonat and TNS) | Order of onset | **Confirmed** | "always lies inside its macroscopic counterpart" (p. 506), and Fig. 4(b) shows the separation |
| N2 | Theory note L40 | Onset of failure is "the first occurrence of a bifurcation in the fundamental solution, for which all cells deform identically", computed "based on a Bloch wave analysis on the deformed unit cell"; "(online 2005)"; 73(3):505–515 | **Confirmed verbatim** (abstract, p. 505). The online date is **not printed** in the PDF | Volume 73, May 2006, pp. 505–515, DOI as cited. The issue number is not printed on the pages |
| N3 | Theory note L71 ("The established workflow (Triantafyllidis et al. 2006; …) has four steps") | 1 Newton on one cell; 2 Bloch eigenproblem, scanning k; 3 rank-one convexity of the homogenized tangent; 4 mode retrieval | **Needs correction** for TNS | Step 1: Newton (Δλ = 10⁻³, p. 509), but on "one-quarter of the unit cell" for θ = 0 and the full cell otherwise (p. 509). Step 2: TNS evaluate "the signature of the matrix K̂ … by performing an LDU decomposition" (p. 510), a sign test rather than an eigen-solve. Step 3: confirmed (Eqs. 2.10–2.12). Step 4: **not in TNS**. The Fig. 5 mode is "courtesy of Dr. J. C. Michel" (p. 511) |
| N4 | Theory note L82–87 | Built surfaces "in strain or stress space by following all loading paths"; "a numerical algorithm which reduces the size of stability matrices by an order of magnitude"; porous and particle-reinforced solids; wavelength dependence | **Confirmed** (abstract quotes verbatim, p. 505), **with two nuances** | (i) Stress space appears only as a possibility in the abstract; all computations are radial paths in principal **log-strain** space (Eq. 3.2; Figs. 4–9). (ii) "Porous" means soft inclusions, μ_f/μ_m = 0.02 (p. 510) |
| N5 | Theory note L90 (Zhang: two condensation methods, one following TNS) | TNS has a condensation | **Confirmed** | Eqs. 3.7–3.12 |
| N6 | Theory note L114 | Exact tangent at the converged nonlinear state, not a linearized prebuckling | **Confirmed** | L = ∂²W/∂F∂F at the principal solution (Eq. 3.5); "standard incremental Newton-Raphson technique" (p. 509); K(λ) is the "full (real and uncondensed) stiffness matrix" (p. 509) |
| N7 | Impl. note L157 (complete in the GMT sense; the pinned center is exactly the k = 0 Bloch matrix) | Completeness; k = 0 pin | **Consistent** | TNS restate GMT's completeness (p. 505). At ω = 0 TNS constrain corner A instead of the center (p. 510); one fixed node suffices for inertia (inference) |
| N8 | Impl. note L167 | The same quotes plus "a numerical algorithm which reduces the size of stability matrices by an order of magnitude"; "abstract via Crossref" | **Confirmed verbatim** against the full text | p. 505. The size reduction is from about 2N²×2N² to about 4N×4N for the square cell (p. 509) |
| N9 | Impl. note L183 (Zhang's description: "the inversion of K_NN … only once at each loading step"; "more efficient … when the sampling grid … is large") | Cost of the condensation | **Confirmed in substance**. The second phrase is Zhang's assessment | "the inversion of K_ii(λ) is done only once (based on an LDU decomposition that takes advantage of the banded nature" (p. 510). TNS say only that the signature calculations are "much quicker" |
| N10 | Impl. note L225–226 (gaps: TNS's algorithms and "tabulated critical values" not retrieved; condensation "not read") | Gap statement | **Superseded** | The full text has now been read. It contains **no tables**, only Figs. 4–9. The condensation is a Schur complement onto V_b = (V_l, V_d, V_A) (Eqs. 3.7, 3.12), the same object as the note's Haynsworth construction, but TNS do not state that K_ii must be positive definite for the sign test to hold |
| N11 | Impl. note L353 ("failure surfaces of porous neo-Hookean solids with circular holes in strain space") | Description | **Needs correction** | They are two-phase compressible neo-Hookean solids with κ/μ = 98 and 9.8 in both phases (Eq. 3.3; p. 511, footnote 4). Circular **inclusions** at f = π/16, with μ_f/μ_m from 0.02 ("voids") to 50 ("rigid") (p. 510; Fig. 7), in square and diagonal arrangements (Fig. 2), at θ = 0 and π/18 (Fig. 9). "Strain space" is correct |
| N12 | Impl. note L377 (TNS "not accessible … or report values only in figures") | Gap statement | **Partly superseded**: accessible now; values appear **only in figures** (confirmed) | Section 3.5 gives digitized values |

### 2.2 Claims that belong to other papers (TNS corroboration where relevant)

- **Triantafyllidis & Maker 1985.**
  - Addendum L56 and theory note L24: "consistently and sometimes considerably overestimates the stable region".
  - Impl. note L168: "the macroscopic stable region consistently envelops the microscopic stable region".
  - Both are correctly attributed to T&M 1985, not to TNS. TNS add that T&M established the link "for the special case of an incompressible, hyperelastic composite under plane strain conditions" (p. 505), which supports B-S1's use of the incompressible limit (addendum L178). TNS's general counterpart is "always lies inside" (p. 506).
- **Geymonat, Müller & Triantafyllidis 1993.** Bloch completeness, and the equivalence of a long-wave mode with loss of ellipticity (addendum L56; theory note L25; impl. note L157, L166).
  - TNS p. 505 restates both, written by a co-author: GMT "proved the completeness of the Bloch wave representation", and a mode of infinite wavelength "can be detected as a loss of ellipticity of the corresponding homogenized tangent moduli".
  - Inequality chains and theorem numbers are in the GMT note ([geymonat_muller_triantafyllidis_1993.md](geymonat_muller_triantafyllidis_1993.md)).
- **Triantafyllidis & Bardenhagen 1996 (TNS ref. [3]).** The onset-of-failure-surface concept itself: "[3] defined the onset-of-failure surfaces in stress and strain space for periodic solids" (p. 505).
  - Theory note L46 cites T&B 1996 for scale effects only. Its bibliographic data (JMPS 44:1891–1928) match TNS ref. [3].
- **TNS ref. [5] (aluminum honeycombs).**
  - The observation that λ_m(ω) often has "a singular point in the neighborhood of ω = 0" is credited to [5] (p. 507).
  - The PDF lists [5] as "Triantafyllidis & Bardenhagen 1996, JMPS 46, 1089–1124". Volume 46 is a 1998 volume, and as recalled the honeycomb paper is Triantafyllidis & Schraad 1998 (unverified; probably a reference-list error in TNS).
- **Nestorović & Triantafyllidis 2004 (TNS ref. [6]).** TNS validated their code against the "analytical calculations of [6]" on "a rectangular unit cell … made of two different layers" (pp. 510–511). This is the precedent for B-S1 and the source to read on compressible laminate micro and macro onsets (Section 5, E6).
- **Michel 2005, "work in preparation" (TNS ref. [11]).**
  - Fig. 5 (the eigenmode) is Michel's computation.
  - Michel repeated "some of the loss of ellipticity calculations … (based on a much finer unit cell f.e.m. mesh)" (p. 511).
  - Michel also supplied the "preliminary calculations" on imperfections and post-bifurcation (p. 514).
  - That this work became Michel, Lopez-Pamies, Ponte Castañeda & Triantafyllidis 2007 is plausible but unverified. If so, Michel et al. 2007 should hold finer-mesh versions of TNS's porous results.
- **Abeyaratne & Triantafyllidis.**
  - TNS cite the 1984 J. Appl. Mech. 51:481–486 paper (ref. [13]) for loss of ellipticity of a porous periodic Mooney–Rivlin-type solid (p. 511). This matches theory note L41.
  - Impl. note L239 attributes the same statement to "Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134".
  - TNS cite only the 1984 paper for porous solids, so L239's citation should be checked. As recalled, the 1981 IJSS paper concerns shear bands in homogeneous plane strain (unverified).
- **Zhang, Feng & Khandelwal 2021.** Zhang's own wording, not TNS's:
  - "The micro-instability … occurs either before or simultaneously";
  - "the homogenization results are only valid before the bifurcation occurs";
  - "no need for adding extra constraints" at k ≠ 0.
- **Grep hits unrelated to TNS:**
  - theory note L45 (Schnaidt 1993), L48 (Elliott 2006), L207 (Combescure 2020);
  - impl. note L287, L327 and L337 (Slesarenko & Rudykh; N&T 2004) and L352 (Michel 2007).
  - Theory note L307 quotes T&M 1985 ("loss of ellipticity in the homogenized incremental equilibrium equations"). TNS's p. 505 wording is consistent with it.
  - Nothing in TNS contradicts any of these lines.

### 2.3 Report statements that TNS support but do not cite

| Location | Statement | TNS support |
|---|---|---|
| Addendum L58 | k = 0 cell-periodic modes differ from long-wave modes as k → 0; the latter belong to the 𝔸_M test | "two physically different types of modes: the strictly periodic modes ω = 0, v = p, and the long-wavelength modes ω → 0" (p. 507), plus Eq. 2.13. TNS predates Zhang 2021 |
| Addendum L58 and L78 (validity of 𝔸_M) | 𝔸_M is meaningful only while the k = 0 problem is unique | "assuming of course that L^H exists, which requires uniqueness in the solution of the unit cell boundary value problem, thus excluding a first bifurcation with ω_c = 0" (p. 508) |
| Addendum L131 | Grids; condensation | Full-zone grid with "grid spacing is π/36 in each direction" in ω_iL_i, swept at every λ (p. 510); Eq. 3.12 |
| Addendum L139 | "The literature offers no shortcut here" for half-cell stability | TNS solve the principal problem on "one-quarter of the unit cell" for θ = 0, with symmetry conditions, and run the Bloch analysis on the full cell (pp. 509–510). This is the same split as NDR-0051, without a discussion of what the quarter cell misses |
| Addendum L178 (B-S1) | Laminates as oracle | TNS's code check on two-layer cells (pp. 510–511) |
| NDR-0045 | Validity rule | Onset of failure as the first bifurcation of the principal branch on which "all cells deform identically" (p. 505) |

---

## 3. Benchmark data (Task 2)

### 3.1 Microstructures (verified unless marked)

| Item | Square (Fig. 2(a)) | Diagonal (Fig. 2(b)) |
|---|---|---|
| Unit cell | 2L × 2L, with L₁ = L₂ = L (half-widths) | 2L × 4L, with L₁ = L, L₂ = 2L |
| Inclusions per cell | 1, centered | 2 (center and corners; a centered rectangular lattice) |
| Radius | R = 0.5L (p. 510) | R = 0.5L |
| Volume fraction | π/16 = 0.19635, "approximately 20%" (p. 510) | π/16 (stated equal, p. 508) |
| Mesh | "N × N rectangular (in parameter coordinate space) grid" of bilinear quadrilaterals, two DOFs per node (p. 509). N is **not reported**. The Fig. 3 schematic shows about 20 elements per side (digitized count), body-fitted with a mapped circle | "N × 2N" grid (p. 509) |
| Integration rule | **Not reported**, so volumetric locking at κ/μ = 98 cannot be excluded | — |
| Rigid-inclusion contact (dot-dash lines) | ε₁ = ln 0.5 = −0.6931 or ε₂ = ln 0.5 (inference, derived; matches Figs. 6 and 9(b) to within 0.002) | ε₁ = ln 0.5; ε₂ = ln 0.25 = −1.3863; and λ₁² + 4λ₂² = 1 between (−0.693, −0.837) and (−0.144, −1.386) (inference, derived; matches Fig. 8(b)) |

### 3.2 Material models and constants (verified)

- Each phase follows W = μ/2(I₁ − ln I₂ − 2) + κ/2(√I₂ − 1)² (Eq. 3.3), which is "polyconvex in the sense of Ball", so each phase is rank-one convex at any deformation (p. 509).
- **Compressibility.** κ/μ = 98 ("nearly incompressible") and 9.8, the same in both phases (footnote 4).
- **Stiffness contrast μ_f/μ_m.**
  - 0.02 ("voids"; Figs. 4, 8(a), 9(a));
  - 0.1 and 0.5 (Fig. 7(a));
  - 2 and 10 (Fig. 7(b));
  - 50 ("rigid"; Figs. 6, 8(b), 9(b)).
  - For μ_f/μ_m = 1 "the solid is uniform and strongly elliptic at all strains" (p. 513).
- **Comparison with NEMOpt (inference, derived).**
  - Shear parts are identical.
  - The volumetric stress terms are κJ(J − 1)F⁻ᵀ (TNS) and λ ln J F⁻ᵀ (NEMOpt, [NEMOpt.jl L281](../../Full_Domain/NEMOpt.jl#L281)). Both have the same small-strain limit when λ = κ.
  - The ratio J(J − 1)/ln J is 0.985 at J = 0.99, 0.926 at J = 0.95 and 0.854 at J = 0.90. Matching only ν₀ therefore reproduces TNS's law to about 1.5% in the volumetric term inside the nearly incompressible matrix, but poorly in the soft inclusion, whose J departs much further from 1.
  - TNS's tangent is NEMOpt's L284 with λ replaced by κJ(2J − 1) in the F⁻ᵀ⊗F⁻ᵀ term and λ ln J replaced by κJ(J − 1) in the F⁻ᵀ_il F⁻ᵀ_kj term.

### 3.3 Loading paths (verified)

- Symmetric macroscopic stretch, with no rotation: ⟨F⟩ = ⟨U⟩ = R(θ) diag(λ₁, λ₂) R(θ)ᵀ (Eq. 3.1).
- Radial paths in principal log strain: ln λ₁ = ε₁ = λ cos φ and ln λ₂ = ε₂ = λ sin φ (Eq. 3.2), so the load parameter λ equals |ε|.
- Principal-axes angle θ = 0 ("orthotropic") or π/18. Paths φ ∈ [0, 2π) are marched "starting from λ = 0" (p. 509). The φ resolution is not reported.
- Key rays:
  - φ = π: uniaxial strain compression along X₁, ε₂ = 0;
  - φ = 5π/4: balanced (equibiaxial) compression;
  - φ = 3π/4: isochoric pure shear;
  - φ ∈ [0, π/2]: biaxial tension.

### 3.4 Tabulated values

None. Every result is a curve in Figs. 4 and 6–9. Fig. 5 is a mode shape.

### 3.5 Critical strains (digitized; λ = |ε|; all θ = 0 unless noted)

Uncertainties combine pixel scale, line width, JPEG blur and overlaps:
- ±0.002 in λ for Fig. 4(b) (2.5e−4 per pixel);
- ±0.005 for Fig. 6(b) (1.0e−3 per pixel), ±0.01 at the corner near φ = 5π/4;
- ±0.005–0.01 for Fig. 7 (3.7e−3 per pixel);
- ±0.03–0.04 on pure-shear rays, which graze the curves.

TNS's own load step adds a one-sided resolution of Δλ = 10⁻³ (p. 509). At φ = π and 3π/2 the values come from rays at 1–2° off the axis and are extrapolated to the axis, because the frame and axis lines hide the curve ends.

**Table A. Square array, soft inclusions (μ_f/μ_m = 0.02), Fig. 4(b).**

| φ | κ/μ = 98, Bloch onset | κ/μ = 98, loss of ellipticity | κ/μ = 9.8, Bloch onset | κ/μ = 9.8, loss of ellipticity |
|---|---|---|---|---|
| 180° (uniaxial strain) | 0.137 ± 0.003 (coincident) | 0.137 ± 0.003 | 0.178 ± 0.003 | 0.194 ± 0.003 |
| 190° | — (hidden by legend) | — | 0.154 | 0.170 |
| 200° | 0.106 | 0.108 | 0.138 | 0.156 |
| 205° | 0.102 | 0.104 | 0.134 | 0.152 |
| 215° | 0.098 | 0.100 | 0.128 | 0.146 |
| 225° (equibiaxial) | **0.096** (ε₁ = ε₂ = −0.068) | **0.098** (−0.069) | **0.126** (−0.089) | **0.146** (−0.103) |
| 240° | 0.100 | 0.102 | 0.130 | 0.148 |
| 250° | 0.106 | 0.108 | 0.138 | 0.156 |
| 260° | 0.118 (coincident) | 0.118 | 0.154 | 0.170 |
| 270° | 0.137 (by symmetry) | 0.137 | 0.178 (0.175 on the 269° ray) | 0.194 (0.191 on the 269° ray) |

Notes on Table A:
- The curves are symmetric about ε₁ = ε₂, as TNS state (p. 511), and the digitized pairs at 200°/250° and 215°/240° agree.
- At κ/μ = 98 the two curves separate only for φ ≈ 196°–250° (the dashed segments run from about (−0.105, −0.03) to (−0.04, −0.105)), by about 0.002 in λ, two of TNS's load steps. This separation is at the limit of digitization.
- In compression the κ/μ = 98 surface lies close to ε₁ + ε₂ ≈ −0.137 to −0.139, a nearly constant macroscopic ln J ≈ −0.137 (inference, observation).
- In stretches, κ/μ = 9.8 at 225° gives a Bloch onset at λ₁ = λ₂ = 0.9148 and loss of ellipticity at 0.9019. At 180°, it gives 0.8369 and 0.8237. For κ/μ = 98: 0.9343 and 0.9330 at 225°, and 0.8720 at 180°.
- Fig. 4(a), κ/μ = 98, full plane (±0.03): the pure-shear rays give λ ≈ 0.78 (ε = ±0.55) at both 135° and 315°. There is "no instability under biaxial tension (0 ≤ φ ≤ π/2)" for |ε_i| ≤ 1.5 (p. 511).
- Check on TNS's text: "critical strains increase by approximately 40%" from κ/μ = 98 to 9.8 (p. 511). Digitized, the increase is 42% (loss of ellipticity at 180°), 49% (loss of ellipticity at 225°), 30% (Bloch onset at 180°) and 31% (Bloch onset at 225°), which is consistent.

**Table B. Square array, stiff inclusions (μ_f/μ_m = 50), Fig. 6(b).** The Bloch onset and loss of ellipticity coincide everywhere: "always coincident" (p. 512).

| φ | κ/μ = 98 | κ/μ = 9.8 | Rigid contact (exact) |
|---|---|---|---|
| 180° | 0.454 ± 0.005 (stretch 0.635) | 0.545 ± 0.005 (0.580) | 0.693 |
| 200° | 0.542 | 0.626 | — |
| 225° | 0.895 ± 0.01 (ε = −0.633; stretch 0.531) | 0.925 ± 0.01 (ε = −0.654; 0.520) | 0.980 |
| 250° | 0.542 | 0.626 | — |
| 270° | 0.454 (by symmetry) | 0.545 | 0.693 |

Fig. 6(a), κ/μ = 98: pure shear λ ≈ 0.49–0.50 (ε = ±0.35 ± 0.01). In uniaxial tension the curve meets the axes near ε ≈ 1.47 ± 0.03, at the edge of the plot.

**Table C. Square array, κ/μ = 98, stiffness contrast, Fig. 7.** Loss of ellipticity only. For contrasts other than 0.02, TNS state that micro and macro coincide "in all cases except for the voided solid … near balanced biaxial compression" (p. 513).

| μ_f/μ_m | φ = 180° | φ = 225° | φ = 135° |
|---|---|---|---|
| 0.02 | 0.137 (Table A) | 0.098 | 0.78 ± 0.04 |
| 0.1 | 0.20 ± 0.005 | 0.145 ± 0.005 | 0.93 ± 0.04 |
| 0.5 | 0.79 ± 0.01 | 0.561 ± 0.005 | 1.55 ± 0.03 |
| 2 | 1.08 ± 0.01 | none in the plotted range (λ > 2.12) | 1.42 ± 0.02 |
| 10 | 0.56 ± 0.008 | 0.99 ± 0.01 (beyond the rigid-contact strain, which does not apply at this contrast) | 0.635 ± 0.01 |
| 50 | 0.455 ± 0.008 | 0.895 ± 0.01 | 0.49 ± 0.01 |

This reproduces TNS's statement that "the first fivefold change … produces considerably smaller changes … than the next fivefold change" (p. 513). At φ = 180°, the voided side goes 0.137 → 0.20 → 0.79 and the rigid side 0.455 → 0.56 → 1.08.

**Table D. Other cases (lower precision).**

| Case | Values (λ) | Source |
|---|---|---|
| Square, soft, θ = π/18 | Unchanged from θ = 0 in biaxial compression (180°–270°, within ±0.003); larger critical strains where one strain is tensile | Fig. 9(a) |
| Square, rigid, θ = π/18 | 0.49 at 180° (against 0.455 at θ = 0), 0.569 at 250° (against 0.544), about 0.895 at 225° (unchanged) | Fig. 9(b) |
| Diagonal, soft | Crosses ε₂ = 0 near ε₁ = −0.165 and ε₁ = 0 near ε₂ = −0.15 (±0.01); about 0.111 ± 0.01 at 225° (ε = −0.079); separate micro curve visible from about (−0.05, −0.12) to (0.5, −0.6), i.e., also in the tension–compression quadrant | Fig. 8(a) |
| Diagonal, rigid | 0.455 at 180° (as for the square); no onset before contact on the 225° and 235° rays; at 270° the onset lies essentially at contact (ε₂ ≈ −1.38) | Fig. 8(b) |

### 3.6 Mode types and which instability comes first (verified unless marked)

- **The k = 0 cell-periodic mode is never first.** "In all the numerical calculations reported here … ω_c ≠ 0" (p. 508).
- **Square, soft (0.02), κ/μ = 98.**
  - For "most loading path angles" the long-wave mode (ω_c → 0) comes first and the micro and macro onsets coincide.
  - The exception is "a region in the biaxial compression zone (about φ = 5π/4)" where the mode is at (ω₁L₁)_c = (ω₂L₂)_c = π/2 (p. 511), i.e., the M point, a 2×2 antisymmetric pattern.
  - Digitized, that region spans φ ≈ 196°–250°.
- **Square, soft, κ/μ = 9.8.**
  - The Bloch onset precedes loss of ellipticity over the whole third quadrant shown in Fig. 4(b), by 8% (φ = π) to 14% (φ = 5π/4) in λ (digitized).
  - "The critical wavelength … remains the same as in the nearly incompressible case" (p. 511). The φ = 5π/4 mode is plotted in Fig. 5: "antisymmetric in the unit cell", with a 2×2 repeating unit (pp. 511–512).
  - The mode at φ = π is not shown separately (no explicit statement).
- **Square, stiff (50), both κ/μ, θ = 0 and π/18.** The long-wave mode is always first: "always coincident" (p. 512; p. 513).
- **Diagonal, soft.** There is "a much broader range of the loading path angle φ for which a microbuckling instability precedes" loss of ellipticity (p. 513; the text says "microscopic loss of ellipticity", a typo for macroscopic).
- **Diagonal, rigid.** Some biaxial-compression paths show "no microscopic instability (and remains macroscopically elliptic) until inclusion contact" (p. 513).
- **θ = π/18, soft.** In biaxial compression the mode is "local in nature; it is antisymmetric on the unit cell with (ω₁L₁)_c = (ω₂L₂)_c = π/2" (p. 513).
- **Interface instabilities** (Biot type): "For the calculations reported here no such instabilities have been detected" (p. 511).
- **Post-bifurcation (preliminary, from Michel).** For porous solids "the bifurcated solutions are unstable", while "the opposite occurs for the case of inclusions" (p. 514). This is weak evidence, but it favors constraining instability in porous designs (NDR-0059/0063 context).

### 3.7 Errata and inconsistencies found in TNS

1. **Fig. 9 shows the square arrangement.** The text says Fig. 9 "compares the onset-of-failure surfaces … with diagonal microgeometry" (p. 513), but the caption says "square distribution". Three observations support the caption:
   - the insets show aligned 3×3 arrays;
   - the rigid-contact line in Fig. 9(b) sits at ε ≈ −0.69 in both directions, where the diagonal cell would give −1.386 in ε₂;
   - the θ = 0 curves reproduce Figs. 4 and 6 (square), for example 0.096/0.098 at 225° in Fig. 9(a) against 0.111 for the diagonal cell in Fig. 8(a).
2. **"About five times" against "an order of magnitude".** The rigid-inclusion critical strains are "about five times" the porous ones (p. 512) and "an order of magnitude larger" under balanced biaxial compression (p. 514). Digitized, the ratio is 0.895/0.096 ≈ 9.3 at φ = 5π/4 and 0.455/0.137 ≈ 3.3 at φ = π, so p. 514 is right for equibiaxial loading and p. 512 is a loose average.
3. **Sec. 2 notation** (Eqs. 2.3 and 2.6) does not match the Sec. 3 half-width convention (Section 1).
4. **Typos.**
   - "microscopic loss of ellipticity" (p. 513, should read macroscopic);
   - "macroscopic and macroscopic" (p. 513);
   - "orthographically loaded" (p. 511);
   - Eq. 2.3 repeats "X₁ + n₁L₁".
   - The size of the full matrix for the diagonal cell is written "2N² × 4N²" (p. 509), presumably 4N² × 4N².
5. **Reference list.**
   - [5] carries authors and year that do not match volume 46 (Section 2.2; unverified);
   - [14] gives the journal as "J. Exp. Zool." with "Lopez-Pamiez", presumably J. Elasticity 76 (unverified).

### 3.8 Can a TNS case serve as a NEMOpt benchmark?

**Not as NEMOpt stands.** Five obstacles, all verified in code unless marked:

1. **Material.** NEMOpt's volumetric term is λ/2(ln J)² with λ/μ = 1.5 ([NEMOpt.jl L250, L281–284](../../Full_Domain/NEMOpt.jl#L250)); TNS use κ/2(J − 1)² with κ/μ = 98 or 9.8.
2. **Phase law.** A soft phase at g = 0.02 falls on NEMOpt's linear branch.
   - The blend switches at ρ^p = ω₂ = 0.05 with β₂ = 200 ([NEMOpt_main.jl L20](../../Full_Domain/NEMOpt_main.jl#L20); [NEMOpt.jl L56, L297–298](../../Full_Domain/NEMOpt.jl#L56)), and H_Steep(0.02) = 6.1e−6 (computed). TNS's soft phase is neo-Hookean.
   - A stiffness ratio of 50 needs the matrix at g = 0.02, which hits the same switch.
   - Lowering ω₂ only approximates the fix: H_Steep(0.02) = 0.982 at ω₂ = 0.01 and 0.997 at 0.005 (computed).
3. **True voids.** NEMOpt's void is the linear law εℂ with ε = 5e−4, so it cannot represent TNS's nearly incompressible soft filler, and TNS do not model an empty hole.
4. **Elements.** NEMOpt integrates Q4 elements with 2×2 Gauss points ([NEMOpt_main.jl L22](../../Full_Domain/NEMOpt_main.jl#L22)). At κ/μ = 98 (ν₀ = 0.495) volumetric locking is expected (inference, standard). TNS do not state their integration rule.
5. **Geometry.** NEMOpt's structured grid gives a staircase circle (R = 20 elements at 80×80); TNS's mesh is body-fitted.

**What does fit.**
- The cell frame is identical: NEMOpt's 2×2 cell is TNS's cell with L = 1.
- TNS prescribe the full macroscopic F, as NEMOpt does, so no mixed control is needed.
- On rays where F_M is diagonal and the log strains keep a fixed ratio, NEMOpt's linear ramp F(t) = I + t(F_M − I) traces exactly TNS's path: φ = π uses F_M = diag(s, 1), and φ = 5π/4 uses F_M = sI.
  - Then ε_i(t) = ln(1 + t(s_i − 1)), and t_c maps to TNS's λ_c (inference).
  - On other rays the paths differ, but the solids are elastic and "the onset-of-failure surface is independent of the load path" for them (p. 506), so the crossing point in (ε₁, ε₂) is what should be compared.
- Rectangular cells for the diagonal case are parameterized: Lx, Ly, Nelx and Nely are separate inputs ([NEMOpt_main.jl L10–11](../../Full_Domain/NEMOpt_main.jl#L10)). Whether every routine supports Lx ≠ Ly is untested.

**Minimum changes for a reproduction.** All are test-only and none touches FEAD2026:
- **(a) Volumetric option.** `volumetric = :lnJ2` (default) or `:Jm1sq`. With U(J) = κ/2(J − 1)², P = μ(F − F⁻ᵀ) + κJ(J − 1)F⁻ᵀ, and 𝔸 = μ δ_ik δ_jl + κJ(2J − 1) F⁻ᵀ_ij F⁻ᵀ_kl + [μ − κJ(J − 1)] F⁻ᵀ_il F⁻ᵀ_kj (inference, derived). It is to be checked by the NDR-0020 AD and central-difference tests. Z is not needed, because the benchmark computes no gradients.
- **(b) Two-phase material map.** Per-element (μ_e, κ_e), with H_Steep ≡ 1, filter and projection bypassed, and a centroid-in-circle phase rule. The discrete volume fraction is reported.
- **(c) Compressive Newton.** NDR-0031 and NDR-0032 are required, because Rep_MU = 0 for F_M = sI and for diag(s, 1) ([NEMOpt.jl L375](../../Full_Domain/NEMOpt.jl#L375)).
- **(d) Bloch grid.** The TRIM M point, plus a full-zone grid of at least 36×36 at the onset load, to match TNS's resolution.

**Recommended case set (B-S2a), in priority order.**

| Priority | Case | Load (NEMOpt F_M; ramp to beyond onset) | Assert | Why |
|---|---|---|---|---|
| 1 | Square, soft 0.02, κ/μ = 9.8 | sI with s ≈ 0.88 | First failure is Tier 2 at the M point, before Tier 0b; Bloch onset at ε = −0.089 and loss of ellipticity at −0.103 (Table A); Tier 0a passes throughout | Largest micro–macro gap (14%); mode documented (Fig. 5); mild locking risk (ν₀ = 0.454) |
| 2 | Same material | diag(s, 1) with s ≈ 0.80 | Bloch onset −0.178 before loss of ellipticity −0.194; record the critical θ (TNS do not state it here) | Second ray; tests ordering at an 8% gap |
| 3 | Square, stiff 50, κ/μ = 9.8 | sI and diag(s, 1) with s ≈ 0.5 | No finite-k onset before loss of ellipticity; loss of ellipticity at −0.654 and −0.545; contact at ln 0.5 is a geometry sanity bound | Exercises Tier 0b, E-S4 and "micro = macro" |
| 4 | Square, soft, κ/μ = 98 | sI | Assert only t_micro ≤ t_macro and an onset near −0.068 to −0.069; the gap is below digitization uncertainty | Locking-prone; qualitative only |
| 5 | Square, κ/μ = 98, contrast 0.1–10 | diag(s, 1) | Loss-of-ellipticity values from Table C | Optional ellipticity regression |

**Acceptance proposal.** Tolerances are pre-registered before any run.
- (i) Hard qualitative criteria:
  - the critical wave vector is the M point for case 1;
  - the order is Bloch before ellipticity, with a gap above 5%, for cases 1 and 2;
  - the Bloch onset coincides with loss of ellipticity, within the load-step tolerance, for case 3;
  - there is no k = 0 failure before any onset.
- (ii) A body-fitted reference with TNS's energy and a locking-free element (for example Q2, or a mixed u–p formulation), run in Ferrite, Gridap or a commercial code as NDR-0042 allows, must match the digitized values to within 7%. That margin is a ±3% digitization band plus an allowance for TNS's unknown mesh.
- (iii) NEMOpt's staircase cell must match the body-fitted reference to within a tolerance fixed from an 80×80 and 160×160 refinement study.

**With NEMOpt's own material** (λ/μ = 1.5, true ε-voids), no TNS case is a reproduction. TNS's trend suggests larger critical strains and a larger micro–macro gap for the more compressible solid: going from κ/μ = 98 to 9.8 raised both and widened the gap (p. 511; Table A). Extrapolating this trend to λ/μ = 1.5 is inference, not a benchmark.

---

## 4. Bloch method in TNS and comparison with the planned module (Task 3)

### 4.1 TNS's algorithm, step by step (verified)

1. **Principal solution.**
   - Newton–Raphson with "step size Δλ = 10⁻³", usually "three iterations per step", and "an accuracy of 10⁻⁶ (based on a Euclidean displacement norm)" (p. 509).
   - For θ = 0 the computation uses "one-quarter of the unit cell" with symmetry displacement conditions. For θ ≠ 0 it uses "the entire unit cell … by imposing appropriate coupling conditions" on opposite faces (p. 509).
2. **Bloch reduction.**
   - Matching nodes on opposite faces (p. 509).
   - The DOFs are partitioned as V = (V_i, V_b, V_t) with V_b = (V_l, V_d, V_A) and V_t = (V_r, V_u, V_B, V_C, V_D) (Eq. 3.7).
   - The phase relations V_r = e^{2iω₁L₁}V_l, V_u = e^{2iω₂L₂}V_d and corners B, C, D from A give V_t = A(ω)V_b (Eqs. 3.9–3.10).
3. **Condensation.**
   - The Hermitian K̂(λ, ω) = K_bb + ĀᵀK_tb + K_bt A + ĀᵀK_tt A − [ĀᵀK_ti + K_bi] K_ii⁻¹ [K_ib + K_it A] (Eq. 3.12).
   - It is about 4N × 4N (square) or 6N × 6N (diagonal), against about 2N² × 2N² uncondensed (p. 509).
   - K_ii is factored once per λ by a banded LDU (p. 510).
4. **Test.**
   - At each λ, "the signature of the matrix K̂ … is evaluated by performing an LDU decomposition", and "Positive definiteness is lost when the lowest entry of the Cholesky diagonal matrix D vanishes" (p. 510).
   - This is a sign test only: no eigenvalues, no margins, no mode shapes.
5. **Wave-vector sampling.**
   - The test is applied "for each λ on a fine grid (grid spacing is π/36 in each direction)" over ω_iL_i ∈ [0, π) (pp. 509–510). That is a uniform full-zone 36 × 36 grid of 1,296 points with phase spacing π/18 (10°). It includes all four TRIM points, and no symmetry reduction is used.
6. **k = 0 and rigid modes.** "To investigate the special case when ω_i = 0 (i.e., when the eigenmode has period 2L_i) the additional constraint V_Ai = 0 is also imposed" (p. 510). At ω = (0, 0) this pins corner A, which equals NEMOpt's center pin for inertia (inference).
7. **Long-wave limit.**
   - Not resolved by the grid: the nearest nonzero grid point is at θ = π/18.
   - It is covered by the macroscopic test, B(λ) = min over unit a and n of a_i n_j L^H_ijkl a_k n_l (Eq. 2.11), with λ_h = inf over ω → 0 of λ_m(ω) (Eq. 2.13). L^H follows from the periodic cell problem (Eq. 2.10).
   - The numerical resolution of the min over a and n is **not reported**.
8. **Failure surface.**
   - λ_c is the first λ at which K̂ loses definiteness at some grid point (p. 510), or at which B(λ) vanishes.
   - "The locus of all the macroscopic stresses or strains evaluated at λ_c for all loading paths constitutes the microscopic onset-of-failure surface" (p. 507); the macroscopic surface is built the same way from B = 0 (Eq. 2.12).
   - That the plotted micro curve is min(grid onset, λ_h) where the long-wave mode governs is implied by Eqs. 2.8 and 2.13, not stated (inference).
9. **Code checks.**
   - Bloch and ellipticity predictions on layered cells were checked against N&T 2004.
   - Michel repeated some loss-of-ellipticity calculations on "a much finer unit cell f.e.m. mesh" (pp. 510–511).

### 4.2 Comparison with the planned NEMOpt module (addendum §§ "Bloch phases ride on the existing master–slave map", "The ellipticity test is exact", NDR-0048 to NDR-0050)

| Aspect | TNS 2006 | NEMOpt plan | Status |
|---|---|---|---|
| Phase map | Faces and corners slaved with e^{2iω_iL_i} (Eq. 3.9) | D1 slave rows × e^{i(n₁θ₁ + n₂θ₂)} from the D2 offsets | **Same**, with θ_i = 2ω_iL_i |
| Operator | Exact tangent at the converged principal state; condensed K̂ | K(θ) = T(θ)ᴴKT(θ) = Σ e^{im·θ}K_m; condensation optional | **Same object** |
| Verdict | LDU sign of K̂ | Three-state shifted Cholesky at ±τλ_ref; LDLᵀ only as a diagnostic | **Compatible**; NEMOpt is stricter (margin band) |
| Sampling | Uniform full-zone 36×36 grid at every load step; no symmetry reduction | Level 1: TRIM points (real arithmetic); Level 2: IBZ contour ("not guaranteed"); Level 3: 16×16 or 32×32 half-zone grid refined near the axes | **Compatible.** All TNS finite-k criticals lie at the M point, a TRIM point. The 16×16 grid (22.5°) is coarser than TNS's 10° and the 32×32 grid (11.25°) slightly coarser, so reproductions need at least 36×36 full zone |
| k = 0 | Grid point; corner A pinned | Center pinned only at θ = (0, 0) (Tier 0a, both parity sectors) | Same at θ = (0, 0); **conflict on the axes** (C1) |
| k → 0 | Ellipticity of L^H (Eq. 2.13); grid points no closer than θ = π/18 | A small ball around θ = 0 excluded; Tier 0b exact quartic test; E-S4 consistency for θ ∈ [1e−3, 1e−1] | **Same principle**; TNS give the earlier primary statement |
| Existence of the homogenized tangent | Needs a unique cell solution, i.e., no ω_c = 0 bifurcation (p. 508) | Tier 0a and Tier 0b both at every converged step | **Consistent**; state the dependence (E2) |
| Domain for equilibrium and for stability | Quarter cell (θ = 0), Bloch on the full cell | Half cell (C2) or full cell; Bloch on the reconstructed full cell (NDR-0051) | **Consistent**; TNS are a precedent (E5) |
| Load path | Radial in principal log strain (Eq. 3.2) | Linear F ramp | **Different** in general; identical on diagonal-F_M rays with a constant strain ratio; compare crossing points otherwise |
| Load resolution | Δλ = 10⁻³ marching, no bisection | 10 legacy steps; t_c bisection in audits | Compatible |
| Eigenvalue metric (for margins) | H¹-seminorm in the theory (Eqs. 2.4, 2.7); none in the numerics | Identity (lumped) for diagnostics; pseudo-mass for optimization | No conflict for signs; see C5 |
| Modes | Not computed by TNS (Fig. 5 by Michel) | Lanczos margins and modes | — |

### 4.3 Conflicts and cautions

- **C1. The V_Ai = 0 rule must not be copied.**
  - Read literally, TNS constrain corner A's component i whenever ω_i = 0 alone, i.e., also along the zone axes ω = (0, ω₂ ≠ 0) and (ω₁ ≠ 0, 0).
  - No rigid mode exists there, because a uniform translation violates the Bloch condition in the other direction. The constraint is therefore spurious.
  - By Cauchy interlacing, a one-DOF constraint can lift a single negative eigenvalue above zero. It can only delay detection, for example of a 1×2 or 2×1 mode at an X point, where the addendum's exp1 already showed that a pin raises the spectrum (inference).
  - TNS's reported criticals (the M point, or ω → 0) are unaffected.
  - NEMOpt keeps its rule: pin only at θ = (0, 0). If the phrase meant "ω = 0" (both components), the two methods agree. The text is ambiguous.
- **C2. TNS's condensation silently needs K_ii to be positive definite.** The sign of K̂ equals the full Bloch verdict only then (Haynsworth). TNS never state this. If NEMOpt adopts condensation, it must check K_ii by Cholesky per state, or track In(K_ii) as exp3 did.
- **C3. Unit conversion.** Translate TNS's wave numbers with θ_i = 2ω_iL_i. A slip would put the benchmark's critical mode at X instead of M and make B-S2a fail for the wrong reason.
- **C4. Grid resolution.** The Level 3 16×16 grid is coarser than TNS's. This does not matter for TNS's M-point criticals, which Level 1 catches exactly. For a claim that "no other k is critical", use at least 36×36 full zone.
- **C5. Margins under a gradient metric (option, not a conflict).**
  - TNS's and GMT's Rayleigh quotient uses the H¹-seminorm. With that metric the long-wave eigenvalue tends to a finite limit, min over n of the smallest eigenvalue of Q(n), instead of vanishing like |θ|². This would make margins comparable across θ and give E-S4 an unscaled form (inference).
  - The identity metric stays the default for verdicts (Sylvester).
  - Whether round-off near θ = 0 improves is untested.
- **No conflict** with the long-wave plan. TNS do exactly what the addendum plans: the ball near k = 0 is left to the 𝔸_M ellipticity test, and k = 0 is checked separately (Tier 0a covers TNS's case "ω_c → 0 is a regular point … still a local mode that is D periodic", p. 507).

---

## 5. Recommended edits (Task 4; not applied)

Old text is quoted from the current `reports/NEMOpt buckling and stability addendum.md`. The GMT note ([geymonat_muller_triantafyllidis_1993.md](geymonat_muller_triantafyllidis_1993.md), its E3–E9) proposes edits to the same sentences at L46, L56, L58 and L178. Where both apply, **apply the GMT edits first and append the TNS additions below**. The "append" lines are written so that they fit after either version. Link targets: from `reports/` the TNS link is `https://doi.org/10.1115/1.2126695`, and this note is `../research_notes/NEMOpt%20supplied%20papers%20verification/triantafyllidis_nestorovic_schraad_2006.md`.

**E1 (L56; after GMT's E3/E4 if applied).**
- Old: "Microscopic instability therefore precedes macroscopic instability or coincides with it."
- New: "Microscopic instability therefore precedes macroscopic instability or coincides with it: the microscopic onset-of-failure surface "always lies inside its macroscopic counterpart" ([Triantafyllidis, Nestorović & Schraad 2006, p. 506 and Eqs. 2.8, 2.13](https://doi.org/10.1115/1.2126695); verified), and a co-author's summary there confirms the two GMT results used here (p. 505)."

**E2 (L58; append after the sentence's citation, in whichever version is current).**
- Old: "so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); verified)."
- New: "so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); [Triantafyllidis et al. 2006, p. 507 and Eq. 2.13](https://doi.org/10.1115/1.2126695), who separate "strictly periodic modes ω = 0" from "long-wavelength modes ω → 0"; verified). The homogenized tangent itself exists only while the one-cell problem has a unique solution, "thus excluding a first bifurcation with ω_c = 0" ([Triantafyllidis et al. 2006, p. 508](https://doi.org/10.1115/1.2126695)), so Tier 0a must pass for Tier 0b to mean anything."

**E3 (L131, grid sentence).**
- Old: "and Level 3, for final designs, scans a uniform 16×16 or 32×32 half-zone grid refined near the axes, in the manner of Zhang et al.'s 100×100 grid with three refined strips ([Zhang et al. 2021](https://arxiv.org/abs/2010.02371); verified)."
- New: "and Level 3, for final designs, scans a uniform 16×16 or 32×32 half-zone grid refined near the axes, in the manner of Zhang et al.'s 100×100 grid with three refined strips ([Zhang et al. 2021](https://arxiv.org/abs/2010.02371); verified) and of Triantafyllidis et al.'s uniform full-zone grid with phase spacing π/18, swept at every load increment ([Triantafyllidis et al. 2006, pp. 509–510](https://doi.org/10.1115/1.2126695); verified; their ω_iL_i uses the half-width L_i, so θ_i = 2ω_iL_i). Reproductions of their results (B-S2a) use at least a 36×36 full-zone grid."

**E4 (L131, condensation sentence).**
- Old: "Boundary condensation by Haynsworth inertia additivity reproduced the direct inertia at 144 of 144 wave vectors ([exp3.jl](../research_notes/NEMOpt%20buckling%20and%20stability%20addendum/evidence/stability_experiments/exp3.jl); [Haynsworth 1968](https://doi.org/10.1016/0024-3795(68)90050-5); verified), but it pays off only above a few hundred k-points per state, since a direct complex Cholesky already costs 0.05–0.12 s, so it stays optional."
- New: "Boundary condensation is Triantafyllidis et al.'s device: interior DOFs are eliminated once per load step, and only a Hermitian Schur complement on the left, bottom and corner DOFs, about 4N×4N for an N×N cell, is sign-tested per wave vector by LDU ([Triantafyllidis et al. 2006, Eqs. 3.7–3.12](https://doi.org/10.1115/1.2126695); verified). That sign alone gives the full verdict only if K_ii is positive definite, which the paper does not state; Haynsworth inertia additivity, In(K) = In(K_ii) + In(S), removes the restriction and reproduced the direct inertia at 144 of 144 wave vectors ([exp3.jl](../research_notes/NEMOpt%20buckling%20and%20stability%20addendum/evidence/stability_experiments/exp3.jl); [Haynsworth 1968](https://doi.org/10.1016/0024-3795(68)90050-5); verified). It pays off only above a few hundred k-points per state, since a direct complex Cholesky already costs 0.05–0.12 s, so it stays optional."

**E5 (L139).**
- Old: "The literature offers no shortcut here: point-symmetric homogenization papers report large speedups and say nothing about bifurcation,"
- New: "The closest precedent is Triantafyllidis et al., who computed the principal solution on one quarter of the cell under loading along the axes of orthotropy and ran the Bloch analysis on the full cell ([Triantafyllidis et al. 2006, p. 509](https://doi.org/10.1115/1.2126695); verified), without discussing which modes the reduced domain misses. Otherwise the literature offers no shortcut: point-symmetric homogenization papers report large speedups and say nothing about bifurcation,"

**E6 (L46, row 7; append to the resolution cell after GMT's E6 if applied).**
- Old (current text): "Whether micro and macro onsets coincide for compressible laminates is an open item, not an acceptance criterion."
- New (append): "… not an acceptance criterion. Triantafyllidis et al. 2006 validated their Bloch and ellipticity code on two-layer cells against [Nestorović & Triantafyllidis 2004](https://doi.org/10.1016/j.jmps.2003.06.001) (pp. 510–511), the further source for compressible laminate onsets."

**E7 (L178, B-S1 source cell; append).**
- Old: "[Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148) |" (the last entry of the B-S1 row)
- New: "[Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148) (incompressible, per Triantafyllidis et al. 2006, p. 505); layered cells as a Bloch-code check: [Triantafyllidis et al. 2006, pp. 510–511](https://doi.org/10.1115/1.2126695) (verified) |"

**E8 (L179, B-S2 row; replace with two rows).**
- Old: "| B-S2 | Square array of circular holes, porosity 0.5–0.6, uniaxial plane-strain compression | First instability at the M point θ = (π, π), a 2×2 pattern found in real arithmetic; ellipticity lost later; below porosity ≈ 0.34, localization first. Quantitative acceptance only against an independent body-fitted computation with NEMOpt's exact material | [Bertoldi et al. 2010](https://doi.org/10.1002/adma.200901956) (verified: "localization occurs prior to the pattern transformation" below 0.34; reading it as macroscopic-first is inference); [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) (non-peer-reviewed: onset at engineering strain −0.03, ellipticity lost at −0.055 for a Rivlin rubber) |"
- New:
  - "| B-S2a | Square array of circular inclusions after Triantafyllidis et al. 2006: 2×2 cell (L = 1), R = 0.5, volume fraction π/16; both phases W = μ/2(I₁ − 2) − μ ln J + κ/2(J − 1)² with the same κ/μ (9.8 primary, 98 secondary), inclusion-to-matrix shear ratio 0.02 (their "voids", soft filler) or 50 ("rigid"); strain control along F_M = sI (φ = 5π/4) and F_M = diag(s, 1) (φ = π) | Soft, κ/μ = 9.8, φ = 5π/4: first instability at the M point θ = (π, π), a 2×2 antisymmetric pattern found in real arithmetic, at principal log strain −0.089, before loss of ellipticity at −0.103; φ = π: Bloch onset −0.178 before loss of ellipticity −0.194. Stiff, κ/μ = 9.8: no finite-k onset before loss of ellipticity at −0.545 (φ = π) and −0.654 (φ = 5π/4). Soft, κ/μ = 98: assert only t_micro ≤ t_macro near −0.068. Tier 0a passes up to every onset. Values are digitized (±0.003 soft, ±0.01 stiff); a body-fitted, locking-free run with the same law must match within 7%, and NEMOpt's staircase cell must match that run within a tolerance pre-registered from an 80×80/160×160 refinement | [Triantafyllidis et al. 2006, Eqs. 3.1–3.3, Figs. 4–6](https://doi.org/10.1115/1.2126695) (peer-reviewed; verified; no tables, values digitized in [the TNS note](../research_notes/NEMOpt%20supplied%20papers%20verification/triantafyllidis_nestorovic_schraad_2006.md)) |"
  - "| B-S2b | Square array of circular holes, porosity 0.5–0.6, uniaxial compression with zero lateral macroscopic traction, which full-F_M control cannot impose without an outer P_M22 = 0 loop | Qualitative only: first instability at the M point, ellipticity lost later; below porosity ≈ 0.34, localization first under this load path (at porosity 0.196, Triantafyllidis et al. 2006 still find the M-point mode first near equibiaxial compression) | [Bertoldi et al. 2010](https://doi.org/10.1002/adma.200901956) (verified: "localization occurs prior to the pattern transformation" below 0.34; reading it as macroscopic-first is inference); [Bertoldi 2008 slides](https://www2.msm.ctw.utwente.nl/sluding/TEACHING/CTW_Mechanics/MSM_Bertoldi_workshop220908.pdf) (non-peer-reviewed: onset at engineering strain −0.03, ellipticity lost at −0.055 for a Rivlin rubber); quantitative targets pending Michel et al. 2007 or Bertoldi et al. 2008 |"

**E9 (L185).**
- Old: "Two benchmark limits deserve emphasis. The hole-array numbers in the literature come from Rivlin, Gent or near-incompressible materials, which differ from NEMOpt's compressible Bonet–Wood-type law, and staircase holes with residual void stiffness shift the critical strain by an amount not known in advance (inference); those numbers are cross-code comparisons, not reproductions."
- New: "Two benchmark limits deserve emphasis. Most hole-array numbers in the literature come from Rivlin, Gent or near-incompressible materials. The exception, Triantafyllidis et al. 2006, shares NEMOpt's shear part μ/2(I₁ − 2) − μ ln J and differs only in the volumetric term, κ/2(J − 1)² with κ/μ = 98 or 9.8 against λ/2(ln J)² with λ/μ = 1.5; their "voids" are soft inclusions with the same κ/μ, not empty holes (verified). With a test-only volumetric option and a two-phase nonlinear material map, B-S2a becomes a reproduction up to digitization and discretization error. Staircase inclusions and, at κ/μ = 98, volumetric locking of the fully integrated Q4 elements shift the critical strain by an amount not known in advance (inference), and B-S2b and the remaining hole-array numbers stay cross-code comparisons."

**E10 (L165).**
- Old: "the two that do (B-S1 and B-S2) run in the Default profile"
- New: "the three that do (B-S1, B-S2a and B-S2b) run in the Default profile"

**E11 (L247, L252, L314; consistency).**
- L247: "Compressive benchmarks (B-S1, B-S2)" → "Compressive benchmarks (B-S1, B-S2a, B-S2b)".
- L252: "Body-fitted references for B-S2 may use" → "Body-fitted references for B-S2a and B-S2b may use".
- L314 (NDR-0058 acceptance): "B-S1 and B-S2 run in Default" → "B-S1, B-S2a and B-S2b run in Default".

**E12 (L306, NDR-0050). Four substring replacements.**
- "no pin for θ ≠ 0;" → "pin only at θ = (0, 0), never on the lines θ₁ = 0 or θ₂ = 0 (Triantafyllidis et al. 2006 constrain a corner DOF whenever one ω_i = 0; not adopted);"
- "optional Haynsworth condensation;" → "optional boundary condensation as in Triantafyllidis et al. 2006, Eq. 3.12, with Haynsworth inertia additivity so that K_ii need not be positive definite; reproductions of published grids use at least their 36×36 full-zone resolution;"
- "| E-S3, E-S4 (θ ∈ [1e−3, 1e−1]), E-S6, B-S1, B-S2 |" → "| E-S3, E-S4 (θ ∈ [1e−3, 1e−1]), E-S6, B-S1, B-S2a, B-S2b |"
- "| Literature: Åberg & Gudmundson 1997; Triantafyllidis et al. 2006; Zhang et al. 2021; Wang & Sigmund 2023 |" → "| Literature: Åberg & Gudmundson 1997; Triantafyllidis et al. 2006 (phase relations, full-zone grid, condensation, long-wave limit via L^H); Zhang et al. 2021 (no pin at k ≠ 0); Wang & Sigmund 2023 (IBZ) |"

**E13 (L310, NDR-0054).**
- Old: "B-S1 asserts t_micro ≤ t_macro; B-S2 against an independent body-fitted computation (NDR-0042 unaffected); compressive cases need NDR-0031 and NDR-0032"
- New: "B-S1 asserts t_micro ≤ t_macro; B-S2a against digitized Triantafyllidis et al. 2006 values through a body-fitted, locking-free computation with their energy, which needs two test-only options that leave FEAD2026 untouched: the volumetric term κ/2(J − 1)², and a two-phase material map with both phases on the nonlinear branch, because a phase at g = 0.02 falls below the ω₂ = 0.05 switch of H_Steep; B-S2b qualitative until Michel et al. 2007 or Bertoldi et al. 2008 is read and a P_M22 = 0 outer loop exists (NDR-0042 unaffected); compressive cases need NDR-0031 and NDR-0032"

**E14 (L304, NDR-0048; optional).**
- Old: "margins by shift-invert Lanczos on P·G·P (= K⁺), never from the pinned spectrum"
- New: "margins by shift-invert Lanczos on P·G·P (= K⁺), never from the pinned spectrum; the gradient (H¹-seminorm) metric of Triantafyllidis et al. 2006 (Eqs. 2.4, 2.7) is an option for margins comparable across θ, to be evaluated"

**E15 (L301 and L305, NDR-0045 and NDR-0049; append after GMT's E8/E9 if applied).**
- NDR-0045 origin: append "; Triantafyllidis et al. 2006 (onset of failure as the first bifurcation of the principal branch, p. 505)".
- NDR-0049 literature: append "; Triantafyllidis et al. 2006 (Eq. 2.11; L^H requires a unique cell solution, p. 508)".

**Research notes (for consistency; optional).**
- `instability_theory_detection.md` L71.
  - Old: "The established workflow (Triantafyllidis et al. 2006; Bertoldi & Boyce 2008; Zhang et al. 2021) has four steps:"
  - New: "The established workflow (Triantafyllidis et al. 2006 for steps 1–3, with an LDU sign test instead of an eigen-solve in step 2; Bertoldi & Boyce 2008; Zhang et al. 2021) has four steps:"
- `stability_implementation_verification.md` L353.
  - Old: "(failure surfaces of porous neo-Hookean solids with circular holes in strain space)."
  - New: "(failure surfaces in log-strain space for compressible neo-Hookean two-phase solids, κ/μ = 98 and 9.8, with circular inclusions at volume fraction π/16 and μ_f/μ_m from 0.02 ("voids") to 50; read in full in [the TNS note](../NEMOpt%20supplied%20papers%20verification/triantafyllidis_nestorovic_schraad_2006.md))."
- Same file, L225–226.
  - Delete "and Triantafyllidis et al. 2006" from L225 and add: "Triantafyllidis et al. 2006 was later read in full; it has no tables (see the TNS note)."
  - Replace L226 with: "Triantafyllidis et al. 2006 condense onto the left, bottom and corner DOFs (their Eq. 3.12) and sign-test the Schur complement by LDU without stating that K_ii must be positive definite; the Haynsworth-based version above generalizes it."
- Same file, L377. Remove "Triantafyllidis et al. 2006" from the "not accessible" list and add "Triantafyllidis et al. 2006 reports values only in figures; digitized values are in the TNS note."
- Same file, L239. Check "Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134" against the 1984 J. Appl. Mech. paper that TNS cite for porous solids (Section 2.2).

---

## 6. Digitization record

- **Inputs.** The embedded figure images were extracted with `pdfimages`: Figs. 4–9 at 300 ppi, about 950 × 1,900 pixels per two-panel figure.
- **Axis frames.** Located from rows and columns with more than 400 dark pixels. The tick values on the frames give the data-to-pixel map, with origin lines checked against the frame midpoints to within 1–2 pixels.
- **Sampling.** For each ray angle φ, 4,000–6,000 points along ε = λ(cos φ, sin φ) were classified as black, red, green or blue by RGB thresholds. Run midpoints give λ.
- **Masks.** Legend boxes, insets and arrow regions were masked. Where a dash gap or an arrow hid a curve, the ray was moved by 1–5° or offset perpendicularly by up to ±0.006, and symmetry about ε₁ = ε₂ was checked.
- **Visual checks.** Fig. 4(b) near the diagonal and near ε₂ = 0, and Fig. 8(a), were also checked by eye on 4–5× enlargements with overlaid data grids.
- **Scratch files.** Scripts and images were left in the session scratchpad (temporary, not part of the repository). The procedure above suffices to redo the extraction.

## 7. Gaps

- Not reported by TNS:
  - the mesh size N and the integration rule;
  - the φ resolution of the radial paths;
  - the resolution of the min over a and n in Eq. 2.11;
  - how L^H is computed in practice;
  - eigenvalue magnitudes and margins.
- The mode type at φ = π for κ/μ = 9.8 is not stated explicitly.
- The diagonal-geometry numbers (Table D) are low precision.
- The identity of TNS ref. [11] with Michel et al. 2007, and the reference-list errors in [5] and [14], are unverified.
- The volumetric option and the 𝔸 formula of Section 3.8 are derived here and not yet tested in code.
