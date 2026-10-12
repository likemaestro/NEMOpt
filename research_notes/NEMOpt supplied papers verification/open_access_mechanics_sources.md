# Open-access sources for three mechanics citations: verification notes

Date: 2026-10-12. Scope:
1. Ohno, Matsuda & Wu (2001), the point-symmetric half-cell reduction (addendum NDR-0051).
2. Kazemi, Vaziri & Norato (2020), the Poisson's-ratio definition credited for article Eq. 54 (main report, NDR-0037).
3. The 2D Hashin–Shtrikman (HS) upper bulk bound for anisotropic two-phase cells (main report, NDR-0064).

Reports checked:
- `reports/NEMOpt refactoring strategy literature.md`, called "main report" below.
- `reports/NEMOpt buckling and stability addendum.md`, called "addendum" below.

Neither report was edited.

Source policy:
- Only legal sources were used: arXiv; publisher open-access pages (J-STAGE); institutional repositories (Swansea Cronfa; NSF Public Access Repository, a U.S. federal funder repository comparable to OSTI); the author's own homepage (gwmilton.com); and the official JSPS KAKEN grant database.
- No shadow library was used.
- All downloads stayed in the session scratchpad (`oa_downloads/<name>/`). None was copied into the repository.

Labels used below:
- **(derived)** marks my own algebra.
- **(computed)** marks numbers from short scratch scripts that are not in the repository.
- **(inference)** marks reasoning that no source states.

## Verdicts at a glance

| # | Report claim | Verdict | Decisive source |
|---|---|---|---|
| 1a | Half-cell reduction with odd fluctuations, credited to Ohno et al. 2001 (main L36; addendum L15, L72, NDR-0051) | **Supported, from secondary sources.** The primary 2001 text was not available. | Matsuda et al. 2002 p. 538; KAKEN 09450046; Saavedra Flores 2011 §4.3 |
| 1b | Reduction "only on the symmetric branch"; full/half agreement never certifies stability | **Correct, but derived.** No open source on the reduction discusses stability, uniqueness or bifurcation. Ohno *assumes* a point-symmetric internal state. | Same as 1a |
| 1c | Addendum L15 cites Ohno et al. 2001 next to the parity-split and stability argument | **Over-attributed.** Ohno supports the reduction only. | Same as 1a |
| 2 | Article Eq. 54, ν = 2𝔸₁₁₂₂/(𝔸₁₁₁₁ + 𝔸₂₂₂₂), "after Kazemi et al." (main L240) | **Not supported.** Kazemi et al. Eq. 24 is ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)]: 3D, linear, cubic lattices. | CMAME 363:112895 §6.3, p. 8 |
| 3a | The isotropic 2D HS bulk bound applied to the hydrostatic-strain modulus of an anisotropic cell is "an assumption still to be confirmed" (main L171; NDR-0064 gate) | **Confirmed. The bound holds** for any geometry of two well-ordered isotropic phases, in 2D, with the void as a limit. | Milton 2002 §23.4, Eq. 23.30, p. 466 |
| 3b | "Do not apply the isotropic shear bound to anisotropic cells" (NDR-0064 a) | **Correct.** Trace-type and single-modulus replacements exist. | Milton §23.5, Eqs. 23.40, 23.44, 23.47; Eq. 23.30 |
| 3c | Nonlinear-bounds bullets in NDR-0064 | **Mostly supported.** The convexity wording needs refining, and the Müller (1987) DOI is wrong. | Milton pp. 14, 271, 282–286 |

## Access log

| Target | Route tried | Outcome |
|---|---|---|
| Ohno, Matsuda & Wu 2001, IJSS 38:2867 | Web search; KAKEN; J-STAGE; Nagoya University profile pages; Wayback Machine | No legal open copy found. The Wayback Machine (web.archive.org) is unreachable from this session: the proxy tunnel closes on every request. |
| Saavedra Flores & de Souza Neto 2010, Eng. Comput. 27:551 | Emerald; Swansea Cronfa record 14268 | Full text paywalled. A search-engine snippet of the Cronfa record says no full text is held there; the page itself could not be loaded because the host was blocked. Substitute: the author's 2011 open PhD thesis, Ch. 4. |
| Saavedra Flores 2011 thesis | Cronfa direct download | Plain curl was blocked (both cronfa hosts). Obtained through WebFetch from the same Cronfa URL. |
| Kazemi et al. 2020 | arXiv; NSF Public Access Repository | Both obtained. The NSF file is the Elsevier-typeset version of record. |
| Milton 2002 | gwmilton.com, the link supplied by the user | Obtained. The domain is "Home Page of Graeme Milton", but the homepage itself links only the companion book `ExtTheoryCompositesNOPRINT.pdf`. |

## Item 1: Ohno, Matsuda & Wu (2001), point-symmetric half cell

### Sources obtained (all open; the primary 2001 text was not obtained)

- **S1.1**: Matsuda, Ohno, Tanaka & Shimizu (2002), *JSME Int. J. Ser. A* 45(4):538–544 (English). On p. 544, its ref. (13) is Ohno, Matsuda & Wu 2001.
- **S1.2**: Matsuda, Ohno, Tanaka & Shimizu (2003), *J. Soc. Mater. Sci. Japan* 52(10):1258–1264 (Japanese). Its ref. 13 is Ohno, Matsuda & Wu 2001 (p. 1264).
- **S1.3**: Ohno (2000), "Application of homogenization theories to inelastic analysis", JSME CMD conference 13:753–754 (Japanese). Its ref. (3) is Ohno, Matsuda & Wu, IJSS, "in print".
- **S1.4**: JSPS KAKEN project 09450046, final-report summary (Ohno, Nagoya University, FY1997–1998).
- **S1.5**: Saavedra Flores (2011), PhD thesis, Swansea University, Ch. 4 "Symmetry conditions". On p. 29 the chapter refers to Saavedra Flores & de Souza Neto (2010).
- **S1.6**: Okumura, Ohno & Noguchi (2000), JSME CMD conference 13:769–770 (Japanese). Related to symmetric bifurcation, but not to the half cell.

### What the sources say

**Premise.** The reduction is conditional on a point-symmetric internal state, which is assumed rather than proven:
- S1.1, p. 538: Ohno et al. showed that "if the microscopic distributions of stress and strain are symmetric with respect to the center of each unit cell", the perturbed velocity field is point-symmetric about the facet centers.
- S1.4: "by assuming the point symmetry of internal distributions with respect to a cell center".
- S1.3, §3, p. 753: 「内部状態がセル中心に関して点対称分布する場合」 ("the case where the internal state is point-symmetric about the cell center").

**Conclusion.**
- S1.4: the perturbed displacement rate "is point-symmetric with respect to cell facet centers as well as the cell center" and "vanishes at these centers". This "led to introducing a semiunit cell".
- S1.2, p. 1258, items (1) and (2) say the same and add 「基本セルの半分を解析領域に取り」 ("half of the basic cell is taken as the analysis domain").
- S1.1, p. 541 applies the reduction because "the hexagonal unit cell has the point symmetry with respect to the cell center".

**Setting.** Ohno's theory family is small-strain and rate-type. S1.1, Eqs. 13–19 (p. 540), gives σ̇ = c(ε̇ − β). S1.3, §2 decomposes the strain rate into elastic plus creep parts. The group's finite-deformation homogenization came later (S1.3 §5; S1.6). Ohno et al. 2001 is therefore, on all open evidence, an infinitesimal-strain theory (inference; the primary text was not read).

**General definition (Saavedra Flores).**
- S1.5, §4.3, p. 33: "point symmetry is defined by the invariance of the RVE material properties and kinematical constraints to rotations of 180°".
- Also on p. 33: any imposed macroscopic strain "will have the same representation in both coordinate systems".
- Eq. 4.8 then gives ũ(y⊕) = −ũ(y⊖) for y⊕ = −y⊖ (Eq. 4.9).
- Eq. 4.10, p. 34: the fluctuation vanishes at the center z⁵ and the edge midpoints z⁶–z⁹. The corners are not listed, but they also vanish, because periodicity and oddness together force it (derived). This matches NEMOpt's six prescribed nodes, which the addendum (L72) describes as four C2-center classes counted with their periodic images.
- p. 29: Ohno's use was "limited to the analysis of periodic media". Saavedra Flores adds linear boundary displacements (§4.3.2, p. 36).
- p. 29: the derivation is "restricted to two-dimensional problems under infinitesimal strains". The examples, however, include "Large strain formulation" (§4.5.5, p. 50).

**Stability and bifurcation.**
- None of S1.1–S1.5 connects the half-cell reduction with stability, uniqueness, bifurcation or buckling. A text search of thesis Ch. 4 for "unique", "bifurc", "stabil" and "buckl" found nothing.
- Saavedra Flores's step on p. 33 ("requires the equilibrium problem to have identical representation … As a result, any fluctuation field satisfies" Eq. 4.8) tacitly assumes a unique solution (inference). Ohno avoids that step by assuming the symmetric state outright.
- In small-strain rate problems with a positive-definite elastic stiffness, the rate solution is unique, so the question does not arise there. In finite-strain hyperelasticity, which is NEMOpt's setting, it does (inference).
- Because the primary text was not read, "Ohno et al. 2001 say nothing about stability" is unconfirmed at the primary level.

**Related but distinct: "symmetric bifurcation".**
- The Ohno group's "microscopic symmetric bifurcation" (S1.3 §5; S1.6 §2.3, p. 770) means that reversing the sign of the spontaneous perturbation does not change the macroscopic state: 「この自発的速度場の符号を変えても巨視的状態の変化には影響を及ぼさない」 ("changing the sign of this spontaneous velocity field does not affect the change of the macroscopic state").
- Its conditions, Eqs. 15–18 of S1.6, are: a nontrivial Y-periodic solution of the incremental cell problem (Eq. 15), subject to Eq. 16, with zero macroscopic strain and stress rates (Eqs. 17–18).
- This is not point symmetry and not the half cell. It is the open conference precursor of Ohno, Okumura & Noguchi (2002, JMPS 50:1125), which the addendum research notes list as unread. Because such modes produce zero macroscopic stress rate, the result agrees with the addendum's statement that 𝔸_M cannot see them (inference).

### Verdicts on report claims

- **Main L36**: "parity rules of the point-symmetric half cell (odd fluctuations …)" and "A 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear load case". **Supported** by S1.5, p. 33 (macroscopic strain invariance; Eq. 4.8) and S1.1, p. 538. At finite strain, R F Rᵀ = F for R = −I (derived).
- **Main L36**: "but only on the symmetric branch". **Correct, but derived.** The sources assume (Ohno) or tacitly need (Saavedra Flores) a unique point-symmetric solution.
- **Addendum L15** (group-theoretic specialist row) cites Ohno, Matsuda & Wu 2001 beside the odd/even block split and the stability claims. **Over-attributed.** Ohno supports the odd-fluctuation half cell only.
- **Addendum L72** (six prescribed nodes; −1 cut-edge pairs; +1 side-edge pairs). **Consistent** with S1.5, Eq. 4.10 and Fig. 4.4, and with S1.4 ("vanishes at these centers").
- **Addendum NDR-0051**, "Origin: … (Ohno et al. 2001)". **Supported as the origin.** The record should state the scope (below).

### Recommended edits

- **Main report L36:** "A 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear load case, but only on the symmetric branch." → "A 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear load case, but only on the symmetric branch (Ohno et al. 2001 assume a point-symmetric internal state; the restriction to the symmetric branch is derived, not stated there)."
- **Addendum L15:** "[Ohno, Matsuda & Wu 2001](https://doi.org/10.1016/S0020-7683(00)00188-8); [Palais 1979]" → "[Ohno, Matsuda & Wu 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) (half-cell reduction only; no stability content in open restatements); [Palais 1979]"
- **Addendum NDR-0051 Context:** "Origin: FEAD 2026 domain-reduction choice (Ohno et al. 2001); W&T 2020 consider no stability." → "Origin: FEAD 2026 domain-reduction choice (Ohno et al. 2001, small-strain elastic–viscoplastic: assumes a point-symmetric internal state and derives fluctuations that are odd about the cell center and facet centers and vanish there; no stability or bifurcation content per Matsuda et al. 2002, JSME Int. J. A 45:538, and Saavedra Flores 2011, §4.3); W&T 2020 consider no stability."
- **Research note** `instability_theory_detection.md` L251–252 and L325 (optional): the Ohno–Okumura–Noguchi condition can now be cited from its open precursor (Okumura, Ohno & Noguchi 2000, JSME CMD 13:769, §2.3, Eqs. 15–18). The absence of stability content is now confirmed for Saavedra Flores 2011, Ch. 4, as well as for the 2010 abstract.

## Item 2: Kazemi, Vaziri & Norato (2020), Poisson's-ratio definition

### Sources

- arXiv preprint 1910.13313 v1 (only version, October 2019).
- Version of record, *CMAME* 363 (2020) 112895, from NSF PAR 10206014. The file carries the journal header; it was "Received 3 June 2019 … accepted 30 January 2020".
- Equations 22–24 are numbered identically in both versions (arXiv p. 9; published pp. 7–8).

### What the paper says

- §6.3, p. 7: "In this work, we only consider effective lattices with cubic symmetry."
- Eq. 22: K = (1/9)(sum of the nine C_iijj). Eq. 23: G = (1/3)(C₂₃₂₃ + C₁₃₁₃ + C₁₂₁₂).
- **Eq. 24, p. 8: ν(z) = C₁₁₂₂ / [2(C₁₁₂₂ + C₁₂₁₂)].** No justification and no isotropy statement accompany it.
- Linear, 3D, small-strain homogenization:
  - p. 8: "six applied unit strains ε⁰⁽ᵏˡ⁾, k, l = 1, …, 3", with the bilinear form in Eq. 33.
  - §7.1, p. 10: "All materials are homogeneous, isotropic, and linearly elastic".
  - The mesh is 64 × 64 × 64.
- The paper contains no 2D formula, no (C₁₁₁₁ + C₂₂₂₂) form and no finite strain. A text search finds "1111" only in Eq. 22.
- For an isotropic tensor, Eq. 24 equals λ/[2(λ + μ)] = ν. For a cubic tensor, the uniaxial-stress ratio along a cube axis is C₁₁₂₂/(C₁₁₁₁ + C₁₁₂₂). The two agree only when C₁₁₁₁ − C₁₁₂₂ = 2C₁₂₁₂, i.e., under isotropy. Eq. 24 is therefore exact only for isotropic tensors (derived).

### Comparison with article Eq. 54

- Eq. 54 is 2D. For a square-symmetric linear tensor (C₁₁₁₁ = C₂₂₂₂) it equals C₁₁₂₂/C₁₁₁₁, the 2D uniaxial-stress Poisson's ratio (derived).
- The two formulas differ even for the isotropic FEAD2026 solid in plane strain (λ = 57.69, μ = 38.46):
  - Eq. 54 gives λ/(λ + 2μ) = 0.4286.
  - The Kazemi form restricted to in-plane indices gives λ/[2(λ + μ)] = 0.300 (computed).
- **Verdict: the article's attribution, and thesis Eq. 7.5's "(Kazemi et al., 2020)", is inaccurate.** The main report repeats it. NDR-0037's *decision* to keep Eq. 54 is unaffected; only the credit is wrong.

### Recommended edits

- **Main report L240:** "| Eq. 54 (after Kazemi et al.), an objective toward −1 |" → "| Eq. 54, an objective toward −1; the article credits Kazemi et al. (2020), whose Eq. 24 is a different, 3D linear formula for cubic lattices, ν = C₁₁₂₂/[2(C₁₁₂₂ + C₁₂₁₂)] |"
- **Main report NDR-0037:** "keep ν (Eq. 54, a deliberate deviation from Wallin & Tortorelli's uniaxial-stress ν)" → "keep ν (Eq. 54, a deliberate deviation from Wallin & Tortorelli's uniaxial-stress ν; for a square-symmetric linear tensor it equals the 2D uniaxial-stress ratio C₁₁₂₂/C₁₁₁₁; it does not appear in Kazemi et al. 2020, whose Eq. 24 is 3D and cubic, so the documentation must not cite Kazemi et al. as its source)"
- **Research notes**, the same correction:
  - `wallin_tortorelli_2020.md` L41, L214, L244 ("(Kazemi et al.)").
  - `numerical_method_choices.md` L310–314. Its open item at L351 is now resolved.
- Whether the article needs a corrigendum for the attribution is the authors' decision.

## Item 3: anisotropic HS bounds (Milton 2002), shear analogs, nonlinear bounds

### Source

- G. W. Milton, *The Theory of Composites*, author-hosted no-print PDF.
- The edition is from Milton & Patton Publishing: "xxviii + 719 p.", ISBN-13 978-1-4835-6919-2.
- Page numbers below are the printed numbers (PDF page = printed page + 25). The page count matches the 2002 Cambridge edition, but its page-by-page pagination was not checked.

### Findings

**F3.1: an energy bound valid for any geometry in 2D (the decisive result).**
- §23.4, p. 465 concerns bounds that "do not depend on the set of constants c(ξ)", i.e., independent of the microgeometry.
- p. 465 also says Eq. 23.28 completely characterizes "the set of all possible (W, E₀, f₁) triplets".
- p. 465: "In two dimensions the formulas for the bounds have been explicitly calculated and depend on the eigenvalues ε₁ and ε₂". The stated conditions are "κ₁ > κ₂ > 0 and μ₁ > μ₂ > 0".
- **Eq. 23.30, p. 466**, upper bound, second regime:

  ε₀·C*ε₀ ≤ μ₁(ε₁ − ε₂)² + [κ₁κ₂ + μ₁(f₁κ₁ + f₂κ₂)]/(μ₁ + f₁κ₂ + f₂κ₁)·(ε₁ + ε₂)²,

  valid if (μ₁ + f₁κ₂ + f₂κ₁)|ε₁ − ε₂| ≤ f₁(κ₁ − κ₂)|ε₁ + ε₂|.
- p. 466: these bounds were "first obtained by Gibiansky and Cherkaev (1984)" and later by Allaire and Kohn (1993a).
- At ε₀ = I (ε₁ = ε₂ = 1) the regime condition holds trivially, and I:C*:I = C₁₁₁₁ + C₂₂₂₂ + 2C₁₁₂₂ = 4κ_NEMOpt. Hence κ_NEMOpt ≤ [κ₁κ₂ + μ₁(f₁κ₁ + f₂κ₂)]/(μ₁ + f₁κ₂ + f₂κ₁), which is algebraically the 2D HS upper bulk bound κ₁ + f₂/[1/(κ₂ − κ₁) + f₁/(κ₁ + μ₁)] (derived).
- Numerically it is identical: 15.3846, 21.3675 and 28.8462 with voids, and 15.4411, 21.4257 and 28.9050 with ε = 5e−4 (computed). These match NDR-0064's figures.
- The same equation at ε₀ = I also gives the HS *lower* bulk bound for this modulus through Eq. 23.29, third regime (derived). With voids that lower bound is zero.

**F3.2: the trace-bound route (the same conclusion).**
- §23.5, pp. 468–473. The 3D upper "bulk modulus type bound", Eq. 23.43 (p. 471): f₂ Tr[Λ_h(C₁ − C*)⁻¹] ≤ 1/[3(κ₁ − κ₂)] − f₁/(3κ₁ + 4μ₁).
- p. 471: the bounds "for elastically anisotropic composites were obtained independently by Zhikov (1988, 1991) and Milton and Kohn (1988)".
- For cubic symmetry, p. 473: "The bounds on κ* remain the same as those given by Hashin and Shtrikman (1963b)".
- 2D analog (derived from the general form, Eq. 23.33, with Tr[Λ_hΓ] = 1/[2(κ + μ)]): f₂ Tr[Λ_h(C₁ − C*)⁻¹] ≤ 1/[2(κ₁ − κ₂)] − f₁/[2(κ₁ + μ₁)].
- By Cauchy–Schwarz, (I:I)² ≤ (I:A⁻¹:I)(I:A:I) with A = C₁ − C*, which yields the same HS bound on ¼ I:C*:I (derived). For square symmetry, I is an eigenvector of C*, so the bound follows directly (derived).
- §24.2, p. 501 derives the 2D bound only for "a two-dimensional isotropic composite", so do not cite §24.2 for anisotropy.

**F3.3: conditions of validity for NDR-0064 (a).**
- The cell has two isotropic, well-ordered phases (κ₁ > κ₂, μ₁ > μ₂). The FEAD2026 solid and ε-material (C_ε = εC_s) qualify.
- The void case is a limit:
  - Milton states Eq. 23.30 for κ₂, μ₂ > 0.
  - He treats voids as limits (p. 467: "By considering the limiting case where one phase is void").
  - C* increases with the phase tensors, so the void bound follows by continuity (derived).
- The bound applies to the energy under hydrostatic average strain, ¼ I:C*:I, which is exactly NEMOpt's κ. It covers any geometry and any effective anisotropy (square, orthotropic or general).
- It is a linear, 2D (plane-elasticity) result with the 2D bulk modulus κ = λ + μ. Milton links 2D results to "fiber reinforced materials" (p. 474).
- It does not apply to gray SIMP designs, which are not two-phase. This agrees with NDR-0064.

**F3.4: shear.**
- The isotropic HS shear bound, Eq. 23.45, holds "when the composite is isotropic" (p. 471).
- Anisotropic replacements are trace bounds:
  - Eqs. 23.40 and 23.44 (3D).
  - Eq. 23.47 for cubic symmetry (p. 473): "a measurement μ*₁ can give information about μ*₂, and vice versa".
  - Eqs. 23.48 and 23.49 ("trace bounds", p. 473).
  - The general Eq. 23.52 (Milton & Kohn 1988).
- **2D square-symmetric trace bound** (derived from Eq. 23.33 with M = Λ_s and Tr[Λ_sΓ] = (κ + 2μ)/[2μ(κ + μ)]):

  f₂{1/[2(μ₁ − μ*₁)] + 1/[2(μ₁ − μ*₂)]} ≤ 1/(μ₁ − μ₂) − f₁(κ₁ + 2μ₁)/[2μ₁(κ₁ + μ₁)],

  with μ*₁ = (C₁₁₁₁ − C₁₁₂₂)/2 and μ*₂ = C₁₂₁₂. At μ*₁ = μ*₂ it reduces to the 2D HS shear bound (13.417 at f₁ = 0.6 with voids; computed).
- **A single shear modulus** (e.g., G₁₂ = C₁₂₁₂): use Eq. 23.30 at pure-shear strain (eigenvalues ±1):

  μ ≤ [μ₁μ₂ + κ₁(f₁μ₁ + f₂μ₂)]/(κ₁ + f₁μ₂ + f₂μ₁)

  This gives 12.41, 16.03 and 19.89 with voids at f₁ = 0.4, 0.5 and 0.6, against isotropic HS values of 7.40, 10.12 and 13.42 (computed).
- An anisotropic cell's G₁₂ can therefore legitimately exceed the isotropic HS shear bound by about 48%. At F_M = I, G_Σ = 𝔸₁₂₁₂ + 𝔸₂₁₂₁ = 2C₁₂₁₂ (stress-free reference; derived).

**F3.5: nonlinear bounds.**
- p. 271: "Talbot and Willis (1985) extended the Hashin and Shtrikman variational principles to nonlinear media".
- Also p. 271: Ponte Castañeda (1991) obtained inequalities "based on comparisons with linear inhomogeneous media".
- §13.7 (pp. 282–286) bounds the effective potential W*(E₀) (Eq. 13.40).
- p. 283 places the setting "assuming that the displacements are still small, so the problem is geometrically linear".
- Also p. 283: "it is often mathematically convenient to assume that W(x, A) is a strictly convex function". This assumption serves the dual principle, Eqs. 13.47–13.49.
- The Talbot–Willis lower bound needs the growth condition Eq. 13.52 (p. 285). The upper bound, Eq. 13.54, needs W to grow no faster than quadratically (p. 286).
- p. 286: "any upper bound on the effective tensor of a linear comparison composite" gives an upper bound on the nonlinear potential.
- Also p. 286: specific nonlinear bounding problems are not discussed.
- §25.1, p. 531 gives a nonlinear translation bound, Eq. 25.8.
- Ponte Castañeda (1989, finite elasticity) is not cited in the book.
- §1.7, p. 14: for nonlinear periodic elasticity the homogenized energy "has to be computed from solutions in increasingly large blocks of unit cells" (Braides 1985; Müller 1987).

### Verdicts on report claims

- **Main L171**, "an assumption still to be confirmed": **resolved, confirmed** (F3.1–F3.3).
- **NDR-0064 (a)** (definition of κ; F_M = I only; thresholded designs with the ε-material; never gray designs): **supported**. The precondition for its hard gate is met.
- **NDR-0064 (a)**, "Do not apply the isotropic shear bound to anisotropic cells": **correct** (F3.4). A replacement should be added.
- **NDR-0064**, "Hashin–Shtrikman bounds are theorems for linear elastic moduli": **supported** (Milton Ch. 23 is linear; nonlinear bounds in §13.7 concern the energy).
- **Talbot & Willis bullet**: **supported** (p. 271; §13.7). Page numbers can now be cited.
- **Ponte Castañeda bullet**: the 1991 part is **supported** (p. 271; pp. 285–286). The 1989 part remains secondary only.
- **"These bound the effective energy, not tangent moduli"**: **supported** (§13.7).
- **"They require convexity that finite-strain hyperelastic energies lack"**: **partly supported.**
  - Milton's §13.7 is geometrically linear.
  - Convexity is assumed for the dual principles, and the primal Talbot–Willis inequality needs growth conditions rather than convexity.
  - That finite-strain W(F) is non-convex is standard, but Milton does not state it.
- **"multi-cell infimum (Müller 1987 …)"**: **supported** by §1.7, p. 14. However, **the report's link is wrong**: doi:10.1007/s00205-016-0992-0 resolves on Crossref to Duerinckx & Gloria (2016), *ARMA* 221:1511–1584. Müller (1987), *ARMA* 99:189–212, is doi:10.1007/BF00284506 (Crossref; Milton's bibliography, p. 16).

### Recommended edits

- **Main report L171:** "- applying the isotropic bulk bound to a square-symmetric cell's hydrostatic-strain modulus is an assumption still to be confirmed." → "- the bound applies to the hydrostatic-strain modulus κ = ¼(𝔸₁₁₁₁ + 𝔸₂₂₂₂ + 𝔸₁₁₂₂ + 𝔸₂₂₁₁) of any anisotropic two-phase cell, square-symmetric or not: Milton (2002, §23.4, Eq. 23.30, p. 466) gives the optimal 2D upper energy bound for any geometry of two well-ordered isotropic phases, and at hydrostatic strain it reduces to the 2D Hashin–Shtrikman bulk bound (it reproduces 28.85 and 28.905; the void is the limit of the ε-material)."
- **NDR-0064, acceptance-evidence cell:** "before (a) becomes a hard gate, confirm in the composites literature (e.g., Milton 2002) that the bulk bound holds for the hydrostatic-strain modulus of anisotropic cells" → "the bulk bound's validity for the hydrostatic-strain modulus of anisotropic cells is confirmed (Milton 2002, Eq. 23.30, p. 466: 2D, two isotropic phases with κ₁ > κ₂ > 0 and μ₁ > μ₂ > 0, the void as a limit), so (a) can be a hard gate for thresholded two-phase designs at F_M = I"
- **NDR-0064 (a):** "Do not apply the isotropic shear bound to anisotropic cells." → "Do not apply the isotropic shear bound to anisotropic cells; for a single shear modulus (G₁₂ = 𝔸₁₂₁₂ at F_M = I) use Milton's Eq. 23.30 at pure-shear strain (12.41, 16.03 and 19.89 at solid fraction 0.4, 0.5 and 0.6 with voids; computed), and for square-symmetric cells optionally the 2D trace bound on the pair of shear moduli (derived from Milton's Eqs. 23.33 and 23.44)."
- **Main report, Talbot bullet:** "- Talbot & Willis (1985) extended the variational principles to nonlinear media ([Milton, *Theory of Composites*](https://resolve.cambridge.org/core/books/theory-of-composites/variational-principles-and-inequalities/F55B22AA23653A981F56C1D0C9F22D78))." → "- Talbot & Willis (1985) extended the Hashin–Shtrikman variational principles to nonlinear media, in a geometrically linear setting (Milton 2002, p. 271 and §13.7, pp. 282–286; [author-hosted copy](https://gwmilton.com/TheoryCompositesNOPRINT.pdf))."
- **Main report, convexity bullet:** "- They require convexity that finite-strain hyperelastic energies lack. Once instabilities appear, the homogenized energy may need a multi-cell infimum ([Müller 1987, ARMA 99:189](https://link.springer.com/doi/10.1007/s00205-016-0992-0), cited from a reference list)." → "- They are posed for geometrically linear problems and rely on growth conditions, and on convexity for the dual principles (Milton 2002, §13.7, pp. 283–286); finite-strain hyperelastic energies are not convex in F. Once instabilities appear, the homogenized energy may need a multi-cell infimum, computed "from solutions in increasingly large blocks of unit cells" (Milton 2002, §1.7, p. 14, citing Braides 1985 and [Müller 1987, ARMA 99:189–212](https://doi.org/10.1007/BF00284506))."
- **Optional, Ponte Castañeda bullet:** append "(the 1991 comparison bounds: Milton 2002, p. 271 and pp. 285–286)".

### Side finding: the shear formula in the bounds-check script

`research_notes/NEMOpt elastic bounds check/hs_simp_check.py` computes the 2D HS shear bound with `z = m1 * (k1 + 2 * m1) / (2 * (k1 + m1))`. The 2D form is `z = k1 * m1 / (k1 + 2 * m1)`, i.e., 1/(μ₁ + ζ₁) = (κ₁ + 2μ₁)/[2μ₁(κ₁ + μ₁)] (derived).

Evidence that the derivation method is right: the same Green-operator calculation in 3D reproduces Milton's Eqs. 23.44–23.45 exactly.

Effect of the error:
- The printed bounds with voids are 7.96, 10.82 and 14.23; the correct values are 7.40, 10.12 and 13.42.
- The SIMP shear check is unchanged: there is no excess at p = 3 with either formula (maximum ratio 0.9996 with the script's formula and 0.9998 with the corrected one).
- The report text quotes no shear-bound values.

## What was obtained, and what still needs the user

| Item | Obtained? | Source URL |
|---|---|---|
| Milton, *The Theory of Composites* (full text) | Yes, author-hosted | https://gwmilton.com/TheoryCompositesNOPRINT.pdf (author homepage: https://gwmilton.com/) |
| Kazemi, Vaziri & Norato 2020, preprint | Yes, arXiv v1 | https://arxiv.org/abs/1910.13313 |
| Kazemi, Vaziri & Norato 2020, version of record | Yes, NSF PAR | https://par.nsf.gov/servlets/purl/10206014 (record: https://par.nsf.gov/biblio/10206014) |
| Ohno, Matsuda & Wu 2001, IJSS 38:2867 | **No legal open copy found** | https://doi.org/10.1016/S0020-7683(00)00188-8 |
| Matsuda et al. 2002, JSME Int. J. A 45:538 (restates Ohno 2001) | Yes, J-STAGE | https://www.jstage.jst.go.jp/article/jsmea/45/4/45_4_538/_pdf |
| Matsuda et al. 2003, J. Soc. Mater. Sci. Japan 52:1258 | Yes, J-STAGE | https://www.jstage.jst.go.jp/article/jsms1963/52/10/52_10_1258/_pdf |
| Ohno 2000, JSME CMD 13:753 (review) | Yes, J-STAGE | https://www.jstage.jst.go.jp/article/jsmecmd/2000.13/0/2000.13_753/_pdf |
| Okumura, Ohno & Noguchi 2000, JSME CMD 13:769 | Yes, J-STAGE | https://www.jstage.jst.go.jp/article/jsmecmd/2000.13/0/2000.13_769/_pdf |
| KAKEN 09450046 final-report summary | Yes, web page | https://kaken.nii.ac.jp/en/grant/KAKENHI-PROJECT-09450046/ |
| Saavedra Flores & de Souza Neto 2010, Eng. Comput. 27:551 | **No** (paywalled) | https://doi.org/10.1108/02644401011044612 |
| Saavedra Flores 2011, PhD thesis, Ch. 4 | Yes, Swansea Cronfa | https://cronfa.swan.ac.uk/Record/cronfa43160 |
| Müller 1987 DOI check | Yes, Crossref metadata | https://api.crossref.org/works/10.1007/BF00284506 |

Still needs the user:
1. **Ohno, Matsuda & Wu 2001, full text**, through institutional ScienceDirect access or a request to the authors. It is needed to confirm the primary wording and equation numbers, the infinitesimal-strain setting, and the absence of any stability statement. Until then, item 1 rests on the authors' own restatements (S1.1–S1.4).
2. **Saavedra Flores & de Souza Neto 2010, full text** (optional). Thesis Ch. 4 appears to reproduce it.
3. **Wayback Machine and Cronfa direct access** are blocked by this session's egress proxy. A check of archived Ohno-lab pages for an author copy would need another environment.
4. **Milton pagination.** If the citation must match the 2002 Cambridge print or the SIAM Classics reprint page for page, someone should spot-check pp. 465–473 against that print.
5. **Decision:** whether to report the article's incorrect Kazemi attribution, for Eq. 54 and thesis Eq. 7.5, as an erratum item. This is the authors' decision.

---
*Evidence note (coordinator).* The bound computations above are preserved in [`evidence/milton_bounds/`](evidence/milton_bounds/) (Python, run with `python3 -I <script>` from that folder): `milton_2330_check.py` (Milton Eq. 23.30 at hydrostatic strain equals the 2D Hashin–Shtrikman bulk bound: 15.38 / 21.37 / 28.85 for voids, 15.44 / 21.43 / 28.905 with the ε-material), `shear_bounds_2d.py` (single-shear bound from Eq. 23.30, e.g. 19.894 at 60% solid, against the isotropic shear bound 13.417), and `shear_hs_2d_check.py` (the corrected 2D shear formula versus the one in `hs_simp_check.py`). The coordinator re-ran all three and reproduced these values. The coordinator also confirmed Kazemi, Vaziri & Norato (2020), Eq. 24 (published version, p. 8): ν = C₁₁₂₂ / [2(C₁₁₂₂ + C₁₂₁₂)], for 3D cubic lattices.

*Update after the user supplied the papers (coordinator, 2026-10-12).* Ohno, Matsuda & Wu 2001 was read in full ([`lazarov_sigmund_2011_kubo_matsuda_2016.md`](lazarov_sigmund_2011_kubo_matsuda_2016.md), §3). Verdicts 1a–1c hold at the primary level: infinitesimal strain (p. 2868); the point-symmetric state is assumed (abstract; p. 2869); fluctuations are odd about C₀ and the facet centers and vanish there (Eqs. 8–9); they also vanish at parallelogram vertices (p. 2870), which replaces the derived corner argument; and there is no uniqueness, stability or bifurcation content. In the table above, read "**Supplied by the user** (version of record)" for Ohno et al. 2001; "Still needs the user" item 1 is resolved. Kubo & Matsuda 2016 (*Mech. Eng. Lett.* 2:16-00248), also supplied, is irrelevant to the half-cell claims: one full cell with misaligned periodicity, elastic constituents with damage, and no symmetry or stability content.
