# Lazarov & Sigmund (2011), Kubo & Matsuda (2016) and Ohno, Matsuda & Wu (2001): verification notes

Date: 2026-10-12.

Scope:
- Check every claim that the two reports and the research notes attribute to Lazarov & Sigmund (2011) and to Ohno, Matsuda & Wu (2001) against the primary texts.
- Assess Kubo & Matsuda (2016), which the user supplied in place of Ohno et al. before Ohno et al. itself arrived.
- Propose exact edits.

No report, decision record or other research note was edited.

## Sources and conventions

User-supplied PDFs, all read in full with the Read tool:
- **LS2011**: Lazarov, B.S. & Sigmund, O. (2011), "Filters in topology optimization based on Helmholtz-type differential equations", *IJNME* 86(6):765–781, doi:10.1002/nme.3072. Wiley version of record, published online 28 December 2010, 17 pp.
- **KM2016**: Kubo, G. & Matsuda, T. (2016), "Effects of laminate misalignment on macroscopic strength and microscopic damage development of plain-woven laminates", *Mech. Eng. Lett.* 2:16-00248, doi:10.1299/mel.16-00248, 8 pp.
- **OMW2001**: Ohno, N., Matsuda, T. & Wu, X. (2001), "A homogenization theory for elastic–viscoplastic composites with point symmetry of internal distributions", *IJSS* 38:2867–2878, doi:10.1016/S0020-7683(00)00188-8, 12 pp.

Comparison texts:
- **FEAD 2026**: Güven & Özdemir, *Finite Elem. Anal. Des.* 260:104599, §3–4 (Eqs. 41–50) and §6, from the user-supplied PDF.
- **Code**: `Domain_Reduction/NEMOpt-DR.jl` L242–299 and `Full_Domain/NEMOpt.jl` L700–752 (working tree).

Text was extracted to the session scratchpad only to confirm absences by search. No PDF was copied into the repository.

Conventions:
- Pages are journal pages: LS2011 pp. 765–781, OMW2001 pp. 2867–2878, KM2016 pp. 1–8, FEAD 2026 by article page.
- Equation and figure numbers are each paper's own.
- **(derived)** marks my own algebra. **(computed)** marks a short scratch script that is not in the repository. **(inference)** marks reasoning that no source states.
- Line numbers refer to the files as of 2026-10-12 and may shift, because another editor is working on the reports. Each edit quotes the old text so it can be found.

## Verdicts at a glance

| # | Claim (where) | Verdict | Decisive evidence |
|---|---|---|---|
| L1 | r = R/(2√3) links the PDE length to a classic filter radius (strategy NDR-0039; notes) | **Confirmed, with a scope correction** | LS2011 Eq. 34, p. 775: r_h = 2√3 r. It is derived from the **1D** hat function and the **1D infinite-domain** Green's function (Eqs. 30–33). The same match in 2D gives R ≈ 3.65 r (derived). |
| L2 | One note writes the relation with R and r swapped (`numerical_method_choices.md` L431) | **A notation clash, not a math error** | top82 calls the PDE length Rmin. LS2011 calls it r and calls the hat radius r_h. |
| L3 | Neumann BCs preserve volume | **Confirmed** | Eqs. 12–14, p. 769 |
| L4 | Periodic BCs also preserve volume | **Not in the paper** (true: derived, and stated by W&T 2020) | LS2011 never mentions periodic BCs. It says any BC other than Eq. 12 makes integral (14) differ from one and needs normalization (p. 769). Periodic BCs are an exception it does not discuss. |
| L5 | Neumann BCs distort or attract material at boundaries; padding | **Not in the paper** | No discussion of boundary artifacts, padding or Dirichlet BCs (full reading and text search) |
| L6 | The filter matrix is constant (design-independent) | **Confirmed (implicit)** | Eq. 19 has no design variable. "the condition number of K_f does not depend on the design" (p. 771). The preconditioner "is not modified during the whole optimization process" (p. 771). |
| L7 | Factor once per run (NDR-0015; strategy L268) | **Confirmed, for small problems** | "if the optimization problem is relatively small, the matrix K_f can be factorized and the factorization can be reused" (p. 771). The primary recommendation is PCG with a preconditioner built once. |
| L8 | Gradients use the same filter matrix, one solve per scalar function | **Confirmed** | Eqs. 26–27, p. 772: solve "K_h a_f = s_f", then apply T. K_f is symmetric (Eq. 19), so one factorization serves both directions (derived). |
| K1 | KM2016 restates the point-symmetric half cell | **No** | It cites Ohno et al. 2000/2001 only for the rate-form homogenization and its incremental procedure (pp. 2 and 4). It uses one full cell with a new "misaligned" periodicity. |
| K2 | KM2016 on uniqueness, stability or bifurcation | **Nothing** | Text search: no hits |
| K3 | KM2016 adds to NDR-0051 | **No. It is irrelevant to our claims.** | — |
| O1 | Ohno et al. assume a point-symmetric internal state (strategy L36; NDR-0051) | **Confirmed (primary)** | Abstract; p. 2869; p. 2877 |
| O2 | Small-strain and elastic–viscoplastic (NDR-0051) | **Confirmed (primary)** | "Infinitesimal strain is assumed" (p. 2868); Eqs. 17 and 32 |
| O3 | Fluctuations are odd about the cell center and the facet centers and vanish there (NDR-0051; addendum L72) | **Confirmed (primary)**; the paper adds the vertices | Eqs. 8–9, p. 2870. Vertices of parallelogram cells: p. 2870, stated without proof. |
| O4 | The reduction holds for any macroscopic load, including shear (strategy L36) | **Confirmed at infinitesimal strain**; the finite-strain case is derived | "the transformation does not affect macroscopic loading" (p. 2870); θ = 22.5° and 45° and off-axial examples (Figs. 6–9) |
| O5 | No stability, uniqueness or bifurcation content (addendum L15; NDR-0051) | **Confirmed (primary)** | No such statement. The symmetry of the rate solution is argued by invariance ("ought to satisfy", p. 2870), which tacitly needs a unique rate solution (inference). |
| O6 | FEAD 2026 §3: the symmetric stress state is "deduced" | **Stronger than its source** | Ohno et al. assume it. At finite strain the deduction holds only on the symmetric branch (derived). |

---

## 1. Lazarov & Sigmund (2011)

### 1.1 What the paper says

**Filter and boundary conditions (§3.1, p. 769).**
- The density filter is "defined implicitly as a solution of the Helmholtz PDE" −r²∇²ρ̃ + ρ̃ = ρ (Eq. 11).
- It is solved "with homogeneous Neumann boundary conditions", ∂ρ̃/∂n = 0 (Eq. 12), "imposed on the boundary of the design domain".
- r "plays a similar role as R in the original filter", where R is the support radius of the convolution filter (Eq. 5).
- The filtered field is the convolution (5) with the Green's function of (11)–(12) as its kernel. G ≥ 0 (Eq. 13) and ∫G dξ = 1 (Eq. 14), "and therefore the filter defined by (11) and (12) preserves the volume".
- "If the specified boundary conditions are different from (12), then the integral (14) differs from one and the solution of the PDE has to be normalized."
- For first-order triangles, tetrahedra, quadrilaterals and hexahedra, "the Neumann boundary condition (12) is automatically fulfilled" and needs no "special treatment around the boundaries" (p. 771). It is the natural BC of the weak form.

**Not in the paper** (full reading; the text search finds no hits):
- periodic BCs;
- padding;
- Dirichlet BCs;
- unit cells;
- any boundary artifact of the Neumann BC.

**A notation slip to avoid copying.** The anisotropic form (§3.2, Eqs. 15–16, K_d = Σ r_i² v_i v_iᵀ) is printed as "∇ᵀK_d∇ρ̃ + ρ̃ = ρ". Eq. 29 has the same sign. Eq. 11 and the weak form (Eq. 19) imply a minus sign. Documentation should copy Eqs. 11 and 19, not Eq. 15.

**Length-scale equivalence (§5.1, pp. 774–775).**
- "A relation between the length parameters can be obtained by equating the second moments of the filter weighting functions."
- Eq. 30 gives the 1D hat function, w_h(x) = (1/r_h)(1 − |x|/r_h) on [−r_h, r_h], "where r_h is the radius of the support domain". Eq. 31 gives its second moment, I_h = r_h²/6.
- "As the Green's functions for a finite domain depend on the position of the point, the equivalent length scale is obtained using the Green's functions of the PDE for an infinite domain. In 1D they are given as" w_g(x) = (1/2r) exp(−|x|/r) (Eq. 32), with second moment I_g = 2r² (Eq. 33).
- Equating Eqs. 31 and 33 gives **r_h = 2√3 r (Eq. 34)**.
- The authors apply the 1D result to 2D runs. Fig. 5 uses the classical density filter "with linearly decaying weighting function" at r_h = 0.069 and 0.139 for r = 0.02 and 0.04.
- Helmholtz designs show more contrast "for equivalent length scales", because the Green's function weights nearby points more heavily (p. 775).

**Assumptions of the derivation:**
- an unbounded domain, so no boundary effects;
- 1D kernels;
- a linear hat (cone) weight;
- equality of second moments only.

The paper gives no 2D or 3D derivation.

**Discretization and solution (§3.3, p. 771).**
- The weak form is K_f ρ̃ = p_f (Eq. 18), with K_f,i = ∫[∇N_eᵀK_d∇N_e + N_eᵀN_e]dΩ (Eq. 19) and p_f,i = ρ_i∫N_eᵀdΩ (Eq. 20).
- NEMOpt's `r^2*B'*B + N'*N` and `𝜌_Unfiltered[el]*N'` (NEMOpt.jl L711–712) match Eqs. 19–20 with K_d = r²I.
- On the solver: "As the tangent matrix K_f is always positive definite, the Conjugate Gradient (CG) method is the most effective iterative algorithm".
- Preconditioners: FSAI or block factorization. Because the preconditioner "is not modified during the whole optimization process", it is worth building a good one once and reusing it "in all optimization steps".
- "Alternatively, if the optimization problem is relatively small, the matrix K_f can be factorized and the factorization can be reused in following iterations."
- "even a simple diagonal preconditioner results in very fast convergence".
- In 3D the filter system is 1/3 the size of the state system, so its cost is "less than 1/9 ≈ 11.1%" of the physical problem. Because "the condition number of K_f does not depend on the design", the filter is in practice much cheaper than that.

**Sensitivities (pp. 771–772).**
- ∂c/∂ρ̃ is computed by the adjoint method (Eqs. 21–22).
- With K_f ρ̃ = Tρ (Eq. 24), the paper writes ∂ρ̃/∂ρ = TᵀK_h⁻¹ (Eq. 26) and ∂c/∂ρ = TᵀK_h⁻¹s_f (Eq. 27).
- The inverse "should never be computed explicitly": the sensitivities come from "first solving K_h a_f = s_f, and then applying the averaging operator T on a_f".
- For first-order elements, this is the average of a_f over the element's nodes.
- K_h in Eqs. 26–27 is evidently K_f; the paper never defines K_h.
- K_f is symmetric (Eq. 19), so the gradient solve uses the same matrix and factorization as the forward filter (derived). The paper writes no transpose and does not use the word "symmetric".

**Cost and parallel performance (§5.4–5.5, pp. 768–769, 777–780).**
- **Serial 3D.** The number of CG iterations is "proportional to a⁻¹", and so depends "linearly on the ratio r_h/a". At a = 0.02 with tolerance 10⁻⁷, it takes 14, 22 and 28 iterations for r = 0.02, 0.03 and 0.04.
- **Table II.** The Helmholtz filter takes 0.16–0.34 s. The classical filter takes 0.46–3.8 s with stored neighbors and 3.59–16.31 s with an octree search.
- **Memory.** The classical filter's neighbor storage grows from 12 MB to 1.6 GB for R/a = 2 to 10 at 100 000 cells (pp. 768–769).
- **Recommendation.** "for general applications the PDE filter should be the preferred choice" (p. 779).
- **Multigrid** appears only once: setting up "a multigrid framework for the PDE solution" costs about as much as the neighbor search (p. 778).
- **Parallel.** With few processes and r_h/a < 2, "the Helmholtz sensitivity filter is 20–30% slower" than the traditional implementation. With appropriate preconditioners, "the solution time does not depend on the length parameter" (p. 779).
- **Solvers and scaling.** PCG for the filter and MinRES with an FSAI preconditioner for elasticity. Scaling is near-linear up to 728 processes on a Cray XT4/XT5 (Fig. 15).

### 1.2 Derived: the second-moment match in 2D and 3D

**Derivation (derived).**
- The free-space Green's function of −r²∇²G + G = δ has the Fourier transform 1/(1 + r²|k|²), so its per-axis second moment is 2r² in any dimension.
- For the normalized d-dimensional cone w ∝ (1 − |x|/R)₊, the per-axis second moment is (d + 1)R²/[(d + 2)(d + 3)].
- Equating the two gives R = r·√[2(d + 2)(d + 3)/(d + 1)].

| d | R / r | Difference from 2√3 |
|---|---|---|
| 1 | 3.464 = 2√3 | none: LS2011 Eq. 34 |
| 2 | 3.651 = √(40/3) | about 5% larger |
| 3 | 3.873 = √15 | about 12% larger |

A scratch quadrature of the 1D, 2D (K₀) and 3D Green's functions reproduced the 2r² moment and these ratios (computed).

What this means:
- The 2√3 factor is a 1D result. DTU codes (top82) and the 2016 review apply it in 2D and 3D as a convention.
- The convention is fine, but "equivalent" should not be read as exact in 2D.
- For NEMOpt, r = 1 element corresponds to a cone radius of 3.46 elements by the convention and 3.65 elements by the 2D moment match.
- The discrete filter (LS2011 Eq. 7) on a coarse grid differs from both.

### 1.3 Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| L-a | Strategy NDR-0039 (L467); `numerical_method_choices.md` L28, L420, L439 and log item for L420; `efficient_fe_adjoint.md` L208 | r = R/(2√3), "a second-moment match on an unbounded domain", cited via the 2016 review | **Confirmed in the primary text, with a scope note** | Eq. 34, p. 775. The match is 1D (Eqs. 30–33). Cite LS2011 Eq. 34 directly. |
| L-b | `numerical_method_choices.md` L431 | "R_min = r_min/(2√3)" | **Notation clash** (correct in top82's letters, reversed against L420) | LS2011 writes r_h = 2√3 r. The correction in `open_access_optimization_papers.md` stands. |
| L-c | `numerical_method_choices.md` L420, L439 | r = 1 element ↔ R ≈ 3.46 elements | **Confirmed as the convention** | Eq. 34. The 2D moment match gives 3.65 (derived). |
| L-d | `efficient_fe_adjoint.md` L221 (via the 2016 review) | "volume preserving for Neumann boundary conditions" | **Confirmed (primary)** | Eqs. 12–14, p. 769 |
| L-e | `numerical_method_choices.md` L437; `efficient_fe_adjoint.md` L227 | Periodic BCs also preserve volume; the real difference is boundary distortion | **Not in LS2011.** Both statements are correct but rest on other sources. | LS2011 discusses neither periodic BCs nor boundary artifacts. Its p. 769 sentence ("integral (14) differs from one") holds for Dirichlet BCs, for example, but not for periodic BCs, under which ρ = 1 gives ρ̃ = 1 (derived). Keep W&T 2020 (p. 9) and Dalklint et al. 2021 (p. 1753) as the sources. |
| L-f | `efficient_fe_adjoint.md` L74, L205, L226 | The filter matrix is design-independent: assemble and factor it once | **Confirmed** | Eq. 19; p. 771 |
| L-g | Strategy L268; NDR-0015 (L443) | Factor once per run and share it across objective, constraint and gradients | **Confirmed.** Factorization is for "relatively small" problems; PCG otherwise. | p. 771; Eqs. 26–27 |
| L-h | Strategy NDR-0023 (L451) | One filter back-substitution per scalar function | **Consistent** | Eq. 27: one solve with K_f for each gradient s_f |
| L-i | `efficient_fe_adjoint.md` L205 | The advantages are memory and scaling with filter radius, especially in 3D with iterative solvers | **Confirmed** | Memory: pp. 768–769. CG iterations ∝ r_h/a: p. 778. Table II; p. 780. "Slower in small 2D meshes with direct solves" is top88's finding, not LS2011's. |
| L-j | `efficient_fe_adjoint.md` L226 | Cholesky applies because the operator with periodic tying is SPD | **Partly confirmed** | "K_f is always positive definite" (p. 771, with Neumann BCs). Master–slave tying PᵀK_fP stays SPD (derived). |
| L-k | `numerical_method_choices.md` L431 | The convolution filter regularizes more strongly for large radii; the PDE filter's singular Green's function makes it less robust | **Not in LS2011** (secondary claims, not checked here) | LS2011 reports more contrast for the Helmholtz filter at equal length scale (p. 775). It does not discuss singularity; its 1D kernel is bounded. |
| L-l | `efficient_fe_adjoint.md` L208, L231; `numerical_method_choices.md` L446; `open_access_optimization_papers.md` §6 (L-f) | Full text not accessible; derivation not read | **Resolved** | The user supplied the version of record. |

The reports themselves attribute no boundary-condition or volume-preservation claim to LS2011. Those claims appear only in the research notes (L-d, L-e).

### 1.4 Recommended edits

**Strategy report** (`reports/NEMOpt refactoring strategy literature.md`)

1. **NDR-0039 (L467).**
   - Old: "(Helmholtz r versus cone R ≈ 2√3·r, a second-moment match on an unbounded domain: [Lazarov, Wang & Sigmund 2016, Eq. 10](https://backend.orbit.dtu.dk/ws/files/140201797/main.pdf), attributing it to Lazarov & Sigmund 2011)"
   - New: "(Helmholtz r versus cone R = 2√3·r: [Lazarov & Sigmund 2011, Eq. 34](https://doi.org/10.1002/nme.3072), p. 775, from equal second moments of the 1D hat function and the 1D infinite-domain Green's function (Eqs. 30–33), restated in [Lazarov, Wang & Sigmund 2016, Eq. 10](https://backend.orbit.dtu.dk/ws/files/140201797/main.pdf); the same match in 2D gives R ≈ 3.65·r (derived), so the factor is a convention, not an exact 2D equivalence)"
2. **Speedup table, filter row (L268), evidence cell.**
   - Old: "[Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072) (filter definition; full text not read)"
   - New: "[Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072), p. 771: K_f does not depend on the design (Eq. 19), and for "relatively small" problems "the matrix K_f can be factorized and the factorization can be reused"; the gradient solves with the same K_f (Eqs. 26–27). For large problems they recommend PCG with a preconditioner built once."
3. **NDR-0015 (L443).**
   - Old: "| NDR-0015 | Assemble and factor the filter once per run; share it across objective, constraint and gradients | B |"
   - New: "| NDR-0015 | Assemble and factor the filter once per run; share it across objective, constraint and gradients (K_f is design-independent and its factorization reusable: Lazarov & Sigmund 2011, Eq. 19 and p. 771; the gradient solves with the same symmetric K_f, Eqs. 26–27) | B |"
4. **NDR-0023 (L451), optional.** After "and the step is exact by linearity.", insert: " Lazarov & Sigmund 2011 already compute each filtered gradient with one solve K_f a_f = s_f (Eq. 27, p. 772)."
5. **"Revisions after full-text verification", "Applied: decision log NDR-0001 to NDR-0044".** Add a bullet: "- NDR-0015, NDR-0039 and the filter row of the speedup table now cite Lazarov & Sigmund 2011 directly (Eqs. 19, 26–27 and 34; p. 771), read in full from the user-supplied version of record. NDR-0039 adds that the 2√3 factor comes from a 1D moment match (about 3.65 in 2D, derived). Source: Lazarov–Kubo–Ohno note, §1.4."

**Research notes.** Append these to each file's "Corrections after full-text verification" log and keep the original text.

6. **`numerical_method_choices.md`, existing log items for L420 and L431.** These are under "Filter length convention (Lazarov & Sigmund 2011)".
   - In the L420 item, old: "(Lazarov, Wang & Sigmund 2016, Eq. 10, citing Lazarov & Sigmund 2011)"
   - New: "(Lazarov & Sigmund 2011, Eqs. 30–34, p. 775: r_h = 2√3 r, from 1D kernels on an infinite domain; restated in Lazarov, Wang & Sigmund 2016, Eq. 10)"
   - In the L431 item, old: "per Lazarov, Wang & Sigmund 2016 (Eq. 10, ref. [93]); the 2011 full text has no legal open copy"
   - New: "Lazarov & Sigmund 2011 write r_h = 2√3 r, with r the PDE length and r_h the hat radius (Eq. 34); the user supplied the version of record"
7. **Same file, new item for L431, second sentence.** Add: "The two robustness statements come from the arXiv preprints, not from Lazarov & Sigmund 2011. That paper reports more contrast for the Helmholtz filter at equal length scale (p. 775) and does not discuss a singular Green's function."
8. **Same file, new item for L437.** Add: "Lazarov & Sigmund 2011 prove volume preservation for Neumann BCs only (Eqs. 12–14, p. 769) and say that other BCs require normalization. They do not consider periodic BCs, for which the w = 1 argument above holds."
9. **Same file, new item for L446 (Gaps).** Add: "Resolved: the derivation was read (Lazarov & Sigmund 2011, Eqs. 30–34). It is one-dimensional; the 2D moment match gives R ≈ 3.65 r (Lazarov–Kubo–Ohno note, §1.2, derived)."
10. **`efficient_fe_adjoint.md`, new log items.**
    - L208: "Lazarov & Sigmund 2011 was read in full (user-supplied). r_h = 2√3 r is their Eq. 34, a 1D moment match."
    - L205 and L226: "Confirmed. K_f is design-independent (Eq. 19) and "the condition number of K_f does not depend on the design" (p. 771). For relatively small problems the factorization "can be reused in following iterations" (p. 771). Gradients solve with the same K_f (Eqs. 26–27)."
    - L231: "Resolved. The paper says nothing about periodic BCs. Its timings are 3D: 0.16–0.34 s against 0.46–3.8 s for a classical filter with stored neighbors (Table II). In parallel, the Helmholtz sensitivity filter is 20–30% slower for r_h/a < 2 with few processes (p. 779)."
11. **`open_access_optimization_papers.md`, §6.**
    - Append: "Update (same day): the user supplied the version of record. L-a and L-d are now verified in the primary text (Eq. 34; Eqs. 12–14). L-e is confirmed (p. 771: factorization reuse for small problems). L-f is resolved. See `lazarov_sigmund_2011_kubo_matsuda_2016.md`, §1."
    - In the final table, row 6, change "**Needs the user**" to "**Supplied by the user** (version of record)".

---

## 2. Kubo & Matsuda (2016)

### 2.1 What the paper does

**Purpose and method.**
- It is a multiscale damage analysis of plain-woven E-glass/vinylester laminates with "laminate misalignment".
- It is "based on the homogenization theory for misaligned internal structures developed by the authors" (Matsuda et al. 2011), using "a novel boundary condition for unit cell analysis" (abstract).

**Kinematics.**
- u̇_i = Ḟ_ij y_j + u̇#_i (Eq. 1, p. 3).
- The cell boundary is split into Γ1–Γ12 according to the misalignments a and b (Fig. 2).
- u̇# is periodic along the inclined α and β directions and along y1 and y2 (Fig. 3, p. 3).
- This is a shifted periodicity on one full cell, not a symmetry reduction.

**Constitutive law and damage.**
- σ̇_ij = c_ijkl ε̇_kl, with c the "elastic stiffness" (Eq. 2, p. 4).
- The fiber bundles are "transversely isotropic elastic" and the matrix is "isotropic elastic" (§4.2, p. 6).
- Damage follows Hoffman's criterion (Eqs. 10–11). Damaged stiffnesses are reduced to 1/100 (§3.2, p. 5).
- The strain measure is never stated. Eq. 2 is a small-strain rate law (inference).

**How it cites Ohno et al.** Three times, none of them about symmetry:
- the framework "for nonlinear time-dependent composites (Ohno et al., 2000, 2001)" (p. 2);
- Eq. 6 "has the same form as that of the original homogenization theory (Suquet, 1987; Ohno et al., 2000, 2001)" (p. 4);
- "for more detail on the computational procedure" (p. 4).

**Other observations.**
- The reference list includes Matsuda, Nimiya, Ohno & Tokuda 2007, "Homogenization using reduced domain of analysis" (*Compos. Struct.* 79:493). The text cites it only as an earlier elastic–viscoplastic analysis (p. 2).
- The only hit for "symmetr" is the title of Ohno et al. 2001 in the reference list.
- There are no hits for "unique", "bifurc", "stabil", "buckl" or "half".

### 2.2 Answers to the questions

- **The point-symmetry reduction or its BCs.** Not stated or restated. The paper uses a full unit cell. Its "novel boundary condition" generalizes periodicity to misaligned layers.
- **Material models.** In this paper, linear elastic constituents with stepwise damage by stiffness reduction. Elastic–viscoplastic models appear only in the cited lineage (Ohno 2000, 2001; Matsuda 2007, 2011).
- **Uniqueness, stability, bifurcation.** Nothing. The stress drops in Fig. 7 come from damage (transverse cracks, then fiber-bundle failure), not from an instability analysis.
- **NDR-0051.** Nothing new. The paper neither reduces a cell by symmetry nor compares reduced and full cells.
- **Bottom line.** It is irrelevant to the reports' claims and should not be cited for the half-cell reduction. Its only link to the topic is bibliographic: it lists Matsuda et al. 2007, the group's reduced-domain paper, which the FEAD article already cites as ref. [19].

### 2.3 Recommended edits

- **Reports.** None.
- **`open_access_mechanics_sources.md`, optional.** Add to the access log: "Kubo & Matsuda 2016, *Mech. Eng. Lett.* 2:16-00248: read in full (user-supplied). Irrelevant to the half-cell claims: one full cell with misaligned periodicity, elastic constituents with damage, no symmetry or stability content."

---

## 3. Ohno, Matsuda & Wu (2001)

### 3.1 What the paper says

**Setting.**
- "Infinitesimal strain is assumed to take place in Ω macroscopically and microscopically" (§2, p. 2868).
- The formulation is in rates: u̇_i = Ḟ_ij y_j + ċ_i + u̇#_i (Eq. 1). Here u̇# is the Y-periodic "perturbed velocity".
- F_ij is called the "macroscopic deformation gradient", but only its symmetric part Ė_ij enters (Eq. 3).
- The constituents "exhibit linear elasticity and nonlinear creep" (Eq. 17), with symmetric elastic constants c_ijkl (Eq. 18).
- The examples use Hookean fibers (Eq. 31) and a matrix with power-law creep "in addition to linear elasticity" (Eq. 32), for a SiC/titanium system (Eq. 36). This is an elastic–viscoplastic composite.
- Only "the macroscopically uniform case" is treated (p. 2871).

**The assumption.**
- Abstract: the theory is "rebuilt on the assumption that constituents, stress and strain distribute point symmetrically with respect to the center of each unit cell".
- p. 2869: "Let us now assume that the internal distributions, responsible for perturbed velocity, satisfy the symmetry with respect to the center of Y, C₀".
- Conclusions (p. 2877): "on the assumption that the internal distributions of constituents, stress, strain, etc. ... are symmetric with respect to the unit cell centers".
- So the symmetry of the internal structure **and** of the stress and strain fields is assumed, not derived.
- Orthotropy is not required. The examples include hexagonal arrays and off-axial loading.

**What the paper derives.**
1. **Geometric step** (Eqs. 5–7, p. 2869). Symmetry about C₀ plus Y-periodicity implies symmetry about every facet center C_i, with N = 4 for squares and N = 6 for hexagons. The combined symmetry is called "C-symmetry" (p. 2870).
2. **Invariance step** (p. 2870). A 180° transformation about any C_k leaves the internal distributions unchanged, and "the transformation does not affect macroscopic loading". "Therefore, the field of perturbed velocity ... ought to satisfy the C-symmetry": u̇#(P_k) = −u̇#(P) (Eq. 8).
3. **Zeros.** Hence "u̇# vanishes at C_k (k = 0, 1, …, N)" (Eq. 9). The strain and stress rates are C-symmetric (Eqs. 10–11).
4. **Vertices.** "Incidentally, it can be shown that for parallelogram unit cells, ... Eqs. (8)–(11) also hold at the cell vertices". For parallelepipeds they also hold at edge midpoints (p. 2870). This is stated, not proven.

**Semiunit cell and boundary conditions (§3, pp. 2870–2871; Figs. 2 and 4).**
- The semiunit cell Ỹ is "half a unit cell". On each part of its boundary, "either the C-symmetry or the Y-periodicity can be imposed".
- Square cell (Fig. 2a): C-symmetry on AB (the cut through C₀) and on CD (the facet through C₁); Y-periodicity on AD and BC.
- Hexagonal cell (Fig. 2b): C-symmetry on all sides.
- The boundary integrals of σ̇_ij n̄_j u̇#_i (Eq. 12) and of u̇#_i n̄_j (Eq. 13) vanish. The weak form (Eq. 16) and the cell problems for χ^{kl} and φ (Eqs. 21–22) therefore keep their unit-cell form.
- χ and φ obey the same C-symmetric or Y-periodic conditions (p. 2872). The equations "have the same forms, except for the boundary conditions" (p. 2873).

**FE implementation (pp. 2874–2875).**
- Fig. 4 marks the nodes at which χ and φ "are set to be zero due to the C-symmetric boundary condition". In the square semiunit cell (Fig. 4a) these are only C₀ and C₁.
- The C-symmetric pairs are imposed "in the same way as the Y-periodicity using the penalty method".
- The ends of AB and the top corners still vanish, because each carries two pair constraints: u(A) = −u(B) from C-symmetry and u(A) = u(B) from periodicity (derived).
- The analysis is semi-3D: 2D meshes with three displacement components, uniform in the fiber direction (p. 2875).

**Validation and cost.**
- Half-cell results reproduce earlier full-cell results. The Fig. 7 meshes "are just the upper halves of the deformed meshes of unit cells computed in the previous work". The off-axial curve "is completely the same as the previous one", computed on the full unit cell (p. 2876).
- "computation time can be reduced at most to about one eighth by adopting semiunit cells" (p. 2878).

### 3.2 Uniqueness, stability and bifurcation

- **Nothing in the paper.** The full reading and a text search for "unique", "bifurc", "stabil" and "buckl" both found nothing.
- **The invariance step needs uniqueness.** "ought to satisfy the C-symmetry" (p. 2870) presumes that the internal state and the loading determine u̇# uniquely (inference).
- **In Ohno's setting, uniqueness holds.** At a known state, Eq. 20 is a linear elastic cell problem with a known creep right-hand side. With positive-definite elastic constants, u̇# is unique up to a translation, which Eq. 9 removes (inference).
- **The symmetric state persists.** Eq. 11 makes σ̇ C-symmetric, so a symmetric state stays symmetric from step to step (inference; the paper does not state this induction).
- **Finite strain is different.** With a hyperelastic tangent, uniqueness can fail, and the argument then holds only on the symmetric branch (derived). This is NEMOpt's case, and it is already the reports' position.

### 3.3 Comparison with the open restatements

The restatements used in `open_access_mechanics_sources.md` are faithful: Matsuda et al. 2002, p. 538; KAKEN 09450046; Ohno 2000; Saavedra Flores 2011. The primary text adds:
- the explicit "Infinitesimal strain is assumed" (p. 2868), which settles that note's inference;
- the **constituents** as part of the assumption (Matsuda et al. 2002 mention only stress and strain);
- the vertex statement for parallelogram cells (p. 2870), which that note had derived;
- the penalty implementation, with only C₀ and C₁ zeroed explicitly in Fig. 4a;
- the half-versus-full validation (p. 2876) and the "about one eighth" cost remark (p. 2878);
- confirmation at the primary level that the paper does not discuss stability.

### 3.4 Comparison with FEAD 2026 (Eqs. 41–50) and NEMOpt

| Point | Ohno et al. 2001 | FEAD 2026, §3–4 | `NEMOpt-DR.jl` L242–299 |
|---|---|---|---|
| Kinematics | Infinitesimal strain, rate form (p. 2868) | Finite strain, total displacement, Newton | Same as FEAD |
| Symmetry required | Point symmetry of constituents, stress and strain about C₀ (assumed) | "for an orthotropic microstructure", with the symmetric stress state "deduced" (p. 8) | Quadrant-mirrored (orthotropic) density |
| Fluctuation parity | u̇#(P_k) = −u̇#(P) (Eq. 8) | ũ_m = −ũ*_m (Eq. 42) | D1 = −1 on top and cut-edge pairs |
| Zeros | C₀ and facet centers (Eq. 9); parallelogram vertices (p. 2870) | Center, "corner and mid-edge points" (p. 8) | Six prescribed nodes: two top corners, top midpoint C₁, two cut-edge ends, center C₀ |
| Half-cell BCs | C-symmetry on AB and CD; Y-periodicity on AD and BC (Fig. 2a) | TL ↔ TR and ML ↔ MR antisymmetric (Eqs. 43–44, 48–49); left ↔ right periodic | −1 on top and cut-edge pairs; +1 on side pairs |
| Constraint method | Penalty (p. 2875) | Master–slave elimination (Eq. 50) | Master–slave |
| Loading | Any uniform loading (p. 2870); off-axial examples | Any F_M, including load case 2 with shear | Same as FEAD |
| Stability | Not discussed | Not discussed for the reduction | NDR-0051 adds diagnostics |

Findings:
- **The constraint sets are equivalent (derived).** NEMOpt prescribes all six zero-fluctuation nodes because master–slave elimination cannot give one node two pair relations. Ohno's penalty pairs force the four extra zeros implicitly.
- **FEAD's "deduced" is stronger than Ohno's "assume".** FEAD §3 (p. 8): "When combined with point symmetry and Bloch type boundary conditions, for an orthotropic microstructure, it is deduced that the stress states at Y and Y* shall be identical." Ohno et al. assume this. The deduction needs a unique solution, which at finite strain holds only on the symmetric branch (derived).
- **The reduction itself does not need orthotropy (inference).**
  - Ohno's argument uses only 180° (C2) symmetry of the constituents.
  - FEAD itself notes that "not all point symmetric microstructures result in orthotropic response".
  - So the half-cell constraints in `NEMOpt-DR.jl` would also serve C2-only designs. The restriction to orthotropic designs comes from quadrant mirroring.
  - This bears on the article's conclusion that the method "is applicable to orthotropic micro-structures" (§6), not on any report claim.
- **The FEAD introduction credits refs. [18–20] for "(non-)linear homogenization".** Ohno et al. cover only material nonlinearity at infinitesimal strain. The finite-strain extension is the FEAD authors' own. The strategy report derives that a 180° rotation leaves every F_M invariant.
- **Half-versus-full agreement is how Ohno et al. validate their theory** (p. 2876), in a setting where the solution is unique. NDR-0051's point therefore stands: in NEMOpt's finite-strain setting, agreement is not evidence of stability.

### 3.5 Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| O-a | Strategy L36 | Parity rules: odd fluctuations, χ fields and adjoints; even densities | **Confirmed for fluctuations, χ (and φ) and densities; derived for adjoints** | Eq. 8; p. 2872 (χ and φ BCs); constituents assumed point-symmetric (abstract). Ohno et al. do no optimization. |
| O-b | Strategy L36 | A 180° rotation leaves every F_M invariant, so the reduction holds for shear | **Confirmed at infinitesimal strain; derived at finite strain** | "the transformation does not affect macroscopic loading" (p. 2870); Figs. 6–9 (θ = 22.5° and 45°; off-axial φ) |
| O-c | Strategy L36 | "Ohno et al. 2001 assume a point-symmetric internal state; the restriction to the symmetric branch is derived, not stated there" | **Confirmed (primary)** | Abstract; p. 2869; no stability content |
| O-d | Strategy L36 | "Full-versus-half agreement ends at the first symmetry-breaking bifurcation" | **Not from Ohno et al.** (derived; already amended at strategy L490) | Ohno et al. validate by agreement only where the solution is unique (p. 2876). |
| O-e | Addendum L15 | Ohno et al.: "half-cell reduction only; no stability content in open restatements" | **Confirmed; the primary text can now be cited** | As O5 |
| O-f | Addendum L72 | Six prescribed nodes are four C2-center classes counted with periodic images; −1 cut-edge pairs, +1 side pairs | **Consistent with Ohno et al.** | Eq. 9; vertex remark, p. 2870; Fig. 2a. Ohno et al. zero C₀ and C₁ explicitly and get the rest from the pair constraints (Fig. 4a; derived). |
| O-g | NDR-0051 Context (addendum L275–278) | Small-strain elastic–viscoplastic; assumes a point-symmetric state; fluctuations odd about the center and facet centers and zero there; no stability content | **Confirmed (primary)**; replace the secondary citations | p. 2868; p. 2869; Eqs. 8–9 and 17 |
| O-h | NDR table (addendum L314) | Origin: "Ohno et al. 2001; Palais 1979; derived" | **Acceptable** | Ohno et al. support the reduction and its BCs only. |
| O-i | `instability_theory_detection.md` L189 | Ohno et al. are "the basis for half-cell analysis of point-symmetric cells" | **Confirmed** | Title; §3 |
| O-j | `instability_theory_detection.md` L252 | No source analyzes the stability limits of the half cell | **Confirmed for the primary text** | — |
| O-k | `numerical_method_choices.md` L383 | Point symmetry allows half-cell analysis under arbitrary macroscopic loading "(not re-verified)" | **Confirmed, with conditions** | p. 2870. The conditions are infinitesimal strain and an assumed point-symmetric state. |
| O-l | `numerical_method_choices.md` L28 (table) | Valid for any F_M "as long as the solution stays on the symmetric branch" | **Consistent** | The condition is derived; Ohno et al. do not state it. |
| O-m | `reproduction_baseline.md` L244 | Thesis: zero fluctuation at edge midpoints and center | **Consistent but incomplete** | Ohno et al. also have zeros at parallelogram vertices (p. 2870). FEAD Eq. 50 prescribes the corners. |
| O-n | `open_access_mechanics_sources.md`, Item 1 | Points pending the primary text: infinitesimal strain (inference), corners (derived), no stability content ("unconfirmed at the primary level") | **All now confirmed** | p. 2868; p. 2870; full reading |

### 3.6 Recommended edits

**Addendum** (`reports/NEMOpt buckling and stability addendum.md`)

1. **Panel table, group-theoretic specialist row (L15).**
   - Old: "[Ohno, Matsuda & Wu 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) (half-cell reduction only; no stability content in open restatements);"
   - New: "[Ohno, Matsuda & Wu 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) (half-cell reduction only: it assumes a point-symmetric internal state at infinitesimal strain and contains no uniqueness, stability or bifurcation statement; read in full);"
2. **NDR-0051, Context (worked template, L275–278).** Keep the code-block line breaks.
   - Old:
     ```
     Origin: FEAD 2026 domain-reduction choice (Ohno et al. 2001, small-strain elastic-viscoplastic:
     assumes a point-symmetric internal state and derives fluctuations that are odd about the cell center
     and facet centers and vanish there; no stability or bifurcation content per Matsuda et al. 2002,
     JSME Int. J. A 45:538, and Saavedra Flores 2011, Sec. 4.3); W&T 2020 consider no stability.
     ```
   - New:
     ```
     Origin: FEAD 2026 domain-reduction choice (Ohno, Matsuda & Wu 2001, IJSS 38:2867, read in full:
     infinitesimal strain (p. 2868), linear elasticity plus nonlinear creep (Eq. 17); assumes that
     constituents, stress and strain are point-symmetric about the cell center (p. 2869) and shows that
     the perturbed velocity is then odd about the center and every facet center and vanishes there
     (Eqs. 8-9; also at parallelogram vertices, p. 2870); semiunit cell with C-symmetric cut and top
     edges and periodic sides (Fig. 2a); validated against full-cell results (p. 2876); no uniqueness,
     stability or bifurcation content). FEAD 2026 Sec. 3 calls the symmetric stress state "deduced"; at
     finite strain that holds only on the symmetric branch (derived). W&T 2020 consider no stability.
     ```
3. **Half-cell paragraph (L72), optional.**
   - After "([NEMOpt-DR.jl L242–299](../Domain_Reduction/NEMOpt-DR.jl#L242); verified)", insert: "; Ohno et al. impose the same pairs, zero the fluctuation at the center and the facet centers (Eqs. 8–9) and note that it also vanishes at the vertices of parallelogram cells (p. 2870) ([Ohno, Matsuda & Wu 2001](https://doi.org/10.1016/S0020-7683(00)00188-8))".
4. **NDR table (L314), Origin column, optional.**
   - Old: "Ohno et al. 2001; Palais 1979; derived"
   - New: "Ohno et al. 2001 (reduction and BCs only); Palais 1979; derived"
5. **Revisions log.** Add a bullet: "- Panel row (L15), NDR-0051 Context and NDR table: Ohno, Matsuda & Wu 2001 is now cited from the primary text (user-supplied): infinitesimal strain, an assumed point-symmetric state, Eqs. 8–9, no stability content. FEAD's "deduced" is flagged as valid only on the symmetric branch. Source: Lazarov–Kubo–Ohno note, §3.6."

**Strategy report** (`reports/NEMOpt refactoring strategy literature.md`)

6. **Panel table, mechanician row (L36).**
   - Old: "but only on the symmetric branch (Ohno et al. 2001 assume a point-symmetric internal state; the restriction to the symmetric branch is derived, not stated there)."
   - New: "but only on the symmetric branch (Ohno et al. 2001, at infinitesimal strain, assume that "constituents, stress and strain distribute point symmetrically" and note that the 180° transformation "does not affect macroscopic loading", p. 2870; the restriction to the symmetric branch is derived, not stated there)."
   - Optional, to match the amendment at L490: change "Full-versus-half agreement ends at the first symmetry-breaking bifurcation." to "The guarantee of full-versus-half agreement ends at the first symmetry-breaking bifurcation; agreement itself never certifies stability (addendum NDR-0051)."
7. **Revisions log.** Add a bullet: "- Panel table, mechanician row: Ohno et al. 2001 is now quoted from the primary text (the assumption, and invariance of the loading, p. 2870). Source: Lazarov–Kubo–Ohno note, §3.6."

**Research notes.** Append to the correction logs.

8. **`open_access_mechanics_sources.md`, Item 1.**
   - Append: "Update (same day; the user supplied the primary text): Ohno, Matsuda & Wu 2001 was read in full. Verdicts 1a–1c hold at the primary level. Infinitesimal strain (p. 2868). The point-symmetric state is assumed (abstract; p. 2869). Fluctuations are odd about C₀ and the facet centers and vanish there (Eqs. 8–9). They also vanish at parallelogram vertices (p. 2870), which replaces the derived corner argument. The paper has no uniqueness, stability or bifurcation content. 'Still needs the user' item 1 is resolved. See `lazarov_sigmund_2011_kubo_matsuda_2016.md`, §3."
   - In the source table, change "**No legal open copy found**" for Ohno, Matsuda & Wu 2001 to "**Supplied by the user** (version of record)".
9. **`instability_theory_detection.md`, log item for L252.** Add: "Confirmed at the primary level: Ohno, Matsuda & Wu 2001 makes no uniqueness, stability or bifurcation statement. Its half-cell validation is agreement with full-cell results (p. 2876), in an infinitesimal-strain rate setting where the rate solution is unique (inference)."
10. **`numerical_method_choices.md`, new log item for L383.**
    - Old: "[PR-bib] Ohno et al. show that point symmetry allows half-cell analysis under arbitrary macroscopic loading (not re-verified)."
    - New: "[PR] Ohno et al. assume point-symmetric constituents, stress and strain at infinitesimal strain, and show that the half cell then holds under any uniform macroscopic loading ("the transformation does not affect macroscopic loading", p. 2870). They verify it against full-cell results under transverse and off-axial loading (p. 2876)."
11. **`reproduction_baseline.md` L244, optional.** Append: "(Ohno et al. 2001, Eq. 9 and p. 2870: the fluctuation also vanishes at the corners of parallelogram cells; FEAD Eq. 50 prescribes them)".

---

## Access and working files

- **Read in full with the Read tool:** LS2011, KM2016 and OMW2001.
- **Read in part:** FEAD 2026 §1, §3–4 and §6, by text search of the user-supplied PDF. `NEMOpt-DR.jl` L242–299 and `NEMOpt.jl` L700–752.
- **Scratchpad only, not in the repository:**
  - text extractions used for absence searches;
  - a 200 dpi crop of OMW2001 Fig. 4, used to identify the zeroed nodes;
  - `moment_match.py`, used for §1.2.
- **Not edited:** the reports, the decision records and the other research notes.
