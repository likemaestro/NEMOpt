# Geymonat, Müller & Triantafyllidis (1993): full-text verification

**Paper.** G. Geymonat, S. Müller & N. Triantafyllidis, "Homogenization of nonlinearly elastic materials, microscopic bifurcation and macroscopic loss of rank-one convexity," *Arch. Ration. Mech. Anal.* 122(3):231–290 (1993), [doi:10.1007/BF00380256](https://doi.org/10.1007/BF00380256). Called "GMT" below.

**How it was read (2026-10-12).**
- Source: the user-supplied PDF, 60 pages, journal pp. 231–290. PDF page n is journal page 230 + n.
- All 60 pages were read as page images. The introduction (§1), Sections 2, 3 and 5, and the laminate results (§6.5–6.6) were read closely. The proofs (§4) and the algebra of §6.2–6.4 were skimmed.
- Pages 237 and 282 were re-rendered at 200 dpi to check formulas. The OCR text layer was used only to find passages; every quote below was checked against the page image.
- The issue number (3) was confirmed on Crossref.

**Bottom line.**
- The reports' core claims are correct, and the full text supports them, but they need hypotheses and theorem numbers added:
  - the inequality chain;
  - the long-wave ⇔ loss-of-ellipticity link;
  - micro before or with macro.
- "Homogenized moduli are meaningful only before the first bifurcation" is not GMT's wording. It should be restated as what they prove.
- The discontinuity at k = 0 is supported:
  - in general only as an inequality (the k → 0 limit ≤ the k = 0 value);
  - explicitly, as a "singular point", in their laminate example.
- NDR-0064's multi-cell sentence is supported by GMT §2, which restates Müller 1987, but its Müller link points to the wrong paper.
- New and useful for NEMOpt:
  - a reproducible compressible-laminate oracle, λ = 0.927, which I reproduced as 0.9267;
  - a theorem that the k-cell tangent equals the one-cell tangent while Bloch-stable;
  - the bound Λ₆ ≤ α_L, which explains why void modes cap every margin;
  - the finding that at a cell-periodic (k = 0) critical point 𝔸_M gives no warning.

---

## Notation map

| GMT symbol | Meaning in GMT | Reports / Zhang et al. | NEMOpt object |
|---|---|---|---|
| 𝕃(x); 𝕃_F₀(x) = ∂²W/∂F²(x, F₀ + ∇φ_F₀) (Eq. 5.4) | Periodic incremental (first-elasticity) tensor at the cell state | 𝔸 | Nominal tangent at the converged one-cell state |
| Λ (Eq. 3.11) = Λ₁ (Bloch, Eq. 4.4) = Λ₃ (all k-cells, Eq. 4.6), by Lemma 4.2 | Global coercivity: inf of ∫(𝕃∇v, ∇v)/∫\|∇v\|² | β(λ) | Bloch infimum over all θ (Tier 2 plus 0a plus 0b) |
| Λ₄ (Eq. 3.20) = Λ₅ (Eq. 3.21) | Shearing a⊗b + ∇q; equivalently liminf as ω → 0 of the Bloch quotient | β(λ; k → 0) | Long-wave limit (Tier 0b region) |
| Λ₆ (Eq. 3.22) | Coercivity on Y-periodic fields (ω = 0) | k = 0 value | K_red positive definiteness (Tier 0a) |
| α_𝕃 (Eq. 3.7) | Pointwise (constituent) ellipticity constant, ess inf over x | — | Weakest constituent: the ε-void |
| 𝕄 (Eqs. 3.16–3.18); 𝕄_F₀ (Eq. 5.6) | One-cell homogenized tensor, defined as an infimum | 𝔸̄ | 𝔸_M from the χ problems |
| α_𝕄 (Eq. 3.29) = min over \|a\| = \|b\| = 1, a, b ∈ ℝᴺ of (𝕄 a⊗b, a⊗b) | Macroscopic ellipticity constant | B(λ) | Tier 0b margin |
| W̄ (Eq. 2.6); Ŵᵏ (Eq. 2.7); Ŵ¹ | True homogenized energy; k-cell energy; one-cell energy | — | NEMOpt computes Ŵ¹ (and 𝔸_M from it) |

The reports' β and B are Zhang et al.'s notation. GMT never use "β", "Hill's second variation" or "Brillouin zone". Their ω ∈ [0, 2π[ᴺ is the Brillouin zone of the unit lattice.

---

## Q1. Claim-by-claim verdicts

Locations:
- "Add." is `reports/NEMOpt buckling and stability addendum.md`.
- "Ref." is `reports/NEMOpt refactoring strategy literature.md`.
- "ITD" is `research_notes/NEMOpt buckling and stability addendum/instability_theory_detection.md`.
- "SIV" is `.../stability_implementation_verification.md` in the same folder.

| # | Claim (location) | Verdict | Evidence in GMT (verbatim where quoted) |
|---|---|---|---|
| 1 | β(λ) ≤ β(λ; k → 0) ≤ B(λ), "proved" by GMT (Add. L56; ITD L29–31) | **Confirmed, with hypotheses** | Thm 3.3(i), p. 241: "Λ ≤ Λ₄ = Λ₅ ≤ Λ₆ ≤ α_𝕃". Cor. 3.5, p. 242: "Let Λ ≥ 0. Then α_𝕄 ≥ Λ₄ and α_𝕄 = 0 if Λ₄ = 0." Lemma 4.2 identifies Λ with the Bloch infimum. The right-hand inequality is proved only for Λ₄ ≥ 0 (Thm 3.3(v)), so beyond the first instability the chain is not asserted. The chain is a theorem about a fixed periodic linear tensor satisfying (3.3)–(3.6). It is checked numerically on p. 284: "the conditions Λ(λ) ≤ Λ₅ˢ(λ) ≤ α_𝕄ˢ(λ) are obviously satisfied". |
| 2 | Loss of rank-one convexity is "exactly" long-wavelength microscopic buckling (Add. L56; ITD L22, L33) | **Confirmed, with qualification** | **Linear level** (the tensor 𝕄 = one-cell 𝔸_M), at a state with Λ ≥ 0 and Λ₆ > 0: if Λ = Λ₄ = 0, 𝕄 "loses strict strong ellipticity" (Thm 3.4(iii)). If Λ = 0 and Λ₄ > 0, 𝕄 "remains strictly strongly elliptic" and a Bloch wave with ω ≠ 0 exists (Thm 3.4(ii)). With Cor. 3.5 this gives α_𝕄 = 0 ⇔ Λ₄ = 0. **Nonlinear W̄**: only long-wave ⇒ loss of uniform rank-one convexity is proved (Thm 5.1(ii); intro p. 234: "long wavelength instabilities for the linearized problem lead to loss of rank-one convexity of the homogenized energy density"). The converse is only "plausible" and proved only under (H1) (remark, pp. 263–264). If the critical mode is cell-periodic (Λ = Λ₆ = 0), the case "is less complete" (p. 243; Prop. 3.6), and 𝕄 can stay elliptic. |
| 3 | Microscopic instability precedes macroscopic instability or coincides with it (Add. L13, L56; ITD L22, L34) | **Confirmed** | Follows from claims 1 and 2: α_𝕄 ≥ Λ₄ ≥ Λ, and α_𝕄 = 0 forces Λ = 0. Laminate example, p. 283: "1 > λ_c⁻ ≥ λ̂_c⁻(0⁺) = λ_h⁻ = λ_H⁻" (first bifurcation no later than loss of ellipticity; the analogous tension inequality is also given). |
| 4 | Stability surface discontinuous at k = 0; cell-periodic mode attained at k = 0, long-wave mode only as k → 0 (Add. L58; ITD L35–38) | **Confirmed in substance; wording needs "can be"** | General result: only Λ₅ ≤ Λ₆ (Thm 3.3(i)). §3.3, p. 243: "Λ₆ measures coercivity with respect to Y-periodic, possibly highly localized deformations"; Λ₄ = Λ₅ "with respect to long-wavelength (ω → 0) perturbations". Explicit jump only in the laminate, §6.2 remark, p. 278: "ω₁ = 0 is a singular point in Λ̂(λ, ω₁)", because near ω₁ = 0 there are "two physically different types of modes". In compression, "Λ₅ˢ < Λ₆ˢ" in all their runs (p. 287). In tension, Λ = Λ₅ˢ = Λ₆ˢ, so there is no jump (p. 286). GMT do not say the surface is "continuous at the origin" when a k = 0 mode governs. That statement (ITD L36) is Zhang's, but it follows from Thm 3.3(i): if Λ = Λ₆, then Λ₅ = Λ₆ (inference). |
| 5 | Homogenized moduli are meaningful only before the first bifurcation of the principal solution (Add. L13; ITD L22, L27) | **Needs correction** (not GMT's wording; restate what is proved) | Proved: (a) Thm 3.3(iv): if Λ ≥ 0, the k-cell tensor 𝕄_k equals the one-cell 𝕄 for every k. (b) Prop. 3.2: if Λ < 0, then for fine enough microstructure (ε < ε₀) "I^ε is not bounded from below", and the linearized problem has no solution. (c) Prop. 5.3 and the remark on p. 263: under hypothesis (H2), one-cell minimizers equal k-cell minimizers, and the expansion with 𝕄_F "holds around F_s as long as Λ(𝕃_F_s) > 0". GMT call the zeros of Λ "the onset of the first bifurcation instability" (p. 281). At Λ = 0, Γ-convergence still holds (Thm 3.4(i)), so "up to" is more accurate than "before". "Only valid before the bifurcation occurs" is Zhang et al.'s paraphrase. |
| 6 | Trichotomy in linearized elasticity: energy unbounded below; strongly elliptic limit; rank-one-nonnegative but not strongly elliptic limit (not found in the reports or notes) | **Confirmed, with a refinement** | Λ < 0: unbounded below (Prop. 3.2). Λ > 0: usual homogenization limit, with (𝕄 a⊗b, a⊗b) ≥ Λ\|a\|²\|b\|² (Thm 3.1, Eq. 3.19). Λ = 0 (Γ-convergence, Thm 3.4(i)) splits in two: Λ₄ > 0 gives a strongly elliptic limit (Bloch-wave case); Λ₄ = 0 gives a rank-one-nonnegative limit with (𝕄 a⊗b, a⊗b) = 0 (shear-band case). Both need Λ₆ > 0. Intro, p. 234: the limit "retains or loses strong ellipticity" by the "character of that solution (Bloch wave or the long-wave-length limit)". The intro cites "Theorem 3.5", a misprint for Corollary 3.5. |
| 7 | GMT give "a rigorous proof that a Bloch-wave representation can be used for the microscale stability of general 3D periodic solids" (ITD L25) | **Needs correction** (dimension) | Lemma 4.2 holds for any N ≥ 2 (Y = [0, 1]ᴺ): Λ (compactly supported test fields on ℝᴺ) = Λ₁ (Bloch waves e^{iω·x}p(x)) = Λ₃ (all k-cell periodic fields). Intro, p. 234: Λ "can equivalently be characterized in terms of Bloch waves". It concerns the linearized operator's coercivity. |
| 8 | If β > 0, the homogenized material is rank-one convex (ITD L32) | **Confirmed (stronger)** | Thm 3.1, Eq. 3.19: Λ > 0 ⇒ (𝕄 a⊗b, a⊗b) ≥ Λ\|a⊗b\|², i.e. strict strong ellipticity with constant Λ. |
| 9 | B(λ) = 0 if β(λ; k → 0) = 0 (Zhang Eq. 56; ITD L33) | **Confirmed** | Cor. 3.5. It is derived from Thm 3.4, whose hypotheses include Λ₆ > 0. |
| 10 | Passing the 𝔸_M rank-one check is necessary but not sufficient (ITD "Ordering in the optimization") | **Confirmed** | Thm 3.4(ii): a Bloch instability with ω ≠ 0 can occur while 𝕄 stays strictly strongly elliptic. Laminate, p. 288: for stiff fraction < 0.15, the first instability is "fiber buckling with ω₁c > 0", before loss of ellipticity. |
| 11 | Elliptic or polyconvex constituents do not guarantee a rank-one-convex homogenized material (ITD L41, via Zhang on Abeyaratne & Triantafyllidis 1984) | **Confirmed**, with a second example | Intro, p. 233, on [AT 84]: the matrix is "polyconvex … and hence always rank-one convex", yet the "homogenized incremental (or linearized) moduli … lose their rank-one convexity at adequately high macroscopic strains". GMT's own laminate, p. 288: with constituents strongly elliptic "at any level of deformation", the composite can "lose ellipticity at fairly low strains, e.g., λ = 0.927" for stiff fraction 0.5. |
| 12 | Porous finding attributed to "Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134" (SIV L239) | **Needs correction** (wrong paper) | GMT attribute the porous homogenization result to [AT 84] = J. Appl. Mech. 51:481–486 (1984), doi:10.1115/1.3167661. [AT 81], "The emergence of shear bands in plane strain", is listed in GMT's references but is not cited for this point. ITD L41 already cites 1984 correctly. |
| 13 | Framing: macro instability = loss of rank-one convexity of the homogenized incremental moduli (ITD L307) | **Confirmed** | Intro, pp. 233–234: macroscopic instabilities are "the failure of strong ellipticity for the homogenized incremental moduli". It is not only the title. |
| 14 | Microscopic bifurcation and macroscopic loss of rank-one convexity are the two stability notions; the long-wave limit connects to rank-one convexity (SIV L166) | **Confirmed** | Thm 3.3(i), Thm 3.4, Cor. 3.5, Thm 5.1(ii). |
| 15 | "Three tests together are complete in the sense of GMT"; "the region near θ = 0 is covered by the ellipticity test" (SIV L157, L263) | **Needs correction** (qualify) | Necessity is confirmed: Λ > 0 requires Λ₆ > 0, Bloch positivity at every ω, and α_𝕄 ≥ Λ₄ > 0. "Covered" holds only in the limit sense of Cor. 3.5 (α_𝕄 = 0 ⇔ Λ₄ = Λ₅ = 0 at Λ ≥ 0), along a path that starts stable. GMT give no radius below which α_𝕄 > 0 guarantees Bloch positivity. The margins are not comparable either: Λ is normalized by ∫\|∇v\|², α_𝕄 by \|a⊗b\|². |
| 16 | "Geymonat et al. 1993 give general theory and a laminate example; their numerical values were not checked" (SIV L379) | **Now checked** | p. 288: λ = 0.927 at stiff fraction 0.5. Reproduced as 0.9267 from GMT Eqs. 6.6, 6.55, 6.56 and 6.42 (see "Laminate reproduction" below). |
| 17 | NDR-0064: "Once instabilities appear, the homogenized energy may need a multi-cell infimum (Müller 1987 …, cited from a reference list)" (Ref. L503) | **Confirmed in substance; the link is wrong** | See Q2a. The link `10.1007/s00205-016-0992-0` is Duerinckx & Gloria 2016, "Stochastic homogenization of nonconvex unbounded integral functionals with convex growth" (Crossref). Müller 1987 is [10.1007/BF00284506](https://doi.org/10.1007/BF00284506) (Crossref: ARMA 99:189–212, matching GMT's [Mü 87]). |
| 18 | "Periodic homogenization defines the effective energy as a minimum over periodic fluctuations" (Ref. L511; not attributed to GMT) | **Needs correction** | GMT Eqs. 2.6–2.7: W̄(F) = inf over k ∈ ℕ of Ŵᵏ(F), an infimum over fluctuations periodic on kY. It equals the one-cell minimum for convex W ("Marcellini [Ma 78]", p. 237). The Voigt-type bound survives, because zero fluctuation is admissible for every k: W̄ ≤ Ŵ¹ ≤ cell average of W(x, F̄). |
| 19 | Ponte Castañeda 1989 (Proc. R. Soc. A 422:147) bounds the effective energy (Ref. L501–502; "verified only through secondary sources") | **Corroborated** (another secondary source) | GMT, p. 237: "see Ponte-Castaneda [PC 89] for some non-trivial bounds on W̄". Reference list: Proc. Roy. Soc. London 422A (1989) 147–171. |
| 20 | Newton converges to unstable equilibria; track all checks at every increment (Add. L13, L58) | **Not in the paper** (labeled inference in the report, correctly) | Related only: the remark on p. 263 propagates one-cell = k-cell along a load path "as long as Λ > 0". On p. 259, GMT note that excluding discontinuous bifurcation "is implicit in many numerical schemes to compute minimizers (e.g., by path-following techniques)". |
| 21 | "Hill's second variation … at every wave vector of the Brillouin zone" (Add. L56) | **Consistent**; not GMT's wording | GMT work with the coercivity constant Λ over compactly supported fields on ℝᴺ and prove it equals the Bloch infimum over ω ∈ [0, 2π[ᴺ (Lemma 4.2). They cite Hill only for linear homogenization ([H 63], [H 65]), not his 1957 criterion. |

---

## Q2a. When does the one-cell energy equal the true (multi-cell) homogenized energy?

**Yes, this is Müller 1987.** GMT §2 is explicitly a review of "[Mü 87] (cf. also [Br 85])", where [Mü 87] = S. Müller, "Homogenization of nonconvex integral functionals and cellular elastic materials," ARMA 99 (1987) 189–212, and [Br 85] = A. Braides, Rend. Accad. Naz. XL 9 (1985) 313–322.

What GMT state:
1. **General nonconvex formula** (Thm 2.3, "cf. [Mü 87, Thm. 1.3, Cor. 2.3]"; Eqs. 2.6–2.7).
   - W̄(F) = inf over k of Ŵᵏ(F), where Ŵᵏ is the energy per cell minimized over fluctuations periodic on k×…×k supercells.
   - Hypotheses: p-growth (2.1)–(2.3) and a bounded Lipschitz Ω. "Theorem 2.3 does not require any convexity assumption on W" (p. 237).
2. **Convex case.** Marcellini [Ma 78] proved that for convex W, W̄ = Ŵ¹ (p. 237; p. 265 for strictly convex W with quadratic growth).
3. **Strict inequality can occur.** "[Mü 87, Thm. 4.3] gives an example where W̄ < Ŵ" (p. 237; in context Ŵ = Ŵ¹).
   - The intro (p. 233) ties the convex/nonconvex difference to "possible instabilities such as buckling on the microscale (see Sections 1 and 4 in [Mü 87])".
4. **Linearized level.** Thm 3.3(iv): if Λ ≥ 0, then (𝕄_k F, F) = (𝕄 F, F) for all k.
   - The multi-cell incremental tensor equals the one-cell tensor while the cell is Bloch-stable.
   - If Λ < 0, the linearized energy is unbounded below (Prop. 3.2).
5. **Nonlinear level near a state F₀.**
   - Hypothesis (H1) *assumes* W̄ = Ŵ¹ near F₀, with W^{1,∞}-continuous minimizers. GMT: "Equation (5.1) states that near F₀ the homogenized energy is given by the solution of an auxiliary problem on only one periodic cell, while (5.3) requires that no discontinuous bifurcation of minimizers occurs" (p. 259).
   - Prop. 5.3 gives a sufficient condition. Suppose (H2) holds: W̄(F₀) = Ŵ¹(F₀), and k-cell minimizers exist and are W^{1,∞}-continuous for every k. Suppose also Λ(F₀) > 0. Then, for small t, every k-cell minimizer equals the one-cell minimizer, so W̄(F₀ + tH) = Ŵ¹(F₀ + tH) and (H1) holds.
   - The remark on p. 263 extends this along a rank-one load path F_s "as long as Λ(𝕃_F_s) > 0".

**Verdict on NDR-0064's sentence.** It is supported. GMT is a better, now verified, pinpoint than "cited from a reference list". GMT also sharpen it: the one-cell value is exact for convex energies. Starting from a state where it is exact, it stays exact along the load path while the cell is stable against all Bloch waves, under the continuity hypothesis (H2). After the first Bloch instability, nothing guarantees it. GMT do not prove that the one-cell value then fails. They show only that it can (Müller's example).

## Q2b. Hypotheses GMT require

| Part | Hypotheses |
|---|---|
| §2, Thm 2.3 (nonlinear Γ-limit) | W(x, F) Y-periodic, Y = [0, 1]ᴺ, N ≥ 2. Growth and coercivity: c\|F\|ᵖ ≤ W ≤ C(1 + \|F\|ᵖ), \|∂W/∂F\| ≤ C(1 + \|F\|^{p−1}), p > 1 (Eqs. 2.1–2.3). Ω bounded Lipschitz. **No convexity.** Throughout the paper W(x, ·) is also assumed quasiconvex and pointwise strongly elliptic (p. 237), to isolate instabilities due to "the way in which the different base materials … are mixed". |
| §3 (linearized) | 𝕃 measurable and Y-periodic (3.3); major symmetry only, (𝕃A, B) = (A, 𝕃B) (3.4); bounded (3.5); pointwise strong ellipticity with a uniform c > 0 (3.6). Dirichlet data on a bounded Lipschitz Ω. Convergence of solutions needs Λ > 0 (Thm 3.1). Γ-convergence needs Λ ≥ 0 and Λ₆ > 0 (Thm 3.4). **No minor symmetry**, consistent with using the full nominal 𝔸_M. |
| §5.2, Thm 5.1 | (2.1)–(2.3), (H1), W(x, ·) ∈ C³ with locally bounded third derivatives (5.9–5.10). (i) needs Λ₄(F₀) > 0. (ii) needs Λ₄(F₀) = 0 and Λ₆(F₀) > 0; (ii) also holds with only W̄(F₀) = Ŵ¹(F₀). |
| §5.3, Prop. 5.3 | (2.1)–(2.3), (5.9), (5.10), (H2) and Λ(F₀) > 0. |
| §5.4, Thm 5.4 (homogenization and linearization commute) | Strictly (uniformly) convex C² W with quadratic growth (5.28–5.32), described as "an unrealistic assumption in elasticity" (p. 259). |
| §6 (laminate) | Plane strain; isotropic layers; perfect bonding; lateral stress σ = 0; layers strongly elliptic at every strain (assumption (E)); Λ₁ˢ(λ) ≥ 0. Energy (6.52) = μ/2(v₁² + v₂² − 2 − 2 ln v₁v₂) + κ/2(v₁v₂ − 1)², with κ = νE/((1 + ν)(1 − 2ν)) (Lamé's λ, although GMT call it the bulk modulus). |

- **Small loads.** Not assumed globally. Thm 5.1 and Prop. 5.3 are local (t → 0, or t ≤ δ₀) around any F₀ satisfying the hypotheses.
- **Growth caveat.** GMT concede that energies blowing up as det F → 0 violate (2.2)–(2.3). For (6.52): "not compatible with the growth assumptions" (p. 284). They argue it is harmless because det F stays bounded away from 0. NEMOpt's −μ ln J + λ/2 (ln J)² has the same issue. It affects only the nonlinear results (§2, §5), not the linearized theorems (§3) applied to the tangent at a given state.

## Q2c. Rigorous versus formal or conjectural

- **Rigorous theorems:**
  - Thm 2.3 (quoted from Müller and Braides);
  - Thm 3.1, Prop. 3.2, Thm 3.3, Thm 3.4, Cor. 3.5, Prop. 3.6, Thm 3.7;
  - Lemmas 4.1 and 4.2;
  - Thm 5.1 (under (H1)), Prop. 5.2, Prop. 5.3 (under (H2)), Thm 5.4, Props. 5.5 and 5.6;
  - Props. 6.1 and 6.2.
- **Hypotheses with no proof available:** (H1) and (H2). GMT: "There are, however, no rigorous mathematical results which would allow one to deduce (5.3) from reasonable assumptions on W" (p. 259). They also do not prove that W̄ is twice differentiable; they work with second-order expansions along rank-one lines.
- **Formal or conjectural statements:**
  - "A formal calculation suggests" that ∂²W̄/∂F² is given by the linearized homogenization (p. 258). This is proved only under (H1) along rank-one directions, or for convex W.
  - That W̄ stays rank-one convex at a Bloch-wave (ω ≠ 0) instability is "plausible" (p. 263).
  - The case Λ = Λ₆ = 0 "is less complete" (p. 243).
- **Numerical or semi-analytic (§6):**
  - Λ is found by "a numerical search" over ω₁ (intervals (0, 7) and (0, 24)).
  - Existence of the λ_c(ω₁) branches is "tacitly assumed" (p. 281).
  - Simple roots are assumed in the proof of Prop. 6.1.
  - The lim/inf interchange (6.47) is "verified in all the numerical calculations" (p. 282), not proved.
  - The "≈ 0.15" stiff-fraction threshold and the critical stretches are numerical results for one material set.

## Q2d. Anything specific to 2D, plane strain or symmetric cells?

- **Dimension.** The general theory (§§2–5) holds in any dimension N ≥ 2, on the unit cube cell. Only §6 is 2D plane strain. Šverák's counterexample (rank-one convexity ⇏ quasiconvexity) is cited for N ≥ 3 only (p. 237), so it does not settle the 2D question.
- **Cell symmetry.** The general theory needs none. In §6:
  - cells "symmetric about x₁" (Figs. 6.2–6.3) and asymmetric cells (Figs. 6.4–6.5) give "very similar" results (p. 286);
  - the critical Bloch modes are always symmetric or antisymmetric across the cell in x₂ (ω₂ = 0 or π), "with the symmetric case occurring most frequently" (pp. 278, 286);
  - Λ̂(λ, ω₁) = Λ̂(λ, −ω₁) (Eq. 6.19), which is the time-reversal halving the addendum uses.
- **Point symmetry.** GMT say nothing about C2 point symmetry, half cells or symmetric criticality. The addendum's half-cell conclusions rest only on its other sources.

---

## Q3. New facts relevant to NEMOpt's planned stability checks

1. **k = 0 inertia check (Tier 0a).**
   - Λ₆ ≥ Λ₄ = Λ₅ (Thm 3.3(i)), so a positive-definite K_red never bounds the long-wave constant. Tier 0b is indispensable. This is now verified rather than relayed.
   - Conversely, at a cell-periodic critical point with Λ ≥ 0, every null mode q satisfies ∫Re(𝕃 a⊗b, ∇q) = 0 (Eq. 4.15). The right-hand sides of the χ problems are therefore orthogonal to the null mode. 𝕄 stays finite, can stay strongly elliptic (Prop. 3.6), and gives no warning. Only the inertia check sees this instability (inference from Eq. 4.15 and Prop. 3.6).
2. **𝔸_M is a variational object only while K_red is positive definite.**
   - GMT define 𝕄 as an infimum (Eq. 3.18). It is attained when Λ₆ > 0 (Thm 3.3(vi)).
   - If Λ₆ < 0, that infimum is −∞. NEMOpt's χ solve would then return a stationary value with no homogenization meaning (inference). This supports reporting no homogenized property beyond t_c (NDR-0045).
3. **Every margin is capped by the weakest constituent.** Λ ≤ Λ₆ ≤ α_𝕃 = ess inf over x of the pointwise ellipticity (Thm 3.3(i); proved with oscillating test fields concentrated at one point, §4.2(c)).
   - For a design with ε-voids, every GMT-normalized margin is at most about the void's ellipticity, whatever the metric that does not down-weight voids.
   - This is the continuum reason behind the addendum's observations: void modes of order ε dominate, raw margins are meaningless, and the shift τ must stay below ε (inference).
4. **Supercells add nothing while stable.**
   - Thm 3.3(iv): while Λ ≥ 0, the k-cell incremental tensor equals the one-cell 𝔸_M, so homogenizing on a supercell is redundant.
   - Lemma 4.2: Λ is also the infimum over all k-cell periodic fields. Proof part (c) decomposes k-periodic fields into Bloch waves at ω ∈ (2π/k){0, …, k − 1}ᴺ. The 2×2 supercell is therefore exactly the four points θ ∈ {0, π}² (the TRIM points).
   - Any finite set of supercells or θ points gives only an upper bound on Λ. This supports labeling the TRIM and contour scans "not guaranteed".
5. **Macroscopic ellipticity check (Tier 0b, NDR-0049).**
   - GMT's α_𝕄 (Eq. 3.29) uses real unit a, b and needs only major symmetry. That is consistent with testing the full 4×4 𝔸_M without minor-symmetry reduction.
   - For orthotropic 𝔸_M (laminates, or any cell with diagonal F_M and mirror symmetry), Eq. 6.42 gives a closed-form scan over φ ∈ [0, π/2], and Eq. 6.49 a root test. Both are useful as an independent check of the quartic test.
   - Eq. 6.49 as printed has misprints. The squared term should read (M₁₁₂₂ + M₁₂₂₁)², not (M₁₁₂₂ + M₂₂₂₂)². The inf and sup ranges for λ_H± are swapped relative to Eq. 6.48. The correct quartic (derived here) is det Q = M₁₁₁₁M₂₁₂₁b₁⁴ + [M₁₁₁₁M₂₂₂₂ + M₁₂₁₂M₂₁₂₁ − (M₁₁₂₂ + M₁₂₂₁)²]b₁²b₂² + M₁₂₁₂M₂₂₂₂b₂⁴.
6. **Bloch scans (Tier 2) on laminates.**
   - Critical ω₂ (across the layers) was always 0 or π.
   - The critical ω₁ (along the layers) varies continuously and grows as the stiff layers thin (p. 288). The searches used ω₁ ∈ (0, 24), that is, wavelengths down to about 0.26 cell heights.
   - B-S1 therefore needs a scan along θ₁ with θ₂ ∈ {0, π}, not TRIM points alone, and a mesh that resolves stiff-layer buckling when the stiff fraction is below about 0.15.
   - Coincidence of the micro and macro onsets: for stiff fraction above about 0.15, λ_c⁻ = λ_H⁻ (first instability long-wave, ω₁ → 0⁺). Below about 0.15, stiff-layer buckling comes first (pp. 287–288). This partly answers the addendum's open item for one compressible material pair.
7. **Compressible laminate oracle (computed, see below).**
   - GMT's homogenized laminate tangent (Eq. 6.6) gives B-S1 a closed form along the σ₂₂ = 0 path.
   - It reproduces GMT's 0.927 and the addendum's incompressible-limit values (0.9045 and 0.8668). It also supplies ν = 0.3 values that the incompressible formula cannot.
8. **Constituent ellipticity.** GMT assume pointwise strong ellipticity exactly to isolate structural instabilities. They note that in practice "all the important instability phenomena occur well below" loss of ellipticity of the weakest layer (p. 283). This supports the addendum's "every instability is structural" framing.
9. **Continuum versus discrete.** All GMT results are continuum statements. The discrete Q4 analogues (for example Λ_h ≤ Λ₄,h ≤ α_𝕄,h) are not covered by their proofs. The Bloch product e^{iω·x}q is not in the FE space, so the discrete chain must be checked numerically under refinement (inference; consistent with the open item in SIV).

---

## Laminate reproduction (computed; inference)

**Method.** For each layer:
1. Solve σ₂₂ = 0 for λ₂ at stretch λ along the layers (Eq. 6.55 for GMT's law; bisection for the Bonet–Wood law).
2. Form the incremental moduli (Eq. 6.56, or the code's tangent with ln J).
3. Homogenize with Eq. 6.6.
4. Scan the acoustic tensor's minimum eigenvalue over φ ∈ [0, π/2] (1,800 steps), then bisect for the first zero below λ = 1.

A scratch script (not kept in the repository) did this. Sanity check: the moduli at λ = 1 equal isotropic linear elasticity (C₁₁₁₁ = λ + 2μ, C₁₁₂₂ = λ, C₁₂₁₂ = C₁₂₂₁ = μ).

**GMT's material pair** (soft E 0.1, ν 0.33; stiff E 1.5, ν 0.48), loss of ellipticity λ_H⁻:

| Stiff fraction H_b | GMT law (Eq. 6.52) | Same, with NEMOpt's (ln J)² volumetric term |
|---|---|---|
| 0.1 | 0.8410 | 0.8409 |
| 0.2 | 0.8951 | 0.8950 |
| 0.3 | 0.9154 | 0.9154 |
| 0.5 | **0.9267** (paper: "λ = 0.927", p. 288) | 0.9267 |
| 0.7 | 0.9154 | 0.9154 |

The curve's shape matches the dotted λ_H⁻ curve in Fig. 6.6 by eye.

**NEMOpt-like pair** (equal ν in both phases, as in the pure-SIMP B-S1 mode; Young's modulus ratio 10), Bonet–Wood law:

| ν | vf 0.5 | vf 0.2 |
|---|---|---|
| 0.3 (FEAD2026) | 0.8677 | 0.8162 |
| 0.45 | 0.8947 | 0.8530 |
| 0.49 | 0.9025 | 0.8639 |
| 0.499 | 0.9043 | 0.8665 |
| Incompressible formula (1 − μ̆/μ̄)^{1/4} | 0.9045 | 0.8668 |

The convergence to the incompressible formula independently supports the addendum's harmonic-mean reading of Slesarenko & Rudykh Eq. 30.

**Caveats.**
- These values hold on the σ₂₂ = 0 path. NEMOpt prescribes F_M, so B-S1 must prescribe F_M = diag(λ, Σ_p H_p λ₂,p(λ)) to reach the same states.
- They are loss-of-ellipticity values. The Bloch onset equals them only where GMT found coincidence (their material, stiff fraction above about 0.15). For the NEMOpt pair this is an expectation, not a result.

---

## Q4. Recommended edits (do not apply here; old text → new text)

### `reports/NEMOpt buckling and stability addendum.md`

**E1 (L13, panel row, claim).**
- Old: "Homogenized moduli are meaningful only before the first bifurcation of the principal solution."
- New: "One-cell homogenized moduli are certified only up to the first bifurcation of the principal solution: while the cell is stable against Bloch waves of every wave vector, the one-cell tangent equals every multi-cell tangent, and beyond that point the linearized energy of a body with enough cells is unbounded below."

**E2 (L13, source cell).**
- Old: "[Geymonat, Müller & Triantafyllidis 1993](https://doi.org/10.1007/BF00380256) (content relayed by [Zhang, Feng & Khandelwal 2021](https://doi.org/10.1002/nme.6802), preprint [arXiv:2010.02371](https://arxiv.org/abs/2010.02371))"
- New: "[Geymonat, Müller & Triantafyllidis 1993, Thm 3.3, Cor. 3.5, Prop. 3.2, Prop. 5.3](https://doi.org/10.1007/BF00380256) (full text verified)"

**E3 (L56).**
- Old: "Geymonat, Müller and Triantafyllidis proved that the Bloch infimum obeys **β(λ) ≤ β(λ; k → 0) ≤ B(λ)**, where B = min over unit m and N of (m⊗N):𝔸̄:(m⊗N) is the macroscopic ellipticity indicator, and that loss of rank-one convexity is exactly long-wavelength microscopic buckling. Their 1993 full text could not be read, so these statements are **verified as relayed** by Zhang et al.'s restatement ([Geymonat et al. 1993](https://doi.org/10.1007/BF00380256); [Zhang et al. 2021, Eqs. 55–56](https://arxiv.org/abs/2010.02371))."
- New: "For the incremental tensor of a periodic cell with strongly elliptic constituents, Geymonat, Müller and Triantafyllidis proved that, while the cell is microscopically stable (β ≥ 0), the Bloch infimum obeys **β(λ) ≤ β(λ; k → 0) ≤ B(λ)**, where B = min over unit m and N of (m⊗N):𝔸̄:(m⊗N) is the macroscopic ellipticity indicator; in their notation this is Λ ≤ Λ₄ = Λ₅ ≤ α_𝕄 (Thm 3.3(i) and Cor. 3.5, with Λ the Bloch infimum by Lemma 4.2). At the first microscopic instability, with the cell-periodic constant still positive, 𝔸̄ loses strong ellipticity exactly when the critical mode is long-wavelength (Thm 3.4(ii)–(iii)); for the nonlinear homogenized energy only the implication from a long-wave mode to loss of rank-one convexity is proved (Thm 5.1(ii)). These statements are **verified** against the full text ([Geymonat et al. 1993, §§3–5](https://doi.org/10.1007/BF00380256); restated in [Zhang et al. 2021, Eqs. 55–56](https://arxiv.org/abs/2010.02371))."

**E4 (L56, next sentence).**
- Old: "Microscopic instability therefore precedes macroscopic instability or coincides with it."
- New: "Microscopic instability therefore precedes macroscopic instability or coincides with it; Geymonat et al.'s laminate shows both cases ([Geymonat et al. 1993, p. 283 and §6.6](https://doi.org/10.1007/BF00380256))."

**E5 (L58).**
- Old: "and the stability surface is discontinuous at the origin; a computation at k = 0 cannot see it, so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); verified)."
- New: "and the stability surface can be discontinuous at the origin, because the long-wave limit can lie strictly below the k = 0 value; a computation at k = 0 cannot see it, so the long-wave region belongs to the 𝔸_M test ([Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371); [Geymonat et al. 1993, Thm 3.3(i) and the remark on p. 278](https://doi.org/10.1007/BF00380256), where the jump appears at every compressive state of their laminate; verified)."

**E6 (L46, row 7, resolution cell).**
- Old: "Whether micro and macro onsets coincide for compressible laminates is an open item, not an acceptance criterion."
- New: "Coincidence is an expectation, not an acceptance criterion: for one compressible plane-strain laminate (Young's modulus ratio 15), micro and macro onsets coincide when the stiff-layer fraction exceeds about 0.15, and finite-wavelength buckling of the stiff layers comes first below it ([Geymonat et al. 1993, §6.6](https://doi.org/10.1007/BF00380256))."
- Also append to that row's source cell: "; [Geymonat et al. 1993, Thm 3.3(i), Cor. 3.5, §6.6](https://doi.org/10.1007/BF00380256)".

**E7 (L178, B-S1, expectation and source cells).**
- Old: "Bloch onset ≤ macroscopic onset | [Slesarenko & Rudykh 2017, Eqs. 29–30](https://doi.org/10.1016/j.jmps.2016.11.002) (harmonic-mean reading of a garbled Eq. 30 is inference); [Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148) |"
- New: "Bloch onset ≤ macroscopic onset; compressible oracle on the σ₂₂ = 0 path (F_M = diag(λ, mean lateral stretch)): loss of ellipticity at λ = 0.8677 (volume fraction 0.5) and 0.8162 (0.2) for ν = 0.3 and stiffness ratio 10 (computed from Geymonat et al. Eq. 6.6), and at λ = 0.927 for Geymonat et al.'s material pair at stiff fraction 0.5 (their p. 288; reproduced as 0.9267) | [Slesarenko & Rudykh 2017, Eqs. 29–30](https://doi.org/10.1016/j.jmps.2016.11.002) (harmonic-mean reading of a garbled Eq. 30 is inference, cross-checked: the exact compressible laminate tangent gives 0.9043 and 0.8665 at ν = 0.499); [Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148); [Geymonat et al. 1993, Thm 3.3(i), Cor. 3.5, Eqs. 6.6, 6.55, 6.56](https://doi.org/10.1007/BF00380256) |"

**E8 (L301, NDR-0045).**
- Old: "New capability: Geymonat et al. 1993; Zhang et al. 2021"
- New: "New capability: Geymonat et al. 1993 (Thm 3.3(iv), Prop. 3.2, Prop. 5.3); Zhang et al. 2021"

**E9 (L305, NDR-0049).**
- Old: "Literature: Geymonat et al. 1993; Zhang et al. 2021; closed form derived"
- New: "Literature: Geymonat et al. 1993 (Eq. 3.29, Cor. 3.5; orthotropic closed form Eq. 6.42); Zhang et al. 2021; closed form derived"

### `reports/NEMOpt refactoring strategy literature.md` (NDR-0064)

**E10 (L503).**
- Old: "- They require convexity that finite-strain hyperelastic energies lack. Once instabilities appear, the homogenized energy may need a multi-cell infimum ([Müller 1987, ARMA 99:189](https://link.springer.com/doi/10.1007/s00205-016-0992-0), cited from a reference list)."
- New: "- They require convexity that finite-strain hyperelastic energies lack. Without convexity the homogenized energy is an infimum over fluctuations periodic on k×k supercells for all k, which can lie strictly below the one-cell value ([Müller 1987, ARMA 99:189–212](https://doi.org/10.1007/BF00284506), as restated in [Geymonat et al. 1993, Thm 2.3 and Eqs. 2.6–2.7](https://doi.org/10.1007/BF00380256); verified there). The one-cell value is exact for convex energies, and from a state where it is exact it stays exact while the cell is stable against Bloch waves of every wave vector, under a continuity hypothesis on the cell minimizers ([Geymonat et al. 1993, Prop. 5.3](https://doi.org/10.1007/BF00380256))."

**E11 (L511).**
- Old: "Periodic homogenization defines the effective energy as a minimum over periodic fluctuations. Taking zero fluctuation as the trial field gives W̄(F̄) ≤ (1/|Ω|)∫ W_ρ(x, F̄) dV (standard derivation). The bound holds only if the computed state is the energy minimizer,"
- New: "Periodic homogenization defines the effective energy as an infimum over fluctuations periodic on k×k supercells for all k, and NEMOpt computes the one-cell (k = 1) value, which bounds it from above ([Geymonat et al. 1993, Eqs. 2.6–2.7](https://doi.org/10.1007/BF00380256)). Zero fluctuation is admissible for every k, so W̄(F̄) ≤ Ŵ¹(F̄) ≤ (1/|Ω|)∫ W_ρ(x, F̄) dV (standard derivation). The computed one-cell energy obeys the bound only if the computed state is the one-cell energy minimizer,"

### Research notes (for consistency)

**E12 (ITD L24).**
- Old: "The full text was not accessible this session (the University of Michigan Deep Blue repository copy sits behind a Cloudflare challenge). Its content is therefore taken from the detailed secondary summary in Zhang, Feng & Khandelwal 2021 (next items)."
- New: "The full text was later read from a user-supplied copy; the items below are confirmed with the qualifications in [the verification note](../NEMOpt%20supplied%20papers%20verification/geymonat_muller_triantafyllidis_1993.md)."

**E13 (ITD L25).**
- Old: "for the microscale stability of general 3D periodic solids of infinite extent"
- New: "for the microscale stability of periodic solids of infinite extent in any dimension N ≥ 2 (Lemma 4.2)"

**E14 (ITD L36).**
- Old: "and the stability surface β(k) is continuous at the origin."
- New: "and the stability surface β(k) is continuous at the origin (Zhang et al.; not stated by Geymonat et al., but it follows from their Thm 3.3(i): if Λ = Λ₆ then Λ₅ = Λ₆)."

**E15 (ITD L44).** Append to the Müller 1987 item: " — formula verified as restated in Geymonat et al. 1993, Thm 2.3 and Eqs. 2.6–2.7 (Müller's own text not read)."

**E16 (ITD L64).**
- Old: "The Geymonat et al. (1993) full text was not read. The precise hypotheses (regularity, convexity of constituents, and the exact relation between k = 0 and k → 0 infima) are relayed via Zhang et al. (2021)."
- New: "Resolved: the Geymonat et al. (1993) full text was read; for hypotheses and the k = 0 versus k → 0 relation (Λ₅ ≤ Λ₆) see the verification note."

**E17 (SIV L157).**
- Old: "Three tests together are complete in the sense of Geymonat–Müller–Triantafyllidis:"
- New: "Three tests together cover every mode class that Geymonat–Müller–Triantafyllidis distinguish; near θ = 0 the coverage holds only in the limit sense of their Cor. 3.5 and along a load path that starts stable:"

**E18 (SIV L239).**
- Old: "[Abeyaratne & Triantafyllidis 1981, IJSS 17:1113–1134](https://doi.org/10.1016/0020-7683(81)90092-5), as cited in Zhang et al."
- New: "[Abeyaratne & Triantafyllidis 1984, J. Appl. Mech. 51(3):481–486](https://doi.org/10.1115/1.3167661), as cited in Zhang et al. and in Geymonat et al. 1993, p. 233"

**E19 (SIV L263).**
- Old: "This follows the Geymonat et al. definitions; the region near θ = 0 is covered by the ellipticity test."
- New: "These conditions are necessary by Geymonat et al.'s definitions. The ellipticity test covers the region near θ = 0 only in the sense that, at a state with Λ ≥ 0, α_𝕄 = 0 if and only if Λ₄ = Λ₅ = 0 (Thm 3.3(i), Cor. 3.5). They give no radius below which α_𝕄 > 0 guarantees Bloch positivity, so the ball radius is an engineering choice."

**E20 (SIV L379).**
- Old: "Geymonat et al. 1993 give general theory and a laminate example; their numerical values were not checked."
- New: "Geymonat et al. 1993 give general theory and a laminate example; their quoted loss-of-ellipticity stretch λ = 0.927 (stiff fraction 0.5, p. 288) was reproduced as 0.9267 from their Eqs. 6.6, 6.55 and 6.56 (verification note)."

### Evidence labels

- Change "verified as relayed" to **"verified"** for the chain and the long-wave equivalence (Add. L56).
- Change "(content relayed by …)" to "(full text verified)" (Add. L13).
- Keep "verified" on the k = 0 versus k → 0 subtlety (Add. L58), adding the GMT pinpoint.
- Change Müller 1987 from "cited from a reference list" to "verified as restated by Geymonat et al. 1993" (Ref. L503). Müller's own text remains unread.
- "Homogenized moduli are meaningful only before…" should not carry "verified" until reworded per E1.

---

## Gaps

- Müller 1987 itself was not read. Its content is verified only as GMT restate it (Thm 2.3 and the example in its Thm 4.3).
- Zhang et al.'s Eqs. 55–56 were not re-read here. The mapping β ↔ Λ, β(k → 0) ↔ Λ₄ = Λ₅ and B ↔ α_𝕄 rests on the definitions as quoted in the reports.
- GMT's Bloch-onset values λ_c⁻ (the solid curve in Fig. 6.6) are given only graphically, apart from the coincidence statement. Only the loss-of-ellipticity value 0.927 is quoted in the text, and only that value was reproduced.
- The ν = 0.3 laminate values are this note's computation from GMT's formulas. They have not been checked against an FE run.

---
*Evidence note (coordinator).* The laminate calculations above are preserved in [`evidence/gmt_laminate/`](evidence/gmt_laminate/) (NumPy-free Python, stdlib only). The coordinator re-ran both scripts and reproduced every quoted value: 0.9267 for GMT's quoted 0.927 (stiff fraction 0.5); for NEMOpt's Bonet–Wood law at ν = 0.3 and stiffness ratio 10, 0.8677 (volume fraction 0.5) and 0.8162 (0.2); at ν = 0.499, 0.9043 and 0.8665, approaching the incompressible formula's 0.9045 and 0.8668.
