# Verification against Bendsøe & Sigmund (1999), "Material interpolation schemes in topology optimization"

Source: M. P. Bendsøe and O. Sigmund, *Archive of Applied Mechanics* 69 (1999) 635–654, doi:[10.1007/s004190050248](https://doi.org/10.1007/s004190050248). Received 6 April 1999, accepted 19 May 1999. The supplied PDF (20 pages) was read in full on 2026-10-12.

Method:
- Page numbers below are journal page numbers (PDF page 1 = p. 635).
- The PDF text layer garbles symbols (ρ → "q", ν → "m", ≥ dropped). Every formula below was therefore read from the page images, and Eqs. 7–11 and 20–21 from 220-dpi renders. Wording in quotes comes from the text layer, checked against the images.
- Quotes are verbatim, so they keep the paper's British spelling ("grey", "analysed").
- Numbers marked "computed" come from a scratch extension of `research_notes/NEMOpt elastic bounds check/hs_simp_check.py`. The scratch copy (`bs1999_check.py`) is outside the repository; the formulas it uses are listed in §5.
- Tags: **[Paper]** means stated in the paper. **[Derived]** means my algebra from the paper's equations. **[Computed]** means numerical.
- Both reports were modified by another process at 20:28 during this check; it added NDR-0065. Every "old" string quoted in §4 was re-checked afterwards and is still present verbatim at the line cited.

---

## Bottom line

1. **The condition is confirmed verbatim.** 2D, Eq. 11, p. 639: p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} "(in 2-D)". 3D, Eq. 13, p. 640: p ≥ max{15(1 − ν⁰)/(7 − 5ν⁰), (3/2)(1 − ν⁰)/(1 − 2ν⁰)} "(in 3-D)". The hedges "usually attributed" and "exact form … not confirmed" in NDR-0064 can be removed.
2. **ν⁰ is the 2D Poisson's ratio of a generic 2D isotropic base material.** The paper never says plane stress or plane strain. It writes the bounds "for plane elasticity" and defines κ⁰ and μ⁰ through the intrinsic 2D relations (Eq. 9). Evaluating the condition with the in-plane ν⁰ = ν/(1 − ν) = 0.4286 is therefore **correct** for NEMOpt's plane-strain model. This gives p\* = 3.5, governed by the bulk term (the shear term gives 2.8).
   - Using the 3D ν = 0.3 would be wrong. It gives 3.08, governed by shear, which is the plane-stress answer.
3. **The NDR-0064 numbers are reproduced.** The p = 3 bulk excess peaks at 1.0163 × the bound at H = 14/15 ≈ 0.933 (analytic and computed), and it vanishes at p = 3.4996 ≈ 3.5 (bisection).
4. **New finding: the evidence script's 2D shear bound is wrong.** In `hs_upper`, `z = m1(k1 + 2m1)/(2(k1 + m1))` should be ζ = κμ/(κ + 2μ). The script's void shear bounds (7.96 / 10.82 / 14.23 at 0.4 / 0.5 / 0.6) are 6–9% too high; the correct values are 7.40 / 10.12 / 13.42, matching B&S Eq. 6. The bulk formula is correct, and no report quotes the shear numbers. NDR-0064's conclusions survive, because p = 3 still satisfies the shear bound (max ratio 0.99998).
5. **New support for an open NDR-0064 item.** B&S Eq. 21 is a strain-energy upper bound for *any* planar composite of void and an isotropic material, anisotropic ones included. At hydrostatic strain it reduces exactly to the 2D Hashin–Shtrikman bulk bound [Derived; Computed]. This supports applying the bulk bound to the hydrostatic-strain modulus of square-symmetric cells, an item the report lists as "still to be confirmed". At pure shear the same equation gives a looser bound than the isotropic HS shear bound, which supports "do not apply the isotropic shear bound to anisotropic cells".
6. **Claims checked in the notes.**
   - The "p ≥ 3 for ν⁰ = 1/3 in 2D" realizability claim is correct; it can be upgraded from [PR-bib] to [PR] and made precise.
   - The paper does **not** discuss RAMP.
   - The paper does **not** present modified SIMP as a void regularization. The algebraic form E_min + ρ^p(E⁰ − E_min) appears only as the two-material SIMP (Eq. 23, citing Sigmund & Torquato 1997), which the paper shows "will never satisfy the Hashin-Shtrikman bounds for all densities" (p. 647).
   - Its minimum-stiffness device is a lower bound on the **density variable**, 0 < ρ_min ≤ ρ (p. 638). That is neither an additive floor nor NEMOpt's clip on the stiffness.

---

## 1. The SIMP / Hashin–Shtrikman compatibility condition

### 1.1 The derivation chain in the paper [Paper]

| Item | Statement (formula from the page image) | Eq., page |
|---|---|---|
| SIMP | C_ijkl(ρ) = ρ^p C⁰_ijkl; base material "isotropic", characterized by "the Young' s modulus E⁰ and the Poisson ratio ν⁰" | Eq. 4, p. 638 |
| Premise | The SIMP tensor is isotropic "with a Young's modulus varying with ρ and a constant Poisson ratio, independent of ρ" | p. 638 |
| HS bounds | "written here for plane elasticity and for the limit of one phase being void": 0 ≤ κ ≤ ρκ⁰μ⁰/((1 − ρ)κ⁰ + μ⁰), 0 ≤ μ ≤ ρκ⁰μ⁰/((1 − ρ)(κ⁰ + 2μ⁰) + κ⁰) "(in 2-D)" | Eq. 6, p. 638 |
| Young's-modulus bound | 0 ≤ E ≤ E\* = ρE⁰/(3 − 2ρ) "(in 2-D)", citing Torquato et al. 1998 | Eq. 7, p. 639 |
| E-only condition | ρ^p E⁰ ≤ ρE⁰/(3 − 2ρ) for all ρ, "true if and only if p ≥ 3" | Eq. 8, p. 639 |
| Stronger condition | "the SIMP model presumes that the Poisson's ratio is independent of the density, and this leads to a stronger condition" | p. 639 |
| 2D moduli | κ⁰ = E⁰/(2(1 − ν⁰)), μ = E⁰/(2(1 + ν⁰)) "(in 2-D)" (the superscript 0 on μ is missing in print) | Eq. 9, p. 639 |
| Bound in E, ν form | ρ^pE⁰/(2(1 − ν⁰)) ≤ ρE⁰/(4 − 2(1 + ν⁰)ρ); ρ^pE⁰/(2(1 + ν⁰)) ≤ ρE⁰/(2(1 − ρ)(3 − ν⁰) + 2(1 + ν⁰)) | Eq. 10, p. 639 |
| **2D condition** | **p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} (in 2-D)**, "which in itself implies p ≥ 3" | **Eq. 11, p. 639** |
| Term origin | "The inequality p ≥ 2/(1 − ν⁰) comes from the bulk modulus bound, while the inequality p ≥ 4/(1 + ν⁰) is due to the shear modulus bound." | p. 639 |
| Examples | p\*(1/3) = 3; p\*(1/2) = 4; p\*(0) = 4; p\*(1) = ∞; p\*(−1) = ∞ "(in 2-D)"; "p\* = 3 holds only for ν⁰ = 1/3" | Eq. 12, p. 639 |
| Meaning | condition (11) "implies that the SIMP model can be made to satisfy the Hashin-Shtrikman bounds" | p. 639 |
| Caveat | "It is still an open problem if all material parameters satisfying the bounds also can be realized as composites" | p. 639 |
| Scope | "The discussion above holds for planar problems." | p. 640 |
| **3D condition** | **p ≥ max{15(1 − ν⁰)/(7 − 5ν⁰), (3/2)(1 − ν⁰)/(1 − 2ν⁰)} (in 3-D)**; derivation "omitted" | **Eq. 13, p. 640** |
| 3D examples | p ≥ 3 for ν⁰ = 1/3; p ≥ 2 for 1/5; p ≥ 15/7 for 0; p ≥ 5/2 for ν⁰ → −1; p → ∞ for ν⁰ → 1/2; "for ν = 1/3 we have the same bounds in 2-D and in 3-D" | Eq. 14, p. 640 |

Further notes on the 3D form:
- The paper does not say which 3D term is bulk and which is shear. [Derived] The (3/2)(1 − ν⁰)/(1 − 2ν⁰) term is the bulk condition: it equals 1 + 3κ⁰/(4μ⁰) from the 3D void HS bulk bound. The other term is therefore the shear condition.
- All transcribed examples in Eqs. 12 and 14 were re-evaluated and match [Computed].

### 1.2 Which ν⁰, and plane stress vs plane strain

- **What the paper does not say.** It never mentions plane stress or plane strain. It says "plane elasticity" (p. 638), "(in 2-D)" (Eqs. 6, 7, 9, 11, 12) and "planar problems" (p. 640).
- **ν⁰ is a 2D Poisson's ratio.** Eq. 9 is the intrinsic 2D isotropic relation between the 2D bulk modulus, the 2D Young's modulus and the 2D Poisson's ratio. The 2D examples run to ν⁰ = ±1 (Eq. 12), the limits of the 2D isotropic range, whereas the 3D examples stop at ν⁰ → 1/2 (Eq. 14). So in Eq. 11, ν⁰ is the base material's **2D (in-plane) Poisson's ratio**, defined through its 2D bulk and shear moduli: ν⁰ = (κ⁰ − μ⁰)/(κ⁰ + μ⁰) [Derived from Eq. 9].
- **The model is a generic 2D isotropic material.** Plane stress is the special case in which the 2D constants equal the 3D E and ν. Under plane strain the in-plane tensor is a 2D isotropic tensor with κ_2D = λ + μ and μ, so:
  - E_2D = E/(1 − ν²) = 109.89;
  - ν_2D = λ/(λ + 2μ) = ν/(1 − ν) = 0.428571;
  - Eq. 9 inverts exactly: E_2D/(2(1 − ν_2D)) = 96.1538 = λ + μ and E_2D/(2(1 + ν_2D)) = 38.4615 = μ [Computed].
- **The paper warns against mixing dimensions.** An interpolation built on 3D spherical-inclusion properties and used for planar problems violates the 2D HS bounds, "a feature stemming from the incompatibility of planar elasticity and spherical inclusions" (p. 642).
- **Verdict: ν⁰ = 0.4286 is the correct input for NEMOpt's plane-strain model.** It is consistent with notes §5, which already establish 0.4286 as the in-plane plane-strain ratio. Using the 3D ν = 0.3 would apply the plane-stress reading to a plane-strain model.
- **NEMOpt's density variable.** NEMOpt evaluates both the stiffness, max(H^p, ε), and the volume, ∫H(ρ̃), at Gauss points on the projected density H(ρ̃) ([NEMOpt.jl L59–60, L75–88](../../Full_Domain/NEMOpt.jl#L59)). The paper's ρ therefore corresponds to H, not to ρ̃.

### 1.3 Checks against the FEAD2026 constants (E = 100, ν = 0.3, plane strain) [Computed]

| Quantity | Value |
|---|---|
| Eq. 6 (void limit) vs the script's bulk formula | max relative difference 3.5e−8: **agrees** |
| Eq. 6 vs the script's shear formula | max relative difference **8.7e−2: disagrees** (see §1.4) |
| Eq. 11 with ν⁰ = 0.4286 | bulk term 3.5000, shear term 2.8000 → **p\* = 3.5** |
| Eq. 11 with ν⁰ = 0.3 (wrong for plane strain) | bulk term 2.857, shear term 3.077 → p\* = 3.08 (shear-governed) |
| p = 3, bulk SIMP/HS ratio | max **1.01630** at H = 0.9333. Analytic peak at H\* = 2(1 + a)/(3a) = 14/15, where a = κ⁰/μ⁰ = 2.5, giving the ratio H\*²(1 + a)/3 = 1.01630 |
| p = 3, shear SIMP/HS ratio | max 0.99998 (no excess) |
| Smallest p with no excess (bisection) | bulk 3.4996, shear 2.7997 |
| p = 3.5 | bulk max 1.00000 (touches at H → 1) |
| Eq. 8 (E-only) at p = 3 | max ratio 1.00000. **An E-only check would miss the plane-strain bulk excess**, as the paper's "stronger condition" remark predicts |
| Plane-stress contrast (κ_2D = 71.43), p = 3 | bulk fine; shear max 1.00047 at H = 0.988 |
| Void HS+ bulk at 0.4 / 0.5 / 0.6 | 15.3846 / 21.3675 / 28.8462; with the ε-phase: 15.4411 / 21.4257 / 28.9050 (NDR-0064 values **confirmed**) |
| κ = 28.4368 / HS+(0.6, void) | 0.9858 (the report's 98.6% **confirmed**) |

[Derived] General form of the check: with SIMP, the bulk ratio is ρ^(p−1)(1 + (κ⁰/μ⁰)(1 − ρ)). It stays ≤ 1 iff p ≥ 1 + κ⁰/μ⁰ = 2/(1 − ν⁰), which reproduces the bulk term of Eq. 11. The shear ratio gives p ≥ 2 + 2μ⁰/κ⁰ = 4/(1 + ν⁰).

### 1.4 Defect in the evidence script (`hs_simp_check.py`)

- **Where.** Line 17 uses `z = m1 * (k1 + 2 * m1) / (2 * (k1 + m1))`.
- **The correct form.** The 2D HS shear bound is μ_U = μ₁ + f₂ / [1/(μ₂ − μ₁) + f₁/(μ₁ + ζ₁)] with ζ₁ = κ₁μ₁/(κ₁ + 2μ₁). Equivalently, f₁/(μ₁ + ζ₁) = f₁(κ₁ + 2μ₁)/(2μ₁(κ₁ + μ₁)). The script's `z` is not this ζ.
- **Corrected values.** With ζ fixed, the general formula matches B&S Eq. 6 to 2.8e−8.

| Solid fraction | Script shear HS+ (void) | Correct (= B&S Eq. 6) | Correct, ε-phase |
|---|---|---|---|
| 0.4 | 7.9576 | 7.3964 | 7.4174 |
| 0.5 | 10.8173 | 10.1215 | 10.1423 |
| 0.6 | 14.2255 | 13.4168 | 13.4372 |

- **Impact.** Reports and notes do not quote the shear numbers, and the p = 3 shear conclusion is unchanged (0.9996 → 0.99998, still ≤ 1). The script is, however, linked as NDR-0064 evidence, and its shear output is non-conservative. The fix is listed in §4.6.
- **The same `z` line was copied** into `plane_stress_comparison.py` (NDR-0065 evidence, created concurrently at 20:28). That script uses only the bulk output, so its printed numbers are unaffected.

---

## 2. Claims attributed to the paper in the reports and notes

Grep for "Bendsøe" / "Bendsoe" / "s004190050248" / "Sigmund 1999" in `reports/` and `research_notes/`.

| # | Location | Claim | Verdict | Evidence |
|---|---|---|---|---|
| C1 | Main report L508 (NDR-0064 "Computed" bullet) | p = 3 exceeds the bulk bound by ≤ 1.6% near 0.93 and the excess vanishes at 3.5; this "matches the 2D compatibility condition usually attributed to Bendsøe & Sigmund (1999), p ≥ max{2/(1 − ν₀), 4/(1 + ν₀)}, evaluated with the in-plane Poisson's ratio ν₀ = 0.4286 … exact form … not confirmed" | **Confirmed, with stronger wording warranted.** The formula is Eq. 11, p. 639, verbatim. The in-plane ν⁰ is correct (§1.2). Numbers reproduced (§1.3). The attribution is direct, not "usually attributed" | Eq. 11; "(in 2-D)"; "The inequality p ≥ 2/(1 − ν⁰) comes from the bulk modulus bound" (p. 639) |
| C2 | Main report L506–507, NDR-0064 row L515 | 2D HS upper bulk bound 15.38 / 21.37 / 28.85; ε-phase 15.44 / 21.43 / 28.905 | **Confirmed.** The void values equal B&S Eq. 6. The ε-phase values use the general two-phase HS form, which the paper does not print (it gives only the void limit) | Eq. 6, p. 638 |
| C3 | Script L29 | "attribution to Bendsoe & Sigmund 1999 unconfirmed" | **Now confirmed.** The script's shear formula is wrong (§1.4) | Eq. 11, p. 639 |
| C4 | Notes L48 | Bibliographic entry: *Arch. Appl. Mech.* 69(9–10):635–654 | **Volume, year, pages and title confirmed** from the PDF header ("Archive of Applied Mechanics 69 (1999) 635–654"). The issue number (9–10) is not printed in the PDF, so it is neither confirmed nor contradicted | p. 635 |
| C5 | Notes L54 ([PR-bib]) | "discuss physical realizability of SIMP: p must exceed a Poisson-ratio-dependent bound (p ≥ 3 for ν0 = 1/3 in 2D)" | **Correct but imprecise.** (i) It is "p ≥ p\*(ν⁰)", not "exceed". (ii) In 2D, p\* ≥ 3 for every ν⁰, with equality "only for ν⁰ = 1/3". (iii) The condition is compatibility with the HS bounds, which is necessary for realization. Full realizability "is still an open problem" (p. 639), but microstructures were built by inverse homogenization for ν⁰ = 1/3 with p = 3 (Fig. 3, p. 640) and for ν⁰ = 0 and 0.5 with p = 4 (Fig. 4, p. 641). Upgrade the tag to [PR] | Sec. 3.2, Eqs. 6–12, pp. 638–639 |
| C6 | Notes L36 Takeaway | "Classic SIMP keeps a lower bound on the design variable." (no citation) | **Supported** by the paper | "a small lower bound, 0 < ρmin ≤ ρ, is usually imposed, in order to avoid a singular FEM problem" (p. 638) |
| C7 | Notes L53 ([PR-bib]) | Sigmund 2007 and Andreassen et al. 2011 "introduced and popularized" modified SIMP E_min + x^pΔE, which "replaced the older approach of a lower bound x_min" | **Partly supported. "Introduced" is too strong.** The older approach is documented in B&S p. 638 (C6). The additive algebraic form appears in B&S as the **two-material** SIMP, C(ρ) = ρ^p C¹ + (1 − ρ^p)C² (Eq. 23, p. 646, "as suggested in [41]", i.e., Sigmund & Torquato 1997). It is presented for two real materials, not as a void regularization. Sigmund 2007 / Andreassen 2011 themselves were not checked here | Eq. 23, p. 646 |
| C8 | Notes L55, L71; julia_architecture.md L96 | RAMP E = E⁰ρ/(1 + q(1 − ρ)) | **Not attributed to B&S, and B&S does not discuss RAMP** (the paper predates Stolpe & Svanberg 2001). [Derived] With q = 2 the RAMP form equals the paper's 2D HS Young's-modulus bound ρE⁰/(3 − 2ρ) (Eqs. 7 and 17) | Eqs. 7, 17, pp. 639, 642 |
| C9 | Main report NDR-0033 (L461), L88, L212, L229 | Modified SIMP ε + (1 − ε)H^p cited to Ferrari & Sigmund 2020; the clip is Wallin & Tortorelli's Eq. 71 | **Not attributed to B&S; no conflict.** The paper is relevant context, though: its floor is a density lower bound (C6); its two-material SIMP (Eq. 23) has the modified-SIMP form; and Eq. 23 cannot satisfy the HS bounds for all densities (p. 647). See §3.1 and §4.4 | pp. 638, 646–647 |

Other Bendsøe hits are different works and out of scope; they were not verified here:
- Bendsøe & Sigmund (2003, book), in the buckling notes;
- Neves, Sigmund & Bendsøe (2002);
- Amir, Bendsøe & Sigmund (2009);
- Bendsøe (1989), notes L47.

The buckling addendum report does not cite B&S 1999.

---

## 3. What the paper says on the three topics

### 3.1 Minimum stiffness: clipped vs additive floor

- **The only void device for single-material SIMP is a density lower bound** (p. 638): "a small lower bound, 0 < ρmin ≤ ρ, is usually imposed, in order to avoid a singular FEM problem". The stiffness floor ρ_min^p C⁰ follows from the box constraint, and the volume counts ρ ≥ ρ_min.
- **Neither of NEMOpt's laws appears in the paper.**
  - **Clip.** It contains no clip on the stiffness (max(ρ^p, ε)).
  - [Derived] NEMOpt's clip satisfies max(H^p, ε) = max(H, ε^{1/p})^p. That is classic SIMP evaluated at the floored density ε^{1/p}, applied inside the stiffness law only. Unlike the paper's ρ_min, H and the volume may fall below ε^{1/p} at constant stiffness, with zero stiffness sensitivity. The effective ρ_min also moves with p: 0.0005 at p = 1, 0.0224 at p = 2, 0.0794 at p = 3.
  - **Additive floor.** The paper does not present E_min + ρ^p(E⁰ − E_min) for material and void. The same algebra is its two-material SIMP (Eq. 23, p. 646).
- **HS behavior of a soft second phase** (p. 647): "The derivative at zero density of the lower bound in (25) is positive", so "a SIMP model in the form (23) will never satisfy the Hashin-Shtrikman bounds for all densities", although "it is possible to keep the SIMP model fairly close". Its Eq. 25 example uses ν = 1/3 for both phases. The slope argument is general [Derived].
- **Applied to NEMOpt's ε = 5e−4** [Computed]:
  - Against the **void** HS bounds, the floor εC⁰ is too stiff for H < 1.75e−3 (bulk) and H < 1.40e−3 (shear). The paper's density floor has no such excess when p ≥ p\*.
  - Against a two-phase composite with the **ε-material** as phase 2:
    - the clip falls to 0.892 × the HS lower bulk bound at H = 0.0794;
    - modified SIMP falls to 0.986 × at H ≈ 0.015;
    - near H = 0.93 both still exceed the upper bulk bound by about 1.6% (clip 1.0159, modified 1.0160).
  - Bound compatibility therefore does not separate the two laws.
- **Why the void matters** (Sec. 7, p. 650): low-density areas "are also included in the analysis for each feasible design", which causes difficulties "both in the formulation of the problem as well as in the numerical treatment" (stress singularity; localized modes, §3.4).

### 3.2 Poisson's ratio of the interpolated material

- **SIMP keeps the base Poisson's ratio at every density**: "a constant Poisson ratio, independent of ρ" (p. 638). This is exactly why the bulk and shear bounds give "a stronger condition" than the Young's-modulus bound (p. 639).
  - [Derived] In NEMOpt's linear limit every gray Gauss point carries f(H)·C⁰, so its in-plane ν is 0.4286 at all H. Under plane strain the bulk bound is the binding one (p\* = 3.5).
- **The HS-upper-bound interpolation** (Eq. 17, p. 642): E(ρ) = ρE⁰/(3 − 2ρ) and ν(ρ) = (1 − ρ(1 − ν⁰))/(3 − 2ρ). Here "independent of the Poisson ratio of the base material, the low volume fraction limit has a Poisson ratio equal to 1/3". It is realizable by "an isotropic rank-3 lamination".
  - [Computed] With ν⁰ = 0.4286: ν(0) = 0.333, ν(0.5) = 0.357, ν(1) = 0.4286.
- **Reuss–Voigt interpolation** (Eqs. 18–19, p. 642): HS-compatible iff α ≤ min{(1 − ν⁰)/2, (1 + ν⁰)/4}, at most 1/3, which "is only possible if ν⁰ = 1/3".
- **Multi-material sections** (Sec. 5, p. 646) restrict the analysis to the 2D case "where both base materials as well as the interpolations have Poisson's ratio equal to 1/3".
- **Dimensional consistency** (p. 642): the 3D spherical-inclusion model violates the 2D HS bounds (see §1.2).

### 3.3 Must intermediate densities respect HS for designs to be meaningful?

The paper's answer is **no for black-and-white designs, yes for interpreting gray**.

- p. 636: "one can in essence choose to ignore the physical relevance of intermediate steps which may include 'grey'". However, the question "is often raised", since most schemes give designs "not completely clear of 'grey'". Also, "the physical realization of all feasible designs plays a role when interpreting results from a premature termination of an optimization algorithm."
- p. 638: a 0/1 design is one "for which the performance has been evaluated with a correct physical model".
- p. 643, on p-continuation: "This procedure is thus a compromise, since initial designs will be analysed using an interpolation which is not realizable as a composite structure." Continuation runs "until the final design is arrived at for a power satisfying (11) or (13)".
  - [Derived] NEMOpt stops at p = 3, below the plane-strain p\* = 3.5, so even its terminal interpolation is not HS-compatible in gray. At p = 1 it is the Voigt interpolation (Eq. 16), which "cannot be realized by composites of material and void" (p. 641).
- p. 651 (Conclusions): "a better computational scheme can be obtained if one allows for a violation of the bounds on properties of composites", especially "where the bounds do not allow for a high enough penalization".
- p. 650: "micromechanical considerations should always play a role in the development of interpolation schemes".

**Consequence for NEMOpt.** This supports NDR-0064's rule: apply HS checks to thresholded designs (or to the F_M = I evaluation of them), never to raw gray designs. A 1.6% excess in gray elements is a property of p = 3 SIMP under plane strain, not a bug.

### 3.4 Other passages relevant to NEMOpt

- **HS bulk bound for anisotropic cells** (Sec. 4, p. 644). "For planar problems, any composite, constructed from void and an isotropic, linearly elastic material" satisfies the complementary-energy bound Eq. 20 (from Allaire & Kohn, ref. [57]) and the strain-energy bound Eq. 21 (from Jog, Haber & Bendsøe 1994, ref. [20]). Eq. 21 "holds for any strain tensor ε with principal strains ε_I, ε_II ordered such that |ε_I| ≤ |ε_II|".
  - [Derived] At hydrostatic strain ε = e·I, only Eq. 21's second case applies. It gives ½C:ε:ε ≤ ρE⁰e²/(2 − ρ(1 + ν⁰)) = 2κ_HS e², where κ_HS is Eq. 6. Since ½C:ε:ε = 2κ_hydro e² with κ_hydro = ¼(C₁₁₁₁ + C₂₂₂₂ + C₁₁₂₂ + C₂₂₁₁), it follows that **κ_hydro ≤ κ_HS for any anisotropic void/solid cell** in the linear limit.
  - [Computed] Eq. 21 equals 2κ_HS to machine precision at ρ = 0.4, 0.5, 0.6 and 0.93, using E_2D and ν_2D. This answers the main report's "assumption still to be confirmed".
- **Shear for anisotropic cells.** [Derived] At pure shear, Eq. 21's first case bounds 𝔸₁₂₁₂ + 𝔸₂₁₂₁ (= 2C₁₂₁₂ at F = I) by ρE⁰/(2 − ρ(1 − ν⁰)) [Computed]:

  | Solid fraction | Eq. 21 bound on 2C₁₂₁₂ | Isotropic HS (2μ_HS) |
  |---|---|---|
  | 0.4 | 24.81 | 14.79 |
  | 0.5 | 32.05 | 20.24 |
  | 0.6 | 39.79 | 26.83 |

  This confirms that the isotropic shear bound must not be applied to anisotropic cells, and it offers an optional anisotropic replacement check.
- **Nonlinear problems** (Sec. 6.2, p. 650). These pose "an even greater theoretical challenge", owing partly to "difficulties in deriving bounding theorems for such problems". This is consistent with NDR-0064(c).
- **Stability** (Sec. 7, p. 651). In intermediate designs it "can become unclear what are the relevant values to consider", with "localized modes which appear in low density regions and which should be filtered out". The paper cites Neves, Rodrigues & Guedes 1995 for buckling. This is early primary support for the addendum's void and soft-band mode separation and for NDR-0060.
- **Penalization trade-off** (p. 643): "a too severe penalization of intermediate density can lead to designs which are local minima, and which are very sensitive to the choice of the initial design".
- **Terminology** (p. 636): the label "'fictitious material model' is actually misleading", because suitably penalized SIMP is realizable.

---

## 4. Recommended edits (old → new)

These are recommendations only; the reports were not edited. Quotation marks inside the new text are the paper's words.

### 4.1 Main report, NDR-0064 "Computed" bullet (L508)

Old:
> - **p = 3 SIMP exceeds the bulk bound** by up to 1.6% at a density of about 0.93. The excess vanishes at p = 3.5, which matches the 2D compatibility condition usually attributed to Bendsøe & Sigmund (1999), p ≥ max{2/(1 − ν₀), 4/(1 + ν₀)}, evaluated with the in-plane Poisson's ratio ν₀ = 0.4286 ([Bendsøe & Sigmund 1999](https://link.springer.com/doi/10.1007/s004190050248); the formula's exact form in the paper was not confirmed).

New:
> - **p = 3 SIMP exceeds the bulk bound** by up to 1.6% (1.0163 times the bound at projected density H = 14/15 ≈ 0.93). The excess vanishes at p = 3.5. This is Bendsøe & Sigmund's 2D condition p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} (Eq. 11, p. 639), whose bulk term gives 3.5 here and whose shear term gives 2.8 ([Bendsøe & Sigmund 1999](https://link.springer.com/doi/10.1007/s004190050248); verified in the full text). The paper derives it for a generic 2D isotropic material from the 2D Hashin–Shtrikman bounds (Eq. 6) and κ⁰ = E⁰/(2(1 − ν⁰)), μ⁰ = E⁰/(2(1 + ν⁰)) (Eq. 9). So ν⁰ is the in-plane Poisson's ratio, which under plane strain is ν/(1 − ν) = 0.4286; the 3D ν = 0.3 would give the plane-stress value 3.08. In 2D, p\* = 3 "holds only for ν⁰ = 1/3".

Optionally append to the "Consequence" bullet (L509):
> Bendsøe & Sigmund (1999, pp. 636, 651) likewise treat bound compatibility of gray as optional when designs end black-and-white.

### 4.2 Main report, NDR-0064 table row (L515)

(a) Old:
> Never apply it to raw gray designs (p = 3 SIMP exceeds the bound by up to 1.6%) or to finite-strain tangents.

New:
> Never apply it to raw gray designs or to finite-strain tangents. Gray elements can exceed the bound by up to 1.6%, because p = 3 is below the plane-strain threshold p\* = 3.5 of Bendsøe & Sigmund 1999, Eq. 11.

(b) Old:
> Do not apply the isotropic shear bound to anisotropic cells.

New (the second sentence is optional):
> Do not apply the isotropic shear bound to anisotropic cells. Bendsøe & Sigmund's (1999) strain-energy bound (Eq. 21) allows 𝔸₁₂₁₂ + 𝔸₂₁₂₁ up to 39.8 at solid fraction 0.6, against 26.8 from the isotropic bound (derived), and may serve as an anisotropic shear check.

(c) Acceptance evidence. Old:
> before (a) becomes a hard gate, confirm in the composites literature (e.g., Milton 2002) that the bulk bound holds for the hydrostatic-strain modulus of anisotropic cells

New:
> the bulk bound's validity for the hydrostatic-strain modulus of anisotropic void/solid cells follows from Bendsøe & Sigmund's (1999) Eq. 21. That bound covers any planar void/isotropic composite and reduces at hydrostatic strain to the 2D Hashin–Shtrikman bulk bound (derived; checked numerically at solid fractions 0.4–0.93). Before (a) becomes a hard gate, cross-check this against the original statement (Allaire & Kohn 1993, Quart. Appl. Math. 51:643; or Milton 2002).

(d) Acceptance evidence. Old:
> Bound formulas reproduce the computed values (15.38 / 21.37 / 28.85; 28.905 with ε)

New:
> Bound formulas reproduce the computed values (15.38 / 21.37 / 28.85, equal to Bendsøe & Sigmund 1999 Eq. 6; 28.905 with ε), after the evidence script's 2D shear formula is corrected to ζ = κμ/(κ + 2μ)

### 4.3 Main report, κ plausibility caveat (L171–174)

Old (L173):
> - p = 3 SIMP itself exceeds the bulk bound by up to 1.6% in gray elements under plane strain (computed);

New:
> - p = 3 SIMP itself exceeds the bulk bound by up to 1.6% in gray elements under plane strain (computed). This is expected, because p = 3 is below the 2D threshold p\* = 3.5 of Bendsøe & Sigmund (1999, Eq. 11) for the in-plane Poisson's ratio 0.4286;

Old (L174):
> - applying the isotropic bulk bound to a square-symmetric cell's hydrostatic-strain modulus is an assumption still to be confirmed.

New:
> - applying the bulk bound to a square-symmetric cell's hydrostatic-strain modulus is supported by Bendsøe & Sigmund's (1999) Eq. 21. That strain-energy bound covers any planar void/isotropic composite and reduces at hydrostatic strain to the 2D Hashin–Shtrikman bulk bound (derived; to be cross-checked against the original before any hard gate).

Optional, in L171, after "upper bound of 28.85 for a porous solid at 60% solid", add:
> (Bendsøe & Sigmund 1999, Eq. 6)

### 4.4 Main report, NDR-0033 (L461)

Old (decision column):
> Stiffness interpolation. The clip `max(H^p, ε)` is Wallin & Tortorelli's reference law (Eq. 71, ε = 1e−5 there).

New:
> Stiffness interpolation. The clip `max(H^p, ε)` is Wallin & Tortorelli's reference law (Eq. 71, ε = 1e−5 there).
> - It equals classic SIMP evaluated at max(H, ε^{1/p}). Classic SIMP, however, bounds the design variable itself, 0 < ρ_min ≤ ρ, "to avoid a singular FEM problem" (Bendsøe & Sigmund 1999, p. 638). The clip instead lets H and the volume fall below ε^{1/p} at constant stiffness.
> - Modified SIMP ε + (1 − ε)H^p has the form of Bendsøe & Sigmund's two-material SIMP (Eq. 23), which "will never satisfy the Hashin-Shtrikman bounds for all densities" (p. 647).
> - At p = 3 under plane strain, neither law is bound-compatible (both exceed the bulk bound by about 1.6% near H = 0.93). Bound compatibility therefore does not decide between them.

The other columns are unchanged. Optionally add to the L212 row's evidence column:
> classic SIMP floor on the density: Bendsøe & Sigmund 1999, p. 638

### 4.5 Notes: `numerical_method_choices.md`

(a) L36. Old:
> Classic SIMP keeps a lower bound on the design variable.

New:
> Classic SIMP keeps a lower bound on the design variable ("0 < ρmin ≤ ρ, is usually imposed, in order to avoid a singular FEM problem", Bendsøe & Sigmund 1999, p. 638). The clip equals classic SIMP evaluated at max(H, ε^{1/p}), with the floor applied to the stiffness only, not to the variable or the volume.

Separately, L36's claim "neither classic nor modified SIMP" and its "SWITCH … new default = modified SIMP" are already superseded by the main report (Wallin & Tortorelli Eq. 71; NDR-0033 keeps the clip as default). Consider a pointer to that.

(b) L53. Old:
> Sigmund 2007 and Andreassen et al. 2011 introduced and popularized "modified SIMP" (E_min + x^p ΔE, E_min ≈ 10⁻⁹E0). It replaced the older approach of a lower bound x_min on the design variable with E = x^p E0, and it decouples void stiffness from the penalization.

New:
> Sigmund 2007 and Andreassen et al. 2011 popularized "modified SIMP" (E_min + x^p ΔE, E_min ≈ 10⁻⁹E0) as a void regularization (not re-verified). The same algebraic form appears earlier as the two-material SIMP C(ρ) = ρ^p C¹ + (1 − ρ^p)C² ([PR] Bendsøe & Sigmund 1999, Eq. 23, citing Sigmund & Torquato 1997). It replaced the older approach of a lower bound x_min on the design variable with E = x^p E0 ([PR] Bendsøe & Sigmund 1999, p. 638), and it decouples void stiffness from the penalization.

(c) L54: move out of the [PR-bib] list. Old:
> Bendsøe & Sigmund 1999 discuss physical realizability of SIMP: p must exceed a Poisson-ratio-dependent bound (p ≥ 3 for ν0 = 1/3 in 2D).

New, as a [PR] bullet under Cited Findings:
> [PR] Bendsøe & Sigmund 1999 (supplied PDF, read in full). SIMP satisfies the Hashin–Shtrikman bounds for a void/solid composite iff:
> - in 2D: p ≥ p\*(ν⁰) = max{2/(1 − ν⁰), 4/(1 + ν⁰)} (Eq. 11, p. 639; the bulk and shear terms respectively);
> - in 3D: p ≥ max{15(1 − ν⁰)/(7 − 5ν⁰), (3/2)(1 − ν⁰)/(1 − 2ν⁰)} (Eq. 13, p. 640).
>
> In 2D, p\* ≥ 3, with equality "only for ν⁰ = 1/3". ν⁰ is the 2D Poisson's ratio (Eq. 9), so NEMOpt's plane-strain solid (ν⁰ = 0.4286) needs p ≥ 3.5.
>
> Compatibility is necessary for realization. Whether every compatible tensor is realizable "is still an open problem" (p. 639), but microstructures were designed for ν⁰ = 1/3 with p = 3 (Fig. 3) and for ν⁰ = 0 and 0.5 with p = 4 (Fig. 4).

(d) L55. Old:
> RAMP, E = E0 ρ/(1+q(1−ρ)), has non-zero slope at ρ = 0.

New:
> RAMP, E = E0 ρ/(1+q(1−ρ)), has non-zero slope at ρ = 0 (Stolpe & Svanberg 2001; Bendsøe & Sigmund 1999 predates and does not discuss RAMP). With q = 2 it coincides with their 2D Hashin–Shtrikman Young's-modulus bound ρE⁰/(3 − 2ρ) (Eqs. 7, 17; algebra).

(e) §1 Inferences. Add after L68:
> - **Bound compatibility of the floors** (computed, plane strain):
>   - Against void, εC⁰ exceeds the 2D HS bulk bound for H < 1.75e−3.
>   - Against an ε-phase, the clip falls to 0.89 × the HS lower bulk bound at H = 0.0794, and modified SIMP to 0.99 × near H = 0.015.
>   - Near H = 0.93 both exceed the upper bulk bound by about 1.6% at p = 3.
>   - This is consistent with Bendsøe & Sigmund 1999, p. 647, and is irrelevant for black-and-white designs (p. 651).

(f) §3 Cited Findings. Add after L207:
> - [PR] Bendsøe & Sigmund 1999, p. 643: p-continuation "is often advisable", raising p "until the final design is arrived at for a power satisfying (11) or (13)". It is "a compromise, since initial designs will be analysed using an interpolation which is not realizable as a composite structure." NEMOpt's terminal p = 3 is below the plane-strain p\* = 3.5 (computed).

(g) §5. At L332, append:
> This in-plane ratio is the ν⁰ that enters Bendsøe & Sigmund's 2D SIMP condition (1999, Eqs. 9 and 11).

At L353 (Gaps), append:
> Bendsøe & Sigmund 1999 evaluate their 2D condition at ν⁰ = ±1 "(in 2-D)" (Eq. 12), consistent with the 2D range −1 ≤ ν ≤ 1. This is indirect support only.

(h) L48. Optionally append:
> Volume, year and pages confirmed from the PDF; the issue number is not printed in the PDF.

### 4.6 Evidence script `research_notes/NEMOpt elastic bounds check/hs_simp_check.py`

This is outside the four requested targets, but it is the linked NDR-0064 evidence.

- L17. Old:
  > `z = m1 * (k1 + 2 * m1) / (2 * (k1 + m1))`

  New:
  > `z = k1 * m1 / (k1 + 2 * m1)  # 2D HS shear; matches Bendsoe & Sigmund 1999 Eq. 6 (void limit)`
- L29. Old:
  > `(attribution to Bendsoe & Sigmund 1999 unconfirmed)`

  New:
  > `(Bendsoe & Sigmund 1999, Eq. 11, p. 639; nu0 = 2D Poisson's ratio)`
- Optional. Add a check that the bulk/shear ratios equal B&S Eq. 6 in closed form, and print the analytic peak H = 14/15 with ratio 1.0163.
- Apply the same L12 fix to `plane_stress_comparison.py`, and have it also print the shear ratio (plane stress: 1.00047).

### 4.7 Main report, NDR-0065 linear-constants table (L527–532)

NDR-0065 was added concurrently and is outside the requested targets, but the paper bears on it directly.

- **The issue.** The plane-stress cell "p = 3 SIMP against that bound: stays at or just below it (ratio 0.9999)" is true for the **bulk** bound only. With plane-stress constants, ν⁰ = 0.3, and B&S Eq. 11 is shear-governed: p\* = 4/1.3 = 3.08. So p = 3 SIMP exceeds the 2D HS **shear** bound by up to 0.05%: ratio 1.00047 at H = 0.988, analytic and computed with the corrected shear formula.
- **Edit.** Add one row to the table:

  | | Plane strain (FEAD2026) | Plane stress |
  |---|---|---|
  | SIMP threshold p\* (Bendsøe & Sigmund 1999, Eq. 11) | 3.5 (bulk-governed) | 3.08 (shear-governed; p = 3 exceeds the shear bound by up to 0.05%) |

- **Generalized plane strain** [Derived; not covered by B&S]. Freeing one uniform out-of-plane stretch can only lower the energy relative to plane strain, and it is more constrained than pointwise plane stress. So the cell's in-plane moduli satisfy plane stress ≤ GPS ≤ plane strain. The plane-strain HS bound (κ_s = λ + μ) therefore stays a valid, possibly loose, upper bound in GPS mode, whereas the plane-stress HS bound is not guaranteed to bound GPS cells. Mention this in the "recompute the NDR-0064 bound values per mode" requirement.

### 4.8 Optional: buckling addendum

- **"Soft-band and void modes are separated before any verdict" (L143).** Add after the first sentence:
  > Bendsøe & Sigmund (1999, p. 651) already flag "localized modes which appear in low density regions and which should be filtered out" for buckling and vibration, citing Neves, Rodrigues & Guedes (1995).
- **NDR-0060 (L316).** Add "Bendsøe & Sigmund 1999, §7" to its sources.

---

## 5. Gaps and caveats

- **Two strong conclusions are mine, not the paper's.**
  - The plane-strain mapping (ν⁰ = ν/(1 − ν)) is inferred from Eq. 9 and the 2D examples; the paper does not state it.
  - The use of Eq. 21 for the hydrostatic-strain modulus of anisotropic cells is derived from the equation as printed.
  - Both are straightforward, but cross-check Eq. 21 against Allaire & Kohn (1993) or Jog, Haber & Bendsøe (1994) before NDR-0064(a) becomes a hard gate.
- **Eq. 21's case ordering.** It is printed as |ε_I| ≤ |ε_II|. The hydrostatic and pure-shear cases used here do not depend on the ordering.
- **Scope of the HS checks.** All checks are linear-limit (F_M = I) statements. As NDR-0064 already says, they do not apply to finite-strain tangents.
- **Not checked.** The issue number "9–10" (not printed in the PDF), and the Sigmund 2007 / Andreassen et al. 2011 history claims.
- **Formulas used in the scratch check** (`bs1999_check.py`):
  - plane strain: κ = λ + μ, ν_2D = (κ − μ)/(κ + μ), E_2D = 4κμ/(κ + μ);
  - B&S Eq. 6 (void limit);
  - two-phase 2D HS: κ_U = κ₁ + f₂/[1/(κ₂ − κ₁) + f₁/(κ₁ + μ₁)] and μ_U = μ₁ + f₂/[1/(μ₂ − μ₁) + f₁/(μ₁ + κ₁μ₁/(κ₁ + 2μ₁))], with the lower bounds by swapping the phases;
  - B&S Eq. 21 evaluated with E_2D and ν_2D;
  - bisection on p for the smallest exponent with SIMP/HS ≤ 1 on a 99,999-point grid in ρ ∈ [1e−4, 1 − 1e−4].
