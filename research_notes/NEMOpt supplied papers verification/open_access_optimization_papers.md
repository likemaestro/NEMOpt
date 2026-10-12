# Open-access verification of six optimization papers cited by the NEMOpt reports

Date of all retrievals: 2026-10-12. Scope: find legal open-access full texts of six papers, check every claim the two reports and the research notes attribute to them, and propose exact edits. The reports and decision records were not edited.

## Method and access rules

- Sources used: publisher open-access or free-access pages, DTU Orbit, Lund University Publications (LUP), the Internet Archive availability API, arXiv API, Crossref and Semantic Scholar metadata. No shadow library was contacted.
- Each download went into its own empty directory under the session scratchpad (`oa_downloads/<short-name>/`). PDFs were read with the Read tool only. Nothing was executed from those directories, and no PDF was copied into the repository.
- Page citations: Thomsen et al. are cited by the page numbers printed on the DTU accepted author manuscript (AAM p. 1 to 23); equation and figure numbers are the paper's own. Dalklint et al. and Haftka et al. are cited by journal page.
- Verdicts: **confirmed**; **needs correction** (the paper says something different or more specific); **not in the paper** (the paper does not contain the claim); **not verifiable** (no legal full text obtained).

## Summary of access

| # | Paper | Obtained? | Version and legal source |
|---|---|---|---|
| 1 | Thomsen, Wang & Sigmund 2018, CMAME 339:115–136 | Yes | Accepted author manuscript (peer-reviewed version), DTU Orbit: [backend.orbit.dtu.dk/ws/files/147927301/MARAC_1_s2.0_S004578251830210X_main.pdf](https://backend.orbit.dtu.dk/ws/files/147927301/MARAC_1_s2.0_S004578251830210X_main.pdf) (record: [orbit.dtu.dk](https://orbit.dtu.dk/en/publications/buckling-strength-topology-optimization-of-2d-periodic-materials-/)). The `orbit.dtu.dk/files/...` link returns a Cloudflare 403; the `backend.orbit.dtu.dk/ws/files/...` link works. |
| 2 | Dalklint, Wallin & Tortorelli 2021, SMO 64:1751–1763 | Yes | Version of record, CC BY 4.0 ("Open access funding provided by Lund University"): [link.springer.com/content/pdf/10.1007/s00158-021-03012-z.pdf](https://link.springer.com/content/pdf/10.1007/s00158-021-03012-z.pdf). LUP record [dd94e13b](https://lup.lub.lu.se/search/publication/dd94e13b-814e-4993-a67a-42f622e1fd58) has no file. |
| 3 | Sigmund 2022, SMO 65(11):315 | **No** (abstract only) | Paywalled ("preview of subscription content" on [Springer](https://link.springer.com/article/10.1007/s00158-022-03427-2)); the [DTU Orbit record](https://orbit.dtu.dk/en/publications/on-benchmarking-and-good-scientific-practise-in-topology-optimiza/) has no file; Semantic Scholar lists no open copy; no arXiv version. |
| 4 | Haftka, Zhou & Queipo 2019, SMO 60(2):405–409 (editorial) | Yes | Publisher free access (page metadata: access type "permanently-free"): [link.springer.com/content/pdf/10.1007/s00158-019-02298-4.pdf](https://link.springer.com/content/pdf/10.1007/s00158-019-02298-4.pdf). The same session was refused the paywalled Sigmund 2022 PDF, so this was not subscription access. |
| 5 | Svanberg 2007, "MMA and GCMMA – two methods for nonlinear optimization" (KTH notes) | **No** (copies exist in the Wayback Machine but could not be reached) | `people.kth.se/~krille/` now returns 404 with "This page is no longer in use"; `smoptit.se` shows "Account Suspended". The Internet Archive lists snapshots of `mmagcmma.pdf` at 20210204193223, 20220120004729 and 20250316062505, and of `gcmma07.pdf` at 20230603101120. Every connection to `web.archive.org` (also `wayback.archive.org`) was reset right after the TLS handshake, so none could be downloaded. |
| 6 | Lazarov & Sigmund 2011, IJNME 86(6):765–781 | **No** (abstract only) | The [DTU Orbit record](https://orbit.dtu.dk/en/publications/filters-in-topology-optimization-based-on-helmholtztype-different/) has no file (no postprint); Wiley returned a Cloudflare 403; Semantic Scholar marks it closed. Claims were checked against a secondary source by the same group: the Lazarov, Wang & Sigmund (2016) postprint on DTU Orbit ([main.pdf](https://backend.orbit.dtu.dk/ws/files/140201797/main.pdf)). |

**Note on arXiv:1809.05300.** The task brief suggested that arXiv:1809.05300 might be Thomsen et al. It is not. The arXiv API lists 1809.05300 as Ferrari & Sigmund, "Revisiting topology optimization with buckling constraints" (submitted 14 Sep 2018), and the research notes already cite it correctly for Ferrari & Sigmund. arXiv searches by author and title found no version of Thomsen et al. 2018, Dalklint et al. 2021, Sigmund 2022, Haftka et al. 2019 or Lazarov & Sigmund 2011.

---

## 1. Thomsen, Wang & Sigmund (2018): read in full (AAM, 23 pp.)

### What the paper actually does

- **Model.** "The employed model is based on linear elasticity, linear homogenization theory, a linearized stability condition and Bloch-wave analysis" (AAM p. 2). A linearized criterion "imposes a restriction of small strain fields prior to instability" (p. 3). Macroscale instability "is beyond the scope of the present paper" (p. 3).
- **Load and prestress.** A macroscopic **stress** σ̄ is prescribed; ε̄ = Ē⁻¹σ̄ and the element prestress is σ_e = E_e(I − B_eX_e)ε̄ (Eq. 12).
- **Eigenproblem, as written.** (K − PG)φ = 0 (Eq. 10). After the Bloch transformation it is (K_k − P G_k)φ_k = 0 with K_k = T(k)ᴴKT(k) (Eq. 16). The optimization uses the inverse load factor λ = 1/P in (G_k(σ̄) − λ_i K_k)φ_ik = 0 (Eq. 19b). In this paper K₀ denotes the stiffness reduced at **k = 0** (Eq. 15), not a constitutive part of K(k).
- **Objective.** A min–max of λ_i(k) (Eq. 19a), implemented "as a bound formulation [25]" with K = 25 eigenvalues and I = 10 bands (p. 8). It is **not** KS aggregation. Multiple eigenvalues are handled with a simplified Seyranian-type generalized gradient that "is only approximative" (Appendix A, p. 20).
- **Wave-vector sampling.** In the biaxial case the problem is "discretized into 7 points along the boundary of the irreducible Brillouin zone": the symmetry points "and one near Γ to capture the long wavelength failure mode" (p. 11, Fig. 6). Uniaxial and mixed loading use "an extended number of optimization points" (p. 13, Fig. 11), roughly a dozen including several near Γ. Boundary-only sampling is justified by studies that "have shown (but not proven)" it suffices; final designs get "a complete analysis" over the zone (p. 5).
- **Interpolation.** Modified SIMP (ρ_min + (1 − ρ_min)ρ^p)E₀ for K and fᴶ, pure ρ^pE₀ for G (Eq. 17). The purpose is "to eliminate the notoriously known problem of artificial modes"; ρ_min = 10⁻⁴ (p. 7).
- **Settings.** 200 × 200 Q4 mesh, E₀ = 1, ν = 0.3, plane stress, V* = 0.3, p = 3, convolution **density filter** with R = 0.03d, η = 0.5, β doubled every 100 iterations to 256, MMA (p. 10). No wall-clock or CPU times are reported anywhere.
- **Symmetry.** Design symmetry is imposed through a mapping ρ = Lμ (Eq. 18). It enforces material symmetry and also "substantially lower[s] the computational cost of the stability analysis as it reduces the unique part of the Brillouin zone" (p. 7). Conclusion: "few geometric symmetries yields the highest strength potential at the design load but also reduces overall load robustness" (p. 20).
- **Results.** Strength gains of 344% (square) and 379% (triangle) in the biaxial case (p. 11; Table 1 ratios 4.44, 2.35 and 4.79), and bulk-modulus losses of 17–26% (p. 12). Optimized bands are "tightly stacked" (p. 11), and multiple critical modes make "a catastrophic post-buckling response" likely (p. 20).

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| T-a | Addendum L191; notes `stability_in_topology_optimization.md` L112 | Linearized; [K₀(k) + λK_σ(k)]φ = 0; a few wave vectors; after small-strain homogenization | **Needs correction** (substance confirmed, equation not Thomsen's) | Substance: linear homogenization and a linearized criterion (p. 2), 7 points in the biaxial case (p. 11). The equation is (K_k − P G_k)φ_k = 0 (Eq. 16), or (G_k − λK_k)φ = 0 with λ = 1/P (Eq. 19b), and K₀ means the k = 0 reduction (Eq. 15). |
| T-b | Notes L112 | "KS aggregation over (k, mode) gives the objective" | **Needs correction** for Thomsen | A min–max bound formulation (Eq. 19, p. 8). KS belongs to Wang & Sigmund 2021 and 2023. |
| T-c | Notes L112 | Sampling = "IBZ vertices plus one small-\|k\| point" | **Needs correction** (more specific) | IBZ vertices, edge midpoints and one point near Γ: 7 points in the biaxial case (Fig. 6); about a dozen for uniaxial and mixed loading (Fig. 11). |
| T-d | Addendum L211; notes L364; `instability_theory_detection.md` L144, L220 | Imposed geometric symmetry shapes failure surfaces | **Confirmed** (now from the full text) | Abstract: symmetry constraints "can be exploited to control the shape of the failure surfaces"; trade-off on p. 20. |
| T-e | `instability_theory_detection.md` L107 | Quote "homogenization theory, a linearized stability criterion and Floquet–Bloch theory" | **Confirmed** | Abstract (p. 1). |
| T-f | `stability_implementation_verification.md` L171 | Quote "homogenization, a linearized stability criterion, and Floquet-Bloch theory" | **Needs correction** (misquote) | The abstract reads "homogenization theory, a linearized stability criterion and Floquet-Bloch theory". |
| T-g | `instability_theory_detection.md` L144; notes L119 | First-order hierarchical designs, strength gains, stiffness loss | **Confirmed** | Abstract; Table 1; p. 17. |
| T-h | `instability_theory_detection.md` L114 | Thomsen et al. use a linear prebuckling approximation | **Confirmed** | p. 3 ("restriction of small strain fields prior to instability"). |
| T-i | Notes L180; `stability_implementation_verification.md` L225 (gaps) | k-point count, failure-surface resolution, mesh, CPU cost and gains unverified | **Now resolved, except two items** | Mesh, k-points and gains are given above. **CPU cost is not in the paper.** The failure-surface resolution (number of stress directions) is **not in the paper**. |
| T-j | Notes L54 | Ferrari & Sigmund cite Thomsen et al. for an interpolation that works against artificial modes | **Consistent** | Eq. 17 and the sentence on p. 7. |
| T-k | `stability_implementation_verification.md` L146 | (K + t*dK/dt)φ = 0 is "the finite-strain analogue of the linearized buckling of Thomsen et al." | **Needs correction** (loose analogy) | Thomsen's G is a linear prestress stiffness at small strain (Eq. 11). The finite-strain one-point decomposition K ≈ K_o + πK_g is in Dalklint et al. 2021 (Eqs. 18 and 21), which is the closer precedent. |
| T-l | Addendum L131 | IBZ contour sampling is "customary but not guaranteed" | **Supported by Thomsen** (not cited there) | "(but not proven)" (p. 5). |
| T-m | Addendum L207 | Every buckling reference read uses a smooth floor | **Consistent** | Eq. 17 with ρ_min = 10⁻⁴ in K. |
| T-n | NDR-0058 and NDR-0062 origins | Thomsen et al. as an origin | **Acceptable**, with a qualifier for NDR-0058 | NDR-0062: design symmetry through ρ = Lμ is tied to IBZ reduction (p. 7). NDR-0058: the paper gives a linear prebuckling analogue, not a finite-strain predictor. |
| T-o | Notes Q2 inference | Cost per iteration ≈ n_k × n_modes eigensolves; "no wall-clock data" | **Confirmed: not in the paper** | No timings in the AAM. In 2D, Thomsen et al. use 7 or more k-points and K = 25 eigenvalues. |

### Recommended edits

1. `reports/NEMOpt buckling and stability addendum.md` L191.
   - Old: "The DTU line of buckling-strength microstructure design is linearized, solving [K₀(k) + λK_σ(k)]φ = 0 at a few wave vectors after small-strain homogenization ([Thomsen, Wang & Sigmund 2018](https://doi.org/10.1016/j.cma.2018.04.031); [Wang & Sigmund 2021](https://doi.org/10.1016/j.jmps.2021.104415); verified)"
   - New: "The DTU line of buckling-strength microstructure design is linearized: after linear homogenization under a prescribed macroscopic stress it solves a Bloch-reduced linear buckling problem, (K_k − P G_k)φ_k = 0 in Thomsen et al.'s notation (their Eq. 16), at a few sampled wave vectors, seven on the irreducible-zone boundary in their biaxial case ([Thomsen, Wang & Sigmund 2018](https://doi.org/10.1016/j.cma.2018.04.031), Eqs. 10–19 and Fig. 6, verified from the DTU Orbit accepted manuscript; [Wang & Sigmund 2021](https://doi.org/10.1016/j.jmps.2021.104415); verified)"
2. Addendum L211.
   - Old: "because imposed geometric symmetry is known to shape failure surfaces ([Thomsen et al. 2018](https://doi.org/10.1016/j.cma.2018.04.031); verified at abstract level)"
   - New: "because imposed geometric symmetry is known to shape failure surfaces: fewer symmetries give more strength at the design load and less robustness to other loads ([Thomsen et al. 2018](https://doi.org/10.1016/j.cma.2018.04.031), §5.4 and §6; verified)"
3. Addendum L131, insertion after the Harrison and Maurin citation, before "Sampling therefore runs in three levels".
   - New sentence: "Thomsen et al. likewise rest boundary-only sampling on numerical studies "(but not proven)" and re-check final designs over the whole zone ([Thomsen et al. 2018](https://doi.org/10.1016/j.cma.2018.04.031), §2.1; verified)."
4. Addendum L207 (optional).
   - Old: "Every buckling reference read uses a smooth floor ([Ferrari, Sigmund & Guest 2021](https://doi.org/10.1007/s00158-021-02854-x); [Wang & Sigmund 2023](https://doi.org/10.1017/pma.2023.5); verified)"
   - New: "Every buckling reference read uses a smooth floor ([Ferrari, Sigmund & Guest 2021](https://doi.org/10.1007/s00158-021-02854-x); [Wang & Sigmund 2023](https://doi.org/10.1017/pma.2023.5); [Thomsen et al. 2018](https://doi.org/10.1016/j.cma.2018.04.031), Eq. 17; [Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z), Eq. 8; verified)"
5. Addendum NDR-0058, Origin column (L314).
   - Old: "Lindgaard & Dahl 2013; Thomsen et al. 2018; code fact (L375)"
   - New: "Lindgaard & Dahl 2013; Dalklint et al. 2021 (one-point approach, Eqs. 18 and 21); Thomsen et al. 2018 (linear prebuckling analogue); code fact (L375)"
6. `research_notes/NEMOpt buckling and stability addendum/stability_in_topology_optimization.md` L112.
   - Old: "then a linearized eigenproblem [K₀(k) + λK_σ(k)]φ = 0 with Floquet–Bloch boundary conditions, sampled at a handful of wave vectors (IBZ vertices plus one small-|k| point for long-wave modes). KS aggregation over (k, mode) gives the objective."
   - New: "then a linearized eigenproblem with Floquet–Bloch boundary conditions, (K_k − P G_k)φ_k = 0 in Thomsen et al. (Eq. 16), sampled at a handful of wave vectors (IBZ vertices, in Thomsen et al. also edge midpoints, plus one point near Γ for long-wave modes). Wang & Sigmund aggregate over (k, mode) with KS; Thomsen et al. use a min–max bound formulation (Eq. 19)."
7. Same file, L117: replace "[PR-abs] — [DTU Orbit](https://orbit.dtu.dk/en/publications/buckling-strength-topology-optimization-of-2d-periodic-materials-/)" with "[PR, accepted manuscript read in full] — [DTU Orbit AAM](https://backend.orbit.dtu.dk/ws/files/147927301/MARAC_1_s2.0_S004578251830210X_main.pdf)".
8. Same file, L180 (gap).
   - Old: "- Thomsen et al. (2018) full text could not be retrieved: the DTU Orbit PDF returned 403 or a bot page. The number and placement of k-points during optimization, the failure-surface resolution (number of stress directions), mesh sizes, CPU cost and the quantitative strength gains are therefore not verified here."
   - New: "- Thomsen et al. (2018), read from the DTU Orbit AAM: 7 k-points on the IBZ boundary for biaxial loading (Fig. 6) and about a dozen for uniaxial and mixed loading (Fig. 11); K = 25 eigenvalues and I = 10 bands; 200 × 200 Q4 elements; strength gains of 344% (square) and 379% (triangle) at bulk-modulus losses of 17–26%. The paper reports no CPU times and does not state how many stress directions trace each failure surface."
9. `research_notes/NEMOpt buckling and stability addendum/stability_implementation_verification.md` L146.
   - Old: "the finite-strain analogue of the linearized buckling of Thomsen et al. 2018"
   - New: "a finite-strain analogue of linear buckling (Thomsen et al. 2018, Eq. 10), closest to the "one point approach" of Dalklint et al. 2021, K ≈ K_o + π_jK_g (Eqs. 18 and 21)"
10. Same file, L171.
    - Old: "use "homogenization, a linearized stability criterion, and Floquet-Bloch theory""
    - New: "use "homogenization theory, a linearized stability criterion and Floquet-Bloch theory""
    - Also change "(abstract via the [DTU Orbit record](https://orbit.dtu.dk/en/publications/buckling-strength-topology-optimization-of-2d-periodic-materials-/); full text not accessible here)" to "(DTU Orbit accepted manuscript)".
11. Same file, L225: delete "Thomsen et al. 2018 (k-sampling details, verification examples, cost) and of".

---

## 2. Dalklint, Wallin & Tortorelli (2021): read in full (CC BY, 13 pp.)

### What the paper actually does

- **Energy transition, exact form.** "the deformation gradient is redefined as F = 1 + γ∇u" (p. 1754), with
  - **Eq. 9:** γ(ν̄) = H_{β̃,η̃}(ν̄^q), "where β̃ = 500 and η̃ = 0.01 replaces β and η in (6)". The argument is the **projected** filtered density ν̄ raised to the SIMP exponent q, and H is the tanh projection of Eq. 6.
  - **Eq. 10:** W(∇u; γ) = W^PH(γ∇u) + W^L(∇u) − W^L(γ∇u) = W^PH(γ∇u) + (1 − γ²)W^L(∇u), so that "γ = 0, we therefore obtain W = W^L".
  - **Eq. 8:** stiffness interpolation χ(ν̄) = ν̄^q(1 − δ_o) + δ_o with δ_o = 10⁻⁹ (a smooth floor).
- **Eigenproblem ("one point approach").** K ≈ K_o + π_jK_g near the critical point (Eq. 18); K_g = Σ∫γ²GᵀYG dV (Eq. 20); K_oφ_j = −π_jK_gφ_j with λ^c = πλ^op, at a stable operating load and for π > 1 (Eq. 21, p. 1755).
- **Why the artificial modes disappear.** "the denominator in the Rayleigh quotient (40) tends to zero as γ → 0, thereby removing the artificial modes from the interval I" (p. 1757). Portal-frame check (Table 1): conforming mesh π₁ = 6.0077 with 46 modes; SIMP ersatz π₁ = 2.3705 with 3325 modes; SIMP with energy transition π₁ = 6.0078 with 46 modes.
- **Aggregation.** A p-norm of 1/λ_i^c (Eqs. 27–28): "The nondifferentiability of degenerate critical points is handled by p-norm aggregation" (p. 1761).
- **Rejected alternatives.** On Lindgaard & Dahl: "different interpolations for different stiffness matrix contributions result in loss of convergence" (p. 1752). Matrix penalizations in Eq. 21 "deteriorate the convergence of the equilibrium analysis" (p. 1757).
- **Other facts relevant to NEMOpt.**
  - The Helmholtz filter uses a Neumann BC, which "is well-known to cause an artificial "attraction" of material to the domain boundaries" (p. 1753).
  - Filter length l = 40/(2√3) mm (§9.1, p. 1758).
  - "MMA-scheme with default parameters, cf. Svanberg (1987)" (p. 1758).
  - q is continued from 2 to 6 in steps of 0.05 every 5th iteration; β then rises from 1 to about 12 by β → 1.3β.
  - Replication-of-results statement: data in the manuscript and code "shared by contacting the corresponding author" (p. 1762).

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| D-a | Addendum L143 | Void elements without geometric stiffness cannot buckle; "the same effect that energy-interpolation schemes exploit"; "verified at abstract level" | **Confirmed**; upgrade the label | K_g carries γ² (Eq. 20); modes leave the interval as γ → 0 (p. 1757). The removal works in the load-factor form (Eq. 21), not in an identity-metric spectrum, which agrees with the addendum's own caution. |
| D-b | Addendum L17 (panel) | Artificial low-density modes need a consistent device, never separate geometric-stiffness interpolation at finite strain | **Confirmed** for the second half | p. 1752 and p. 1757 quotes above. Dalklint et al. use energy transition, not a pseudo-mass; the row's pseudo-mass rests on Zhang et al. |
| D-c | Notes L224–227; `instability_theory_detection.md` L276 | Strain-energy transition; Newton along the fundamental path; eigenvalue estimate; minimize displacement with buckling and volume constraints; Helmholtz filter; MMA | **Confirmed** | Abstract; §§4–6 and 9. |
| D-d | Notes L227 | Extends the Lindgaard–Dahl precritical approach with the Wang et al. 2014 interpolation and an aggregation function (via Zhang et al.) | **Confirmed** from the primary text | One-point approach, citing Lindgaard & Dahl (p. 1752); Wang et al. 2014 energy transition (p. 1752, Eq. 9); p-norm aggregation (Eqs. 27–28). |
| D-e | Notes L276 (gap) | Exact transition formula, eigenproblem and aggregation unknown | **Now resolved** | Eqs. 9–10, Eq. 21, Eq. 28. |
| D-f | Notes L191, L415, L433(b) | Precritical extrapolation and load-factor form; elements in void cannot enter the denominator | **Confirmed** | Eqs. 18, 21 and 40. |
| D-g | Addendum L211 | "Dalklint et al. (a conforming mesh)" as final-design verification practice | **Needs correction** (ambiguous year) | The notes attribute conforming-mesh verification of designs to Dalklint et al. 2022. The 2021 paper uses a conforming mesh only as the reference in the artificial-mode test (Table 1), not to re-analyze its optimized designs. |
| D-h | Strategy report L213 and NDR-0034 | The option Ĥ(H(ρ̃)^p) with (500, 0.01), attributed to Kumar et al. 2021 | **Confirmed, with an extra primary precedent** | Dalklint et al. 2021, Eq. 9 uses exactly (500, 0.01) on the projected density raised to q, in the Wallin–Tortorelli group. |
| D-i | NDR-0052, NDR-0060 origins | Dalklint et al. 2021 as an origin of the artificial-mode policy | **Confirmed** for NDR-0052 | For NDR-0060, Dalklint et al. 2021 is a further primary source for rejecting separate interpolation and for the load-factor form. |

### Recommended edits

1. Addendum L143.
   - Old: "([Wang, Lazarov, Sigmund & Jensen 2014](https://doi.org/10.1016/j.cma.2014.03.021); [Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z); verified at abstract level)"
   - New: "([Wang, Lazarov, Sigmund & Jensen 2014](https://doi.org/10.1016/j.cma.2014.03.021), verified at abstract level; [Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z), Eqs. 9, 10 and 20, verified: their geometric stiffness carries γ², so void modes drop out of the load-factor problem K_oφ = −πK_gφ, p. 1757)"
2. Addendum L205: after "([Lindgaard & Dahl 2013](https://doi.org/10.1007/s00158-012-0832-2); verified)", insert "; Dalklint et al. reach the same conclusion and exclude such penalizations ([Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z), pp. 1752 and 1757; verified)".
3. Addendum L211.
   - Old: "Dalklint et al. (a conforming mesh)"
   - New: "Dalklint et al. 2022 (a conforming mesh)"
4. Strategy report L213 and NDR-0034 (L462).
   - In L213, old: "Ĥ(H(ρ̃)^p) with (500, 0.01) ([Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4))"
   - L213 new: "Ĥ(H(ρ̃)^p) with (500, 0.01) ([Kumar et al. 2021](https://doi.org/10.1007/s00158-020-02764-4); [Dalklint et al. 2021, Eq. 9](https://doi.org/10.1007/s00158-021-03012-z), from Wallin and Tortorelli's group)"
   - In NDR-0034, append to the option text: "(precedents: Kumar et al. 2021; Dalklint et al. 2021, Eq. 9)".
5. Addendum NDR-0060, Origin column (L316).
   - Old: "Zhang et al. 2023; Du & Olhoff 2007; Lindgaard & Dahl 2013"
   - New: "Zhang et al. 2023; Du & Olhoff 2007; Lindgaard & Dahl 2013; Dalklint et al. 2021 (load-factor form, Eq. 21; penalizations rejected, p. 1757)"
6. `stability_in_topology_optimization.md` L223: change "[PR-abs]" to "[PR, CC BY full text]". After L224, add one bullet: "Transition: F = 1 + γ∇u, γ = H_{500,0.01}(ν̄^q) on the projected density (Eq. 9), W = W^PH(γ∇u) + (1 − γ²)W^L(∇u) (Eq. 10); eigenproblem K_oφ = −πK_gφ with K_g ∝ γ² (Eqs. 20–21); p-norm aggregation (Eqs. 27–28); Table 1: 3325 modes with plain SIMP against 46 with the transition."
7. Same file, L276: replace the gap with "Resolved from the CC BY version of record (see L223 entry)."

---

## 3. Sigmund (2022): not obtained (abstract only)

Bibliographic facts confirmed from Crossref and the Springer page: *SMO* 65(11), article 315, published online 2022-10-25, article type "Educational Paper".

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| S-a | `verification_of_refactoring.md` L134, L297; `numerical_method_choices.md` L218 | "the responsibility of fair evaluation of contributions is up to authors"; "an open question"; "lists recommendations, tools and concepts" | **Confirmed (abstract)** | Abstract: "In the lack or impracticality of common benchmarks, the responsibility of fair evaluation of contributions is up to authors". |
| S-b | Strategy report L178 | "Following Sigmund's point that fair evaluation is the authors' responsibility, final designs should be re-analyzed with one frozen evaluator, thresholded and compared on M_nd and IoU" | **Partly not verifiable** | The first clause is confirmed, but it is conditional ("In the lack or impracticality of common benchmarks"). The frozen evaluator, thresholding, M_nd and IoU cannot be checked in Sigmund 2022; M_nd itself is Sigmund 2007. |
| S-c | Strategy report L37 | Sigmund 2022 supports counting continuation in design iterations and comparing thresholded designs with a common evaluator | **Not verifiable** | Full text not obtained. |
| S-d | Strategy report L120 | Sigmund 2022 among the L4 references for the thresholded re-analysis, M_nd and IoU | **Not verifiable** | As S-b. |
| S-e | `verification_of_refactoring.md` L135 | Converting designs to binary before comparison (via secondary sources) | **Not verifiable** | Still secondary. |

### Recommended edits

1. Strategy report L178.
   - Old: "Following Sigmund's point that fair evaluation is the authors' responsibility, final designs should be re-analyzed with one frozen evaluator, thresholded and compared on M_nd and IoU ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2))."
   - New: "Sigmund notes that, without practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors" ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), abstract; full text not read). This report therefore recommends re-analyzing final designs with one frozen evaluator, thresholding them, and comparing them on M_nd ([Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x)) and IoU."
2. Strategy report L37 and L120: after each "[Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2)", add " (abstract only)".
3. `numerical_method_choices.md` L218.
   - Old: "Per a search summary of its abstract, it puts "the responsibility of fair evaluation of contributions" on the authors."
   - New: "Its abstract (Springer page, read 2026-10-12) states that, in the absence of practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors"."
4. Keep the gap statements at `numerical_method_choices.md` L243, `verification_of_refactoring.md` L160 and L395, adding: "No legal open copy exists as of 2026-10-12 (Springer paywall; DTU Orbit record without file; no arXiv)."

---

## 4. Haftka, Zhou & Queipo (2019), "Replication of results": read in full (5 pp.)

### What SMO requires, per the editorial

- **The requirement.** SMO is "recently requiring that papers have a section called Replication of Results (RR)" (p. 405). The phrase means "a detailed enough description of the proposed model/method so that the results can be replicated to a degree consistent with the claims and conclusions drawn in the paper (Plesser 2018)" (p. 405).
- **No fixed content.** The editorial prescribes no mandatory checklist. It presents "options available to authors" (p. 405). Its section 3 lists:
  - source code, called "the best option", ideally in a version-controlled repository with a tagged version and a DOI "through, for example, Zenodo.org or figshare.com" (pp. 407–408);
  - executable code;
  - pseudocode that mentions "all tuning and convergence parameters";
  - input files and scripts for commercial codes;
  - digital data for plots;
  - detailed intermediate results (p. 409).
- **Related guidance.**
  - "detailed settings for results need to be documented" (p. 406).
  - Fix random seeds, and document the operating system, configuration and dependencies, or use virtual machines or containers (p. 406).
  - "authors should make repeated runs in order to determine the scatter in the results" (p. 407).
  - The randomness of operating systems "may often be imitated by small perturbation in loading, finite element meshes, or material properties" (p. 407).
  - Put tuning and convergence parameters at the start of the code, as in the 99-line code (p. 408).
- **Terminology.** "There is no consensus on terminology" (p. 405), citing Plesser 2018. Also: "Almost all the papers in SMO are about computational algorithms" (p. 405).

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| H-a | Strategy report L374 | A mandatory "Replication of results" section in SMO | **Confirmed** | p. 405. Its content is optional and illustrative, not prescribed. |
| H-b | Strategy report L42 (panel) | Pin and archive the environment, give releases a DOI, and define "reproduced" explicitly | **Confirmed** for the DOI, environment and degree-of-replication points | pp. 405–408. Pre-registering decisions is **not in the paper**; that point rests on the other sources in the row. |
| H-c | `verification_of_refactoring.md` L296 | "detailed enough" description consistent with the claims; Plesser on terminology; almost all SMO papers concern algorithms | **Confirmed** | p. 405. |
| H-d | `verification_of_refactoring.md` L221 | The editorial cites Plesser on the lack of consensus | **Confirmed** | p. 405. |
| H-e | `verification_of_refactoring.md` L296 | Quote "It is desirable that papers, especially on new methodology, should help readers reproduce results presented in the paper" | **Not in the paper** | This wording is not in the editorial; it is attributed to the submission guidelines via a search snippet and remains unverified. |
| H-f | `verification_of_refactoring.md` L345 | The RR section can point to the release DOI, parameter table, deviation register and verification report | **Consistent** with the editorial's options | pp. 407–409. |

### Recommended edits

1. Strategy report L374.
   - Old: "the mandatory "Replication of results" section of *Structural and Multidisciplinary Optimization* ([Haftka et al. 2019](https://doi.org/10.1007/s00158-019-02298-4));"
   - New: "the "Replication of results" section that *Structural and Multidisciplinary Optimization* requires ([Haftka et al. 2019](https://doi.org/10.1007/s00158-019-02298-4), p. 405). The editorial fixes no content but recommends source code with a tagged version and DOI, all tuning and convergence parameters, plot data and intermediate results, and repeated runs to show scatter (pp. 406–409);"
2. Strategy report NDR-0010 (L438).
   - Old: "Full-run acceptance: a legacy ensemble (size, perturbations, and metrics J, ν, κ, Vf, M_nd, thresholded ν, IoU) plus the article values at printed precision"
   - New: "Full-run acceptance: a legacy ensemble (size, perturbations, and metrics J, ν, κ, Vf, M_nd, thresholded ν, IoU), following SMO's advice to make repeated, slightly perturbed runs to show scatter (Haftka et al. 2019, p. 407), plus the article values at printed precision"
3. Strategy report NDR-0044 (L472): append "(SMO: "detailed settings for results need to be documented", Haftka et al. 2019, p. 406)".
4. `verification_of_refactoring.md` L296: keep the guideline quote, flagged as unverified. Change "(wording via search snippet; direct fetch blocked by Springer authentication)" to "(guideline wording via search snippet, unverified; the editorial was read in full from the publisher's free PDF on 2026-10-12 and does not contain this sentence)".

---

## 5. Svanberg (2007), "MMA and GCMMA – two methods for nonlinear optimization": not obtained

### Access

- `https://people.kth.se/~krille/mmagcmma.pdf` returns 404 with "This page is no longer in use" (confirms `numerical_method_choices.md` L281).
- The KTH profile (`kth.se/profile/krille`) lists contact details only. `smoptit.se` shows "Account Suspended".
- Wayback snapshots exist, per the archive.org availability API:
  - `https://web.archive.org/web/20210204193223/https://people.kth.se/~krille/mmagcmma.pdf`
  - `https://web.archive.org/web/20220120004729/https://people.kth.se/~krille/mmagcmma.pdf`
  - `https://web.archive.org/web/20250316062505/https://people.kth.se/~krille/mmagcmma.pdf`
  - `https://web.archive.org/web/20230603101120/https://people.kth.se/~krille/gcmma07.pdf`, titled "MMA and GCMMA, versions September 2007" per a search index.
- Every TLS connection to `web.archive.org` and `wayback.archive.org` was reset right after the handshake, and the fetch tool refuses that host. The user needs to open one of these snapshots or allow `web.archive.org` in the environment's network settings.

### What can be checked instead: the GPL Python port the report actually cites

This source is `GCMMA-MMA-Python`, `src/mmapy/mma.py` and `examples/gcmma_beam.py` on master, fetched as text and not executed. Its README says the MATLAB original was "taken from http://www.smoptit.se/ under the GNU General Public License" and cites `mmagcmma.pdf` as retrieved in 2017.

- **mmasub defaults (L37–39):** `move 0.5`, `asyinit 0.5`, `asydecr 0.7`, `asyincr 1.2`, `asymin 0.01`, `asymax 10`, `raa0 0.00001`, `albefa 0.1`.
- **epsimin (L96):** mmasub hard-codes `epsimin = 0.0000001`. The same value ends the interior-point subproblem solver (`while epsi > epsimin`, L351).
- **concheck (L677–685):** conservative only if `f_app + epsimin >= f_new` for the objective and all constraints.
- **raaupdate (L637–643):** uses the slack `0.5*epsimin`; raa0 becomes `min(1.1*(raa0 + Δ), 10*raa0)`.
- **GCMMA example driver:** sets `epsimin = 1e-7` and `raa0eps = 1e-6`, and caps inner iterations with `while conserv == 0 and innerit <= 15`.

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| V-a | Strategy report L208; `numerical_method_choices.md` L25, L252, L262, L285 | Svanberg's GCMMA conservativeness tolerance is 1e−7 | **Confirmed for the Python port only**; not verifiable in Svanberg 2007 | In the port, 1e−7 is `epsimin`, primarily the subproblem tolerance, reused as the concheck slack; it is set by the GCMMA driver, not fixed inside `concheck`. |
| V-b | `numerical_method_choices.md` L260–263 | Asymptote rules and defaults: asyinit 0.5, asyincr 1.2, asydecr 0.7, bounds [0.01, 10]·(xmax − xmin), move 0.5, albefa 0.1; raa update rule | **Confirmed for the port**; not verifiable in Svanberg 2007 | mma.py L37–39, L103–115, L637–643. |
| V-c | Task brief | Whether Svanberg allows changing parameters between iterations | **Not verifiable** | The port's API accepts the parameters on every call (L37–39), which shows that it is technically possible, not what Svanberg advises. The notes' "no source verifies resetting asymptotes after continuation" (L244) stands. |
| V-d | `numerical_method_choices.md` L298 | 2002 and 2007 statements on differentiability and the inner-loop tolerance not re-read | **Still open** | — |

### Recommended edits

1. Strategy report L208.
   - Old: "far above the 1e−7 conservativeness tolerance in Svanberg's GCMMA code; NonconvexMMA's check has no tolerance at all ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python))"
   - New: "far above the slack epsimin = 1e−7 that the GPL Python port of Svanberg's MATLAB code adds in its GCMMA conservativeness check (also the port's subproblem tolerance; Svanberg's 2007 notes not re-read); NonconvexMMA's check has no tolerance at all ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python), `concheck` and `examples/gcmma_beam.py`)"
2. `numerical_method_choices.md` L25.
   - Old: "far larger than GCMMA's conservativeness tolerance (1e-7)."
   - New: "far larger than the 1e-7 slack (epsimin) in the conservativeness check of the Python port of Svanberg's GCMMA."
3. `numerical_method_choices.md` L252.
   - Old: "within ~1e-7 in Svanberg's code"
   - New: "within epsimin = 1e-7 in the Python port of Svanberg's code (not checked against Svanberg 2007)"
4. `numerical_method_choices.md` L281: append "Wayback snapshots exist (20210204193223, 20220120004729, 20250316062505; gcmma07.pdf 20230603101120) but `web.archive.org` was unreachable from the research container on 2026-10-12."
5. Do not cite Svanberg 2007 for any parameter value until a snapshot has been read.

---

## 6. Lazarov & Sigmund (2011): not obtained (abstract only); secondary source used

Abstract (DTU Orbit record): the filtered variable is defined "implicitly as a solution of a Helmholtz‐type differential equation with homogeneous Neumann boundary conditions".

Secondary source by the same group, read in full where relevant: Lazarov, Wang & Sigmund 2016, DTU Orbit postprint (cited by printed page and equation):
- PDE filter (Eq. 9). "The filter is volume preserving for Neumann boundary conditions" (p. 5).
- r_f = R_f/(2√3) (Eq. 10). It is obtained "by matching the second central moments of the Hat function and the Green's function for the PDE filter on unbounded domain [93]" (pp. 5–6); ref. [93] is Lazarov & Sigmund 2011 (p. 40).

Further corroboration: Dalklint et al. 2021 use l = 40/(2√3) mm for the Helmholtz length (p. 1758) and note the Neumann BC's "attraction" of material to boundaries (p. 1753). Periodic volume conservation is stated by Wallin & Tortorelli 2020 (p. 9; quoted and verified in `wallin_tortorelli_2020.md` L280).

### Claims and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| L-a | Strategy report NDR-0039 (L467); `numerical_method_choices.md` L28, L420, L442; `efficient_fe_adjoint.md` L208 | r = R/(2√3) (cone R ≈ 2√3·r) | **Confirmed via the 2016 review** (not read in the 2011 paper) | Lazarov et al. 2016, Eq. 10. It is a second-moment match on an unbounded domain, so it is approximate near boundaries. |
| L-b | `numerical_method_choices.md` L431 | "R_min = r_min/(2√3)" | **Needs correction** (letters reversed against L420) | top82 writes `Rmin = rmin/2/sqrt(3)` with Rmin as the PDE length. Use r_PDE = R_cone/(2√3) throughout. |
| L-c | `numerical_method_choices.md` L420 | r = 1 element corresponds to R ≈ 3.5 elements | **Confirmed** (2√3 = 3.46) | — |
| L-d | `numerical_method_choices.md` L437; `efficient_fe_adjoint.md` L221, L227 | Volume preserved with Neumann BCs, and also with periodic BCs; the real difference is boundary distortion | **Confirmed** | Neumann: Lazarov et al. 2016, p. 5. Periodic: Wallin & Tortorelli 2020, p. 9. Boundary artifact: Dalklint et al. 2021, p. 1753. L227's "my inference" can now cite Wallin & Tortorelli. |
| L-e | Strategy report L268 | Lazarov & Sigmund 2011 support "assemble and factor the Helmholtz filter once per run" | **Not verifiable** in the 2011 paper | top88 (Andreassen et al. 2011) is the verified source. Lazarov et al. 2016 (p. 5) mention "re-utilization of solvers". |
| L-f | `efficient_fe_adjoint.md` L208, L231 | Full text not accessible | **Still true** | No Orbit file; Wiley 403; closed per Semantic Scholar. |

### Recommended edits

1. Strategy report NDR-0039 (L467).
   - Old: "(Helmholtz r versus cone R ≈ 2√3·r)"
   - New: "(Helmholtz r versus cone R ≈ 2√3·r, a second-moment match on an unbounded domain: [Lazarov, Wang & Sigmund 2016, Eq. 10](https://backend.orbit.dtu.dk/ws/files/140201797/main.pdf), attributing it to Lazarov & Sigmund 2011)"
2. Strategy report L268.
   - Old: "[Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072)"
   - New: "[Lazarov & Sigmund 2011](https://doi.org/10.1002/nme.3072) (filter definition; full text not read)"
3. `numerical_method_choices.md` L420.
   - Old: "The DTU convention links the PDE length r to the classic cone or density filter radius R by r = R/(2√3). NEMOpt's r = 1 element therefore corresponds to R ≈ 3.5 elements."
   - New: "The DTU convention links the PDE length r to the classic cone or density filter radius R by r = R/(2√3), from matching second central moments on an unbounded domain (Lazarov, Wang & Sigmund 2016, Eq. 10, citing Lazarov & Sigmund 2011). NEMOpt's r = 1 element therefore corresponds to R ≈ 3.46 elements."
4. `numerical_method_choices.md` L431.
   - Old: "The source of the approximate R_min = r_min/(2√3) relation, per secondary sources."
   - New: "The source of the approximate relation r = R/(2√3) (PDE length r, cone radius R; top82 writes `Rmin = rmin/2/sqrt(3)` with the letters reversed), per Lazarov, Wang & Sigmund 2016 (Eq. 10, ref. [93]); the 2011 full text has no legal open copy."
5. `efficient_fe_adjoint.md` L227.
   - Old: "This is my inference from the Neumann volume-preservation statement plus periodicity."
   - New: "Wallin & Tortorelli 2020 state this explicitly (p. 9: periodicity "also conserves volume"); Dalklint et al. 2021 (p. 1753) note that Neumann BCs attract material to the boundaries."

---

## Access log (failures worth knowing)

| Host or URL | Result | Note |
|---|---|---|
| `orbit.dtu.dk/files/...` and Orbit search pages | 403 (Cloudflare challenge) | Record pages load; use `backend.orbit.dtu.dk/ws/files/...` for files. |
| `link.springer.com` with a browser user agent | JavaScript "Client Challenge" | A plain client that keeps cookies through the `idp.springer.com` redirect is served normally. Open-access and free PDFs downloaded; paywalled PDFs redirect to the abstract page. |
| `portal.research.lu.se` search | 403 | LUP (`lup.lub.lu.se`) works; the Dalklint record has no file. The thesis file on LUP is "utan papers" (without the papers). |
| `onlinelibrary.wiley.com` | 403 (Cloudflare) | — |
| `web.archive.org`, `wayback.archive.org`, `web-beta.archive.org` | Connection reset after the TLS ClientHello | The egress proxy logs "tunnel closed"; `archive.org/wayback/available` works (with occasional 429s). |
| `kth.diva-portal.org` | Connection reset | Same pattern. |
| api.unpaywall.org; api.openalex.org | Rejected a placeholder email; daily budget exhausted | Not needed in the end. The user's email was not sent. |

## Final table

| # | Paper | Status | Source URL |
|---|---|---|---|
| 1 | Thomsen, Wang & Sigmund 2018 | **Obtained** (accepted manuscript) | https://backend.orbit.dtu.dk/ws/files/147927301/MARAC_1_s2.0_S004578251830210X_main.pdf |
| 2 | Dalklint, Wallin & Tortorelli 2021 | **Obtained** (CC BY version of record) | https://link.springer.com/content/pdf/10.1007/s00158-021-03012-z.pdf |
| 3 | Sigmund 2022 | **Needs the user** (no legal open copy; abstract only) | https://link.springer.com/article/10.1007/s00158-022-03427-2 (paywalled) |
| 4 | Haftka, Zhou & Queipo 2019 | **Obtained** (publisher free access) | https://link.springer.com/content/pdf/10.1007/s00158-019-02298-4.pdf |
| 5 | Svanberg 2007 notes | **Needs the user** (a Wayback copy exists; host unreachable here) | https://web.archive.org/web/20220120004729/https://people.kth.se/~krille/mmagcmma.pdf |
| 6 | Lazarov & Sigmund 2011 | **Needs the user** (no legal open copy; abstract only; checked via the 2016 review) | https://orbit.dtu.dk/en/publications/filters-in-topology-optimization-based-on-helmholtztype-different/ (record without file) |

---
*Update after the user supplied the papers (coordinator, 2026-10-12).* The text above is kept as first written.
- **§6, Lazarov & Sigmund 2011.** The user supplied the version of record, and it was read in full ([`lazarov_sigmund_2011_kubo_matsuda_2016.md`](lazarov_sigmund_2011_kubo_matsuda_2016.md), §1). L-a and L-d are now verified in the primary text (Eq. 34; Eqs. 12–14). L-e is confirmed (p. 771: factorization reuse for relatively small problems). L-f is resolved. Final table, row 6: read "**Supplied by the user** (version of record)". The report edits proposed in §6 were superseded by primary citations.
