# Instability theory and detection in finitely strained periodic microstructures (first-order periodic homogenization), with implications for NEMOpt's unit-cell and half-cell models

Scope note. These notes cover theory and detection: which instabilities exist, how each is detected, and what each implies for NEMOpt's one-cell model, its point-symmetric half-cell model and its hinge-like auxetic designs. Notation used throughout:

- F(t) = I + t(F_M − I) is the load ramp NEMOpt already applies (`A(t) = (t/tFinal)*(F_M - I)` in `Full_Domain/NEMOpt.jl` L373). Here t plays the role of the load parameter λ in the literature.
- K_red is NEMOpt's reduced tangent: periodic constraints, center node fixed (`K_Reduced` in `Full_Domain/NEMOpt.jl` L340–368, currently solved with `\`).
- 𝔸_M = ∂P_M/∂F_M is the homogenized nominal (first) elasticity tensor from the four χ problems.
- k is the Bloch wave vector; N and m are unit vectors.

Every bibliographic record below (title, authors, venue, volume, pages, DOI) was checked against OpenAlex or Crossref in this session. Where a source's content could only be checked at title level, or comes from a search-engine summary or a secondary source, the entry is flagged **[content not re-verified]**. Statements marked "derived" are short derivations of mine, not literature claims.

## Q1. Microscopic vs. macroscopic instability: definitions, relation, and where the one-cell "principal solution" stops being meaningful

### Takeaway
Take an infinite periodic solid under prescribed F_M. The principal (fundamental) solution, in which all cells deform identically, is stable only while Hill's second-variation functional stays positive for Bloch-wave perturbations at every wave vector k in the Brillouin zone. That set covers three cases:

- k = 0: one-cell-periodic modes.
- k → 0: long wavelength.
- k ≠ 0: modes that repeat over several cells.

Loss of rank-one convexity, i.e. strong ellipticity, of the homogenized tangent is exactly the k → 0 limit. Microscopic instability therefore precedes macroscopic instability or coincides with it. Past the first instability of any type, one-cell homogenized quantities (P_M, 𝔸_M, ν, κ, μ) describe an unstable equilibrium and lose their physical meaning.

### Cited Findings
- **Earliest comparison (fiber composites).** Triantafyllidis & Maker studied an axially compressed, fiber-reinforced, plane-strain composite. The microscopic instability is bifurcation buckling of the fibers. The macroscopic instability is loss of ellipticity of the homogenized incremental equilibrium equations. The macroscopic approach "consistently and sometimes considerably overestimates the stable region of the composite." — [Triantafyllidis & Maker 1985, J. Appl. Mech. 52(4):794–800](https://doi.org/10.1115/1.3169148)
- **Rigorous general theory.** Geymonat, Müller & Triantafyllidis, "Homogenization of nonlinearly elastic materials, microscopic bifurcation and macroscopic loss of rank-one convexity," Arch. Ration. Mech. Anal. 122(3):231–290 (1993) — [doi:10.1007/BF00380256](https://doi.org/10.1007/BF00380256). The full text was not accessible this session (the University of Michigan Deep Blue repository copy sits behind a Cloudflare challenge). Its content is therefore taken from the detailed secondary summary in Zhang, Feng & Khandelwal 2021 (next items).
  - It gives a rigorous proof that a Bloch-wave representation can be used for the microscale stability of general 3D periodic solids of infinite extent.
  - It shows that "the buckling mode of infinite wavelength at microscale can be equivalently detected as a loss of rank-1 convexity of the homogenized tangent moduli at macroscale."
  - It shows that "the homogenization results are only valid before the bifurcation occurs." — [Zhang, Feng & Khandelwal 2021, IJNME 122:6527–6575](https://doi.org/10.1002/nme.6802) (preprint read in full: [arXiv:2010.02371](https://arxiv.org/abs/2010.02371))
- **Inequality chain from Geymonat et al.** As restated by Zhang et al. (Eq. 55): β(λ) ≤ β(λ; k → 0) ≤ B(λ).
  - β is the infimum over k of the Bloch-wave Rayleigh quotient ∫∇v:𝔸:∇v / ∫∇v:∇v.
  - B(λ) = min over |m| = |N| = 1 of (m⊗N):𝔸̄:(m⊗N) is the macroscopic ellipticity indicator.
  - Consequence: if microscopic stability holds (β > 0), the homogenized material is rank-one convex.
  - Eq. 56: B(λ) = 0 if β(λ; k → 0) = 0, i.e. loss of rank-one convexity *is* long-wavelength microscale buckling. — [Zhang et al. 2021 (arXiv:2010.02371), §4.2](https://arxiv.org/abs/2010.02371)
- **Order of onset.** "The micro-instability (short wavelength buckling) occurs either before or simultaneously with the macro-instability (long wavelength buckling)." — [Zhang et al. 2021, §4](https://arxiv.org/abs/2010.02371) (citing Geymonat et al. 1993 and Triantafyllidis et al. 2006)
- **Two different modes near k = 0.**
  - A short-wavelength mode with k = 0 has the same periodicity as the cell. The infimum is attained at k = 0, and the stability surface β(k) is continuous at the origin.
  - A long-wavelength mode is approached only as k → 0. The surface is discontinuous at the origin there, so a computation at exactly k = 0 does not see it.
  - Source: [Zhang et al. 2021, §4.1](https://arxiv.org/abs/2010.02371)
- **Implicit assumption of homogenization.** Homogenization assumes that one cell stays representative: "the principal loading path is stable, and no other bifurcated solution besides the principal solution exists." — [Zhang et al. 2021, §1](https://arxiv.org/abs/2010.02371)
- **Onset of failure in periodic solids.** For periodic solids of infinite extent, the onset of failure can be defined as "the first occurrence of a bifurcation in the fundamental solution, for which all cells deform identically." The calculations are "based on a Bloch wave analysis on the deformed unit cell." — [Triantafyllidis, Nestorović & Schraad 2006 (online 2005), J. Appl. Mech. 73(3):505–515](https://doi.org/10.1115/1.2126695)
- **Matrix ellipticity does not carry over.** A homogenized periodic porous solid can lose strong ellipticity at large enough loads even when the matrix material is strongly elliptic. Polyconvex constituents likewise do not guarantee a rank-one convex homogenized material. — [Zhang et al. 2021](https://arxiv.org/abs/2010.02371), summarizing [Abeyaratne & Triantafyllidis 1984, J. Appl. Mech. 51(3):481–486](https://doi.org/10.1115/1.3167661)
- **Loss of rank-one convexity need not localize.** Whether a localized (discontinuous) pattern appears on the bifurcated branch depends on the stability of that branch. — [Zhang et al. 2021, §4.2](https://arxiv.org/abs/2010.02371) (their ref. [56])
- **Related theory, bibliographic only.** The following were verified at the bibliographic level only **[content not re-verified]**:
  - Müller's multi-cell homogenization formula for nonconvex periodic energies (the one-cell formula is not the true homogenized energy once nonconvexity matters): [Müller 1987, Arch. Ration. Mech. Anal. 99:189–212](https://doi.org/10.1007/BF00284506)
  - Triantafyllidis & Schnaidt's 2D comparison of microscopic and macroscopic instabilities: [JMPS 41:1533–1565 (1993)](https://doi.org/10.1016/0022-5096(93)90039-I)
  - Scale-size effects on finite periodic solids: [Triantafyllidis & Bardenhagen 1996, JMPS 44:1891–1928](https://doi.org/10.1016/0022-5096(96)00047-6)
  - Hill's uniqueness and stability criterion: [Hill 1957, JMPS 5:229–241](https://doi.org/10.1016/0022-5096(57)90016-9)
  - The analogous "phonon vs. Cauchy–Born" micro/macro split in crystals: [Elliott, Triantafyllidis & Shaw 2006, JMPS 54:161–192](https://doi.org/10.1016/j.jmps.2005.07.009)
- **NEMOpt's own framing.** The NEMOpt article computes 𝔸_M as "the linear response (first order linear approximation of the response at a deformed state) of the unit cell at the state imposed by F_M." It names "instabilities that might arise particularly under compressive macroscopic loading" as open work. — [Güven & Özdemir 2026, FEAD 260:104599, §4 and §6](https://doi.org/10.1016/j.finel.2026.104599)

### Inferences
- **Three nested checks.** For NEMOpt there are three distinct checks, and passing one does not imply passing another:
  1. **k = 0, one-cell periodic.** This is exactly the positive definiteness of NEMOpt's K_red (periodic fluctuations, center node fixed to remove translation).
  2. **k → 0, macroscopic.** Rank-one convexity of 𝔸_M.
  3. **k ≠ 0, multi-cell.** A Bloch-wave check.

  K_red being positive definite proves only check 1.
- **When the one-cell results stop being physical.** The one-cell principal solution stops being physically meaningful at the first load factor t_c ∈ (0, 1] along the ramp F(t) where any of the three checks fails. That one-cell 𝔸_M, ν, κ and μ still converge with Newton beyond t_c (Newton happily converges to unstable equilibria) does not mean they are valid.
- **Instability before t = 1.** An instability can occur at t_c < 1 on the ramp even if the state at t = 1 looks benign. Stability should therefore be tracked at every converged increment, not only at F_M.
- **Ordering in the optimization.** Because the macroscopic check is an upper bound on the microscopic one (β ≤ B), passing the cheap rank-one check of 𝔸_M is necessary but not sufficient. An optimizer that only enforces B > 0 can still produce designs that buckle microscopically with k ≠ 0.
- **Constituent vs. structural instability.** The compressive neo-Hookean constituent is far from losing ellipticity at the published strain levels (≤ 7.5% stretch, ~5% principal compression). Any instability in NEMOpt designs would therefore be structural or geometric (slender ligaments, near-mechanisms), not a constituent-level material instability. **[unverified for NEMOpt's specific volumetric function]**

### Gaps
- The Geymonat et al. (1993) full text was not read. The precise hypotheses (regularity, convexity of constituents, and the exact relation between k = 0 and k → 0 infima) are relayed via Zhang et al. (2021).
- The exact statement of Müller's (1987) multi-cell formula and its applicability to compressible neo-Hookean cells with voids was not re-verified.
- I found no source quantifying how far beyond t_c the one-cell 𝔸_M of a hinge-type auxetic deviates from the true (multi-cell or relaxed) homogenized response.

## Q2. Computational approaches to micro and macro instability in finite-strain homogenization

### Takeaway
The established workflow (Triantafyllidis et al. 2006; Bertoldi & Boyce 2008; Zhang et al. 2021) has four steps:

1. Follow the principal one-cell solution with Newton.
2. At each converged state, solve the Bloch-wave eigenproblem of the exact tangent on the deformed cell, scanning wave vectors.
3. Separately, check rank-one convexity of the homogenized nominal tangent at each step.
4. Retrieve real buckling modes from complex Bloch eigenvectors.

Wave-vector sampling and the treatment of the k → 0 singularity dominate cost and reliability.

### Cited Findings
- **Bloch analysis in standard FE codes.** Åberg & Gudmundson proposed a method by which "standard finite element programs can be used to compute dispersion relations in periodic composites." This is the basis of Bloch-type constraints in standard FE. — [Åberg & Gudmundson 1997, J. Acoust. Soc. Am. 102(4):2007–2013](https://doi.org/10.1121/1.419652)
- **Failure surfaces.** Triantafyllidis, Nestorović & Schraad:
  - built onset-of-failure surfaces in strain or stress space by following all loading paths, with a Bloch analysis on the deformed cell;
  - described "a numerical algorithm which reduces the size of stability matrices by an order of magnitude" for continuum unit cells;
  - applied it to porous and particle-reinforced hyperelastic solids, reporting how the wavelength of the instability depends on microgeometry and macroscopic loading.

  Source: [Triantafyllidis et al. 2006, J. Appl. Mech. 73:505–515](https://doi.org/10.1115/1.2126695)
- **Zhang, Feng & Khandelwal's framework.** Source: [Zhang et al. 2021, IJNME](https://doi.org/10.1002/nme.6802); [arXiv:2010.02371 §4–5](https://arxiv.org/abs/2010.02371)
  - Micro (Bloch) and macro (rank-one) checks within one consistent first-order finite-strain homogenization framework.
  - Three solution methods for the constrained Bloch eigenproblem: two condensation methods (one following Triantafyllidis et al.) and a null-space projection.
  - Retrieval of real-valued buckling modes from complex Bloch waves.
  - Wave-vector space: any primitive cell of reciprocal space suffices; they use k₁, k₂ ∈ [0, 1).
  - For the first bifurcation, it suffices to search for a zero or negative eigenvalue.
  - In their examples:
    - the k-mesh is 100 × 100 over [0, 1)², plus 100 × 100 refinements in the strips (0, 0.01) × (0.01, 1] and (0.01, 1] × (0, 0.01) and the corner (0, 0.01)²;
    - B(λ) is checked "at each loading step";
    - the results are invariant to the choice of RVE.
  - Their examples include a one-cell (1 × 1) critical mode at λ_c = 0.9934. They also include a long-wavelength case under constrained compression in which the Bloch indicator β and the rank-one indicator B vanish together between stretches 0.98244 and 0.98243. Both serve as verification benchmarks.
  - The micro check "is much more computationally demanding, since the length scale of buckling mode is not a priori known." Long-wavelength buckling "can be more efficiently" detected by the rank-one check than by Bloch analysis.
- **Bertoldi & Boyce: band structures under deformation (PRB 77).** They computed full band structures of periodic elastomeric solids at different deformation levels. Under axial compression, the structures "undergo a transformation in their patterned structure upon reaching a critical value of applied load," which closes some band gaps and opens others. — [Bertoldi & Boyce 2008a, PRB 77:052105](https://doi.org/10.1103/PhysRevB.77.052105)
- **Bertoldi & Boyce: waves and instabilities (PRB 78).** They gave analytical and numerical techniques for wave propagation and instabilities in finitely deformed monolithic and periodic elastomers. "Both elastic instabilities and propagation of acoustic waves are strongly influenced by geometric pattern, material properties and loading conditions." — [Bertoldi & Boyce 2008b, PRB 78:184107](https://doi.org/10.1103/PhysRevB.78.184107)
  - **[content not re-verified]** That paper also relates vanishing long-wavelength wave speed to macroscopic loss of ellipticity, and zero-frequency Bloch modes at finite k to microscopic instability.
- **Pattern transformation as an instability.** Pattern transformation in periodic elastomeric cellular solids under uniaxial compression "is triggered by a reversible elastic instability." — [Mullin, Deschanel, Bertoldi & Boyce 2007, PRL 99:084301](https://doi.org/10.1103/PhysRevLett.99.084301); see also [Bertoldi, Boyce, Deschanel, Prange & Mullin 2008, JMPS 56:2642–2668](https://doi.org/10.1016/j.jmps.2008.03.006) **[content not re-verified]**
- **Miehe, Schröder & Becker (2002).** "Computational homogenization analysis in finite elasticity: material and structural instabilities on the micro- and macro-scales of periodic composites and their interaction," CMAME 191(44):4971–5005 — [doi:10.1016/S0045-7825(02)00391-2](https://doi.org/10.1016/S0045-7825(02)00391-2). The abstract and full text were not retrievable (ScienceDirect returned 403; no abstract in Crossref, OpenAlex or Semantic Scholar) **[content not re-verified]**.
  - A later record of Miehe's group frames the dichotomy as follows. Macroscopic *material* instabilities relate to rank-one convexity of the homogenized energy. Microscopic *structural* instabilities relate to coercivity of the microscopic boundary-value problem. — [Univ. Stuttgart PUMA record (later Miehe paper, as summarized by search)](https://puma.ub.uni-stuttgart.de/bibtex/24129e0def09499c1e33f4df6c7edfdda/simtech) **[content not re-verified]**
- **Michel et al. (2007).** "Microscopic and macroscopic instabilities in finitely strained porous elastomers," JMPS 55(5):900–938 — [doi:10.1016/j.jmps.2006.11.006](https://doi.org/10.1016/j.jmps.2006.11.006). The abstract was not retrievable **[content not re-verified]**. As I recall it, the paper compares periodic FE Bloch analyses with second-order homogenization estimates of macroscopic loss of ellipticity.
- **Linearized bifurcation in topology optimization.** Topology optimization of 2D periodic materials for buckling strength has used "homogenization theory, a linearized stability criterion and Floquet–Bloch theory." — [Thomsen, Wang & Sigmund 2018, CMAME 339:115–136](https://doi.org/10.1016/j.cma.2018.04.031)
- **Wave-vector sampling heuristic.** "The critical buckling mode can be captured by sweeping k-vectors along the boundaries of the irreducible Brillouin zone (IBZ) ... determined by the shared symmetries between the microstructure geometry and the macroscopic stress state." — [Wang & Sigmund 2023, Programmable Materials 1 (doi:10.1017/pma.2023.5)](https://doi.org/10.1017/pma.2023.5); text read at [arXiv:2210.00003](https://arxiv.org/abs/2210.00003). The same authors use IBZ-boundary sweeps for 3D designs ([arXiv:2012.01359](https://arxiv.org/abs/2012.01359)).
- **Caveat on IBZ contours.** Maurin et al. studied the "probability that a band-gap extremum is located on the irreducible Brillouin-zone contour" for the 17 plane crystallographic lattices. This implies that extrema need not lie on the contour. — [Maurin, Claeys, Deckers & Desmet 2018, IJSS 135:26–36](https://doi.org/10.1016/j.ijsolstr.2017.11.006) **[content not re-verified; title-level]**
- **Finite strain in topology optimization.** Hyperelastic phononic topology optimization has applied a finite macroscopic deformation to a unit cell with Floquet–Bloch conditions, i.e. Bloch analysis on the pre-deformed cell inside topology optimization. — [Dalklint, Wallin, Bertoldi & Tortorelli 2022, JMPS 163:104849](https://doi.org/10.1016/j.jmps.2022.104849) ([OSTI accepted manuscript](https://www.osti.gov/pages/servlets/purl/1881619))

### Inferences
- **Bloch operator from NEMOpt's existing machinery.** The Bloch stability operator at k uses the same element tangent as K_red at the converged state. Only the dependent-node constraint changes: in NEMOpt's `D1`-type dependency map, each periodic image is multiplied by the phase e^{i k·L}, and the operator is assembled as D(k)ᴴ K D(k). This acts on perturbations of the *increment*, not on the principal solution. A Bloch check can therefore reuse NEMOpt's assembly and partitioning, with complex Hermitian arithmetic.
- **Exact tangent, not a linearized prebuckling.** For finite-strain neo-Hookean cells, the relevant operator is the exact tangent at the converged nonlinear state (Zhang et al.; Triantafyllidis et al.). The linear prebuckling ("linearized bifurcation") approximation used in linear-elastic buckling topology optimization (Thomsen et al. 2018) is unnecessary because NEMOpt already computes the nonlinear state.
- **Resolve k → 0 separately.** The k → 0 singularity means a uniform k-grid can miss long-wavelength modes. The rank-one check of 𝔸_M covers that limit exactly and cheaply, so the k-grid only needs to resolve finite k.

### Gaps
- The abstracts and full texts of Miehe et al. (2002), Michel et al. (2007), Bertoldi et al. (2008, JMPS) and Overvelde & Bertoldi (2014) could not be retrieved (Elsevier paywall and 403s). Their specific algorithmic choices (e.g., n × n "ensemble" cells vs. Bloch analysis in Miehe et al.) are unverified here.
- I found no published wall-clock cost comparison of Bloch scans vs. Newton solves for topology-optimization-scale 2D cells (~10⁴–10⁵ DOF).

## Q3. Instabilities in auxetic, hinge-dominated and rotating-unit metamaterials, including whether tensile or tension–shear loading can trigger them

### Takeaway
In soft periodic porous solids, the canonical negative Poisson's ratio under compression is itself a post-bifurcation, pattern-transformation phenomenon (Mullin et al. 2007; Bertoldi et al. 2010). A one-cell, principal-branch topology optimizer can therefore neither represent it nor be trusted near it.

Hinge-dominated designs with ν ≈ −0.9 are near-mechanisms: in an isotropic idealization, κ/μ ≈ 0.05. That makes them prone to soft modes, including k ≠ 0 modes, and to size and boundary effects.

The published equibiaxial-tension case is low-risk. The tension–shear and bulk "shear" cases contain genuine compressive principal directions or compressive lateral stresses, so instability cannot be ruled out without a check.

### Cited Findings
- **Instability-induced negative Poisson's ratio.** Negative Poisson's ratio in cellular solids made of a matrix with a square array of circular voids "results from an elastic instability, which induces a pattern transformation," with "excellent quantitative agreement" between calculation and experiment. — [Bertoldi, Reis, Willshaw & Mullin 2010, Adv. Mater. 22(3):361–366](https://doi.org/10.1002/adma.200901956)
- **Pore shape tunes the response.** In soft cellular structures with a square array of holes, "by simply changing the shape of the holes the response of porous structure can be easily tuned" (compaction through buckling). — [Overvelde, Shan & Bertoldi 2012, Adv. Mater. 24(17):2337–2342](https://doi.org/10.1002/adma.201104395)
- **Pore shape and nonlinear response.** [Overvelde & Bertoldi 2014, JMPS 64:351–366](https://doi.org/10.1016/j.jmps.2013.11.014) **[content not re-verified]**. As I recall it, the paper relates pore shape to the onset of instability and to the post-buckling negative Poisson's ratio, using periodic cells larger than the primitive cell for the critical modes.
- **Review.** Kochmann & Bertoldi describe the field's "paradigm shift, away from avoiding instability and toward taking advantage thereof." They trace instabilities to "a nonconvex potential energy landscape," responsible for "localization, pattern formation, or structural buckling and snapping." — [Kochmann & Bertoldi 2017, Appl. Mech. Rev. 69(5):050801](https://doi.org/10.1115/1.4037966)
- **Rotating units.** [Grima & Evans 2000, J. Mater. Sci. Lett. 19:1563–1565](https://doi.org/10.1023/A:1006781224002) is the archetypal rotating-squares auxetic mechanism **[content not re-verified beyond title]**.
- **Size effects in mechanism-based metamaterials.** Mechanism-based metamaterials, including auxetics, show "a non-monotonic variation of their stiffness with system size" and sensitivity to textured boundaries. These "stem from the competition between rotation-based deformations ... and ordinary elasticity" and are "controlled by a characteristic length scale." — [Coulais, Kettenis & van Hecke 2018, Nat. Phys. 14:40–44](https://doi.org/10.1038/nphys4269). See also [Coulais et al. 2015, PRL 115:044301 (discontinuous buckling of metabeams)](https://doi.org/10.1103/PhysRevLett.115.044301) **[content not re-verified]**.
- **Isostatic lattices.** When the number of bonds balances the degrees of freedom, a frame is "on the verge of mechanical instability" and has zero-frequency floppy modes. — [Kane & Lubensky 2014, Nat. Phys. 10:39–45](https://doi.org/10.1038/nphys2835)
- **Bistable auxetics.** Networks of "rotating units connected with compliant hinges" can be both auxetic and structurally bistable, i.e. snapping instabilities in auxetic hinge networks. — [Rafsanjani & Pasini 2016, Extreme Mech. Lett. 9:291–296](https://doi.org/10.1016/j.eml.2016.09.001) ([arXiv:1612.05988](https://arxiv.org/abs/1612.05988))
  - **[content not re-verified]** A search-engine summary reported serrated snap-through force–displacement curves in *tensile* tests of such hinged auxetics.
- **Periodic frameworks and prestress (bibliographic only).** **[content not re-verified]**
  - [Guest & Hutchinson 2003, JMPS 51:383–391](https://doi.org/10.1016/S0022-5096(02)00107-2): periodic frameworks cannot be both statically and kinematically determinate.
  - [Hutchinson & Fleck 2006, JMPS 54:756–782](https://doi.org/10.1016/j.jmps.2005.10.008): periodic mechanisms of Kagome-type trusses described by Bloch wave vectors.
  - [Guest 2006, IJSS 43:842–854](https://doi.org/10.1016/j.ijsolstr.2005.03.008): prestress stiffness; tensile prestress stiffens mechanisms, compressive prestress softens them.
- **Optimized buckling designs and stiffness.** Buckling-strength-optimized periodic designs "turn out to be first-order hierarchical type microstructures," and "the gains come at the cost of reductions in stiffness." "Imposing geometric symmetry constraints can be exploited to control the shape of the failure surfaces." — [Thomsen, Wang & Sigmund 2018](https://doi.org/10.1016/j.cma.2018.04.031)
- **Coalescing critical modes.** In optimized microstructures, "the first three buckling bands along the IBZ boundaries ... are tightly clustered," i.e. "buckling for short- and long-wavelengths in different directions coincides." — [Wang & Sigmund, arXiv:2012.01359](https://arxiv.org/abs/2012.01359)
- **Penalizing localized modes.** Penalizing highly localized buckling modes in periodic-microstructure topology optimization (double-scale asymptotic stability condition) produced designs "less prone to localized instability." — [Neves, Sigmund & Bendsøe 2002, IJNME 54(6):809–834](https://doi.org/10.1002/nme.449)
- **Printed large-strain designs.** Topology-optimized architectures "with programmable Poisson's ratios ranging from −0.8 to 0.8 over large deformations of 20% or more" were designed and printed. — [Clausen, Wang, Jensen, Sigmund & Lewis 2015, Adv. Mater. 27:5523–5527](https://doi.org/10.1002/adma.201502485). Related design papers: [Wang, Sigmund & Jensen 2014, JMPS 69:156–174](https://doi.org/10.1016/j.jmps.2014.05.003); [Wang 2018, JMPS 114:303–318](https://doi.org/10.1016/j.jmps.2018.01.013). **[stability-check details in these not re-verified]**
- **Wallin & Tortorelli.** They designed neo-Hookean microstructures "for maximum load carrying capacity subject to negative Poisson's ratio constraints" with Bloch-type displacement boundary conditions, and discuss "domain symmetry of the unit cell and its effect on uniaxial loading conditions." — [Wallin & Tortorelli 2020, Mech. Mater. 145:103324](https://doi.org/10.1016/j.mechmat.2020.103324)
- **NEMOpt facts** (all from [Güven & Özdemir 2026](https://doi.org/10.1016/j.finel.2026.104599), §5.2 and §6):
  - ν = 2(𝔸_M)₁₁₂₂ / ((𝔸_M)₁₁₁₁ + (𝔸_M)₂₂₂₂), targeted at ν₀ = −1.
  - The optimum is −0.9008 for case 1, F_M = 1.075·I.
  - The topology has "internal flexible 'hinge-like' regions that allow rotation."
  - Nonlinear designs have a greater number of hinge-like regions.
  - −0.9218 and −0.9132 are reported in the comparisons.
  - Case 2 is "bi-axial stretching (2.5% in each direction) ... combined with shear."
  - "Excessive element distortion ... in very low density elements at large strains is a limiting barrier."

### Inferences
- **Derived: near-mechanism character.** For a 2D isotropic linear idealization, NEMOpt's ν definition gives ν = (κ − μ)/(κ + μ), so ν = −0.9 ⇒ κ/μ = 0.1/1.9 ≈ 0.053. The ν ≈ −0.9 designs are therefore close to a dilational mechanism (the rotating-units limit ν → −1). Near-mechanisms have low-energy modes, possibly at k ≠ 0 (e.g., alternating rotation senses with doubled period). Hence:
  - a positive-definite K_red (k = 0 only) is weak evidence of stability for these designs;
  - Coulais-type size and boundary effects mean that first-order homogenization itself becomes questionable as κ/μ → 0.
- **Derived: case 1, F_M = 1.075·I (equibiaxial tension).** There is no macroscopic compression. By prestress-stability arguments (Guest 2006, unverified here), tensile prestress stiffens mechanism modes. The risk is low but not zero: hinge regions bend, and members transverse to load paths can carry local compression.
- **Derived: case 2, F_M = [1.025 0.025; 0.025 1.025].**
  - F is symmetric, with principal stretches 1.05 and 1.00 along the ±45° diagonals. This is a plane-strain *uniaxial* stretch along one diagonal with zero lateral strain.
  - For a 2D-isotropic auxetic with lateral strain held at zero, σ_lat = ν σ_axial. With ν ≈ −0.9, the lateral stress is *compressive* and nearly equal in magnitude to the tensile axial stress.
  - So this "tensile" case loads the microstructure in genuine biaxial tension–compression. Microscopic buckling of ligaments aligned with the compressive direction is possible.
- **Derived: bulk case F_M = [1.0 0.05; 0.05 1.0].** The principal stretches are 1.05 and 0.95 (det F = 0.9975): a near-isochoric pure shear with 5% principal *compression*. Bulk-modulus-optimal designs are typically stocky, which lowers the risk, but this is the published case closest to the "compressive loading" the authors flag as open.
- **Unrepresentable post-bifurcation auxeticity.** For any future compressive load case, the Bertoldi-type mechanism (negative Poisson's ratio *by* instability) needs a cell larger than the primitive one and a post-bifurcation branch. A one-cell principal-branch optimizer either cannot represent it, or (worse) reports ν < 0 on an unstable principal branch.
- **Optimizer pressure toward instability.** By analogy with Thomsen et al. and Wang & Sigmund, an optimizer that pushes ν → −1 drives the design toward mechanism-like, near-critical states with clustered low eigenvalues. Stability checks on optimized auxetics should therefore expect near-degenerate critical modes and scan k broadly.

### Gaps
- I found no study that specifically tests the stability of *topology-optimized* finite-strain auxetic unit cells under *tensile* or *tension–shear* macroscopic loading with Bloch analysis. Whether NEMOpt's published designs are stable at their F_M is an open, testable question.
- The prestress-stability statements (Guest 2006) and periodic-mechanism statements (Guest & Hutchinson 2003; Hutchinson & Fleck 2006) are from prior knowledge, verified at title level only.
- The claim about snap-through of hinged auxetics under tension rests on a search-engine summary of Rafsanjani & Pasini (2016), not on the full text.

## Q4. Symmetry-breaking bifurcations and symmetry-reduced (half-cell) domains

### Takeaway
On NEMOpt's point-symmetric (C₂) cell, the reduced tangent block-diagonalizes exactly into two blocks:

- an "odd-field" block: C₂-invariant fluctuations, which are the ones the half cell models;
- an "even-field" block: C₂-anti-invariant fluctuations.

The half cell sees only the odd block. Its Cholesky or eigen check *cannot* detect symmetry-breaking (even-field) modes, and neither can 𝔸_M computed on either domain. A companion even-parity half-cell factorization recovers the k = 0 check at half-cell cost. Modes with k ≠ 0 need a Bloch check, of which only the four time-reversal-invariant points (Γ, X, Y, M) inherit the parity split.

The fact that full- and half-cell results agree does **not** demonstrate stability.

### Cited Findings
- **Point-symmetric homogenization (Ohno et al.).** Ohno, Matsuda & Wu developed homogenization for composites "with point symmetry of internal distributions," the basis for half-cell analysis of point-symmetric cells. — [Ohno, Matsuda & Wu 2001, IJSS 38(16):2867–2878](https://doi.org/10.1016/S0020-7683(00)00188-8)
- **Symmetry in computational homogenization (Saavedra Flores & de Souza Neto).** Their abstract states:
  - They extend Ohno et al. to staggered-translational and point symmetry under three kinematical constraints: periodic, linear boundary displacement and minimal kinematical constraint (uniform traction).
  - "Speed-up factors in excess of seven" arise when both symmetries are combined in fully coupled two-scale analyses.
  - The abstract contains no statement about bifurcation or stability.

  Source: [Saavedra Flores & de Souza Neto 2010, Eng. Comput. 27(4):551–575](https://doi.org/10.1108/02644401011044612)
- **NEMOpt's half cell.** The article states that "by employing the domain reduction technique, admissible design space is restricted since certain symmetry conditions are imposed." The method applies to orthotropic microstructures. Densities are mirrored about the vertical axis, and a 4 × 4 example has four design variables "due to orthotropy." — [Güven & Özdemir 2026, §3, §4 and §6](https://doi.org/10.1016/j.finel.2026.104599)
- **Existing NEMOpt strategy report.** It states that a 180° rotation leaves every F_M invariant, so the half-cell reduction holds even for the shear case, "but only on the symmetric branch. Full-versus-half agreement ends at the first symmetry-breaking bifurcation." — [NEMOpt refactoring strategy literature report, L36 (local file)](../../reports/NEMOpt%20refactoring%20strategy%20literature.md) (which cites [Ohno et al. 2001](https://doi.org/10.1016/S0020-7683(00)00188-8) and [Wallin & Tortorelli 2020](https://doi.org/10.1016/j.mechmat.2020.103324))
- **Group-theoretic tools.** Group-theoretic bifurcation theory and exact block diagonalization of symmetric structural eigenproblems are the standard tools:
  - [Healey 1988, CMAME 67:257–295](https://doi.org/10.1016/0045-7825(88)90049-7)
  - [Healey & Treacy 1991, IJNME 31:265–285](https://doi.org/10.1002/nme.1620310205)
  - [Ikeda & Murota 2019, *Imperfect Bifurcation in Structures and Materials*, 3rd ed., ch. "Group-Theoretic Bifurcation Theory"](https://doi.org/10.1007/978-3-030-21473-9_8)

  **[content not re-verified beyond titles; the block-diagonalization principle is standard]**
- **Honeycomb applications (bibliographic only).** **[content not re-verified]** These apply group theory, Bloch analysis and post-bifurcation stability to finitely strained honeycombs:
  - [Saiki, Ikeda & Murota 2005, Int. J. Bifurc. Chaos 15:497–515 ("flower patterns" on honeycombs)](https://doi.org/10.1142/S021812740501217X)
  - [Combescure, Henry & Elliott 2016, IJSS 88–89:296–318](https://doi.org/10.1016/j.ijsolstr.2016.02.016)
  - [Combescure, Elliott & Triantafyllidis 2020, JMPS 142:103976](https://doi.org/10.1016/j.jmps.2020.103976)
  - [Combescure, Healey & Treacy 2023, J. Nonlinear Sci. 33](https://doi.org/10.1007/s00332-022-09878-7)

  A search summary of the 2020 paper reports that equibiaxial loading of circular-cell honeycombs gives "either a simple or a double bifurcation point depending on the orientation of loading" **[search-summary level]**.
- **Ohno–Okumura–Noguchi line.** The bibliographic records are verified:
  - [Ohno, Okumura & Noguchi 2002, JMPS 50(5):1125–1153](https://doi.org/10.1016/S0022-5096(01)00106-5), "Microscopic symmetric bifurcation condition of cellular solids based on a homogenization theory of finite deformation";
  - [Okumura, Ohno & Noguchi 2002, IJSS 39:3487–3503](https://doi.org/10.1016/S0020-7683(02)00165-8).

  A search summary (sourced to [Okumura's publication list](https://profs.provost.nagoya-u.ac.jp/html/100002069_ronbn_8_en.html)) describes their 2004 JMPS 52:641–666 work **[search-summary level]**:
  - the eigenmode problem with an orthogonality condition;
  - under uniaxial compression of hexagonal honeycombs, a long-wave mode localizing in a cell row;
  - under equibiaxial compression, a "flower-like mode" changing into "an asymmetric, long-wave mode through a sextuple bifurcation."
- **Symmetry and the IBZ.** The irreducible Brillouin zone used for Bloch sweeps is "determined by the shared symmetries between the microstructure geometry and the macroscopic stress state." — [Wang & Sigmund 2023 (arXiv:2210.00003)](https://arxiv.org/abs/2210.00003)
- **Symmetry constraints shape failure surfaces.** Geometric symmetry constraints shape the failure surfaces of optimized periodic materials. — [Thomsen et al. 2018](https://doi.org/10.1016/j.cma.2018.04.031)

### Inferences
All items below are derived from standard representation theory applied to NEMOpt's setting. They are logically checkable and should be confirmed numerically on NEMOpt.

- **Parity decomposition.** Let g act on fluctuation fields by (g·w)(y) = R w(R⁻¹y) with R = −I, i.e. (g·w)(y) = −w(−y). The cell geometry, the even density field, the isotropic constituent and every F_M are invariant under the 180° rotation, so K_red commutes with g. The fluctuation space splits into:
  - V_odd = {w : w(−y) = −w(y)} (g-invariant), which is the half-cell space;
  - V_even = {w : w(−y) = w(y)}.

  Then K_red = K_odd ⊕ K_even exactly. By Sylvester's law, n_neg(K_red) = n_neg(K_odd) + n_neg(K_even).
- **What the half-cell check sees.** The half-cell Cholesky or eigenvalue check tests only K_odd. A failure there implies a full-cell failure. Success says nothing about K_even, i.e. about symmetry-breaking one-cell modes.
- **𝔸_M is blind to K_even on both domains.** The χ right-hand sides come from uniform macroscopic increments E(ij), and every constant tensor is C₂-invariant. The χ fields are therefore odd, and 𝔸_M depends only on K_odd⁻¹. Even on the *full* cell, the 𝔸_M obtained with LU is identical whether or not K_even is indefinite. Neither ν nor the rank-one check of 𝔸_M signals symmetry-breaking one-cell instability.
- **Agreement is not stability.** Starting from a symmetric state, Newton on the full cell stays on the symmetric (possibly unstable) branch except for round-off, and an LU solve does not complain about an indefinite K_red. Full/half agreement past a symmetry-breaking bifurcation is therefore expected, not reassuring. Only a definiteness test on the full K_red (or on K_even) reveals it.
- **Companion even-parity check.** K_even can be assembled on the *same* half-cell mesh by flipping the sign of the C₂-pair constraints (w(−y) = +w(y) instead of −w(y)) and releasing the six points that the odd parity forces to zero (center, corners, mid-edges). Rigid translation must still be removed (e.g., fix the center node as the full cell does), because constants are even fields.

  Cost: one extra factorization of half-cell size. Together, odd + even = the full-cell k = 0 check at about the same total cost as one full-cell factorization. Implementation details must be checked against `Domain_Reduction/NEMOpt-DR.jl`'s constraint construction.
- **Bloch modes at k ≠ 0.**
  - C₂ maps k → −k, which time reversal (complex conjugation of the real tangent) already does. C₂ alone therefore adds no reduction of the k-scan beyond half of the Brillouin zone.
  - At the four time-reversal-invariant points Γ, X = (π, 0), Y = (0, π) and M = (π, π), with phases ±1, i.e. periodic or anti-periodic, which equals 2 × 1, 1 × 2 and 2 × 2 supercells, the Bloch operator is real and commutes with g. It splits into odd and even blocks that can be done on the half cell in real arithmetic.
  - General k requires full-cell complex Bloch analysis.
- **Extra mirror symmetry under case 1.** NEMOpt's designs are orthotropic, i.e. mirror-symmetric about both axes (D₂ ⊃ C₂). Under F_M = 1.075·I the load keeps D₂, so K_red splits into four 1D-irrep blocks and the IBZ is a quarter of the Brillouin zone.

  Under the shear-type cases ([1.025 0.025; 0.025 1.025], [1.0 0.05; 0.05 1.0]), the mirrors map F₁₂ → −F₁₂ and are broken by the load; only C₂ survives. The IBZ is then half of the Brillouin zone (interior included, not just contours).
- **What the half cell cannot follow.** Generic theory (Ikeda & Murota) says:
  - a simple bifurcation in the C₂-breaking (even) sector is a symmetric pitchfork, whose two branches ±φ are related by the symmetry and may be sub- or supercritical;
  - in the C₂-preserving (odd) sector, criticality is generically a limit point or an asymmetric (transcritical) bifurcation.

  The half cell can follow odd-sector post-critical branches but none of the C₂-breaking ones.
- **Validity statement for the half cell.** Half-cell results equal the full one-cell principal branch while K_red is nonsingular, which can be certified with odd + even factorizations. They are *physically* valid only while the Bloch and rank-one checks also pass (Q1).

### Gaps
- The Ohno–Okumura–Noguchi (2002) "microscopic symmetric bifurcation condition" could not be read. Its precise relation between point-symmetric (half-cell-compatible) bifurcation and the homogenized tangent may refine the statements above, and should be read before implementation.
- No reviewed source was found that explicitly analyzes the stability limits of the Ohno-type half-cell reduction. Saavedra Flores & de Souza Neto's abstract does not mention bifurcation. The parity argument above is derived, not cited.
- I did not verify the precise six-node constraint set of NEMOpt-DR in code. The companion-check construction needs confirmation there.

## Q5. Cheap detection indicators and their reliability

### Takeaway
A ladder of indicators, from cheapest to most expensive:

1. Cholesky success or failure of the converged K_red. This is binary and free once the factorization is Cholesky.
2. Inertia from a symmetric-indefinite LDLᵀ, via Sylvester's law. This counts the number of negative modes.
3. The smallest eigenvalues by shift-invert Lanczos reusing the factorization. This gives a continuous margin.
4. An acoustic-tensor scan of 𝔸_M. In 2D this is a 2 × 2 eigenproblem over one angle and costs essentially nothing.
5. Bloch scans. Each k needs a complex factorization.

The first three only cover k = 0 (and, on the half cell, only the odd sector). Low-density regions can generate spurious modes, so expensive checks belong on the thresholded design.

### Cited Findings
- **Inertia by symmetric-indefinite factorization.** Stable symmetric-indefinite factorizations give the inertia of a symmetric matrix: [Bunch & Kaufman 1977, Math. Comp. 31:163–179, "Some stable methods for calculating inertia and solving symmetric linear systems"](https://doi.org/10.1090/S0025-5718-1977-0428694-0). **[title-level]**
- **Shift-invert Lanczos.** The spectral-transformation (shift-invert) Lanczos method is the standard way to compute the eigenvalues of a large sparse symmetric problem closest to a shift: [Ericsson & Ruhe 1980, Math. Comp. 35:1251](https://doi.org/10.2307/2006390). **[title-level]**
- **Search for a sign change.** For the first bifurcation, it suffices to search for a zero or negative eigenvalue of the condensed Bloch matrix. Positive definiteness of the auxiliary matrix lets the generalized problem be solved as a standard one. B(λ) is checked "at each loading step." — [Zhang et al. 2021, §4.3 and §5](https://arxiv.org/abs/2010.02371)
- **Cholesky failure as a diagnostic.** The existing NEMOpt strategy proposes explicit CHOLMOD Cholesky with a logged LDLᵀ/LU fallback (NDR-0022). It notes that "a failed Cholesky of a converged tangent signals a loss of positive definiteness, which has physical meaning as a microstructural instability (inference)." — local `reports/NEMOpt refactoring strategy literature.md` L279, L445; CHOLMOD: [Chen et al. 2008, ACM TOMS](https://doi.org/10.1145/1391989.1391995)
- **NEMOpt today.** NEMOpt currently solves Newton, χ and adjoint systems with Julia's `\` on `K_Reduced` (`Full_Domain/NEMOpt.jl` L367, L448; `Domain_Reduction/NEMOpt-DR.jl` L584, L676). It applies F_M through the ramp A(t) = (t/tFinal)(F_M − I) with Newton at each increment (L373–388). Code inspected locally; see also [Güven & Özdemir 2026, §4](https://doi.org/10.1016/j.finel.2026.104599).
- **Spurious modes from low-density regions.** These are a known failure mode of stability analysis in density-based topology optimization:
  - Nonlinear finite-strain topology optimization with stability constraints notes that "spurious buckling modes in these [low-density] elements can exist in the stability analysis." It suppresses them with a pseudo-mass matrix that assigns small masses to DOFs surrounded only by low-density elements. It also reviews earlier remedies: Neves et al. ignore void contributions; Gao & Ma use an eigenvalue shift. — [Zhang, Khandelwal & Guo 2023, CMAME 413:116119](https://doi.org/10.1016/j.cma.2023.116119) ([arXiv:2303.14911](https://arxiv.org/abs/2303.14911)); [Neves, Rodrigues & Guedes 1995, Struct. Optim. 10:71–78](https://doi.org/10.1007/BF01743533); [Gao & Ma 2015, Comput. Struct. 157:142–152](https://doi.org/10.1016/j.compstruc.2015.05.020)
  - A "strain energy transition approach" removes "artificial buckling modes that often occur in stability constrained topology optimization." There, the fundamental path is traversed with Newton and critical loads are estimated by eigenvalue analysis. — [Dalklint, Wallin & Tortorelli 2021, SMO 64:1751–1763](https://doi.org/10.1007/s00158-021-03012-z)
- **NEMOpt's low-density treatment.** NEMOpt uses a stiffness floor ε and, following Wallin & Tortorelli, a stress that mixes the hyperelastic response with a linear elastic response to fight element distortion in low-density elements. — [Güven & Özdemir 2026, Eqs. 29–30](https://doi.org/10.1016/j.finel.2026.104599)
- **Near-degenerate critical eigenvalues.** Optimized microstructures show tightly clustered critical buckling bands, i.e. near-degenerate critical eigenvalues. — [Wang & Sigmund, arXiv:2012.01359](https://arxiv.org/abs/2012.01359). Degenerate Bloch eigenvalues in finite-strain topology optimization have been handled with symmetric polynomials: [Dalklint et al. 2022](https://doi.org/10.1016/j.jmps.2022.104849) **[search-summary level]**.

### Inferences
- **Cholesky pass/fail.** In exact arithmetic, Cholesky success ⇔ K_red ≻ 0. Near a critical point the pivots become tiny and pass/fail is round-off-sensitive. Use a margin: attempt Cholesky of K_red − σ·M_s, with a small σ > 0 scaled to a reference eigenvalue (e.g., σ = 10⁻³ × λ_min at t = 0) and M_s a mass-like or pseudo-mass matrix that is zero or small on void DOFs. Success then certifies λ_min > σ. On the half cell, this covers only the odd sector (Q4).
- **Inertia counts.** LDLᵀ with Bunch–Kaufman (or another symmetric-indefinite) pivoting gives the number of negative eigenvalues. With a sequence of shifts this becomes spectrum slicing, which is useful to count how many modes crossed between increments.
- **Smallest eigenvalue for a continuous margin.** Shift-invert Lanczos at σ = 0 reuses the factorization and needs only a few dozen triangular solves for the lowest few eigenpairs. This gives a continuous stability margin (for plots, continuation control or optimization constraints) and the mode shape (to separate genuine ligament modes from void-localized spurious modes, e.g., by the mode's energy fraction in solid elements).
- **Acoustic-tensor scan (derived).** Q(N)_ik = (𝔸_M)_iJkL N_J N_L, with N = (cos θ, sin θ), and B = min_θ λ_min(Q(θ)).
  - det Q(θ) is a homogeneous quartic in (cos θ, sin θ), so its real roots can be found exactly from a quartic. Alternatively, a fine θ-grid (e.g., 1°) with local refinement suffices.
  - Strong ellipticity requires both eigenvalues of Q positive for all θ (det > 0 *and* trace > 0), not just det ≠ 0.
  - Cost: microseconds per state.
- **Check frequency.**
  - The binary k = 0 check and the rank-one scan are essentially free, so they can run at every converged load increment of every design iteration.
  - Bloch scans cost N_k complex factorizations per state. They belong at the final design, along its whole load ramp, and optionally every ~10–20 iterations in late continuation (large β), when the topology has stabilized.
  - Checks on gray intermediate designs are unreliable because of spurious low-density modes. The decisive checks should use the projected (thresholded) design, or use pseudo-mass or energy-transition treatment.
- **Note on Newton itself.** Newton convergence failure is not a reliable instability indicator. Newton converges to unstable equilibria and may also fail for unrelated reasons (element distortion, as the NEMOpt article notes).

### Gaps
- I did not verify in Julia documentation whether SparseArrays/CHOLMOD `ldlt` performs 2 × 2 Bunch–Kaufman pivoting. I believe it does not, which would make inertia unreliable for indefinite K. Nor did I verify which Julia-accessible sparse solvers (e.g., MUMPS, Pardiso, HSL MA57/MA97) expose inertia. This belongs to the implementation-focused notes.
- I found no quantitative guidance on choosing the shift margin σ or a mode-localization threshold for spurious-mode rejection in *periodic* finite-strain cells.

## Q6. Is 𝔸_M (the homogenized nominal tangent from the four χ problems) the right operator for the macroscopic ellipticity check?

### Takeaway
Yes, provided it is the full, consistent 4 × 4 nominal tangent ∂P_M/∂F_M at a converged, microscopically stable state, including the initial-stress (geometric) contribution and with no minor-symmetry reduction. Rank-one convexity, B = min (m⊗N):𝔸_M:(m⊗N) > 0, is the correct test.

Positive definiteness of 𝔸_M as a 4 × 4 matrix is the wrong test: it fails for any objective material under compressive mean stress because of rigid rotations. The second elasticity tensor ∂S/∂E alone is also wrong, because it omits the stress term.

### Cited Findings
- **Operator used in the literature.** Hill's stability functional for the microscale is built with the first (nominal) elasticity tensor 𝔸. Macroscopic stability is "the rank-1 convexity of the homogenized tangent moduli," tested with B(λ) = min over |m| = |N| = 1 of (m⊗N):𝔸̄:(m⊗N). When B = 0, N gives the normal of the discontinuity curves and m its nature: simple shear if m ⟂ N, pure splitting if m ∥ N. — [Zhang et al. 2021, Eqs. 44 and 54](https://arxiv.org/abs/2010.02371)
- **Geymonat et al.'s framing.** The 1993 result is stated in terms of loss of rank-one convexity of the *homogenized incremental moduli* (title of [Geymonat et al. 1993](https://doi.org/10.1007/BF00380256)). Macroscopic instability is "loss of ellipticity in the homogenized incremental equilibrium equations." — [Triantafyllidis & Maker 1985](https://doi.org/10.1115/1.3169148)
- **Hill's criterion.** Hill's uniqueness and stability criterion for finite elastic strain is the basis. — [Hill 1957](https://doi.org/10.1016/0022-5096(57)90016-9) **[title-level]**
- **NEMOpt's χ problems.** NEMOpt's χ problems are linear responses of the cell at the state imposed by F_M. The sensitivity contains a term Q that "could be considered as an initial stress like term." The authors caution that their tangent-based moduli deviate "from small strain elasticity tensor as deformation gets larger." — [Güven & Özdemir 2026, §2–4 and §5](https://doi.org/10.1016/j.finel.2026.104599)

### Inferences
- **Derived: why positive definiteness of 𝔸_M is the wrong test.** For any objective energy, W(QF) = W(F). Differentiating twice along Q(s) = exp(sΩ) with Ω skew gives (ΩF):𝔸:(ΩF) = −P:(Ω²F). In 2D, Ω² = −θ²I, so (ΩF):𝔸:(ΩF) = θ² P:F = θ² J tr σ.
  - Hence 𝔸 (and 𝔸_M, since the homogenized energy is also objective) is *indefinite* whenever the mean Cauchy stress is compressive, for any material. A 4 × 4 positive-definiteness test of 𝔸_M would flag every compressive state.
  - ΩF has rank 2, so this does not affect rank-one convexity. That is one reason rank-one convexity, not positive definiteness, is the right macroscopic test.
- **Equivalent spatial form (derived, standard).** The spatial form a_ijkl = J⁻¹ F_jJ F_lL 𝔸_iJkL with n ∝ F⁻ᵀN gives the same sign test. The referential form on 𝔸_M is simplest for NEMOpt.
- **Implementation pitfalls.**
  - 𝔸_M lacks minor symmetries, so a 3 × 3 Voigt reduction must not be used.
  - It should have major symmetry (hyperelastic), which is a cheap sanity check on the four χ solves.
  - The full 4 × 4 must enter Q(N).
- **Order of checks.** 𝔸_M is the Schur complement of the full cell operator. If K_red is indefinite (unstable one-cell state), 𝔸_M is the tangent of an unstable equilibrium and its rank-one test is not meaningful on its own. Check K_red (both parity sectors) first, then rank-one, then Bloch.
- **Macroscopic control type.** The criteria above assume an infinite medium under macroscopic deformation control. Stress-controlled or mixed macroscopic loading, or finite specimens, add structural modes and different criteria; these are outside this note.

### Gaps
- I did not find an explicit statement, in a source I could read, comparing rank-one tests on 𝔸 vs. ℂ = ∂S/∂E for homogenized tangents. The argument above is derived.
- The "symmetric bifurcation" criterion of Ohno et al. (2002), which may relate point-symmetric one-cell bifurcation to singularity of a homogenized tangent under stress control, was not accessible.

## Q7. Recommendations for NEMOpt: what to check, when, and at what cost

### Takeaway
Add a three-tier stability protocol:

- **Tier 0 (every converged load increment, about free):** a k = 0 definiteness test of K_red plus an acoustic-tensor scan of 𝔸_M.
- **Tier 1 (half cell only, about +1 half-cell factorization):** a companion even-parity test.
- **Tier 2 (final designs along the whole load ramp; occasional late iterations):** a Bloch scan over the IBZ.

Report the critical load factor t_c for every published design. Treat any result at t ≥ t_c as invalid for both the full-cell and half-cell models.

### Cited Findings
- **Full analysis practice.** A full multiscale stability analysis checks rank-one convexity at every load step and Bloch over a dense k-mesh with refinement near k = 0. Their verification examples include a one-cell critical mode and a long-wavelength mode. — [Zhang et al. 2021](https://arxiv.org/abs/2010.02371)
- **IBZ-boundary sweeps.** These are used in topology-optimization practice; the IBZ is determined by the symmetries shared by the microstructure and the macroscopic stress — [Wang & Sigmund 2023](https://arxiv.org/abs/2210.00003). Extrema off the IBZ contour are a documented concern — [Maurin et al. 2018](https://doi.org/10.1016/j.ijsolstr.2017.11.006) (title-level).
- **Spurious modes.** Spurious low-density modes must be handled in stability analyses inside topology optimization — [Zhang, Khandelwal & Guo 2023](https://doi.org/10.1016/j.cma.2023.116119); [Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z).
- **Benchmark geometry.** Pattern-transformation benchmarks (a neo-Hookean matrix with a square array of circular holes under compression) have well-documented instability-driven responses. — [Mullin et al. 2007](https://doi.org/10.1103/PhysRevLett.99.084301); [Bertoldi et al. 2010](https://doi.org/10.1002/adma.200901956)
- **Open item named by the authors.** NEMOpt's authors name instabilities under compressive loading as open work. — [Güven & Özdemir 2026, §6](https://doi.org/10.1016/j.finel.2026.104599)

### Inferences
1. **Tier 0: always on, both domains.**
   - **What:**
     - At every converged increment t_j of the ramp F(t), attempt Cholesky of K_red − σM_s (margin σ, pseudo-mass M_s; see Q5).
     - Compute B(t_j) = min_θ λ_min(N(θ)·𝔸_M·N(θ)) on the full 4 × 4 𝔸_M.
     - Log pass/fail, the minimizing θ, and (if available) the LDLᵀ negative-eigenvalue count.
   - **Cost:** zero extra factorizations if NDR-0022 (explicit Cholesky) is adopted. Otherwise, one Cholesky per increment (less than one Newton iteration). The ellipticity scan is negligible.
   - **Use:** flag the design iteration, and stop treating 𝔸_M-derived objectives (ν, κ, μ) as valid when either check fails.
2. **Tier 1: half-cell runs only.**
   - **What:** at least at the final increment of each design iteration (ideally at every increment), factorize the even-parity companion K_even on the half-cell mesh (Q4).
   - **Cost:** about one extra half-cell factorization, so the half-cell route remains cheaper than the full cell. Optionally run a full-cell Tier 0 check on the final half-cell design as a cross-check.
   - **Why:** without this, the half cell (and 𝔸_M on *either* domain) is blind to symmetry-breaking one-cell modes. Full/half agreement does not certify stability.
3. **Tier 2: Bloch scan.**
   - **What:** for every published or final design, at every ramp increment up to and somewhat beyond F_M (e.g., continue to t = 1.25 to report a margin t_c):
     - check the four time-reversal-invariant points first (Γ, X, Y, M). On the half cell these are real, parity-split problems.
     - then scan a k-grid over the IBZ with refinement near k = 0. The IBZ is a quarter of the Brillouin zone for D₂-invariant loads such as 1.075·I, and half of it for the shear-type F_M.
     - use LDLᵀ inertia at shift σ for pass/fail, and shift-invert Lanczos for margins and modes.
   - **Cost:** N_k complex Hermitian factorizations per state. Complex arithmetic is roughly 2–4× real per factorization (my estimate). With N_k ~ 10²–10³, this is feasible for final designs but too costly for every iteration. During late continuation, check only the time-reversal-invariant points plus the IBZ contour every ~10–20 iterations.
   - **Handling spurious modes:** run Tier 2 on the thresholded design, or with pseudo-mass or energy-transition treatment of voids, and inspect the mode shapes.
4. **Verification tests before trusting the checks.**
   - Reproduce Zhang et al.'s examples: a one-cell mode, and a long-wavelength mode in which β and B vanish together.
   - Reproduce the circular-hole square-array benchmark. The critical mode there is expected to need a cell larger than the primitive one **[expected from prior knowledge, not re-verified]**: Tier 0 should *pass* while Tier 2 *fails*, which demonstrates why Bloch is needed.
   - Build a synthetic C₂-symmetric cell whose first mode is even. Tier 0 on the half cell should pass and Tier 1 should fail.
5. **Validity rules to document.**
   - The **one-cell model** is valid on [0, t_c), where t_c is the first failure of Tier 0, Tier 1 or Tier 2.
   - The **half-cell model** equals the one-cell model on [0, t_even) ∩ [0, t_odd), and is physically valid only on [0, t_c).
   - Beyond t_c, report "principal branch unstable." A post-bifurcation response needs supercell or branch-following analysis, ideally with group-theoretic branch switching (Ikeda & Murota; Healey), and is outside the one-cell and half-cell framework.
6. **Priorities for the published cases (inference).**
   - Highest priority: case 2 ([1.025 0.025; 0.025 1.025]), where the auxetic coupling induces lateral compression; and the bulk "shear" case ([1.0 0.05; 0.05 1.0]), which has a 5% principal compression.
   - Lower priority: case 1 (1.075·I).
   - Any future compressive load case requires Tier 2 from the outset, and probably a stability constraint (λ_min(K(k)) ≥ δ over sampled k, and B ≥ δ). Such a constraint has a nonsmooth max–min structure with clustered eigenvalues, as optimized designs show (Wang & Sigmund).

### Gaps
- There are no measured NEMOpt timings for factorizations and Bloch scans. The cost statements are relative estimates.
- Whether the published NEMOpt designs pass Tier 0, Tier 1 or Tier 2 at their F_M is unknown. Running Tier 0 and Tier 1 on the published case-1 and case-2 designs is the cheapest decisive experiment.
- I did not assess the sensitivity analysis (adjoint of eigenvalues along a nonlinear path) needed to turn the checks into optimization constraints. It is covered by finite-strain topology-optimization literature ([Dalklint et al. 2021](https://doi.org/10.1007/s00158-021-03012-z); [Zhang et al. 2023](https://doi.org/10.1016/j.cma.2023.116119)) and is out of this note's scope.

## Corrections after full-text verification (2026-10-12)

The text above is kept as first written. The items below correct it after the supplied papers and open-access sources were read in full. Each item gives the corrected statement and the verification note that supports it. All notes are in [NEMOpt supplied papers verification](../NEMOpt%20supplied%20papers%20verification/). Line numbers refer to this file as of 2026-10-12.

- **L21 (summary).** Two statements need qualification.
  - "Exactly the k → 0 limit" holds at the linear level, at a state with β ≥ 0 and a positive cell-periodic constant: 𝔸̄ loses strong ellipticity if and only if the critical mode is long-wave (Geymonat et al. 1993, Thm 3.4(ii)–(iii), Cor. 3.5). For the nonlinear homogenized energy, only long-wave ⇒ loss of rank-one convexity is proved (Thm 5.1(ii)).
  - "Lose their physical meaning" should read: one-cell homogenized moduli are certified only up to the first bifurcation. While the cell is Bloch-stable, the one-cell tangent equals every multi-cell tangent (Thm 3.3(iv)); beyond it the linearized energy of a large enough body is unbounded below (Prop. 3.2).
  - Source: [GMT note](../NEMOpt%20supplied%20papers%20verification/geymonat_muller_triantafyllidis_1993.md), Q1 claims 2 and 5.
- **L25 ("The full text was not accessible this session …").** Corrected: the full text was later read from a user-supplied copy, and the items below are confirmed with the qualifications in the GMT note. Source: GMT note E12.
- **L26 ("general 3D periodic solids").** Corrected: "for the microscale stability of periodic solids of infinite extent in any dimension N ≥ 2 (Lemma 4.2)". Source: GMT note E13.
- **L33 ("loss of rank-one convexity *is* long-wavelength microscale buckling").** Qualify as in the L21 item: at a state with β ≥ 0, B = 0 if and only if β(k → 0) = 0 (Geymonat et al. 1993, Cor. 3.5, which needs a positive cell-periodic constant). Source: GMT note, Q1 claims 2 and 9.
- **L36 ("the stability surface β(k) is continuous at the origin").** Corrected: "and the stability surface β(k) is continuous at the origin (Zhang et al.; not stated by Geymonat et al., but it follows from their Thm 3.3(i): if Λ = Λ₆ then Λ₅ = Λ₆)". Source: GMT note E14.
- **L37 ("The surface is discontinuous at the origin there").** Read "can be discontinuous". The general result is only Λ₅ ≤ Λ₆ (Thm 3.3(i)); the jump is explicit in Geymonat et al.'s laminate (p. 278). Source: GMT note, Q1 claim 4.
- **L44 (Müller 1987).** Append: "formula verified as restated in Geymonat et al. 1993, Thm 2.3 and Eqs. 2.6–2.7 (Müller's own text not read)". Source: GMT note E15.
- **L64 (Gaps: "The Geymonat et al. (1993) full text was not read …").** Resolved: the full text was read. For hypotheses and the k = 0 versus k → 0 relation (Λ₅ ≤ Λ₆), see the GMT note. Source: GMT note E16.
- **L71 ("The established workflow (Triantafyllidis et al. 2006; …) has four steps").** Corrected: "The established workflow (Triantafyllidis et al. 2006 for steps 1–3, with an LDU sign test instead of an eigen-solve in step 2; Bertoldi & Boyce 2008; Zhang et al. 2021) has four steps". Triantafyllidis et al. do not retrieve modes; their Fig. 5 mode is Michel's computation. Source: [TNS note](../NEMOpt%20supplied%20papers%20verification/triantafyllidis_nestorovic_schraad_2006.md), §5 research-note edits (claim N3).
- **L106 (Michel et al. 2007, "The abstract was not retrievable … As I recall it …").** Corrected: the full text was read from the HAL deposit hal-00111501 (submitted version).
  - The paper compares second-order homogenization estimates of macroscopic loss of ellipticity (random, square and hexagonal arrays of circular holes) with the finite-element one-cell homogenized moduli.
  - First-bifurcation curves come from 2×2 aggregates, after a single coarse Bloch scan of one case (square, porosity π/16).
  - The PDF supplied as "Michel 2007" was in fact Lopez-Pamies & Ponte Castañeda (2007), JMPS 55:1677–1701.
  - Source: [Michel note](../NEMOpt%20supplied%20papers%20verification/michel_et_al_2007.md), E12 and §0.
- **L118 (Gaps: full texts not retrieved).** Corrected: Miehe et al. (2002), Bertoldi et al. (2008, JMPS) and Overvelde & Bertoldi (2014) remain unretrieved; Michel et al. (2007) was later read from HAL. Source: Michel note E12.
- **L251–252 and L325 (Ohno–Okumura–Noguchi condition not read; no source on the stability limits of the half cell).**
  - Update: the condition can now be cited from its open precursor, Okumura, Ohno & Noguchi 2000, JSME CMD 13:769, §2.3, Eqs. 15–18. It requires a nontrivial cell-periodic solution of the incremental cell problem with zero macroscopic strain and stress rates.
  - "Symmetric bifurcation" there means that reversing the mode's sign leaves the macroscopic state unchanged. It is not point symmetry and not the half cell.
  - The absence of stability content is now confirmed for Saavedra Flores 2011, Ch. 4, as well as for the 2010 abstract.
  - The 2002 JMPS paper itself remains unread.
  - Source: [mechanics note](../NEMOpt%20supplied%20papers%20verification/open_access_mechanics_sources.md), Item 1.
