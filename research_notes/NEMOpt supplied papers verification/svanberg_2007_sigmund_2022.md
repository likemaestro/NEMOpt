# Svanberg (2007) and Sigmund (2022): full-text check of the NEMOpt reports' claims

Date: 2026-10-12. Both PDFs were supplied by the user and read in full with the Read tool (Svanberg pp. 1–15, Sigmund pp. 1–10). Their text was also extracted into the session scratchpad, only to search it. No PDF was copied into the repository. The reports, decision records and other notes were **not** edited. Section C gives every change as "old → new".

Line numbers are as read on 2026-10-12, while another editor was working on the reports. Find each edit by its quoted old text, not by its line number.

Verdicts used here:
- **confirmed**: the document says it;
- **confirmed with refinement**: correct, but the document says more, or says it more precisely;
- **needs correction**: the document says something different;
- **not in the document**: the document does not contain the claim.

Quotes keep the original spelling. Sigmund writes "practise" and "visualise", and Svanberg writes "is obtain".

Side check. Section A.3 compares Svanberg's notes with the NonconvexMMA source at v1.1.2 (commit `fa9352c`). This uses a clone in the session scratchpad. NonconvexMMA is not one of the two documents. Its findings are labeled "source reading".

---

## 0. Bibliographic records (from the PDFs)

| | Svanberg 2007 | Sigmund 2022 |
|---|---|---|
| Title | "MMA and GCMMA – two methods for nonlinear optimization" | "On benchmarking and good scientific practise in topology optimization" |
| Type and length | Lecture notes, Optimization and Systems Theory, KTH, Stockholm; 15 pages; no date or version printed | "Educational Paper"; 10 pages ("Page 1 of 10" to "Page 10 of 10"), **not about 24** |
| What the PDF prints | Opening line: "the author's 2007 implementations of MMA and GCMMA in Matlab" (p. 1). PDF metadata: built 2014-10-14 (pdfTeX), so this copy is a later build of the 2007 notes. | Header: *Structural and Multidisciplinary Optimization* (2022) **65:315**, i.e., volume 65, article number 315. The issue number is not printed; Crossref lists issue 11 (earlier note). doi:10.1007/s00158-022-03427-2. Received 7 May 2022, revised 4 August 2022, accepted 3 October 2022, published online 25 October 2022. Responsible Editor: Gregoire Allaire. |
| References | [1] Svanberg 1987, IJNME 24:359–373; [2] Svanberg 2002, SIAM J. Optim. 12:555–573 (p. 15) | 16 references; §2.17 cites Haftka et al. 2019 |
| Access | The former host people.kth.se/~krille/mmagcmma.pdf is dead (earlier note). Read from the user's copy. | Paywalled version of record. Read from the user's copy. |

Suggested reference strings:
- Svanberg K (2007) MMA and GCMMA – two methods for nonlinear optimization. Lecture notes, KTH, Stockholm, 15 pp. Formerly at people.kth.se/~krille/mmagcmma.pdf; read from a user-supplied copy.
- Sigmund O (2022) On benchmarking and good scientific practise in topology optimization. Struct Multidisc Optim 65:315 (10 pp.). https://doi.org/10.1007/s00158-022-03427-2

---

## Part A. Svanberg (2007)

### A.1 What the notes contain

- **§1, problem (Eq. 1.1), p. 1.** The functions f₀, …, f_m are "given, continuously differentiable, real-valued functions on X". The variables y and z are artificial variables, used to model least-squares and min–max problems.
- **§2, practical considerations, p. 3.**
  - Scale the objective so that "1 ≤ f₀(x) ≤ 100" for reasonable x. Scale the constraints so that 1 ≤ g_i^max ≤ 100, and the variables so that 0.1 ≤ x_j^max − x_j^min ≤ 100.
  - If some y_i ≠ 0 at the optimum, "increase the corresponding values of c_i by e.g. a factor 100 and solve the problem again".
  - If a variable reaches an artificial bound, "change this bound and solve the problem again (starting from the recently obtained solution)".
- **§3, ordinary MMA, pp. 3–5.**
  - Subproblem: Eq. 3.1.
  - Approximations: Eqs. 3.2–3.5. Eqs. 3.3–3.4 include the 1.001/0.001 weights and a 10⁻⁵/(x^max − x^min) term.
  - Move bounds: Eqs. 3.6–3.10.
  - Asymptote rules: Eqs. 3.11–3.14.
  - A list of default values (p. 5; see A.2).
- **§4, GCMMA, pp. 6–8.**
  - Outer iterations k and inner iterations ν.
  - The conservativeness test (p. 6).
  - The parameters ρ_i^(k,ν) (Eqs. 4.3–4.9).
  - The first-order property, Eq. 4.10: f̃ and its gradient equal f and ∇f "at the current iteration point".
- **§5, primal–dual interior-point subproblem solver, pp. 9–15.** The KKT conditions (Eq. 5.7), the relaxed conditions (Eq. 5.9), the Newton system (Eqs. 5.11–5.22), the line search (§5.4) and the algorithm (§5.5). The solver stops at "say ε^(ℓ) ≤ 10⁻⁷" (p. 15).
- **What the notes do not contain:**
  - a stopping criterion for the outer iterations;
  - a maximum number of inner iterations;
  - a tolerance in the conservativeness test;
  - continuation, or any change of the f_i during a run;
  - asymptote resets;
  - the names "MMA87" or "MMA02";
  - the 1987 paper's parameter values.

### A.2 Answers to the five questions

**(1) The GCMMA conservativeness test and its tolerance.** An outer iteration ends when f̃_i^(k,ν)(x̂^(k,ν)) ≥ f_i(x̂^(k,ν)) "for all i = 0, 1, …, m" (§4, p. 6).
- Otherwise a new subproblem "is generated and solved at x^(k)". Its approximations are more conservative "for those indices i for which the above inequality was violated" (p. 6).
- **No tolerance appears.** The only 10⁻⁷ in the notes is the stop of the interior-point subproblem solver (§5.5, p. 15).
- The other small constants are:
  - the ρ floor 10⁻⁶ (after Eq. 4.6);
  - the 10⁻⁵ default raa0 of ordinary MMA (Eqs. 3.3–3.4).
- The `f_app + epsimin ≥ f_new` slack of the GPL Python port is therefore the port's own choice: it reuses the subproblem ε. Svanberg does not prescribe it.
- NonconvexMMA's test, `fg[i] > approxfg[i]` (`mma_algorithm.jl` L479), has no tolerance, like the notes.

Inference: the slack probably absorbs the error of the inexact subproblem solution (ε ≈ 10⁻⁷).

**(2) Default values.** The notes do call them defaults.
- On the update rules: "The default rules for updating the lower asymptotes … and the upper asymptotes" (p. 5).
- On the parameters: "most of the explicit numbers in the above expressions are just default values of different parameters in the fortran code" (p. 5). The notes say "fortran" here, although they describe Matlab implementations.
- On changing them: "All these values can be carefully changed by the user." A smaller asyinit is the conservative direction: "a more conservative method is obtain by decreasing move and/or asyinit and/or asyincr" (p. 5).

| Notes' parameter | Default | Where | Called a default? | NonconvexMMA v1.1.2 (source reading) | NEMOpt |
|---|---|---|---|---|---|
| raa0 | 10⁻⁵ | Eqs. 3.3–3.4 | Yes, named (p. 5) | No equivalent in MMA87 (ρ = 0) | — |
| albefa | 0.1 | Eqs. 3.6–3.7 | Yes, named | 0.1, hard-coded (`mma_model.jl` L115–116) | — |
| move | 0.5 | Eqs. 3.6–3.7 (and 3.10) | Yes, named | **Absent**: no move-limit term | — |
| asyinit | 0.5 | Eq. 3.11 | Yes, named | `s_init` = 0.5 (σ = s_init·(x^max − x^min), L413) | **0.25** |
| asydecr | 0.7 | Eq. 3.13 | Yes, named | `s_decr` = 0.7 | 0.7 |
| asyincr | 1.2 | Eq. 3.13 | Yes, named | `s_incr` = 1.2 | 1.2 |
| Asymptote bounds | 0.01 and 10 × (x^max − x^min) | Eq. 3.14 | Only through "most of the explicit numbers" | diff/100 and 10·diff (L440–441) | Same as NonconvexMMA |
| Weights 1.001 and 0.001 | — | Eqs. 3.3–3.4, 4.3–4.4 | Not named | Not used: p = σ²max(0, ∇f) + ρσ/4 | — |
| Initial ρ (GCMMA) | (0.1/n)·Σ_j \|∂f_i/∂x_j\|(x_j^max − x_j^min), floor 10⁻⁶ | Eq. 4.6 | "the following default values are used" (p. 7) | ρ ← max(ρ/10, 10⁻⁵) at each outer iteration, from ρ = 0 (L456) | — |
| ρ increase | min{1.1(ρ + δ), 10ρ} if δ > 0 | Eqs. 4.8–4.9 | A rule, not a parameter | Same form (L481) | — |
| Inner-iteration cap | None | §4 | — | `maxinner` = 10 | **40** |
| Subproblem stop | ε ≤ 10⁻⁷ ("say") | §5.5 | Offered as an example | The dual is solved with Optim.jl, a different solver | — |

These are the defaults of the author's **2007** implementations: "The first versions of these methods were published in [1] and [2]" (p. 1). The notes do not say that these are the 1987 values.

Consequences for the article:
- The article's "default parameters given by [32]" ([32] = Svanberg 1987) cannot be checked against these notes.
- NEMOpt's `s_init` = 0.25 is not the default either in Svanberg 2007 or in NonconvexMMA. It does lie in the direction Svanberg calls more conservative.
- `maxinner` = 40 has no counterpart in the notes.

**(3) Changing the problem functions or parameters between iterations.** The notes never discuss continuation. The relevant facts are these:
- The f_i are "given" (§1).
- The asymptotes carry history from x^(k−1) and x^(k−2) (Eqs. 3.12–3.13), so they assume one fixed problem across iterations.
- GCMMA's test compares f̃^(k,ν), which is built from f_i(x^(k)) and ∇f_i(x^(k)), with f_i(x̂). It therefore assumes the same f_i at both points.
- Within one outer iteration, "the only differences between two inner iterations are the values of" the ρ parameters (p. 7).
- Between outer iterations, ρ restarts from Eq. 4.6 using the current gradients. In Svanberg's GCMMA, only the asymptotes carry memory across a change of the problem.
- "All these values can be carefully changed by the user" (p. 5) is about choosing values, not about changing them during a run.
- The nearest guidance concerns changed problem data, namely c_i or artificial bounds: "solve the problem again (starting from the recently obtained solution)" (§2, p. 3). A new solve starts again at the k = 1, 2 asymptote rule (Eq. 3.11).

Inference:
- A continuation step taken between outer iterations comes closest to this "solve again" pattern.
- Changing p and β inside the objective makes the GCMMA test compare two different functions, which the notes never contemplate. Under MMA02 this happens between the evaluation at x^(k) and the evaluation at x̂.

**(4) Inner iterations and their cost.**
- An outer iteration may have "zero, one, or several inner iterations" (p. 6).
- Each inner iteration solves a new subproblem at the same x^(k), with the same asymptotes and bounds. It needs f_i(x̂^(k,ν)) for all i, i.e., one evaluation of values.
- Gradients are not needed in inner iterations: "there is no need to recalculate the gradients ∇f_i(x^(k)), since x^(k) has not changed. Gradients of the original functions f_i are calculated only once in each outer iteration" (p. 6).
- An outer iteration with ν inner iterations therefore costs ν + 1 value evaluations and one gradient evaluation. The gradient is computed at the accepted point, whose value is already known.
- The test is met "after a finite (usually small) number of inner iterations" (p. 6). No cap is stated.
- Ordinary MMA needs one value and one gradient per iteration.

**(5) Does a quantized objective with an exact gradient violate a stated assumption?** Yes: the assumption that f₀ is **continuously differentiable** (§1, p. 1).
- A value rounded to 10⁻³ is piecewise constant. It jumps at rounding boundaries and has zero derivative elsewhere.
- Supplying the gradient of the unrounded J with it gives a value and a gradient that belong to no single C¹ function. The approximations are then "first order approximations" (Eq. 4.10) of a function that the solver never evaluates.
- In NEMOpt the gradient is not even the exact gradient of the unrounded J, because ν in 200(ν − ν₀)∇ν is also rounded (NDR-0028).

Effect on each algorithm:
- **Ordinary MMA (MMA87).** The objective value enters the subproblem only through the constant r₀ of Eq. 3.5, and the solver "drop[s] the constant r₀^(k,ν) from the objective function" (§5, p. 9).
  - Inference: in exact arithmetic, rounding J alone cannot change an MMA87 iterate. In floating point it can perturb the dual solve only at round-off level, through the constant term.
  - Inference: the full-domain path changes through the rounded ν in ∇J (NDR-0028), and through any stopping test that uses J. Confirm in NonconvexCore that `KKTCriteria` uses no objective value.
- **GCMMA (MMA02).** The rounded value enters the test (§4) and δ (Eq. 4.8). Late in a run the true change per step falls below the 10⁻³ rounding step. The rounded f(x̂) then often equals the rounded f(x^(k)) while f̃(x̂) predicts a decrease. The test fails, and ρ grows (Eq. 4.9) until the predicted decrease nearly vanishes or a rounding boundary is crossed (inference).
  - Svanberg's "finite (usually small)" holds under his assumptions. For a step function the loop may need the cap (inference).
  - In NonconvexMMA each failed test still moves the design (A.3) and advances NEMOpt's continuation counter.
- **Scaling (advice, not an assumption).** J = 100(ν + 1)² is 1 at ν = −0.90 and 0.25 at ν = −0.95. Late in a run, J falls below Svanberg's range 1 ≤ f₀ ≤ 100 (§2, p. 3), and the fixed 10⁻³ rounding becomes a larger relative perturbation (inference).

### A.3 Side check: NonconvexMMA v1.1.2 against the notes (source reading, not one of the two documents)

| Item | Svanberg 2007 | NonconvexMMA `mma_algorithm.jl` |
|---|---|---|
| Asymptote rule | Eqs. 3.11–3.14: start at 0.5, factors 0.7/1.2/1, bounds 0.01 and 10, update from k ≥ 3 | Same structure: `s_init`, `s_decr`, `s_incr`; `outer_iter > 2 && updateσ!` (L246); bounds L440–441 |
| Where an inner iteration solves | "solved at x^(k)", with fixed asymptotes; only ρ changes (pp. 6–7) | **At the trial point.** `updatex!` (L287, L392) sets x to the subproblem solution, and `updateapprox!(dualmodel, optimalx, fg, ∇fg)` (L310) re-centers the approximation and the asymptotes there, whether or not the test passed |
| Evaluations per inner iteration | Values only (p. 6) | Value and Jacobian (`value_jacobian`, L295) |
| Test | f̃ ≥ f for all i, no tolerance | `fg[i] > approxfg[i]` makes the step not upper-bounded (L479); no tolerance |
| ρ | Eq. 4.6 at each outer start; Eq. 4.9 on failure | max(ρ/10, 10⁻⁵) at each outer start (L456); Eq. 4.9 form on failure (L481) |
| Inner cap | None | `iter < prev_iter + options.maxinner` (L264) |

**Consequence.** NonconvexMMA's MMA02 never rejects a step. It moves to each trial point, raises ρ for the following steps of the same outer iteration, and ends the outer iteration only after a conservative step. Two things follow:
- Every objective evaluation after the setup evaluation is a design update, under MMA87 and MMA02 alike.
- What differs from Svanberg's GCMMA is the grouping of design updates into outer iterations, and the fact that a "rejected" step is still taken.

This sharpens the reports' wording on "rejected iterations" (C.1, edit 8). It does not change any decision.

### A.4 Claims attributed to Svanberg, or bearing on his notes, and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| SV-1 | Strategy report L208 (rounding row) | The ±5e−4 noise is "far above the slack epsimin = 1e−7 that the GPL Python port … adds in its GCMMA conservativeness check (… Svanberg's 2007 notes not re-read)"; "NonconvexMMA's check has no tolerance at all" | **Needs correction** (now checkable) | The test has no tolerance (§4, p. 6). 10⁻⁷ only stops the subproblem solver (§5.5, p. 15). NonconvexMMA agrees with the notes. Compare the noise with the predicted decrease, not with a tolerance. |
| SV-2 | Same row, "Effect if changed" (A, both drivers) | Rounding changes the path | **Confirmed with refinement** | MMA02: the value enters the test (§4) and δ (Eq. 4.8). MMA87: the value is a dropped constant (§5, p. 9), so the full-domain effect runs through NDR-0028 and stopping tests (inference). |
| SV-3 | Report L212 (clip row) | GCMMA's theory assumes continuously differentiable functions ("assumption not re-verified") | **Confirmed** for the problem class | §1, p. 1: "continuously differentiable". The notes give no theorem; they cite [2]. |
| SV-4 | Report L211 and L239; NDR-0035; `reproduction_baseline.md` item 12; `wallin_tortorelli_2020.md` Q5 | MMA87 matches Wallin & Tortorelli's plain MMA with Svanberg (1987) values. `s_init` 0.25 and `maxinner` 40 are undocumented. The library defaults are 0.5 and 10. | **Confirmed with refinement** | Svanberg's 2007 defaults are asyinit 0.5 (= `s_init`), asyincr 1.2 and asydecr 0.7 (p. 5). NEMOpt changes only asyinit, toward "a more conservative method". The notes set no inner cap, and they describe the 2007 implementation, not 1987 (p. 1). |
| SV-5 | Report L37 (panel) | Switching between MMA87 and GCMMA changes the path | **Consistent** | GCMMA adds the ρ terms (Eqs. 4.3–4.4) and inner iterations (§4). |
| SV-6 | Report L84 | NonconvexMMA evaluates value and Jacobian per inner iteration; MMA02 repeats until conservative or `maxinner` | **Confirmed** for NonconvexMMA. It **differs from Svanberg's GCMMA.** | Svanberg: inner iterations at a fixed x^(k), values only (p. 6). NonconvexMMA moves to each trial point (A.3). |
| SV-7 | NDR-0029; `numerical_method_choices.md` §3 | Count continuation per outer (design) iteration | **Supported** | Svanberg's iteration index is the outer k, with one accepted point: "the next iteration point becomes (x^(k+1), …)" (p. 6). |
| SV-8 | `numerical_method_choices.md` L25 | "far larger than GCMMA's conservativeness tolerance (1e-7)" | **Needs correction** | As SV-1. |
| SV-9 | `numerical_method_choices.md` L252 | GCMMA accepts an iterate only if each approximation is ≥ the true value, "within ~1e-7 in Svanberg's code" | **Needs correction** (first half confirmed) | §4, p. 6: "for all i = 0, 1, …, m", with no tolerance. NonconvexMMA moves even when the test fails. |
| SV-10 | `numerical_method_choices.md` L252 | "Each rejected inner iteration costs a full nonlinear solve plus gradients in Nonconvex.jl" | **Confirmed** for NonconvexMMA | Svanberg's GCMMA would need values only (p. 6). |
| SV-11 | `numerical_method_choices.md` L256–258 ([PR-bib]) | The convergence results assume continuously differentiable functions | **Confirmed** for the assumption | §1. The theorem is not in the notes. |
| SV-12 | `numerical_method_choices.md` L259–263 | The port's defaults (asyinit 0.5, asyincr 1.2, asydecr 0.7, bounds [0.01, 10], move 0.5, albefa 0.1) and its ρ update | **Confirmed in Svanberg 2007** | Eqs. 3.6–3.7, 3.11–3.14 and 4.8–4.9; p. 5. |
| SV-13 | `numerical_method_choices.md` L262 | The port's `concheck` uses epsimin = 1e-7 | **Confirmed for the port only; not in the document** | — |
| SV-14 | `numerical_method_choices.md` L285 | The noise is "~5000× the 1e-7 tolerance in Svanberg's check" | **Needs correction** | There is no tolerance in the check. |
| SV-15 | `numerical_method_choices.md` L286–288 | Spurious rejections grow ρ, shrink steps and add inner iterations; spurious acceptances occur; continuation accelerates | **Confirmed in mechanism** for Svanberg's GCMMA (Eq. 4.9) | In NonconvexMMA a "rejected" step is still taken (A.3). |
| SV-16 | `numerical_method_choices.md` L289 | "MMA87 has no acceptance test. Rounding mainly affects the KKT and convergence tests and the ν used in the gradient." | **Confirmed and sharpened** | r₀ is dropped (§5, p. 9). |
| SV-17 | `numerical_method_choices.md` L298 (gap) | The 2007 statements on differentiability and on the inner-loop tolerance were not re-read | **Resolved** | §1; §4 (no tolerance, no cap); §5.5 (10⁻⁷ for the subproblem). |
| SV-18 | `numerical_method_choices.md` L281 | Link rot at KTH | **Still true** | The user's copy was read. |
| SV-19 | `numerical_method_choices.md` L240 and L244 (gap) | Resetting asymptotes or ρ after continuation; "no peer-reviewed source" | **Partly addressed; still not peer-reviewed** | The notes are silent on continuation. For changed problem data they advise solving again from the latest solution (§2, p. 3), which restarts the asymptotes (Eq. 3.11). GCMMA's ρ restarts at every outer iteration anyway (Eq. 4.6). |
| SV-20 | `julia_architecture.md` L196 | "Continuation should be an outer loop with explicit restarts" | **Consistent** | §2, p. 3. |
| SV-21 | `julia_architecture.md` L259, item 2 | Rounding matters for MMA02's checks and for the KKT/f tolerances | **Consistent; add the MMA87 refinement** | §5, p. 9. |
| SV-22 | `open_access_optimization_papers.md` §5 (V-a to V-d, edit 5) | 1e−7 "not verifiable in Svanberg 2007"; defaults confirmed for the port; parameter changes not verifiable; a standing rule not to cite 2007 values | **Superseded** | V-a: not in the document (10⁻⁷ is the subproblem stop). V-b: confirmed. V-c: partly answered (§2; p. 5). V-d: resolved. Edit 5 can be lifted. |
| SV-23 | `reproduction_baseline.md` L125; `wallin_tortorelli_2020.md` L341 | Article and thesis: "default MMA parameters presented by (Svanberg, 1987)" | **Not checkable** with these notes | The 1987 paper was not supplied. The 2007 defaults belong to the 2007 implementation. |
| SV-24 | Report L325; `julia_architecture.md` L264 | A clean-room MMA/GCMMA written from Svanberg 1987 and 2002 | **Consistent** | The 2007 notes are a complete algorithmic specification, including the subproblem solver (§5), and are a useful added source. |

The task's list of port defaults (asyinit 0.5, asyincr 1.2, asydecr 0.7, move 0.5, albefa 0.1, raa0 1e−5) is **confirmed** by the notes (p. 5). The "1e−7 conservativeness tolerance" is **not in the document**.

---

## Part B. Sigmund (2022)

### B.1 Recommendations by topic

**Reporting**
- **Continuation schedules.** No recommendation. Continuation appears once: a failed cross-check "could indicate need for regularization or continuation" (§2.3, p. 4). The general rule covers schedules by implication (inference): results must be reproducible "by meticulously listing all parameters, model and optimization settings, or … by including codes as a supplement (Haftka et al. 2019)" (§2.17, p. 10).
- **Non-discreteness (M_nd).** **Not in the document.** "Discreteness" appears only in a list of settings that "may strongly influence objective values" (§2.1, p. 2). The paper gives no measure and no formula.
- **Convergence and stopping (§2.8, p. 7).**
  - "the only rigorous way to conclude that a code has converged is to monitor the fulfillment of the optimality criterion. Everything else must be categorized as stopping criteria—not convergence."
  - "Sloppy stopping criteria are acceptable in many practical cases if motivated and acknowledged in the text".
  - When codes are compared for efficiency, an exact definition of the convergence criterion is "mandatory".
- **Objective values.**
  - Report thresholded as well as gray values. On the MBB example, C = 212.85 against C01 = 190.74, "more than 10% lower" (§2.1, p. 2).
  - "state the objective value for a uniform distribution of material" as a reference (§2.1, pp. 2–3). Example: C_solid = 128.36.
  - Use volume-corrected values when volumes differ (§2.5, p. 5).
  - Void stiffness and interpolation strongly affect objective values (§2.10, p. 8).
  - Check physical sanity (§2.15, p. 9).
  - "Time is over for making statements like "... results seem to be reasonable …"" (§3, p. 10).

**Fair comparison of designs**
- **Comparable settings.** "it is imperative that all settings are directly comparable" (§2.1, p. 2).
- **Thresholding.** "only approaches that use exactly the same material interpolation functions and filter settings are quantitatively comparable. If this is not the case, the only fair comparison is based on thresholded 0-1 solutions" (§2.1, p. 2).
  - The thresholding is volume-preserving (code snippet, after Sigmund & Maute 2013).
  - "For better reference performance", threshold Heaviside-projected designs from the 110-line code.
- **Re-evaluation and cross-checks.**
  - "authors must perform cross-checks whenever possible" (§2.3, p. 4). Table 1 evaluates each design under both boundary-condition models. An off-diagonal minimum "would indicate that the optimization was unsuccessful".
  - Post-evaluate with established solvers on body-fitted meshes. This "becomes even more important if designs are optimized using one software … and post evaluated using another software or using other parameter settings" (§2.5, p. 5).
- **Common evaluators.** The paper uses no such term. The idea is covered by §2.1 (comparable settings), §2.3 (cross-evaluation) and §2.5 (post-evaluation).
- **Benchmarks.**
  - The community "should define common benchmark examples" (abstract).
  - Use literature examples rather than new ones. "All implementation details must be included or better, reference to "benchmark examples" from the literature". "Any change to existing test cases must be motivated and explained". Three test cases are recommended (§2.9, pp. 7–8).
- **Other checks.**
  - Mesh-independent boundary conditions (§2.2).
  - Visualize the density "that is sent to the (FE-)analyser" (§2.4).
  - "always … include the imposed length-scale in illustrations" (§2.6).
  - Compare up-sampled designs with post-optimized ones (§2.7).
  - Microstructure results "must be carefully evaluated and checked for consistency". Claimed goals that the formulation omits, "(c.f. buckling resistance …)", must be "verified by post-analysis" (§2.11, p. 8).
  - "Large scale" means at least one million elements (§2.12). The reports never use the term.

**Code and data sharing**
- Reproducibility comes from listing all settings or from supplying code (§2.17, p. 10).
- Code papers are encouraged for "well-proven, stable and theoretically well-established codes". "Very specific codes or codes that are not theoretically well-founded or highly problem dependent are better left as supplementary material" (§2.17).
- The paper's own Replication of Results statement: "Matlab code snippets are included in the text" (p. 10).
- **Data sharing: not in the document.**

### B.2 Local minima, initial designs and continuation (for NDR-0010 and NDR-0029)

- **Local minima.**
  - "Non-convexity, premature stopping criteria, too strong penalization, programming errors or multiple other challenges may prevent convergence to good local minima" (§2.3, p. 4).
  - "one cannot check if a global optimum has been reached (without significant effort)" (§2.3).
  - "convergence to tree-like structures is witness of local minima" (§2.9, p. 8).
- **Initial designs.**
  - "it is not a shame to rerun the optimization with alternative starting guesses or neighboring solutions until the cross-check is satisfied as long as this is noted in the discussions" (§2.3, p. 4).
  - Qualified "de-homogenization-based starting guesses" are needed for some problems (§2.9).
  - In the force inverter, "the optimizer in the first step removes material from the input point", which may be "impossible to recover from" (§2.9).
- **Continuation.** Only the §2.3 remark above. There is nothing on schedules, on counting units, or on whether a step comes before or after the gradient.

For **NDR-0010**, the paper supports:
- treating different optima as an expected outcome;
- quantitative judgment;
- cross-evaluation;
- reruns from other starting designs, if they are disclosed.

It says nothing about round-off ensembles or statistical acceptance bands.

For **NDR-0029**, the paper is **neutral**. It supports only the reporting of the schedule and its counting unit (§2.17).

### B.3 Claims attributed to Sigmund 2022, and verdicts

| # | Where | Claim | Verdict | Evidence |
|---|---|---|---|---|
| SG-1 | Strategy report L37 (methodologist) | Sigmund 2022 anchors: "the result is a function of the optimization path"; "Count continuation in design-loop iterations, as every source does"; "compare designs after thresholding with a common evaluator" | **Partly confirmed.** Continuation counting is **not in the document.** | Thresholding: §2.1. Cross-checks and local minima: §2.3. Continuation appears only in §2.3, without a unit. |
| SG-2 | Report L120 (L4 row, Basis) | Basis for thresholded re-analysis, M_nd and IoU | **Thresholded re-analysis confirmed; M_nd and IoU not in the document** | §2.1. M_nd is Sigmund 2007. IoU comes from ML preprints. |
| SG-3 | Report L178 | The quote; "abstract; full text not read"; the frozen evaluator, thresholding, M_nd and IoU are the report's own recommendation | **Quote confirmed.** The label is out of date. Thresholding can now be credited to §2.1. A frozen evaluator is **consistent** with §§2.1, 2.3 and 2.5 but is never named. M_nd and IoU are **not in the document.** | Abstract, p. 1. |
| SG-4 | Task wording | "fair evaluation is the authors' responsibility" | **Confirmed as a paraphrase** | Abstract, conditional: "In the lack or impracticality of common benchmarks, the responsibility of fair evaluation of contributions is up to authors". Conclusions, unconditional: "It is authors' responsibility to make quantitative judgements and fair comparisons to the state of the art" (§3, p. 10). |
| SG-5 | Report L586 (revision log) | Sigmund 2022 marked "(abstract only)" because "no legal full text exists" | **Outdated.** Still true of open access, but the full text has now been read. | Add a new log entry (C.1, edit 16). |
| SG-6 | `verification_of_refactoring.md` L134 | "how to fairly evaluate and compare new concepts and ideas to existing ones is an open question"; the responsibility sentence; "lists recommendations, tools, and concepts"; *SMO* 65(11):315 | **Confirmed** | Abstract, p. 1. The PDF prints 65:315 and no issue number. |
| SG-7 | `verification_of_refactoring.md` L135 | Secondary sources: (1) convert designs to binary before comparing; (2) judge success by whether thresholded designs are close | (1) **Confirmed** from the primary text (§2.1). (2) **Needs correction.** | Sigmund compares thresholded **objective values**, not design geometry. Up-sampling fails if (C01₁ − C01₂)/C01₂ > 5% and succeeds at 0.25% (§2.7, p. 7). See also the cross-check table (§2.3). |
| SG-8 | `verification_of_refactoring.md` L155 | Different local optima inside the legacy spread count as "same result", per "Sigmund's (2022) burden-of-fairness principle" | **Consistent** | §2.3 allows reruns from other starting guesses "as long as this is noted". |
| SG-9 | `verification_of_refactoring.md` L160 and L395 (gaps) | M_nd reporting, post-evaluation and code sharing are known only via the abstract | **Resolved** | No M_nd recommendation. Post-evaluation: §2.5. Code: §2.17. |
| SG-10 | `verification_of_refactoring.md` L245 (d) | "the fair-comparison expectations of Sigmund (2022)" | **Confirmed** | §2.1, §2.3, §3. |
| SG-11 | `verification_of_refactoring.md` L297 | The quotes; secondary sources cite the paper for binary conversion and "closeness of thresholded designs" | **Quotes confirmed; second sentence as SG-7** | — |
| SG-12 | `verification_of_refactoring.md` L345 | "Following Sigmund (2022) … same mesh, volume fraction, and evaluator; thresholded designs re-analyzed; M_nd reported" | **Confirmed except M_nd**, which is **not in the document** | §2.1 (comparable settings, thresholding). §2.5 (volume correction when volumes differ). |
| SG-13 | `verification_of_refactoring.md` L418 | Outcomes are path-dependent and non-unique; fair comparison is the authors' job; compare thresholded designs | **Confirmed** for the Sigmund parts | §2.3, §2.9, the abstract, §3 and §2.1. |
| SG-14 | `numerical_method_choices.md` L213 | *SMO* 65(11) | **Confirmed** | 65:315, 10 pp. The issue number comes from Crossref, not the PDF. |
| SG-15 | `numerical_method_choices.md` L218 | "lists recommendations for reproducible TO"; the abstract sentence, via a search summary | **Confirmed**; upgrade to primary | Abstract; §2.17. |
| SG-16 | `numerical_method_choices.md` L243 (gap) | Its advice on continuation and reporting is unconfirmed | **Resolved** | No continuation advice beyond §2.3. Reporting: §2.17, §2.8, §2.1. |
| SG-17 | `open_access_optimization_papers.md` §3, S-a to S-e, and final table row 3 | S-b, S-c and S-d "not verifiable"; "Needs the user" | **Superseded** | S-a: confirmed. S-b: thresholding confirmed (§2.1); frozen evaluator, M_nd and IoU not in the document. S-c: thresholded comparison confirmed; continuation counting not in the document. S-d: as S-b. S-e: binary conversion confirmed (§2.1). |

### B.4 Decision records that Sigmund 2022 supports or contradicts

**Supports**
- **NDR-0010.** Local minima are expected, and comparisons must be quantitative and fair (§§2.3, 2.9, 3). Add a cross-evaluation and the disclosure of reruns from other starting designs (§2.3).
- **L4 acceptance.** Thresholded re-analysis (§2.1).
- **NDR-0044 and the methods paragraph.** List all settings (§2.17). Call a `maxiter` stop a stopping criterion, not convergence (§2.8).
- **NDR-0008 and NDR-0043.** Use literature test cases, and motivate every change to them (§2.9). This supports the Wallin & Tortorelli golden test.
- **Addendum NDR-0057 and NDR-0059.** A claimed goal such as buckling resistance must be in the formulation or verified by post-analysis (§2.11). Post-evaluation on body-fitted meshes (§2.5).
- **NDR-0064.** Check microstructure results for consistency (§2.11). Sigmund's Hashin–Shtrikman interpolation, Eq. 3, assumes a base Poisson's ratio of 1/3 and is a modulus interpolation. It is not a substitute for NDR-0064's bounds.
- **`verification_of_refactoring.md` L241 ("not a paper about the code").** Supported by §2.17: very specific codes "are better left as supplementary material".

**Neutral**
- NDR-0027 and NDR-0028: rounding is not discussed.
- NDR-0029: there is no advice on counting continuation.
- NDR-0035: only "list all settings".

**Tension (no direct contradiction)**
- **NDR-0033 and NDR-0044 (ε = 5e−4).** §2.10 (p. 8) says the 88-line default is "E_min/E₀ = 10⁻⁹" and "values higher than 10⁻⁶ are strongly discouraged". Fig. 7 shows curved members caused by "too high stiffness in void regions".
- The context is linear compliance. Wallin & Tortorelli's 1e−5 also exceeds 10⁻⁶. In finite strain, larger floors are used for Newton convergence (inference).
- Recommendation: state ε, and post-evaluate thresholded designs at a lower void stiffness, where Newton converges.

**Not contradicted, but over-attributed**
- Continuation counting (L37).
- M_nd and IoU (L120, L178, notes L345).

---

## Part C. Recommended edits

### C.1 Strategy report (`reports/NEMOpt refactoring strategy literature.md`; edit in place)

Suggested Svanberg 2007 link: the dead KTH URL, with the label "(KTH notes; user-supplied copy read)". The Wayback snapshots in the earlier note were not compared with the supplied copy.

1. **L37, panel, methodologist, "What they would say".**
   - Old: "Count continuation in design-loop iterations, as every source does, and compare designs after thresholding with a common evaluator."
   - New: "Count continuation in design-loop iterations, as every source that states a counting unit does, and compare designs after thresholding with a common evaluator."
2. **L37, "Anchoring evidence".**
   - Old: "[Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2) (abstract only)"
   - New: "[Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), §§2.1 and 2.3 (thresholded comparison, cross-checks, local minima; it does not say how continuation is counted)"
3. **L120, L4 row, "Basis".**
   - Old: "[Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2) (abstract only)"
   - New: "[Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), §2.1 (thresholded re-analysis; M_nd and IoU are not in it)"
4. **L178.**
   - Old: "Sigmund notes that, without practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors" ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), abstract; full text not read). This report therefore recommends re-analyzing final designs with one frozen evaluator, thresholding them, and comparing them on M_nd ([Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x)) and IoU."
   - New: "Sigmund notes that, without practical common benchmarks, "the responsibility of fair evaluation of contributions is up to authors" ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), abstract and §3). Unless interpolation and filter settings are identical, "the only fair comparison is based on thresholded 0-1 solutions" (§2.1). Designs obtained under different settings should be evaluated under each other's settings (§2.3, Table 1). This report therefore recommends three steps: re-analyze final designs with one frozen evaluator; threshold them with a volume-preserving threshold; and evaluate legacy and refactored designs under each other's evaluator. The designs are then compared on M_nd ([Sigmund 2007](https://doi.org/10.1007/s00158-006-0087-x)) and IoU. The frozen evaluator, M_nd and IoU are this report's choices; Sigmund 2022 names none of them."
5. **L208, rounding row, "Effect if changed".**
   - Old: "**A.** The returned value carries ±5e−4 noise, far above the slack epsimin = 1e−7 that the GPL Python port of Svanberg's MATLAB code adds in its GCMMA conservativeness check (also the port's subproblem tolerance; Svanberg's 2007 notes not re-read); NonconvexMMA's check has no tolerance at all ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python), `concheck` and `examples/gcmma_beam.py`)."
   - New: "**A.** The returned value carries ±5e−4 noise, which late in a run exceeds the decrease the approximation predicts. Svanberg's GCMMA accepts a step only if f̃_i(x̂) ≥ f_i(x̂) for all i, with no tolerance ([Svanberg 2007](https://people.kth.se/~krille/mmagcmma.pdf), §4, p. 6; KTH notes, user-supplied copy read), and NonconvexMMA's test likewise has none. The 1e−7 in Svanberg's notes only stops the subproblem solver (§5.5, p. 15); the GPL Python port reuses it as a slack in `concheck` ([GCMMA-MMA-Python port](https://github.com/arjendeetman/GCMMA-MMA-Python)). A value rounded to 1e−3 is also not the continuously differentiable function that MMA assumes (§1, p. 1). In plain MMA (full domain) the value enters the subproblem only as a dropped constant (§5, p. 9). There, the path changes only through the rounded ν and the stopping tests (inference)."
   - Keep the text from "The gradient factor error …" onward unchanged.
6. **L211, optimizer row, "Status and origin".**
   - Old: "No source documents `s_init` or `maxinner`."
   - New: "No source documents the values `s_init` = 0.25 or `maxinner` = 40. Svanberg's default for the matching asymptote parameter (asyinit) is 0.5, and lowering it gives "a more conservative method". His notes set no cap on inner iterations, which end when the approximations are conservative ([Svanberg 2007](https://people.kth.se/~krille/mmagcmma.pdf), §3, p. 5; §4, p. 6)."
7. **L212, clip row.**
   - Old: "([Svanberg 2002](https://doi.org/10.1137/S1052623499362822); assumption not re-verified)"
   - New: "([Svanberg 2002](https://doi.org/10.1137/S1052623499362822); the assumption is stated for the problem class in Svanberg 2007, §1, p. 1; the 2002 theorem was not re-read)"
8. **L84, resolved-gap row.** After "under MMA02 it depends on how many iterations are rejected.", insert:
   - "Unlike Svanberg's GCMMA, whose inner iterations re-solve from the unchanged x^(k) and need only function values ([Svanberg 2007](https://people.kth.se/~krille/mmagcmma.pdf), §4, p. 6), NonconvexMMA moves to each trial point. It then rebuilds the approximation there from a new value and Jacobian, so every evaluation after the setup evaluation is a design update. A "rejected" step is still taken; only ρ rises (`updatex!` and `updateapprox!(dualmodel, optimalx, fg, ∇fg)`, `mma_algorithm.jl` L287 and L310 at v1.1.2; source reading)."
9. **L239, parameter table, Optimizer row, "Treat as".**
   - Old: "Per preset; record `s_init` and `maxinner` as undocumented values"
   - New: "Per preset; record `s_init` = 0.25 (Svanberg's default asyinit is 0.5) and `maxinner` = 40 (no cap in Svanberg's notes) as undocumented values"
10. **L369, methods paragraph.**
    - Old: "The table should include the FEAD2026 values no source documents: ε = 5e−4, β̂ = 200, `s_init` and `maxinner`."
    - New: "The table should include the FEAD2026 values that no source documents (ε = 5e−4, β̂ = 200, `s_init` and `maxinner`). It should also give the continuation schedule with its counting unit, the stopping rule, and J for a uniform design as a reference value. The stopping rule is called a stopping criterion unless the KKT test was met ([Sigmund 2022](https://doi.org/10.1007/s00158-022-03427-2), §§2.1, 2.8 and 2.17)."
11. **NDR-0010 (L438).**
    - Old: "following SMO's advice to make repeated, slightly perturbed runs to show scatter (Haftka et al. 2019, p. 407), plus the article values at printed precision"
    - New: "following SMO's advice to make repeated, slightly perturbed runs to show scatter (Haftka et al. 2019, p. 407). Add a cross-check that evaluates legacy and refactored final designs under each other's evaluator, and report any rerun from other starting designs (Sigmund 2022, §2.3). Compare with the article values at printed precision"
12. **NDR-0027 (L455).**
    - Old: "| NDR-0027 | Objective value returned to the optimizer (likely defect) |"
    - New: "| NDR-0027 | Objective value returned to the optimizer (likely defect: a value rounded to 1e−3 is not the continuously differentiable function MMA assumes, Svanberg 2007, §1; under MMA02 it enters the conservativeness test, §4; under MMA87 it is only a dropped constant, §5) |"
13. **NDR-0028 (L456).**
    - Old: "| NDR-0028 | ν used in ∇J (likely defect) |"
    - New: "| NDR-0028 | ν used in ∇J (likely defect: the gradient then matches neither the rounded nor the unrounded J, while MMA's approximations must reproduce f and ∇f at the current point, Svanberg 2007, Eq. 4.10) |"
14. **NDR-0029 (L457).**
    - Old: "The article, thesis and Wallin & Tortorelli all count design iterations."
    - New: "The article, thesis and Wallin & Tortorelli all count design iterations. Svanberg's notes index iterations by the outer iteration, each with one accepted point (Svanberg 2007, §4). They are silent on continuation. For changed problem data they advise solving again from the latest solution (§2), which restarts the asymptotes."
15. **NDR-0035 (L463).**
    - Old: "Full-domain MMA87 matches Wallin & Tortorelli's plain MMA; half-domain MMA02 follows the thesis."
    - New: "Full-domain MMA87 matches Wallin & Tortorelli's plain MMA; half-domain MMA02 follows the thesis. Svanberg's own defaults are asyinit 0.5, asyincr 1.2, asydecr 0.7, move 0.5, albefa 0.1 and raa0 1e−5, with no inner-iteration cap (Svanberg 2007, p. 5 and §4). NonconvexMMA maps the first three to `s_init`, `s_incr` and `s_decr`. It has no `move` or raa0 term, and in MMA02 it moves to every trial point (source reading, v1.1.2). Any "default parameters" statement must say whose defaults it means."
16. **NDR-0044 (L472).** Before "(SMO: "detailed settings …", insert:
    - "State whether each run ended on the KKT test or at `maxiter`, and call the latter a stopping criterion, not convergence (Sigmund 2022, §2.8). Give ε = 5e−4 together with Sigmund's warning that void stiffness above 1e−6 is "strongly discouraged" in compliance problems (§2.10). Post-evaluate thresholded designs at a lower void stiffness wherever Newton converges (inference)."
17. **Revision log.** Keep the existing entries at L586, L589, L593 and L643, because they record what was done then. Add a new subsection:
    - "**Applied after reading Svanberg 2007 and Sigmund 2022 (user-supplied full texts, 2026-10-12).** Sigmund 2022 '(abstract only)' replaced by section citations (§§2.1, 2.3, 2.8, 2.17); M_nd, IoU and continuation counting are not in it. The 1e−7 is now attributed to Svanberg's subproblem stop (§5.5); his conservativeness test has no tolerance (§4). Svanberg's defaults and the continuity assumption are cited from the notes. The earlier rule 'cite no Svanberg 2007 parameter value' is lifted. Source: `svanberg_2007_sigmund_2022.md`."

**Optional addendum edit.** In `reports/NEMOpt buckling and stability addendum.md`, NDR-0057, Origin column:
- Old: "Synthesis; precedent Zhang & Khandelwal 2019"
- New: "Synthesis; precedent Zhang & Khandelwal 2019; Sigmund 2022, §2.11 (goals left out of the formulation must be verified by post-analysis)"

### C.2 Research notes

These notes keep their original text and append a "Corrections after full-text verification (2026-10-12)" section. Add the items below to that section. Each item supersedes the earlier item on the same line.

**`research_notes/NEMOpt refactoring strategy literature/numerical_method_choices.md`.** Add a block "**Svanberg (2007) and Sigmund (2022), read in full (user-supplied copies)**". It replaces the blocks "Svanberg's GCMMA tolerance" (L487–490) and "Sigmund (2022)" (L483–485).

1. **L24.** Append: "Svanberg's own default for the asymptote parameter matching `s_init` is 0.5; lowering it gives 'a more conservative method' (Svanberg 2007, p. 5)."
2. **L25.** Corrected: "The rounding error (≈5e-4) exceeds the late-run predicted decreases. GCMMA's conservativeness test has no tolerance (Svanberg 2007, §4, p. 6), and a value rounded to 1e-3 is not continuously differentiable (§1, p. 1)."
3. **L213.** Bibliographic entry: "*SMO* 65:315 (10 pp.); issue 11 per Crossref, not printed in the PDF".
4. **L218.** Corrected: "Sigmund 2022 (read in full) asks authors for four things:
   - to list "all parameters, model and optimization settings" or to supply code (§2.17);
   - to compare only thresholded 0-1 designs unless interpolation and filter settings are identical (§2.1);
   - to cross-check (§2.3);
   - to call a non-KKT stop a stopping criterion (§2.8).

   It mentions continuation once, as a possible remedy when a cross-check fails (§2.3), and gives no advice on schedules or counting."
5. **L240 and L244 (asymptote reset; gap).** Append: "Svanberg 2007 is silent on continuation. For changed problem data it advises solving again from the latest solution (§2, p. 3), which restarts the asymptotes (Eq. 3.11); GCMMA's ρ restarts at every outer iteration anyway (Eq. 4.6). Still no peer-reviewed source."
6. **L243 (gap).** Resolved: "Read in full: no advice on continuation beyond §2.3; reporting advice in §§2.1, 2.8 and 2.17."
7. **L252.** Corrected: "GCMMA accepts an iterate only if f̃_i(x̂) ≥ f_i(x̂) for every i, with no tolerance (Svanberg 2007, §4, p. 6). The 1e-7 in the notes stops the interior-point subproblem solver (§5.5, p. 15); the Python port reuses it as a slack (epsimin) in `concheck`. In Svanberg's GCMMA an inner iteration needs function values only (p. 6). Nonconvex.jl evaluates values and gradients because it moves to each trial point and re-centers there (`mma_algorithm.jl` L287, L310)."
8. **L256–258 ([PR-bib]).** Upgrade: "The assumption is stated for the problem class: f₀ … f_m are 'continuously differentiable' (Svanberg 2007, §1, p. 1). The notes give no theorem; they refer to Svanberg 2002."
9. **L259–263.** Append: "Confirmed in Svanberg 2007: Eqs. 3.6–3.7 (albefa 0.1, move 0.5); 3.11–3.14 (asyinit 0.5, asydecr 0.7, asyincr 1.2, bounds 0.01 and 10); 3.3–3.4 (raa0 1e-5); 4.6 (ρ start, floor 1e-6); 4.8–4.9 (ρ update). All are called defaults on pp. 5 and 7. The `concheck` slack (L262) is the port's addition. NonconvexMMA has no `move` and no raa0 term, and starts each outer iteration with ρ = max(ρ/10, 1e-5) rather than Eq. 4.6."
10. **L281.** Append: "The notes were read from a user-supplied copy (15 pp., undated; PDF built 2014-10-14)." This supersedes "Do not cite Svanberg 2007 for any parameter value until a snapshot has been read."
11. **L285.** Corrected: "Late in the run, predicted decreases fall below ~1e-3. The ±5e-4 noise is then comparable to or larger than the predicted decrease, and neither Svanberg's test (2007, §4) nor Nonconvex.jl's has a tolerance."
12. **L286.** Append: "In NonconvexMMA a failed test does not undo the step: the design moves to the trial point, and only ρ rises."
13. **L289.** Append: "The value itself enters the MMA subproblem only through the constant r₀, which is dropped (Svanberg 2007, §5, p. 9). In exact arithmetic, MMA87 iterates depend on the rounded ν in the gradient, not on the rounded J (inference; confirm that `KKTCriteria` uses no objective value)."
14. **L298 (gap).** Resolved: "§1 (continuously differentiable); §4 (no tolerance, no inner cap, 'finite (usually small)' count); §5.5 (subproblem stop ε ≤ 1e-7)."

**`research_notes/NEMOpt refactoring strategy literature/verification_of_refactoring.md`.** Append to its Corrections section.

1. **L134.** "Confirmed from the full text (abstract, p. 1); 10 pp.; the PDF prints 65:315 without an issue number."
2. **L135 and L297 (secondary practices).** Replace with primary citations:
   - "binary comparison: 'the only fair comparison is based on thresholded 0-1 solutions' unless interpolation and filter settings are identical (§2.1, p. 2), with volume-preserving thresholding;
   - closeness: Sigmund compares thresholded objective values, not geometry. Up-sampling fails if (C01₁ − C01₂)/C01₂ > 5% (§2.7, p. 7); designs are also cross-evaluated (§2.3, Table 1)."

   This supersedes "remain secondary" in the earlier correction item.
3. **L155.** Append: "Sigmund allows reruns 'with alternative starting guesses … as long as this is noted in the discussions' (§2.3, p. 4)."
4. **L157.** Append: "Sigmund 2022, §2.17: list 'all parameters, model and optimization settings' or supply the code."
5. **L160 and L395 (gaps).** Resolved:
   - "Read in full (user-supplied copy). It contains no M_nd recommendation; M_nd is Sigmund 2007.
   - Post-evaluation: §2.5. Code: §2.17. Data sharing is not discussed.
   - It still has no legal open copy."
6. **L345.**
   - Old: "Following Sigmund (2022), compare against the FEAD 2026 designs fairly: same mesh, volume fraction, and evaluator; thresholded designs re-analyzed; M_nd reported."
   - New: "Following Sigmund (2022, §§2.1, 2.3 and 2.5), compare against the FEAD 2026 designs fairly: same mesh, volume fraction and evaluator; thresholded designs re-analyzed and cross-evaluated; volume-corrected values if volumes differ. Report M_nd as well (Sigmund 2007; not part of Sigmund 2022)."
7. **L418.** Append the section numbers "(Sigmund 2022, abstract, §§2.1, 2.3, 3)".
8. **New item for L241 (the next paper is not a code paper).** "Consistent with Sigmund 2022, §2.17: very specific codes 'are better left as supplementary material'."

**`research_notes/NEMOpt refactoring strategy literature/reproduction_baseline.md`.** This note has no Corrections section; add one.
- **Item 12 (L175).** Append: "Svanberg's own default for the matching asymptote parameter (asyinit) is 0.5, and his notes set no inner cap (Svanberg 2007, p. 5 and §4)."
- **Item 2 (L160).** Append: "In NonconvexMMA every inner iteration moves the design to the trial point (`mma_algorithm.jl` L287, L310). Unlike Svanberg's GCMMA (2007, §4), a failed test does not reject the step."
- **L125.** Append: "The 1987 values are not in the supplied 2007 notes, which describe the 2007 implementation."

**`research_notes/NEMOpt refactoring strategy literature/wallin_tortorelli_2020.md`.** Append to its Corrections section.
- **L360.** "Svanberg's 1987 paper was still not read. His 2007 notes list asyinit 0.5, asyincr 1.2, asydecr 0.7, move 0.5, albefa 0.1 and raa0 1e-5 as the 2007 implementation's defaults (p. 5)."
- **L350.** Append: "NEMOpt's `s_init` = 0.25 halves Svanberg's asyinit default; the notes call that direction more conservative (p. 5)."

**`research_notes/NEMOpt refactoring strategy literature/julia_architecture.md`.** This note has no Corrections section; add one.
- **L196.** Append: "Svanberg 2007 (§2, p. 3) advises solving again from the latest solution when problem data change. The notes also specify a complete subproblem solver (§5), useful for an in-house implementation."
- **L259, item 2.** Append: "For MMA87 the value enters only as a dropped constant (Svanberg 2007, §5); the rounded ν in the gradient is what changes the path."

### C.3 Earlier verification note (`research_notes/NEMOpt supplied papers verification/open_access_optimization_papers.md`)

1. **§3 header and summary row 3.**
   - Old: "**No** (abstract only)"
   - New: "**User-supplied version of record read in full** (no legal open copy); see `svanberg_2007_sigmund_2022.md`"
   - Bibliography: add "10 pp.; the PDF prints 65:315 without an issue number".
2. **§3 table.**
   - S-b: "partly confirmed (thresholding, §2.1)".
   - S-c: "partly confirmed (thresholded comparison, §2.1; continuation counting not in the paper)".
   - S-d: as S-b.
   - S-e: "confirmed (§2.1)".
3. **§5 header and summary row 5.**
   - Old: "**No** (copies exist in the Wayback Machine but could not be reached)"
   - New: "**User-supplied copy read in full** (15 pp.)"
4. **§5 table.**
   - V-a: "not in Svanberg 2007; 1e−7 is the subproblem stop (§5.5)".
   - V-b: "confirmed (pp. 4–5, 7–8)".
   - V-c: "partly answered (§2, p. 3; p. 5)".
   - V-d: "resolved".
5. **§5 edit 5.** Mark it "superseded: the notes have been read".
6. **Final table, rows 3 and 5.** Status "Obtained (user-supplied)".

### C.4 What should not change

- **Thresholding and the frozen evaluator in the reports.** Keep them. Thresholding now has a primary source (§2.1). The frozen evaluator stays labeled as the report's own choice.
- **The decisions in NDR-0027, NDR-0028, NDR-0029 and NDR-0035.** Neither document overturns them; the edits above only add support or precision.
- **Effect class A in NDR-0027.** Keep it, because the report's classes are bitwise and MMA02 is affected. Only its explanation for MMA87 changes.

---

## Summary

**Bibliographic facts.**
- Sigmund 2022: *SMO* 65:315, article 315, **10 pages** (not 24). The issue number is not printed; Crossref gives 11.
- Svanberg 2007: 15 pages, undated, describing the 2007 Matlab implementations.

**Svanberg 2007.**
- The conservativeness test f̃_i(x̂) ≥ f_i(x̂) for all i has **no tolerance** (§4, p. 6). The only 10⁻⁷ stops the subproblem solver (§5.5, p. 15).
- The reports' "Svanberg 1e−7 conservativeness tolerance" is therefore **not in the document**. It belongs to the Python port.
- The six defaults (raa0 1e−5, albefa 0.1, move 0.5, asyinit 0.5, asydecr 0.7, asyincr 1.2) are **confirmed**. The notes call them defaults and say they "can be carefully changed".
- NEMOpt's `s_init` = 0.25 halves asyinit, which is the "more conservative" direction. `maxinner` = 40 has no counterpart: the notes set no cap.
- Inner iterations re-solve at a fixed x^(k) and need **values only**, with gradients once per outer iteration (p. 6).
- The notes are silent on continuation. For changed problem data they advise a restart from the latest solution (§2).
- A rounded J **violates the continuous-differentiability assumption** (§1). In plain MMA, however, the value is a dropped constant (§5). The full-domain effect therefore runs through the rounded ν (NDR-0028) and the stopping tests.
- Side finding (source reading): NonconvexMMA's MMA02 moves to every trial point and never rejects a step, unlike Svanberg's GCMMA.

**Sigmund 2022.**
- The responsibility quote is **confirmed**: conditional in the abstract, unconditional in §3.
- Thresholded comparison is **confirmed** (§2.1).
- Cross-evaluation and post-evaluation are **confirmed** (§§2.3, 2.5).
- Stopping versus convergence: §2.8. Settings or code for reproducibility: §2.17.
- **Not in the document:** M_nd, IoU, a "frozen evaluator", how continuation is counted, and data sharing.
- On local minima and starting designs, §2.3 permits reruns if they are reported. This supports NDR-0010's spirit, not its statistics. It is neutral on NDR-0029.
- Tension: ε = 5e−4 against "values higher than 10⁻⁶ are strongly discouraged" (§2.10), which concerns linear compliance.

**Edits.**
- Update the report's L37, L120, L178, L208, L211, L212, L84, L239 and L369, and NDR-0010, NDR-0027, NDR-0028, NDR-0029, NDR-0035 and NDR-0044. Add a revision-log entry.
- Append the correction items to five research notes.
- Mark the earlier note's Sigmund and Svanberg items as superseded.
