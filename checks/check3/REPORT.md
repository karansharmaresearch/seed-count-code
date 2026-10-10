# Referee check 3: `paperB/average-order.tex`

Scope: every lemma, proposition, theorem and corollary of Sections 4 and 5,
checked line by line, plus the statements in Sections 1–3 and 6–7 that depend
on them. The closed form (1.1) is taken as given, as instructed. I also tested
it numerically (see c1). No existing file was edited; everything is in
`checks/check3/`.

## Verdict

**(b) Correct, with fixable gaps and exposition issues.**

I found no error that breaks either result. Theorem 4.1 and Theorem 1.2 (with
the constant G/12 and the error term O(x log² x)) follow from the arguments as
written, apart from one small gap about support when n = 1 (Issue 1). The
missing step takes one line to supply. Every inequality direction, change of
summation order, constant and edge case that I checked holds, both by hand and
numerically. The independent computations agree exactly with the paper's
numbers: T(x) up to 10⁷ by a method that does not use Section 5, and the exact
identity T = M + E₁ + E₂.

## Issues

**1. Support of the n = 1 terms is not covered by the lemma cited.**
Proposition 5.5 (proof): *"By Lemma 5.4(d) and the definition, a term of E₁ is
nonzero only if u′v′ ≤ x and t ≤ 2x"*. Proposition 5.10 (proof): *"By Lemma
5.4(d) … nonzero only if n ≤ x and t ≤ T_n, where T₁ = x"*.

- **What is wrong:** Lemma 5.4(d) bounds t only for n ≥ 2. For n = 1 the
  inequality it gives, t(n−1) ≤ x−n, reads 0 ≤ x−1 and puts no condition on t.
  So the bound t ≤ T₁ = x (and t ≤ 2x) for u′ = v′ = 1 is not justified by the
  cited lemma. This matters in two places: it is what makes E₁ a finite sum,
  and it is what makes the t-sum in Proposition 5.10 have length x.
- **Is the claim true?** Yes. For n = 1 we have ε = 1, so A = Q(Q−1)/2, and
  A ≠ 0 forces Q ≥ 2. Condition (5.4) with q = 2, n = 1, ρ = 0 is
  2(t+1) ≤ x, so t ≤ x/2 − 1. Also, the main term x/(2t)·1[t ≤ x] vanishes
  for t > x. (Alternatively, Lemma 5.4(b) gives Q ≥ 2 ⇒ R ≥ 1 ⇒ t ≤ x.)
  Numerically, the largest t with A₁,₁(t) ≠ 0 is exactly x/2 − 1 (c3).
- **Severity:** minor fixable gap.
- **Fix:** add to Lemma 5.4(d): *"if n = 1 and A ≠ 0 then Q ≥ 2, hence
  t ≤ x/2 − 1"*, and cite it in Propositions 5.5 and 5.10.

**2. Remark 5.14 overstates what the scripts check.**
Remark 5.14: *"Lemma 5.4 was checked for every relevant (u′,v′,t) at x ≤ 1000
… Lemma 5.12 for all n ≤ x ≤ 5·10⁴"* and *"checked against direct evaluation
… for x ∈ {…, 2000, 3000}"*.

- **What is wrong:**
  - `code/verify_errors.py` checks Lemma 5.4 only at x ∈ {7, 20, 61, 150, 400, 1000}.
  - It checks Lemma 5.12 only at x ∈ {100, 1000, 10000, 50000}, for all n ≤ x at those x.
  - `code/verify_rewrite.py` compares T(2000) and T(3000) with hard-coded
    values, not with a direct evaluation.
  - The scripts also use an older numbering (they cite "Lemma 4.3",
    "Theorem 2.4" and so on, from `docs/average_order.tex`).
  - `verify_errors.py` tests the older, weaker E₁ constant (9√2+7), not
    the paper's (7√2+5). The paper's constant is tested in c4 here.
- **Severity:** exposition. The remark itself says these checks are not used
  in the proofs.
- **Fix:** state the exact x values checked, and update the scripts'
  numbering and constants to match paperB.

**3. log S(B) is undefined at primes (Corollary 4.2).**
Corollary 4.2: *"lim sup_{B→∞} log S(B)/log B = ½"*.

- **What is wrong:** S(B) = 0 for prime B, so log S(B) is undefined there.
- **Severity:** exposition.
- **Fix:** restrict to composite B, or state the result as
  lim sup log(1+S(B))/log B = ½.

**4. Letters are reused for different things.**

- **What is wrong:**
  - t is (B−1)/E in the proof of Theorem 4.1 but the bijection parameter in
    Section 5.
  - q is (k+1)/n = uv/k in Lemma 5.2 but the progression index in Lemma 5.4.
  - γ is Euler's constant, while γ₁ is a fit coefficient in Section 6.
  - D is a gcd in Section 4 but the function D(z) = Σ τ(m)/m from Section 3
    onward.
- **Severity:** typo / exposition.
- **Fix:** rename, for example s in Theorem 4.1 Step 1, 𝔇(z) or D₂(z) for the
  function, and a₁, a₀ for the fit coefficients.

**5. Numerical remark on Lemma 5.4(c).**
Remark 5.14: *"largest observed ratio … is about 0.46"*.

- **What I found:** my exhaustive run (all x ≤ 90, all coprime (u′,v′), all
  t ≤ 3x+5) gives 0.382. `verify_errors.py` at its six x values reports 0.462.
- **Is there a contradiction?** No. The inequality has ample room in every test,
  including 2·10⁵ random large (x, n, t, ρ) with arbitrary ρ, where the worst
  ratio is 0.518.
- **Severity:** none. Recorded only for completeness.

No other issues. In particular, I checked the following points and found them
correct:

- **Theorem 4.1:**
  - Step 1: D | t−j, and j ≥ 1 gives D < t.
  - Step 2: both cases, including D > E ⇒ D ≥ 2 ⇒ k ≥ D/2.
  - Step 3: a nonzero term forces DE > k, and CRT then gives at most one k in [1, DE).
  - The strict inequality when S(B) = 0.
- **Lemma 5.1:**
  - (b): both gcd computations, using B−(k+1) = tuv and (B−1)−k = tuv.
  - (c): both directions, including B ≥ k+2 ≥ 3 for the inverse.
- **Lemma 5.2:**
  - the case n | k+1;
  - the case n ∤ k+1, with frac(q) ≥ 1/n > q/B;
  - n | k+1 ⇔ v′ = 1.
- **Proposition 5.3:** B ≤ x ⇔ (k+1)(tk+n) ≤ xn, which is an exact equivalence.
- **Lemma 5.4:**
  - (a): q ≥ 0 is automatic; the q = 0 term vanishes because ρ ≥ 2 there; the
    (n,q) = (1,1) phantom term k = 0 is harmless.
  - (b): the upper bound; the root formula for s; √(n²+4ntx) ≤ n+2√(ntx);
    s ≥ R − 1/t; the case s < 2.
  - (c): the constants 3.5R+2 and 3.5R+2.5.
  - (d): the case n ≥ 2.
- **Proposition 5.5:** the count of 2^ω(n) ordered coprime factorizations; the
  split into two finite sums is legitimate.
- **Lemma 5.6:** h(p) = −2/p, h(p²) = −(p−2)/p, h(p^a) = 0 for a ≥ 3, and the
  Euler factor 1 − 3/p² + 2/p³ = (1−1/p)²(1+2/p). Checked in exact arithmetic
  in c6.
- **Lemma 5.7** and **Proposition 5.8:** the expansion of D(y/c), the use of
  (1+log c)² ≪ c^{1/4}, and the tail Σ_{c>y} ≪ y^{−1/4}.
- **Corollary 5.9:** the interchange Σ_n ∫_n^x = ∫_1^x Σ_{n≤y}, giving
  G/6·log³x, then ×x/2.
- **Proposition 5.10:** Σ_{t≤T} t^{−1/2} ≤ 2√T; T_n ≤ 2x/n; the n = 1 bound
  9.5x ≤ (7√2+5)x.
- **Lemma 5.11:** the monotonicity of c_m; 0 < c_m − γ ≤ 1/m; the passage to
  real y.
- **Lemma 5.12:** Möbius inversion; H(Z/d) = 0 for d > Z; log Z ≥ 0.
- **Proposition 5.13:**
  - the bijection (n, squarefree d | n) ↔ (d, m) with dm ≤ x, and
    2^ω(dm) ≤ 2^ω(d)2^ω(m);
  - Σ₂: (x/2n)(2/Z) = 1, and d ≤ x/(dm) ⇔ d²m ≤ x ⇒ dm ≤ x;
    Σ_{m≤y}τ(m) ≤ y(1+log y) with y = x/d² ≥ 1;
  - Σ₃: log Z < log d.
- **Section 5.7:** the passage from integer x to real x costs O(log³x).
- **Section 1:** the claim S(B) ≥ √(2B)(1+o(1)) under the Sophie-Germain-type
  hypothesis matches Proposition 3.3.

## Steps I could not fully verify

- **Theorems 3.1, 3.2 and Proposition 3.3.** These are imported from the
  companion paper, whose proofs are not in this file. I tested them numerically
  only:
  - the seed definition versus (1.1) for B ≤ 60;
  - S(pq) ≥ 2 gcd(p−1, q−1) for pq ≤ 30000.
- **T(10⁸) = 21 187 268 282.** I reproduced every value up to 10⁷
  independently, but not 10⁸. The fits in Section 6 were recomputed from
  `code/partial_sums_1e8.txt`.
- **Literature statements in Section 2.** Ingham's constant 6/π², the Hardy–Wright
  theorem number (Thm 315, τ(n) = O(n^δ)), and the attributions were checked
  from memory only, not against the sources.
- **The open-questions claim "exactly one B ≤ 10⁸ with S(B)/B ≥ 4/9".** I
  verified it only for B ≤ 10⁷.

## Scripts run (all in `checks/check3/`, output in the matching `.out` file)

| script | what it tests | result |
|---|---|---|
| `c1_closed_form_and_thm41.py` | seed definition (all n up to 2m+64) vs (1.1) for B ≤ 60; Thm 4.1 Steps 1, 2, 3 and the final bound for every (B,k), 3 ≤ B ≤ 3000 | 0 mismatches, 0 violations; max summand/(2√DE) = 0.495; max S/(2ττ√B) = 0.051; S(16) = 6, S(36) = 16 |
| `c2_bijection_floor.py 1500` | Lemma 5.1(a),(b),(c) in both directions ((b) also for gcd(t,u′v′) > 1); Lemma 5.2; uv·n = k(k+1); the equivalence B ≤ x ⇔ (k+1)(tk+n) ≤ xn; B ≤ 1500 | 1 122 751 pairs = 1 122 751 admissible quadruples, 0 errors |
| `c3_lemma54.py 90` | Lemma 5.4(a)–(d) against a brute-force A, every x ≤ 90, every coprime (u′,v′) with n ≤ x+2, every t ≤ 3x+5; the support t ≤ 2x; the n = 1 support; 2·10⁵ random large cases for (b),(c) | 0 errors; max ratio 0.382 (random: 0.518); n = 1 largest nonzero t = x/2−1 (Issue 1) |
| `c4_decomposition.py` | T(x) computed three ways (closed form, Prop 5.3 grouping, M+E₁+E₂ from the definitions) for x ∈ {3,…,7, 10, 30, 97, 100, 257, 500, 1000}; the Prop 5.10 bound; the Prop 5.13 bound and each Σᵢ bound; Lemma 5.12 for all n ≤ x; Lemma 5.11 on a grid y < 3000 including points just below integers | all equal (e.g. T(1000) = 15901 three ways); all bounds hold with large margin (E₁ ≈ −13.2x vs bound 481x at x = 1000); Lemma 5.11 max ratio 0.27 |
| `c5_direct.c` (`./c5 1e7`) | T(x) by divisor-pair enumeration of (1.1), independent of Section 5; Thm 4.1 in exact integers for every B ≤ 10⁷; S(B)/B records | T(10⁶) = 97 603 377, T(10⁷) = 1 475 926 653 (both match the paper); 0 violations; max S/(ττ√B) = 0.102; only B = 36 has S/B ≥ 4/9; runner-up B = 16 (3/8) |
| `c6_constants_mainterm.py` | Lemma 5.6 local factors in exact arithmetic; G by Euler product to 2·10⁷; Prop 5.8 and Cor 5.9 numerically to 10⁷; the Section 6 table and fits | G = 0.286747428, G/12 = 0.023895619; the normalized errors of Prop 5.8 and Cor 5.9 stay bounded (≈0.85 and ≈0.24); table columns 0.0370/0.0352/0.0339 and 0.181/0.183/0.184 reproduced; fit β = 0.193, γ₁ = −0.145, δ = −0.163; free fit 0.0239; local θ = 2.68 |
| `c7_prop33.py` | Prop 3.3 for pq ≤ 30000; S/√B along the family m+1, 2m+1 prime | 0 failures; ratios 1.5–2.0 (≥ √2 asymptotically, as claimed) |
| existing `code/verify_rewrite.py`, `verify_mainterm.py` (run from `code/`), `verify_errors.py`, `sqrt_upper_bound.py` | as described in Issue 2 | all pass (`existing_scripts.out`). `verify_mainterm.py` fails when run from the repo root, because it opens `partial_sums_1e8.txt` by a relative path. |

## Follow-up (2026-09-29): fixes applied

- **Issues 1, 3 and 4:** fixed in `paperB/average-order.tex`.
- **Issue 2:** Remark 5.14 now states the exact ranges tested, and the repo
  scripts were updated to match paperB's numbering and constants.
- **T(10⁸) and the 4/9 claim up to 10⁸:** confirmed independently with
  `c5_seg.c` (4 ranges, combined by `c5_seg_combine.py`, output
  `c5_seg_1e8.out`).

All checks were rerun after the edits (`rerun/`) and all pass. Details are in
`docs/PROGRESS.md`.
