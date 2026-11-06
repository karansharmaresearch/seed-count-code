# seed-count-code

[![DOI](https://zenodo.org/badge/1412788268.svg)](https://doi.org/10.5281/zenodo.24002843)

Code and data for two papers by Karan Sharma and Dheeraj Sonkhla:

- **Paper A:** *Inserting Zeros Detects Primes*
- **Paper B:** *Numbers that divide their zero-insertions: classification, maximal order and average order*

## The problem

In base B, take a two-digit number with digits a and b and insert zeros between them: ab, a0b, a00b, and so on. Call (a, b) a *seed* if the number ab divides every one of these. This happens exactly when Ba + b divides (B − 1)b.

- S(B) is the number of seeds in base B. This is OEIS [A394974](https://oeis.org/A394974). S(B) = 0 exactly when B is prime.
- Z_B is the set of all numbers (of any length) in base B that divide every zero-insertion after the first digit, and N(B) = |Z_B|. Z_10 is OEIS [A398346](https://oeis.org/A398346).
- T(x) = Σ_{B ≤ x} S(B) and W(x) = Σ_{B ≤ x} N(B).

## Layout

```
code/                  seed count S(B), the sum T(x), and checks for the papers
code/zero_insertion/   Z_B and N(B) (has its own README)
data/seed_count/       S(B) records and T(x) up to 10^8
data/zero_insertion/   N(B), smallest/largest members, max length, W(x) up to 10^7
checks/check3/         independent referee-style check of the average-order proof, with outputs
```

### code/

| File | What it does |
|------|--------------|
| `seed_count.py` | S(B) by the closed form and by brute force |
| `verify.py` | closed form vs brute force; S(B) = 0 iff B is prime |
| `largest_seed_check.c` | checks that the largest seed for composite B is B(B−1)/p (p the smallest prime factor) |
| `sqrt_upper_bound.py` | checks the bound S(B) < 2 τ(B) τ(B−1) √B |
| `three_digit_gaps.py` | three-digit case: which gap the zeros go into |
| `poly_base_seeds.py` | the same question with a polynomial base over F_q |
| `fast_seedcount.py` | fast S(B) and T(x) using divisor pairs of B and B−1 |
| `seed_sum.c` | exact T(x), single process |
| `seed_sum_par.c` + `run_par.py` | exact T(x) in parallel chunks (used
