# Code and data: zero-insertion divisors and OEIS A394974

Programs, verification scripts and data for two papers by Karan Sharma and
Dheeraj Sonkhla:

- **"Inserting Zeros Detects Primes"**, preprint (2026), Zenodo,
  https://doi.org/ZENODO_DOI_A
- **"Numbers that divide their zero-insertions: classification, maximal order
  and average order"**, preprint (2026), Zenodo, https://doi.org/ZENODO_DOI_B

Archived version of this repository: https://doi.org/CODE_DOI

## The objects

For a base `B >= 2`, a pair of digits `(a, b)`, `1 <= a, b <= B-1`, is a
**seed** if `Ba + b` divides every zero-insertion `aB^n + b` (`n >= 1`);
equivalently `(Ba + b) | (B-1) b`. `S(B)` is the number of seeds
([OEIS A394974](https://oeis.org/A394974)); for example, in base ten
45 divides 405, 4005, 40005, ... More generally, `Z_B` is the set of numbers
with at least two base-`B` digits, nonzero second and last digits, that divide
every insertion of zeros after their first digit, and `N(B) = |Z_B|`. In base
ten, `Z_10 = {15, 18, 45, 225, 675, 1125, 3375, 5625, 7875}`
([OEIS A398346](https://oeis.org/A398346)).

The computations are consistency checks: none of the proofs depends on them.
The proofs of the second paper have not yet been checked by a human referee;
an AI referee check of its average-order proof is recorded in
`checks/check3/REPORT.md`.

## Requirements

- Python 3.10 or later; `numpy` for `analyze_avg.py`, `fit_constant.py`,
  `verify_mainterm.py` and `checks/check3/c6_constants_mainterm.py`
  (`pip install numpy`). Everything else uses the standard library.
- A C compiler for the `.c` programs. Any C99 compiler works; we used
  `pip install ziglang` and `python -m ziglang cc -O2 -o X X.c`.

## Layout

| path | contents |
|------|----------|
| `code/seed_count.py` | `S(B)` from the closed form and by brute force |
| `code/verify.py` | closed form vs. brute force, primality characterization, lower bounds |
| `code/fast_seedcount.py` | `S(B)` in `O(tau(B) tau(B-1))` operations |
| `code/seed_sum.c`, `code/seed_sum_par.c`, `code/run_par.py` | exact `T(x) = sum_{B<=x} S(B)` (parallel; `-DRECORDS` with `code/records_S.py` gives the record bases) |
| `code/sqrt_upper_bound.py`, `verify_rewrite.py`, `verify_mainterm.py`, `verify_errors.py` | checks of Sections 4-5 of the second paper (maximal and average order) |
| `code/verify_merged.py` | checks of Sections 3 and 6-8 of the second paper |
| `code/largest_seed_check.c` | brute-force check of the largest-seed theorem (first paper, Theorem 12) |
| `code/three_digit_gaps.py`, `code/poly_base_seeds.py` | three-digit seeds and polynomial bases (first paper, Section 5) |
| `code/fit_constant.py`, `code/analyze_avg.py` | fits of `T(x)` against the main term |
| `code/zero_insertion/` | classification of `Z_B`, `N(B)` and its growth (own README) |
| `data/seed_count/` | `T(x)` checkpoints up to `10^8`, records of `S(B)`, largest-seed check output |
| `data/zero_insertion/` | b-files of `N(B)` and related sequences, `W(x) = sum N(B)` up to `10^7`, `N(B)/B` statistics |
| `oeis/` | the fields and b-file of A394974 |
| `checks/check3/` | independent re-check of the average-order proof: scripts `c1`-`c7`, their outputs, and `REPORT.md` |

## Reproducing the main numbers

Run from the repository root unless stated otherwise.

```
python code/verify.py
python code/verify_rewrite.py
python code/verify_mainterm.py
python code/verify_errors.py
python code/sqrt_upper_bound.py
python code/verify_merged.py
cd code/zero_insertion && python zi.py verify 1000 && python verify_growth.py
```

`T(x)` up to `10^8` (compare with `data/seed_count/partial_sums_1e8.txt`;
`T(10^6) = 97603377`, `T(10^7) = 1475926653`, `T(10^8) = 21187268282`):

```
python -m ziglang cc -O2 -o seed_sum_par code/seed_sum_par.c
python code/run_par.py ./seed_sum_par 100000000 8 partial_sums.txt
```

`checks/check3/c5_seg.c` computes the same sums by a different algorithm;
`c5_seg_combine.py` expects its four range outputs in a folder `seg1e8/`
(see the docstring), which is not included.

`W(x) = sum_{B<=x} N(B)` up to `10^7` (compare with
`data/zero_insertion/growth_1e7.txt`):

```
cd code/zero_insertion
python -m ziglang cc -O2 -o n_fast n_fast.c
./n_fast 10000000 growth.txt
```

## Notes

- `code/zero_insertion/zi.py` uses older letters (base `b`, member `k`,
  multiplier `c`); the papers use `B`, `z`, `k`.
- File paths quoted inside `checks/check3/REPORT.md` refer to the authors'
  working repository at the time of the check, not to this one.

## License and citation

MIT License (see `LICENSE`). If you use this code or data, please cite the
papers above and the archived version, https://doi.org/CODE_DOI (see
`CITATION.cff`).
