# seed-count-code

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
| `seed_sum_par.c` + `run_par.py` | exact T(x) in parallel chunks (used for x up to 10^8) |
| `records_S.py` | record values of S(B) up to X |
| `fit_constant.py`, `analyze_avg.py` | fit T(x) against C x (log x)^3 with C = G/12 |
| `verify_rewrite.py`, `verify_mainterm.py`, `verify_errors.py` | numerical checks of Section 5 of Paper B (change of variables, main term, error terms) |
| `verify_merged.py` | checks of the other new statements in Paper B |

### data/

All files are plain text, one row per line, columns separated by spaces. Lines starting with `#` are comments.

- `seed_count/partial_sums_1e8.txt`: `x T(x) U(x)` at checkpoints up to 10^8 (U is the unfloored sum)
- `seed_count/partial_sums_2e6.txt`: `x T(x)` up to 2·10^6
- `seed_count/records_S_1e8.txt`: bases B ≤ 10^8 where S(B) sets a new record
- `seed_count/largest_seed_check_1e4.txt`: output of `largest_seed_check.c` for B ≤ 10^4
- `zero_insertion/b_count_N.txt`, `b_smallest.txt`, `b_largest.txt`, `b_maxlength_L.txt`: b-files for N(B), the smallest and largest member, and the maximum length
- `zero_insertion/table_members_b2_100.txt`: all members of Z_B for 2 ≤ B ≤ 100
- `zero_insertion/records.txt`: record bases of N(B)
- `zero_insertion/growth_1e7.txt`, `n_ratio_stats_1e7.txt`: W(x) and N(B)/B statistics up to 10^7

## Running

Python 3 is needed, plus NumPy for `fit_constant.py` and `analyze_avg.py`. Everything else uses only the standard library.

```
python code/verify.py
python code/verify_merged.py
python code/verify_mainterm.py
python code/zero_insertion/zi.py verify 1000
python code/zero_insertion/zi.py members 10
```

The C programs can be built with any C compiler. We used zig:

```
python -m ziglang cc -O2 -o seed_sum_par code/seed_sum_par.c -lm
python code/run_par.py ./seed_sum_par 100000000 14 partial_sums_1e8.txt
```

The last argument to `run_par.py` before the output file is the number of processes.

## Citation

Please cite the archived release:

> K. Sharma and D. Sonkhla, *Code and data for: Inserting Zeros Detects Primes / Numbers that divide their zero-insertions*, Zenodo, 2026. https://doi.org/10.5281/zenodo.23281208

GitHub's "Cite this repository" button gives the same citation from `CITATION.cff`.

## License

MIT. See `LICENSE`.
