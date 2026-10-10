# code/zero_insertion — zero-insertion divisors Z_B and N(B) = |Z_B|

Code behind Sections 6-9 of "Numbers that divide their zero-insertions:
classification, maximal order and average order". Data files are in `data/zero_insertion/`.

| file | purpose |
|------|---------|
| `zi.py` | classification (Theorem 6.3), count formula (Corollary 7.1), brute-force verifier, b-file writer |
| `verify_growth.py` | checks behind Section 8 (q-decomposition, Corollary 7.3 bounds, summed form of W(x), Theorem 8.1 bounds) |
| `n_fast.c` | N(B) for all B <= X via Proposition 7.2; output for X = 1e7 in `data/zero_insertion/growth_1e7.txt` |
| `n_ratio_stats.c` | N(B)/B statistics (Remark 8.8, Section 9); output in `data/zero_insertion/n_ratio_stats_1e7.txt` |

Usage (from this folder):
- `python zi.py verify 1000` — classification vs brute force, an (a,k,ell) search using only Lemma 3.5, and the definition of A394974
- `python zi.py members 10` — lists Z_10 (= A398346)
- `python zi.py bfiles 20000` — regenerates the b-files in `data/zero_insertion/`
- `python verify_growth.py`
- `python -m ziglang cc -O2 -o n_fast n_fast.c` then `n_fast 10000000 out.txt` (a few minutes)

Note: `zi.py` uses older letters (base `b`,
member `k`, multiplier `c`); the paper uses `B`, `z`, `k`.
The merged-paper checks of the new statements are in `code/verify_merged.py`.
