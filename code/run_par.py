"""Run seed_sum_par over [2, X] in parallel chunks and stitch the checkpoints.

Usage: python run_par.py EXE X [NPROC] [OUTFILE]
Writes lines "x T(x) U(x)" (T = sum S(B), U = unfloored envelope), one per
checkpoint (x = 1000*2^j, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, and X).
Columns 1-2 have the same format as partial_sums_2e6.txt.
"""
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor


def run_chunk(exe, b0, b1):
    out = subprocess.run([exe, str(b0), str(b1)], capture_output=True, text=True, check=True).stdout
    cps, end, best = [], None, None
    for line in out.split("\n"):
        f = line.split()
        if not f:
            continue
        if f[0] == "C":
            cps.append((int(f[1]), int(f[2]), float(f[3])))
        elif f[0] == "E":
            end = (int(f[1]), int(f[2]), float(f[3]))
        elif f[0] == "M":
            best = (float(f[1]), int(f[2]))
    return b0, b1, cps, end, best


def main():
    exe, X = sys.argv[1], int(float(sys.argv[2]))
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    outfile = sys.argv[4] if len(sys.argv) > 4 else f"partial_sums_{X:.0e}.txt".replace("+0", "")
    nchunks = max(1, min(8 * nproc, X // 20000))
    edges = [2 + (X - 1) * i // nchunks for i in range(nchunks + 1)]
    ranges = [(edges[i], edges[i + 1] - 1) for i in range(nchunks)]
    ranges[-1] = (ranges[-1][0], X)
    t0 = time.time()
    with ThreadPoolExecutor(nproc) as ex:
        res = sorted(ex.map(lambda r: run_chunk(exe, *r), ranges))
    rows, S0, U0, best = [], 0, 0.0, (0.0, 0)
    for b0, b1, cps, end, bst in res:
        for c, s, u in cps:
            rows.append((c, S0 + s, U0 + u))
        # a checkpoint that coincides with a chunk end is reported via E only
        S0 += end[1]
        U0 += end[2]
        rows.append((end[0], S0, U0))
        best = max(best, bst)
    cp = set()
    for c, _, _ in rows:
        q = c // 1000
        if c == X or c in (10**6, 2 * 10**6, 5 * 10**6, 10**7, 2 * 10**7, 5 * 10**7, 10**8, 2 * 10**8) \
                or (c % 1000 == 0 and q & (q - 1) == 0):
            cp.add(c)
    with open(outfile, "w") as fh:
        for c, s, u in sorted(set(rows)):
            if c in cp:
                fh.write(f"{c} {s} {u:.6f}\n")
    print(f"X={X}: {nchunks} chunks on {nproc} processes, {time.time() - t0:.1f} s -> {outfile}")
    print(f"max S(B)/B over B != 36, B <= X: {best[0]:.6f} at B = {best[1]}")


if __name__ == "__main__":
    main()
