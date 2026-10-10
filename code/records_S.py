"""Record bases of S(B) = A394974(B) up to X, in parallel.

Needs seed_sum_par.c compiled with -DRECORDS (it then prints the chunk-local
records "R B S(B)"). Every global record is a chunk-local record of its chunk,
so merging the chunks in order with a running maximum gives the exact list.

Usage: python records_S.py EXE X [NPROC] [OUTFILE]
Build: python -m ziglang cc -O2 -DRECORDS -o seed_sum_rec seed_sum_par.c -lm
"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor


def run_chunk(exe, b0, b1):
    out = subprocess.run([exe, str(b0), str(b1)], capture_output=True, text=True, check=True).stdout
    return b0, [(int(f[1]), int(f[2])) for f in (l.split() for l in out.split("\n")) if f and f[0] == "R"]


def main():
    exe, X = sys.argv[1], int(float(sys.argv[2]))
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    outfile = sys.argv[4] if len(sys.argv) > 4 else "records_S.txt"
    nchunks = max(1, min(8 * nproc, X // 20000))
    edges = [2 + (X - 1) * i // nchunks for i in range(nchunks + 1)]
    ranges = [(edges[i], edges[i + 1] - 1) for i in range(nchunks)]
    ranges[-1] = (ranges[-1][0], X)
    with ThreadPoolExecutor(nproc) as ex:
        res = sorted(ex.map(lambda r: run_chunk(exe, *r), ranges))
    best, rec = 0, []
    for _, rows in res:
        for B, S in rows:
            if S > best:
                best = S
                rec.append((B, S))
    with open(outfile, "w") as fh:
        fh.write(f"# Record bases of S(B) = A394974(B) for B <= {X}: B S(B)\n")
        for B, S in rec:
            fh.write(f"{B} {S}\n")
    print(f"{len(rec)} records up to {X} -> {outfile}")


if __name__ == "__main__":
    main()
