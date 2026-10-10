"""Combine the four c5_seg range outputs in seg1e8/ (run as
./c5_seg 3 26000000, 26000001 51000000, 51000001 76000000, 76000001 100000000)
into T(x) at every multiple of 1000 up to 1e8, and compare with all 24
checkpoints of data/seed_count/partial_sums_1e8.txt."""
import glob, os
here = os.path.dirname(os.path.abspath(__file__))
segs = sorted(glob.glob(os.path.join(here, 'seg1e8', 'seg_*.out')),
              key=lambda f: int(os.path.basename(f).split('_')[1]))
base, P, extra = 0, {}, []
for f in segs:
    for line in open(f):
        if line.startswith('P '):
            _, B, v = line.split(); P[int(B)] = base + int(v)
        elif line.startswith('END'):
            tot = int(line.split()[1].split('=')[1]); extra.append(line.strip())
        elif line[0] in 'RV':
            extra.append(line.strip())
    base += tot
ref = {int(a): int(b) for a, b, _ in
       (l.split() for l in open(os.path.join(here, '..', '..', 'data', 'seed_count', 'partial_sums_1e8.txt')))}
print("T(1e8) =", base)
print("checkpoints compared:", len(ref), " mismatches:", [x for x in ref if P.get(x) != ref[x]])
print("\n".join(extra))
