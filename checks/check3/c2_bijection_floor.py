"""check3 / c2: Lemma 5.1 (both directions), Lemma 5.2, Proposition 5.3.

Written independently of code/verify_rewrite.py.
"""
from math import gcd
import sys

BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 700


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


# forward map on every (B,k), B<=BMAX
fwd = {}
errs = []
for B in range(3, BMAX + 1):
    for k in range(1, B - 1):
        u, v = gcd(k + 1, B), gcd(k, B - 1)
        if gcd(u, v) != 1:
            errs.append(("a-coprime", B, k))
        if (B - k - 1) % (u * v):
            errs.append(("a-cong", B, k))
        t = (B - k - 1) // (u * v)
        if t < 1:
            errs.append(("a-t", B, k))
        up, vp = (k + 1) // u, k // v
        n = up * vp
        if gcd(t, n) != 1:
            errs.append(("c-gcd", B, k))
        # Lemma 5.2 exactly
        lhs = ((B - 1) * u * v) // (B * k)
        rhs = (k + 1) // n - (1 if vp == 1 else 0)
        if lhs != rhs:
            errs.append(("floor", B, k, lhs, rhs))
        # (u'v')*uv = k(k+1)
        if n * u * v != k * (k + 1):
            errs.append(("uv", B, k))
        # B<=x  <=>  (k+1)(tk+n) <= xn, tested at x=B and x=B-1
        if not ((k + 1) * (t * k + n) <= B * n):
            errs.append(("cond=", B, k))
        if (k + 1) * (t * k + n) <= (B - 1) * n:
            errs.append(("cond<", B, k))
        q = (k, up, vp, t)
        if q in fwd:
            errs.append(("not-injective", B, k))
        fwd[q] = (B, k)

# backward: every admissible quadruple with B<=BMAX is hit, and (b) holds for
# ALL quadruples (also gcd(t,u'v')>1)
hit = 0
for k in range(1, BMAX - 1):
    for up in divisors(k + 1):
        for vp in divisors(k):
            u, v = (k + 1) // up, k // vp
            t = 1
            while k + 1 + t * u * v <= BMAX:
                B = k + 1 + t * u * v
                if gcd(k + 1, B) != u * gcd(up, t) or gcd(k, B - 1) != v * gcd(vp, t):
                    errs.append(("b", k, up, vp, t))
                if gcd(t, up * vp) == 1:
                    hit += 1
                    if fwd.get((k, up, vp, t)) != (B, k):
                        errs.append(("inverse", k, up, vp, t))
                t += 1
print("pairs (B,k):", len(fwd), " admissible quadruples with B<=BMAX:", hit)
print("errors:", errs[:10], len(errs))
