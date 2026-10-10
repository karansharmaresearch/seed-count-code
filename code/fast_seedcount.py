"""Fast S(B) and cumulative T(x) = sum_{B<=x} S(B), for settling the average order.

Per-B method: a summand floor((B-1)/L_k) of Theorem 2 is nonzero only when
g1*g2 > k, where g1 = gcd(k+1,B) | B and g2 = gcd(k,B-1) | (B-1). Since
gcd(B,B-1)=1 any such (g1,g2) are coprime, and the congruences k = -1 (mod g1),
k = 0 (mod g2) pin a unique residue k0 in [1, g1*g2); only that k0 can give a
nonzero term for the pair (g1,g2). Iterating over divisor pairs therefore costs
O(d(B) d(B-1)) per B instead of O(B).
"""
from math import gcd


def divisors_from_spf(n, spf):
    """All divisors of n using a smallest-prime-factor table."""
    if n <= 1:
        return [1]
    factors = {}
    m = n
    while m > 1:
        p = spf[m]
        c = 0
        while m % p == 0:
            m //= p
            c += 1
        factors[p] = c
    divs = [1]
    for p, c in factors.items():
        cur = []
        pk = 1
        for _ in range(c + 1):
            for d in divs:
                cur.append(d * pk)
            pk *= p
        divs = cur
    return divs


def S_fast(B, dB, dBm1):
    """S(B) given divisor lists of B and B-1."""
    total = 0
    for g1 in dB:
        if g1 == 1:
            continue
        for g2 in dBm1:
            inv = pow(g2 % g1, -1, g1)          # g2 coprime to g1 (g1>=2)
            t0 = (-inv) % g1
            k0 = g2 * t0                          # unique residue in [1, g1*g2)
            if k0 < 1 or k0 > B - 2:
                continue
            if gcd(k0 + 1, B) != g1 or gcd(k0, B - 1) != g2:
                continue
            L = B * k0 // (g1 * g2)
            if L <= B - 1:
                total += (B - 1) // L
    return total


def build_spf(X):
    spf = list(range(X + 1))
    i = 2
    while i * i <= X:
        if spf[i] == i:
            for j in range(i * i, X + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


if __name__ == "__main__":
    import sys
    from seed_count import S, S_brute

    X = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    spf = build_spf(X)

    # --- validate against brute force (small) and the O(B) formula (to X) ---
    bad = []
    for B in range(2, min(X, 400) + 1):
        if S_fast(B, divisors_from_spf(B, spf), divisors_from_spf(B - 1, spf)) != S_brute(B):
            bad.append(B)
    assert not bad, f"S_fast != brute at {bad}"
    bad = []
    for B in range(2, X + 1):
        if S_fast(B, divisors_from_spf(B, spf), divisors_from_spf(B - 1, spf)) != S(B):
            bad.append(B)
    assert not bad, f"S_fast != formula at {bad}"
    print(f"OK  S_fast validated against brute (<=400) and formula (<= {X})")
