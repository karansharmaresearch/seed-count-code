"""Numerical checks for Sections 5.1, 5.3 and 5.4 of the merged paper (zero-insertion.tex)
(the decomposition and the main term).

g(n) = f(n)/n with f(n) = 2^omega(n) phi(n)/n;  G = prod_p (1-3/p^2+2/p^3).

  Lemma 5.6  (a) f(n) = sum over ordered coprime u'v' = n of phi(u')phi(v')/n
                 (the 2^omega(n) count used in Prop 5.5);
             (b) h = f*mu*mu has h(p) = -2/p, h(p^2) = -(p-2)/p, h(p^a) = 0 (a >= 3).
  Prop 5.8   A(y) = sum_{n<=y} g(n) = (G/2) log^2 y + O(1 + log y)
  Cor  5.9   M(x) := (x/2) sum_{n<=x} g(n) log(x/n) = (G/12) x log^3 x + O(x log^2 x)
  Prop 5.3   T(x) = sum_{(u',v') coprime} sum_{t: (t,n)=1} A_{u',v'}(t), with
             A_{u',v'}(t) = sum of w(k,u',v') over k >= 1, u'|k+1, v'|k,
                            t*k(k+1) <= (x-k-1) n
  Prop 5.5   T(x) = M(x) + E1(x) + E2(x) exactly, where
             E1 = sum_{(u',v')} sum_{(t,n)=1} (A_{u',v'}(t) - [t <= x/n] x/(2nt)),
             E2 = sum_{n<=x} 2^omega(n) (x/(2n)) (sum_{t<=x/n,(t,n)=1} 1/t - (phi(n)/n) log(x/n)).

Usage: python code/verify_mainterm.py [path/to/partial_sums_1e8.txt]  (default: data/seed_count/)
(the default data file is found next to this script, so it runs from any
directory)
"""
import os
import sys
from math import gcd, log
import numpy as np

G = 0.2867474335
EULER_GAMMA = 0.5772156649015329


def spf_table(N):
    spf = np.zeros(N + 1, dtype=np.int64)
    for p in range(2, int(N**0.5) + 1):
        if spf[p] == 0:
            block = spf[p * p::p]
            block[block == 0] = p
            spf[p * p::p] = block
    idx = np.nonzero(spf == 0)[0]
    spf[idx] = idx
    return spf


def primes_of(n, spf):
    ps = []
    while n > 1:
        p = int(spf[n])
        ps.append(p)
        while n % p == 0:
            n //= p
    return ps


def harmonic(m):
    if m <= 60:
        return sum(1.0 / j for j in range(1, m + 1))
    return log(m) + EULER_GAMMA + 1 / (2 * m) - 1 / (12 * m * m) + 1 / (120 * m**4)


def coprime_harmonic(Z, ps):
    """sum_{t <= Z, gcd(t, n) = 1} 1/t, where ps are the primes of n (Mobius)."""
    ds = [(1, 1)]
    for p in ps:
        ds += [(d * p, -s) for d, s in ds if d * p <= Z]
    return sum(s / d * harmonic(Z // d) for d, s in ds)


def E2(x, spf):
    tot = 0.0
    for n in range(1, x + 1):
        ps = primes_of(n, spf)
        ph = n
        for p in ps:
            ph = ph // p * (p - 1)
        tot += 2 ** len(ps) * x / (2 * n) * (coprime_harmonic(x // n, ps) - ph / n * log(x / n))
    return tot


def T_M_E1_E2_small(x, spf):
    """T(x) directly, and sum of A, M, E1, E2 straight from the definitions."""
    T = sum(((B - 1) * gcd(k + 1, B) * gcd(k, B - 1)) // (B * k)
            for B in range(3, x + 1) for k in range(1, B - 1))
    E1, sum_A = 0.0, 0
    for up in range(1, x + 1):
        for vp in range(1, x + 1):
            if gcd(up, vp) != 1:
                continue
            n = up * vp
            ks = [k for k in range(1, x - 1) if (k + 1) % up == 0 and k % vp == 0]
            t = 1
            while True:
                A = sum((k + 1) // n - (vp == 1) for k in ks
                        if t * k * (k + 1) <= (x - k - 1) * n)
                a = x / (2 * n * t) if t * n <= x else 0.0
                if A == 0 and a == 0.0:
                    break
                if gcd(t, n) == 1:
                    E1 += A - a
                    sum_A += A
                t += 1
    M = sum(2 ** len(primes_of(n, spf)) * _phi(n, spf) / n**2 * log(x / n) for n in range(1, x + 1)) * x / 2
    return T, sum_A, M, E1, E2(x, spf)


def _phi(n, spf):
    r = n
    for p in primes_of(n, spf):
        r = r // p * (p - 1)
    return r


def sieve_g(N):
    """g(n) = 2^omega(n) phi(n) / n^2 for 0 <= n <= N (g(0) = 0)."""
    ratio = np.ones(N + 1)            # phi(n)/n
    two = np.ones(N + 1)              # 2^omega(n)
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    for p in np.nonzero(is_p)[0]:
        ratio[p::p] *= 1 - 1 / p
        two[p::p] *= 2
    n = np.arange(N + 1, dtype=np.float64)
    g = np.zeros(N + 1)
    g[1:] = two[1:] * ratio[1:] / n[1:]
    return g


def check_lemma_5_6(spf):
    """Lemma 5.6 (a) for n < 3000 and (b) for n < 1500, in exact arithmetic."""
    from fractions import Fraction as Fr

    def fac(n):
        out = {}
        while n > 1:
            p = int(spf[n])
            out[p] = out.get(p, 0) + 1
            n //= p
        return out

    def mob(n):
        F = fac(n)
        return 0 if any(e > 1 for e in F.values()) else (-1) ** len(F)

    for n in range(1, 3000):
        s = sum(_phi(a, spf) * _phi(n // a, spf) for a in range(1, n + 1)
                if n % a == 0 and gcd(a, n // a) == 1)
        assert s == 2 ** len(fac(n)) * _phi(n, spf), n                     # (a)
    f = [Fr(0)] + [Fr(2 ** len(fac(n)) * _phi(n, spf), n) for n in range(1, 1500)]
    mu2 = [0] + [sum(mob(d) * mob(n // d) for d in range(1, n + 1) if n % d == 0)
                 for n in range(1, 1500)]
    for n in range(1, 1500):
        h = sum(f[d] * mu2[n // d] for d in range(1, n + 1) if n % d == 0)
        pred = Fr(1)
        for p, e in fac(n).items():
            pred *= Fr(-2, p) if e == 1 else (Fr(-(p - 2), p) if e == 2 else 0)
        assert h == pred, n                                                # (b)


if __name__ == "__main__":
    spf = spf_table(10**6)
    check_lemma_5_6(spf)
    print("Lemma 5.6 (a) holds for n < 3000 and (b) for n < 1500 (exact arithmetic)")
    for x in (10, 30, 60, 120):
        T, sA, M, e1, e2 = T_M_E1_E2_small(x, spf)
        assert T == sA, (x, T, sA)                                         # Prop 5.3
        assert abs(T - (M + e1 + e2)) < 1e-7 * max(1, T), (x, T, M + e1 + e2)  # Prop 5.5
    print("Prop 5.3 (reordering) and Prop 5.5 (T = M + E1 + E2) hold exactly for x in {10,30,60,120}")

    fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "data", "seed_count", "partial_sums_1e8.txt")
    data = {int(r[0]): r[1] for r in np.loadtxt(fn)}
    N = 10**7
    g = sieve_g(N)
    A = np.cumsum(g)
    Glog = np.cumsum(g * np.log(np.maximum(np.arange(N + 1), 1)))
    print(f"\n{'x':>9} {'(A-G/2 L^2)/L':>14} {'(M-G/12 xL^3)/xL^2':>19} "
          f"{'(T-M)/xL^2':>11} {'E2/xL^2':>9} {'E1/xL^2':>9}")
    for x in sorted(c for c in data if c <= N):
        L = log(x)
        M = x / 2 * (L * A[x] - Glog[x])
        T = data[x]
        e2 = E2(x, spf) if x <= 10**6 else float("nan")
        print(f"{x:>9} {(A[x] - G / 2 * L * L) / L:14.5f} {(M - G / 12 * x * L**3) / (x * L * L):19.5f} "
              f"{(T - M) / (x * L * L):11.5f} {e2 / (x * L * L):9.5f} {(T - M - e2) / (x * L * L):9.5f}")
