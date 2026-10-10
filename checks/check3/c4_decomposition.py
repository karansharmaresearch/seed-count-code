"""check3 / c4: Propositions 5.3, 5.5, 5.10, 5.13 and Lemmas 5.11, 5.12.

For several integer x: T(x) computed three ways
  (1) closed form (1.1) directly,
  (2) the grouped sum of Prop 5.3 with A from Lemma 5.4(a) (Q by exact
      integer search),
  (3) M(x)+E1(x)+E2(x) with each piece computed from its own definition,
and the explicit bounds |E1| <= (7 sqrt2+5) x D(x),
|E2| <= x C1 D(x) + C2 x (1+log x), and each Sigma_i bound separately.
"""
from math import gcd, isqrt, log, sqrt
import sys

def closed_T(x):
    T = 0
    for B in range(3, x + 1):
        for k in range(1, B - 1):
            T += ((B - 1) * gcd(k + 1, B) * gcd(k, B - 1)) // (B * k)
    return T

def sieve_spf(N):
    spf = list(range(N + 1))
    for p in range(2, isqrt(N) + 1):
        if spf[p] == p:
            for m in range(p * p, N + 1, p):
                if spf[m] == m:
                    spf[m] = p
    return spf

N = 20000
spf = sieve_spf(N)
def factor(m):
    f = {}
    while m > 1:
        p = spf[m]; f[p] = f.get(p, 0) + 1; m //= p
    return f
def omega(m): return len(factor(m))
def phi(m):
    r = m
    for p in factor(m): r = r // p * (p - 1)
    return r
def tau(m):
    r = 1
    for e in factor(m).values(): r *= e + 1
    return r
def mu(m):
    f = factor(m)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)
def divisors(m):
    ds = [1]
    for p, e in factor(m).items():
        ds = [d * p ** i for d in ds for i in range(e + 1)]
    return ds
def H(y):
    return sum(1.0 / s for s in range(1, int(y) + 1))
def D(z):
    return sum(tau(m) / m for m in range(1, int(z) + 1))

def Q_of(x, n, t, rho):
    if (n + rho) * (t * (n + rho - 1) + n) > x * n:
        return 0
    lo, hi = 1, 2 * isqrt(x) + 10
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if (n * mid + rho) * (t * (n * mid + rho - 1) + n) <= x * n:
            lo = mid
        else:
            hi = mid - 1
    return lo

def rho_of(up, vp):
    n = up * vp
    if vp == 1:
        return 0
    # k = 0 mod v', k = -1 mod u'
    if up == 1:
        return 1 % n
    k = vp * ((-pow(vp, -1, up)) % up)
    return (k + 1) % n

C1 = sum(2 ** omega(d) * (1 + log(d)) / d ** 2 for d in range(1, N + 1))
C2 = sum(2 ** omega(d) / d ** 2 for d in range(1, N + 1))
# tail of C1 beyond N is < sum tau(d)(1+log d)/d^2 ~ (log N)^2/N ; negligible here
print("C1 ~ %.4f  C2 ~ %.4f (truncated at %d)" % (C1, C2, N))

xs = [int(a) for a in sys.argv[1:]] or [3, 4, 5, 6, 7, 10, 30, 97, 100, 257, 500, 1000]
for x in xs:
    Tc = closed_T(x) if x <= 1200 else None
    Tg = 0; E1 = 0.0; E1_abs_terms = 0.0
    E1_by_n_ok = True
    for n in range(1, x + 1):
        for up in divisors(n):
            vp = n // up
            if gcd(up, vp) != 1: continue
            rho = rho_of(up, vp); eps = 1 if vp == 1 else 0
            tmax = 3 * x + 5
            for t in range(1, tmax + 1):
                if gcd(t, n) != 1: continue
                Q = Q_of(x, n, t, rho)
                A = Q * (Q + 1) // 2 - eps * Q
                Tg += A
                main = x / (2 * n * t) if t * n <= x else 0.0
                E1 += A - main
                if (A != 0 or main != 0) and t > (x if n == 1 else x / (n - 1)):
                    E1_by_n_ok = False
    M = 0.0; E2 = 0.0; S1 = S2 = S3 = 0.0
    for n in range(1, x + 1):
        w = 2 ** omega(n)
        Z = x / n
        M += x / 2 * (w * phi(n) / n / n) * log(Z)
        Delta = sum(1.0 / t for t in range(1, int(Z) + 1) if gcd(t, n) == 1) - phi(n) / n * log(Z)
        E2 += w * x / (2 * n) * Delta
        sq = [d for d in divisors(n) if mu(d) != 0]
        b1 = 0.5772156649 * sum(1 / d for d in sq) + sum(log(d) / d for d in sq)
        b2 = 2 / Z * sum(1 for d in sq if d <= Z)
        b3 = log(Z) * sum(1 / d for d in sq if d > Z)
        if abs(Delta) > b1 + b2 + b3 + 1e-9:
            print("  Lemma 5.12 FAILS at", x, n)
        S1 += w * x / (2 * n) * b1; S2 += w * x / (2 * n) * b2; S3 += w * x / (2 * n) * b3
    Dx = D(x)
    bE1 = (7 * sqrt(2) + 5) * x * Dx
    bE2 = x * C1 * Dx + C2 * x * (1 + log(x))
    print(f"x={x}: T closed={Tc} grouped={Tg} M+E1+E2={M+E1+E2:.6f} | "
          f"E1={E1:.1f} (bound {bE1:.0f}) E2={E2:.1f} (bound {bE2:.0f}) | "
          f"S1={S1:.0f}<= {x/2*C1*Dx:.0f}, S2={S2:.0f}<= {C2*x*(1+log(x)):.0f}, "
          f"S3={S3:.0f}<= {x/2*C1*Dx:.0f} | support T_n ok={E1_by_n_ok}")
    assert Tc is None or Tc == Tg
    assert abs(Tg - (M + E1 + E2)) < 1e-6 * max(1, Tg)
    assert abs(E1) <= bE1 and abs(E2) <= bE2

# Lemma 5.11 on a grid including points just below integers
worst = 0.0
for m in range(1, 3000):
    for y in (m, m + 1e-9, m + 0.5, m + 1 - 1e-9):
        e = abs(H(y) - log(y) - 0.5772156649015329) * y / 2
        worst = max(worst, e)
print("Lemma 5.11: max |H(y)-log y-gamma| / (2/y) on grid y<3000:", round(worst, 6))
