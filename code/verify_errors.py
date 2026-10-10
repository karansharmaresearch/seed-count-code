"""Numerical checks for Sections 5.2, 5.5 and 5.6 of the merged paper (zero-insertion.tex)
(progression sums and the error terms E1, E2).

Notation: coprime (u',v'), n = u'v', rho = (k+1) mod n on the class K(u',v')
(so rho = 0 iff v' = 1), eps = [rho = 0], R = sqrt(x/(n t)),
Q = largest q >= 1 with (nq+rho)(t(nq+rho-1)+n) <= x n (Q = 0 if none),
a = (x/(2nt)) [t <= x/n].

  Lemma 5.4(a)  A_{u',v'}(t) = Q(Q+1)/2 - eps*Q.
  Lemma 5.4(b)  max(0, R-3) <= Q <= R+1.
  Lemma 5.4(c)  |A - a| <= 3.5 R + 2.5.
  Lemma 5.4(d)  A != 0  =>  2n-1 <= x;  t < x/(n-1) if n >= 2;  t <= x/2-1 if n = 1.
  Lemma 5.11    |H(y) - log y - gamma| <= 2/y for real y >= 1, H(y) = sum_{s<=y} 1/s.
  Lemma 5.12    |Delta(n)| <= gamma*S0 + S1 + 2 N_Z/Z + log Z * S_>.
  Prop 5.10     |E1(x)| <= (7 sqrt2 + 5) x frakD(x),  frakD(z) = sum_{m<=z} tau(m)/m.
  Prop 5.13     |E2(x)| <= x C1 frakD(x) + C2 x (1 + log x).

Lemma 5.4 is checked in two ways: (i) exhaustively for EVERY integer
3 <= x <= 90 (all coprime (u',v') with n <= x+2, all t <= 3x+5), with A
computed by brute force over the progression; (ii) for all (u',v',t) with
n <= 2x+2, t <= 2x+2 at x in {7, 20, 61, 150, 400, 1000}, with A computed by
brute force for x <= 150 and from (a) otherwise.
"""
import math
from math import gcd, log, sqrt

EG = 0.5772156649015329


def primes_of(n):
    ps, p = [], 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def rho_of(up, vp):
    """(k+1) mod n on the class k = -1 (mod u'), k = 0 (mod v')."""
    n = up * vp
    for k in range(1, n + 1):
        if (k + 1) % up == 0 and k % vp == 0:
            return (k + 1) % n


def A_direct(x, up, vp, t):
    """A_{u',v'}(t) by brute force over the progression k in K(u',v')."""
    n = up * vp
    k0 = next(k for k in range(1, n + 1) if (k + 1) % up == 0 and k % vp == 0)
    tot, k = 0, k0
    while (k + 1) * (t * k + n) <= x * n:
        tot += (k + 1) // n - (vp == 1)
        k += n
    return tot


def Q_of(x, n, rho, t):
    q = 0
    while (n * (q + 1) + rho) * (t * (n * (q + 1) + rho - 1) + n) <= x * n:
        q += 1
    return q


def H(y):
    return sum(1.0 / s for s in range(1, int(math.floor(y)) + 1))


def tau(m):
    return sum(2 - (d * d == m) for d in range(1, math.isqrt(m) + 1) if m % d == 0)


def check_54(x, up, vp, t, brute, rho):
    """Check Lemma 5.4 (a)-(d) at one (x, u', v', t); return |A-a|/(3.5R+2.5)."""
    n = up * vp
    eps = 1 if rho == 0 else 0
    Q = Q_of(x, n, rho, t)
    Af = Q * (Q + 1) // 2 - eps * Q
    if brute:
        assert A_direct(x, up, vp, t) == Af, (x, up, vp, t)                # (a)
    R = sqrt(x / (n * t))
    assert max(0.0, R - 3) - 1e-12 <= Q <= R + 1 + 1e-12, (x, n, t, Q, R)  # (b)
    a = x / (2 * n * t) if t * n <= x else 0.0
    assert abs(Af - a) <= 3.5 * R + 2.5, (x, n, t)                         # (c)
    if Af != 0:                                                            # (d)
        assert 2 * n - 1 <= x
        assert (t * (n - 1) < x) if n >= 2 else (2 * (t + 1) <= x), (x, n, t)
    return Af, abs(Af - a) / (3.5 * R + 2.5)


def coprime_pairs(nmax):
    for n in range(1, nmax + 1):
        for up in range(1, n + 1):
            if n % up == 0 and gcd(up, n // up) == 1:
                yield up, n // up


if __name__ == "__main__":
    # (i) Lemma 5.4, exhaustive for every x <= 90
    worst_i = 0.0
    for x in range(3, 91):
        for up, vp in coprime_pairs(x + 2):
            rho = rho_of(up, vp)
            assert (rho == 0) == (vp == 1)
            for t in range(1, 3 * x + 6):
                worst_i = max(worst_i, check_54(x, up, vp, t, True, rho)[1])
    print(f"Lemma 5.4 (a)-(d) hold for every x <= 90, all coprime (u',v') with n <= x+2, "
          f"all t <= 3x+5 (A by brute force); max |A-a|/(3.5R+2.5) = {worst_i:.3f}")

    # (ii) Lemma 5.4 at selected larger x; also E1 exactly, for Prop 5.10
    worst_ii = 0.0
    for x in (7, 20, 61, 150, 400, 1000):
        E1 = 0.0
        for up, vp in coprime_pairs(2 * x + 2):
            n = up * vp
            rho = rho_of(up, vp)
            assert (rho == 0) == (vp == 1)
            for t in range(1, 2 * x + 3):
                Af, r = check_54(x, up, vp, t, x <= 150, rho)
                worst_ii = max(worst_ii, r)
                if gcd(t, n) == 1:
                    E1 += Af - (x / (2 * n * t) if t * n <= x else 0.0)
        fD = sum(tau(m) / m for m in range(1, x + 1))
        assert abs(E1) <= (7 * sqrt(2) + 5) * x * fD                        # Prop 5.10
        print(f"  x={x:>5}: E1 = {E1:10.1f}, Prop 5.10 bound (7sqrt2+5) x frakD(x) = "
              f"{(7 * sqrt(2) + 5) * x * fD:10.1f}")
    print(f"Lemma 5.4 (a)-(d) hold for all (u',v',t) with n, t <= 2x+2 at x in "
          f"{{7,20,61,150,400,1000}} ((a) by brute force for x <= 150); "
          f"max |A-a|/(3.5R+2.5) = {worst_ii:.3f}")

    # Lemma 5.11 on a fine grid of real y
    ys = [1 + i / 997 for i in range(0, 200000)]
    mx = max(abs(H(y) - log(y) - EG) * y for y in ys[::7])
    assert mx <= 2
    print(f"Lemma 5.11: max_y y*|H(y) - log y - gamma| = {mx:.4f} <= 2 on a grid of step 7/997 in [1, 201]")

    # Lemma 5.12 for all n <= x, and Prop 5.13
    N = 20000
    C1 = sum(2 ** len(primes_of(d)) * (1 + log(d)) / d ** 2 for d in range(1, N + 1))
    C2 = sum(2 ** len(primes_of(d)) / d ** 2 for d in range(1, N + 1))
    for x in (100, 1000, 10000, 50000):
        E2, fD = 0.0, 0.0
        for n in range(1, x + 1):
            ps = primes_of(n)
            w = len(ps)
            ph = n
            for p in ps:
                ph = ph // p * (p - 1)
            Z = x / n
            ds = [(1, 1)]
            for p in ps:
                ds += [(d * p, -s) for d, s in ds]
            hs = sum(s / d * H(Z / d) for d, s in ds if d <= Z)
            Delta = hs - ph / n * log(Z)
            b = (EG * sum(1 / d for d, _ in ds) + sum(log(d) / d for d, _ in ds)
                 + 2 * sum(1 for d, _ in ds if d <= Z) / Z
                 + log(Z) * sum(1 / d for d, _ in ds if d > Z))
            assert abs(Delta) <= b + 1e-12, (x, n, Delta, b)                  # Lemma 5.12
            E2 += 2**w * x / (2 * n) * Delta
            fD += tau(n) / n
        bnd = x * C1 * fD + C2 * x * (1 + log(x))
        assert abs(E2) <= bnd                                                  # Prop 5.13
        print(f"  x={x:>6}: Lemma 5.12 holds for all n <= x; E2 = {E2:.1f}, "
              f"Prop 5.13 bound = {bnd:.1f} (C1, C2 truncated at d <= {N}: underestimates)")
