"""Numerical checks for the change of variables in the merged paper (zero-insertion.tex),
Section 5.1 (Lemma 5.1, Lemma 5.2, Proposition 5.3).

Notation: for 1 <= k <= B-2, u = gcd(k+1,B), v = gcd(k,B-1), u' = (k+1)/u,
v' = k/v, n = u'v'.

  Lemma 5.1(a)  gcd(u,v) = 1 and B = k+1 (mod uv).
  Lemma 5.1(c)  (B,k) <-> (k,u',v',t) with B = k+1+t*uv, t >= 1, gcd(t,n) = 1,
                u' | k+1, v' | k, is a bijection (checked both directions,
                which also checks the gcd formulas of Lemma 5.1(b)).
  Lemma 5.2     floor((B-1)/L_k) = floor((k+1)/n) - [v' = 1].
  Prop 5.3      T(x) = sum over the quadruples of floor((k+1)/n) - [v'=1],
                evaluated here with the t-sum done first:
                sum_{k<=x-2} sum_{u'|k+1} sum_{v'|k} w * #{t <= Y : gcd(t,n)=1},
                w = floor((k+1)/n) - [v'=1], Y = floor((x-k-1) n / (k(k+1)))
                (B <= x  <=>  t*uv <= x-k-1), with the count of t done both
                directly and by Mobius inversion, sum_{d|n} mu(d) floor(Y/d).
Every T(x) is compared with a direct evaluation of formula (1.1); no
hard-coded values are used.
"""
from math import gcd


def divisors(n):
    ds, i = [], 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n:
                ds.append(n // i)
        i += 1
    return ds


def mobius(n):
    res, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            res = -res
        p += 1
    return -res if n > 1 else res


def T_direct(x):
    return sum(((B - 1) * gcd(k + 1, B) * gcd(k, B - 1)) // (B * k)
               for B in range(3, x + 1) for k in range(1, B - 1))


def T_rewrite(x, use_mobius):
    tot = 0
    for k in range(1, x - 1):
        for up in divisors(k + 1):
            for vp in divisors(k):
                n = up * vp
                w = (k + 1) // n - (vp == 1)
                if w <= 0:
                    continue
                Y = (x - k - 1) * n // (k * (k + 1))      # t*uv <= x-k-1, uv = k(k+1)/n
                if use_mobius:
                    cnt = sum(mobius(d) * (Y // d) for d in divisors(n))
                else:
                    cnt = sum(1 for t in range(1, Y + 1) if gcd(t, n) == 1)
                tot += w * cnt
    return tot


if __name__ == "__main__":
    X = 400
    fwd = set()
    for B in range(3, X + 1):
        for k in range(1, B - 1):
            u, v = gcd(k + 1, B), gcd(k, B - 1)
            assert gcd(u, v) == 1 and (B - k - 1) % (u * v) == 0          # Lemma 5.1(a)
            up, vp, t = (k + 1) // u, k // v, (B - k - 1) // (u * v)
            assert t >= 1 and gcd(t, up * vp) == 1                         # Lemma 5.1(c) (->)
            n = up * vp
            assert ((B - 1) * u * v) // (B * k) == (k + 1) // n - (vp == 1)  # Lemma 5.2
            fwd.add((B, k))
    bwd = set()
    for k in range(1, X - 1):
        for up in divisors(k + 1):
            for vp in divisors(k):
                u, v, n = (k + 1) // up, k // vp, up * vp
                t = 1
                while k + 1 + t * u * v <= X:
                    if gcd(t, n) == 1:
                        B = k + 1 + t * u * v
                        assert gcd(k + 1, B) == u and gcd(k, B - 1) == v    # Lemma 5.1(b),(c) (<-)
                        assert (B, k) not in bwd
                        bwd.add((B, k))
                    t += 1
    assert fwd == bwd
    print(f"Lemmas 5.1-5.2 and the bijection hold for all 3 <= B <= {X} ({len(fwd)} pairs)")

    xs = (3, 4, 5, 10, 30, 97, 100, 257, 1000, 2000, 3000)
    for x in xs:
        d = T_direct(x)
        r2 = T_rewrite(x, True)
        r1 = T_rewrite(x, False) if x <= 1000 else r2
        assert d == r1 == r2, (x, d, r1, r2)
    print("Proposition 5.3 (t-count direct for x <= 1000, and by Mobius) equals the "
          f"direct evaluation of (1.1) for x in {xs}")
