"""Checks for the growth results on N(b) (Section 8 of the merged paper).

  python verify_growth.py            (from code/zero_insertion)

1. The q-decomposition  N(b) = sum_q C_b(q),  where C_b(q) = sum over g*r = q-1,
   g | b-1, gcd(r, n) = 1 (n = (b-1)/g, r < n) of floor((P_b(q)-1)/r) and
   P_b(q) = gcd(b-q, b^oo); and the bound C_b(q) <= (b-q) sigma(q-1)/(q-1)^2.
2. The summed form  T(x) = sum_{g,r,m>=1, g(m+r)<=x-1, gcd(m,r)=1}
   floor((gcd(m,(gr+1)^oo)-1)/r).
3. The two inequalities behind the bounds for T(x):
   lower  T(x) >= (x-1)/6 * (sum_{s <= sqrt(2(x-1)/3)} 1/rad(s) - 1),
   upper  T(x) <= x * sum_{u<x} tau(u)/u * F(u+1),
          F(v) = #{2 <= s <= x : rad(s) | v}.
"""
from math import gcd, isqrt
from fractions import Fraction
import zi


def smooth_part(m, b):
    s, g = 1, gcd(m, b)
    while g > 1:
        m //= g
        s *= g
        g = gcd(m, b)
    return s


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def contributions(b):
    """dict q -> C_b(q)."""
    out = {}
    for g in range(1, b):
        if (b - 1) % g:
            continue
        n = (b - 1) // g
        for r in range(1, n):
            if gcd(r, n) != 1:
                continue
            q = g * r + 1
            c = (smooth_part(b - q, b) - 1) // r
            if c:
                out[q] = out.get(q, 0) + c
    return out


def check_decomposition(B):
    for b in range(3, B + 1):
        C = contributions(b)
        assert sum(C.values()) == zi.count(b), b
        for q, c in C.items():
            # P_b(q) is the b-smooth part of b-q
            assert c <= Fraction((b - q) * sigma(q - 1), (q - 1) ** 2), (b, q, c)
            assert c <= (b - 4) / 2 if q == 3 else True, (b, c)
            assert 9 * c <= 4 * (b - q) if q >= 4 else True, (b, q, c)
    print(f"q-decomposition, C_b(q) <= (b-q)sigma(q-1)/(q-1)^2, C_b(3) <= (b-4)/2, "
          f"C_b(q) <= 4(b-q)/9 (q >= 4): ok for b <= {B}")
    eq = [b for b in range(4, B + 1) if 2 * contributions(b).get(3, 0) == b - 4]
    assert eq == [3 ** j + 3 for j in range(6) if 3 ** j + 3 <= B], eq
    print(f"C_b(3) = (b-4)/2 exactly for b = 3^j + 3 (j >= 0), b <= {B}: ok")


def T_direct(x):
    return sum(zi.count(b) for b in range(3, x + 1))


def T_summed(x):
    tot = 0
    for g in range(1, x):
        for r in range(1, x):
            if g * (1 + r) > x - 1:
                break
            q = g * r + 1
            for m in range(1, (x - 1) // g - r + 1):
                if gcd(m, r) == 1:
                    tot += (smooth_part(m, q) - 1) // r
    return tot


def rad(n):
    r, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        p += 1
    return r * (n if n > 1 else 1)


def bounds(x):
    T = T_direct(x)
    y = isqrt(2 * (x - 1) // 3)
    lower = Fraction(x - 1, 6) * (sum(Fraction(1, rad(s)) for s in range(1, y + 1)) - 1)
    rads = [rad(s) for s in range(x + 1)]
    tau = [0] * (x + 1)
    for d in range(1, x + 1):
        for k in range(d, x + 1, d):
            tau[k] += 1
    F = [0] * (x + 2)
    for s in range(2, x + 1):
        k = rads[s]
        for v in range(k, x + 1, k):
            F[v] += 1
    upper = x * sum(Fraction(tau[u], u) * F[u + 1] for u in range(1, x))
    assert lower <= T <= upper, (x, lower, T, upper)
    print(f"x={x}: lower {float(lower):.4g} <= T={T} <= upper {float(upper):.4g}")


if __name__ == "__main__":
    check_decomposition(300)
    for x in (50, 200, 600):
        assert T_summed(x) == T_direct(x), x
    print("summed form of T(x): ok for x = 50, 200, 600")
    for x in (100, 1000, 3000):
        bounds(x)
