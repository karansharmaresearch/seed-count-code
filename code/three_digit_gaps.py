"""Three-digit zero-insertion: which gap receives the zeros matters.

For digits 1 <= a, c, b <= B-1 put g = a*B^2 + c*B + b and t = c*B + b.
  - Lower gap: g never divides a*B^3 + c*B^2 + b.
  - Upper gap: g | a*B^(i+1) + t for all i >= 1  <=>  g | (B-1)*t.
  - U(B) := #upper-gap seeds is 0 exactly when B is a prime power, and for
    B = u*v (gcd(u,v)=1, 1<u<v) the number g = (B-1)*v^2 is an upper-gap seed.

Checks all of this (brute force over i up to N for small B; the criterion
and the construction for larger B) and prints the table used in the paper.
"""
from math import gcd


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def triples(B):
    for a in range(1, B):
        for c in range(1, B):
            for b in range(1, B):
                yield a, c, b


def brute(B, N=12):
    upper, lower = [], []
    for a, c, b in triples(B):
        g = a * B * B + c * B + b
        if all((a * B ** (i + 1) + c * B + b) % g == 0 for i in range(1, N)):
            upper.append((a, c, b))
        if (a * B ** 3 + c * B * B + b) % g == 0:
            lower.append((a, c, b))
    return upper, lower


def criterion(B):
    return [(a, c, b) for a, c, b in triples(B)
            if ((B - 1) * (c * B + b)) % (a * B * B + c * B + b) == 0]


def U(B):
    """Fast count via (B-1)t = k g, i.e. k a B^2 = m t with m = B-1-k."""
    P, cnt = (B - 1) * B * B, 0
    for m in range(1, B - 1):
        k = B - 1 - m
        step = m // gcd(m, P)
        for a in range(step, B, step):
            if a * k >= m:
                break
            t = k * a * B * B // m
            if t >= B and t % B:
                cnt += 1
    return cnt


def is_upper_seed(B, g):
    a, t = divmod(g, B * B)
    c, b = divmod(t, B)
    return 1 <= a <= B - 1 and c >= 1 and b >= 1 and ((B - 1) * t) % g == 0


if __name__ == "__main__":
    for B in range(2, 41):
        up, lo = brute(B)
        assert lo == [], B
        assert up == criterion(B) and len(up) == U(B), B
    print("B <= 40: lower gap empty; upper gap = criterion = fast count")

    assert brute(10)[0] == [(2, 2, 5), (6, 7, 5)]  # 225 and 675
    print("base 10 upper-gap seeds: 225, 675")

    LIM = 3000
    for B in range(2, LIM + 1):
        f = factor(B)
        assert (U(B) == 0) == (len(f) == 1), B
        if len(f) > 1:
            u = next(iter(p ** e for p, e in f.items()))
            v = B // u
            u, v = min(u, v), max(u, v)
            assert is_upper_seed(B, (B - 1) * v * v), B
    print(f"B <= {LIM}: U(B)=0 iff prime power; (B-1)v^2 construction valid")

    print("nonzero U(B), B<=40:", [(B, U(B)) for B in range(2, 41) if U(B)])
