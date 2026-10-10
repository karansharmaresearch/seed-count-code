"""
Zero-insertion divisors in base b  —  computation and verification.

Z_b = numbers k with >= 2 base-b digits, second and last digit nonzero,
      such that k divides every number obtained by inserting any number
      of 0's between the first and second digits of k.
(Base 10: OEIS A398346.)

Classification (Theorem "Classification" of the paper): Z_b is in bijection with the pairs
    (a, c),  a, c >= 1,  c(a+1) <= b-2,  rad(s) | b,
where e = b-1-c, s = e / gcd(e, a c).  The member is
    k(a, c) = a (b-1) b^mu / e,   mu = mu_b(s) = max_{p|s} ceil(v_p(s)/v_p(b)),
with leading digit a and exactly mu+1 digits.

Usage:
    python3 zi.py verify [B]     brute-force check of the theorem for 2 <= b <= B
    python3 zi.py members b      list Z_b with base-b digits
    python3 zi.py bfiles N       write b-files for b = 2..N into data/zero_insertion/
"""
import sys, os
from math import gcd


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def is_prime(n):
    return n >= 2 and factor(n) == {n: 1}


def coprime_part(n, b):
    """Largest divisor of n coprime to b."""
    g = gcd(n, b)
    while g > 1:
        n //= g
        g = gcd(n, b)
    return n


def mu(s, fb):
    return max(-(-v // fb[p]) for p, v in factor(s).items())


def digits(k, b):
    d = []
    while k:
        d.append(k % b)
        k //= b
    return d[::-1]


def members(b):
    """Z_b via the classification; returns sorted list of (k, a, c, mu)."""
    fb = factor(b)
    out = []
    for c in range(1, b - 2):
        e = b - 1 - c
        for a in range(1, (b - 2) // c):          # c(a+1) <= b-2
            s = e // gcd(e, c * a)
            if coprime_part(s, b) == 1:              # rad(s) | b  (s > 1 automatically)
                m = mu(s, fb)
                num = a * (b - 1) * b ** m
                assert num % e == 0
                out.append((num // e, a, c, m))
    return sorted(out)


def count(b):
    """N(b) = sum_{c=1}^{b-3} floor(A_c / rho_b(c)) (Corollary "Exact count")."""
    tot = 0
    for c in range(1, b - 2):
        A = (b - 2) // c - 1
        if A > 0:
            R = coprime_part(b - 1 - c, b)
            tot += A // (R // gcd(R, c))
    return tot


def in_Zb_bruteforce(k, b, insertions=8):
    d = digits(k, b)
    if len(d) < 2 or d[1] == 0 or d[-1] == 0:
        return False
    a, m = d[0], len(d) - 1
    t = k - a * b ** m
    return all((a * b ** (m + j) + t) % k == 0 for j in range(insertions))


def verify(B):
    """Independent check: brute force over every k with up to L+1 digits,
    L = floor(log2(b-2)) (digit bound of the paper), for small b; formula checks beyond."""
    for b in range(2, B + 1):
        Z = members(b)
        ks = [x[0] for x in Z]
        assert len(set(ks)) == len(ks)
        assert all(in_Zb_bruteforce(k, b) for k in ks), b
        assert digits(ks[0], b) if ks else True
        assert all(len(digits(k, b)) == m + 1 and digits(k, b)[0] == a for k, a, c, m in Z)
        assert count(b) == len(Z), b
        assert (len(Z) == 0) == is_prime(b), b
        for d in range(2, b):
            if b % d == 0:
                assert d * (b - 1) in ks, (b, d)
        # exhaustive search up to b^(L+1): only for small bases (cost b^(L+1))
        L = (b - 2).bit_length() - 1 if b > 2 else 0
        if b ** (L + 1) <= 3 * 10 ** 6:
            brute = [k for k in range(b, b ** (L + 1)) if in_Zb_bruteforce(k, b)]
            assert brute == ks, (b, brute, ks)
        # search over (a, c, m), m <= 60, using only "one condition suffices"
        if b <= 300:
            assert search_acm(b) == ks, b
        # two-digit members vs the definition of A394974
        if b <= 120:
            assert sum(1 for z in Z if z[3] == 1) == a394974(b), b
    print(f"verified 2 <= b <= {B}")


def search_acm(b, M=60):
    """Z_b restricted to <= M+1 digits, found from (b-1)t = c k alone
    (no use of the range or rigidity lemmas)."""
    out = set()
    for c in range(1, b - 1):
        e = b - 1 - c
        for a in range(1, b):
            for m in range(1, M + 1):
                if (c * a * pow(b, m, e)) % e:
                    continue
                t = c * a * b ** m // e
                if b ** (m - 1) <= t < b ** m and t % b:
                    out.add(a * b ** m + t)
    return sorted(out)


def a394974(n):
    """A394974(n) from its definition: pairs (c, d) with n*c+d | (n-1)*d."""
    return sum(1 for c in range(1, n) for d in range(1, n)
               if ((n - 1) * d) % (n * c + d) == 0)


def write_bfiles(N):
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "..", "..", "data", "zero_insertion")
    os.makedirs(out, exist_ok=True)
    fN = open(os.path.join(out, "b_count_N.txt"), "w")
    fL = open(os.path.join(out, "b_maxlength_L.txt"), "w")
    fM = open(os.path.join(out, "b_largest.txt"), "w")
    fS = open(os.path.join(out, "b_smallest.txt"), "w")
    rec, best = [], -1
    for b in range(2, N + 1):
        n = count(b)
        fN.write(f"{b} {n}\n")
        if n > best:
            best = n
            rec.append((b, n))
        if b <= min(N, 2000):
            Z = members(b)
            fL.write(f"{b} {max([x[3] for x in Z], default=-1) + 1}\n")
            fM.write(f"{b} {Z[-1][0] if Z else 0}\n")
            fS.write(f"{b} {Z[0][0] if Z else 0}\n")
    with open(os.path.join(out, "records.txt"), "w") as f:
        for i, (b, n) in enumerate(rec, 1):
            f.write(f"{i} {b} {n}\n")
    with open(os.path.join(out, "table_members_b2_100.txt"), "w") as f:
        i = 0
        for b in range(2, 101):
            for k, *_ in members(b):
                i += 1
                f.write(f"{i} {k}\n")
    print("records:", rec)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "verify"
    if cmd == "verify":
        verify(int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    elif cmd == "members":
        b = int(sys.argv[2])
        for k, a, c, m in members(b):
            print(k, digits(k, b), f"a={a} c={c} length={m+1}")
    elif cmd == "bfiles":
        write_bfiles(int(sys.argv[2]) if len(sys.argv) > 2 else 10000)
