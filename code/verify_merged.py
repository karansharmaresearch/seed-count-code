"""Numerical checks of the new or adapted statements in the merged paper
zero-insertion.tex (numbering as in that paper). Consistency tests only;
nothing in the proofs depends on them.

Run from the repository root:  python code/verify_merged.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "zero_insertion"))
import zi  # noqa: E402  (classification: members(b), count(b), digits, in_Zb_bruteforce)


def is_prime(n):
    return n >= 2 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def is_prime_power(n):
    if n < 2:
        return False
    p = next(d for d in range(2, n + 1) if n % d == 0)
    while n % p == 0:
        n //= p
    return n == 1


def tau(n):
    return sum(1 for d in range(1, n + 1) if n % d == 0)


def ndigits(z, B):
    return len(zi.digits(z, B))


def bpart(m, B):
    """m^{(B)}: largest divisor of m all of whose primes divide B."""
    g = math.gcd(m, B)
    out = 1
    while g > 1:
        out *= g
        m //= g
        g = math.gcd(m, B)
    return out


# Lemma 3.5: z | I_i(z) for all i  <=>  z | (B-1) w, for every z with 2..4 digits.
def check_lemma_3_5(Bmax=24, imax=12):
    n = 0
    for B in range(2, Bmax + 1):
        for z in range(B, min(B ** 4, 200000)):
            d = zi.digits(z, B)
            l = len(d)
            a, w = d[0], z - d[0] * B ** (l - 1)
            lhs = all((a * B ** (l - 1 + i) + w) % z == 0 for i in range(imax + 1))
            assert lhs == ((B - 1) * w % z == 0), (B, z)
            n += 1
    print(f"Lemma 3.5: {n} numbers z (2 <= B <= {Bmax}, up to 4 digits, {imax} insertions): ok")


# Theorem 6.6 and Proposition 6.7: length bound, and when it is attained.
def check_length(Bmax=2000):
    for B in range(3, Bmax + 1):
        Z = [x[0] for x in zi.members(B)]
        J = (B - 2).bit_length() - 1
        Lstar = max((ndigits(z, B) for z in Z), default=0)
        assert Lstar <= J + 1, B
        assert len(Z) < B * math.log(B), B
        if B >= 5:
            assert (Lstar == J + 1) == (B % 4 == 2), B
    assert max(ndigits(z, 4) for z, *_ in zi.members(4)) == 2
    print(f"Theorem 6.6, Proposition 6.7: 3 <= B <= {Bmax}: ok")


# Theorem 6.8 and Theorem 6.9 (a), (c).
def check_prime_and_pp(Bmax=2000):
    S = {}
    with open(os.path.join(os.path.dirname(__file__), "..", "oeis", "b394974.txt")) as fh:
        for line in fh:
            f = line.split()
            if len(f) == 2 and f[0].isdigit():
                S[int(f[0])] = int(f[1])
    for B in range(2, Bmax + 1):
        Z = [x[0] for x in zi.members(B)]
        assert (len(Z) == 0) == is_prime(B), B
        if not is_prime(B):
            assert len(Z) >= S[B] >= tau(B) - 2, B
        N2 = sum(1 for z in Z if ndigits(z, B) == 2)
        assert N2 == S[B], B
        if is_prime_power(B):
            assert N2 == len(Z), B
        elif not is_prime(B):
            assert any(ndigits(z, B) == 3 for z in Z), B
    print(f"Theorem 6.8, Theorem 6.9(a),(c): 2 <= B <= {Bmax} (N_2 = S against the b-file): ok")


# Theorem 6.9(b): N_3(B) = U(B), with U from the definition in Paper A.
U_PAPER_A = {6: 2, 10: 2, 12: 4, 14: 2, 15: 1, 18: 3, 20: 4, 21: 6, 22: 3, 24: 3,
             26: 4, 28: 5, 30: 13, 33: 2, 34: 3, 35: 2, 36: 7, 38: 2, 39: 6, 40: 4}


def U_def(B):
    return sum(1 for a in range(1, B) for c in range(1, B) for b in range(1, B)
               if ((B - 1) * (c * B + b)) % (a * B * B + c * B + b) == 0)


def check_N3_U(Bmax=60):
    for B in range(2, Bmax + 1):
        N3 = sum(1 for z, *_ in zi.members(B) if ndigits(z, B) == 3)
        U = U_def(B)
        assert N3 == U, (B, N3, U)
        if B <= 40:
            assert U == U_PAPER_A.get(B, 0), B
        assert (U == 0) == (is_prime_power(B) or B < 2), B
    print(f"Theorem 6.9(b): N_3 = U (U from its definition) for 2 <= B <= {Bmax}, "
          "and the table of Paper A for B <= 40: ok")


# Corollary 7.1 and Proposition 7.2 against the listing.
def second_formula(B):
    tot = 0
    for g in range(1, B):
        if (B - 1) % g:
            continue
        nu = (B - 1) // g
        for r in range(1, nu):
            if math.gcd(r, nu) == 1:
                tot += (bpart(nu - r, B) - 1) // r
    return tot


def check_formulas(Bmax=1500):
    for B in range(3, Bmax + 1):
        n = len(zi.members(B))
        assert zi.count(B) == n == second_formula(B), B
    print(f"Corollary 7.1, Proposition 7.2: 3 <= B <= {Bmax}: ok")


# Theorem 8.4: the two families, and the last statement about k = 1.
def check_families():
    for j in range(1, 11):
        B = 2 ** j + 2
        Z = {z: (a, k) for z, a, k, m in zi.members(B)}
        fam = [z for z, (a, k) in Z.items() if k == 1]
        assert len(fam) == B - 3, B
        for z, (a, k) in Z.items():
            if k == 1 and a % 2 == 1:
                assert ndigits(z, B) == j + 1, (B, z)
            if B <= 1100 and k == 1:
                assert zi.in_Zb_bruteforce(z, B), (B, z)
    for j in range(1, 8):
        B = 3 ** j + 3
        Z = {z: (a, k) for z, a, k, m in zi.members(B)}
        fam = sorted(a for z, (a, k) in Z.items() if k == 2)
        assert fam == list(range(1, (3 ** j - 1) // 2 + 1)), B
        assert len(fam) == (B - 4) // 2 and (B - 4) % 2 == 0, B
        for z, (a, k) in Z.items():
            if k == 2:
                assert zi.in_Zb_bruteforce(z, B), (B, z)
    assert sorted(z for z, a, k, m in zi.members(12) if k == 2) == [44, 176, 352, 704]
    for B in range(3, 2001):
        k1 = sum(1 for z, a, k, m in zi.members(B) if k == 1)
        pow2 = (B - 2) & (B - 3) == 0
        assert (k1 == B - 3) == pow2, B
    assert [j for j in range(0, 400) if (3 ** j + 1) & 3 ** j == 0] == [0, 1]
    print("Theorem 8.4 (2^j+2, j <= 10; 3^j+3, j <= 7; members of 3^j+3 also by brute "
          "force), k = 1 statement for B <= 2000, 3^j+1 = 2^i only for j = 0, 1: ok")


if __name__ == "__main__":
    check_lemma_3_5()
    check_N3_U()
    check_families()
    check_length()
    check_prime_and_pp()
    check_formulas()
