"""S(n) = A394974: number of base-n seeds whose zero-insertions stay divisible.

S(n) counts pairs (c, d) with 1 <= c, d <= n-1 such that (n*c + d) divides (n-1)*d.
Equivalently, n*c + d divides every c*n^k + d for k >= 1.

Key facts:
    S(n) = sum_{k=1..n-2} floor((n-1) / L_k),
        with L_k = n*k / (gcd(k+1, n) * gcd(k, n-1))
    S(n) = 0  if and only if  n is prime.
"""
from math import gcd


def S(n):
    """Closed-form value of S(n) for a base n >= 2."""
    total = 0
    for k in range(1, n - 1):
        L = (n * k) // (gcd(k + 1, n) * gcd(k, n - 1))
        if L <= n - 1:
            total += (n - 1) // L
    return total


def S_brute(n):
    """Brute-force S(n) by direct enumeration (ground truth)."""
    return sum(
        1
        for c in range(1, n)
        for d in range(1, n)
        if (n - 1) * d % (n * c + d) == 0
    )


def seeds(n):
    """List the seeds (c, d, value = n*c + d) in base n."""
    return [
        (c, d, n * c + d)
        for c in range(1, n)
        for d in range(1, n)
        if (n - 1) * d % (n * c + d) == 0
    ]


def write_bfile(path="b394974.txt", nmax=10000):
    """Write an OEIS-format b-file for n = 2..nmax."""
    with open(path, "w") as f:
        f.write(f"# Table of n, a(n) for n = 2..{nmax}\n")
        for n in range(2, nmax + 1):
            f.write(f"{n} {S(n)}\n")


if __name__ == "__main__":
    print("S(n) for n = 2..30:")
    print([S(n) for n in range(2, 31)])
    print("\nSeeds in base 10:", seeds(10))
    print("Seeds in base 6: ", seeds(6))
