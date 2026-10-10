"""Correctness checks for A394974 / S(n).

Run:  python3 verify.py
"""
from math import isqrt

from seed_count import S, S_brute


def isprime(n):
    return n > 1 and all(n % i for i in range(2, isqrt(n) + 1))


def check_formula(limit=400):
    """The closed form must match brute-force enumeration."""
    bad = [n for n in range(2, limit + 1) if S(n) != S_brute(n)]
    assert not bad, f"formula disagrees with brute force at: {bad}"
    print(f"OK  closed form == brute force for n = 2..{limit}")


def check_primality(limit=3000):
    """S(n) = 0 must hold exactly for prime n."""
    bad = [n for n in range(2, limit + 1) if (S(n) == 0) != isprime(n)]
    assert not bad, f"primality characterization fails at: {bad}"
    print(f"OK  S(n) = 0 iff n prime for n = 2..{limit}")


if __name__ == "__main__":
    check_formula()
    check_primality()
    print("All checks passed.")
