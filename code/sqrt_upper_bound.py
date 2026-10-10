"""Check Theorem 4.1 of the merged paper (zero-insertion.tex): S(B) < 2 tau(B) tau(B-1) sqrt(B).

For every multiplier k with 1 <= k <= B-2, D = gcd(k+1, B) and E = gcd(k, B-1)
satisfy DE < B-1 (Step 1), and each nonzero summand floor((B-1)/L_k) is
< 2 sqrt(B) (Step 2); a nonzero summand has DE > k and distinct k give
distinct pairs (D, E) (Step 3). Also checks the closed form (1.1) (Theorems 3.1-3.2,
from the companion paper) against brute force for small B.
"""
import math
from math import gcd, isqrt


def tau(n):
    return sum(2 - (d * d == n) for d in range(1, isqrt(n) + 1) if n % d == 0)


def S_checked(B):
    total = 0
    pairs = set()
    for k in range(1, B - 1):
        D, E = gcd(k + 1, B), gcd(k, B - 1)
        assert D * E < B - 1, (B, k)                    # Thm 4.1, Step 1
        c = (B - 1) * D * E // (B * k)                  # floor((B-1)/L_k)
        if c:
            assert c <= D * E / k < 2 * math.sqrt(B), (B, k)   # Step 2
            assert D * E > k and (D, E) not in pairs, (B, k)   # Step 3
            pairs.add((D, E))
            total += c
    return total


if __name__ == "__main__":
    for B in range(2, 200):
        brute = sum(1 for a in range(1, B) for b in range(1, B)
                    if ((B - 1) * b) % (B * a + b) == 0)
        assert brute == S_checked(B), B
    print("closed form = brute force, B < 200")

    worst = (0.0, 0)
    for B in range(3, 3001):
        s = S_checked(B)
        bound = 2 * tau(B) * tau(B - 1) * math.sqrt(B)
        assert s < bound, B
        worst = max(worst, (s / bound, B))
    print("DE < B-1, per-k < 2 sqrt(B), and S(B) < 2 tau(B) tau(B-1) sqrt(B) "
          f"(Theorem 4.1) hold for 3 <= B <= 3000 (largest S/bound = {worst[0]:.3f} at B={worst[1]})")
