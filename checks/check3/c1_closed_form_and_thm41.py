"""check3 / c1: sanity of the input closed form, and Theorem 4.1 step by step.

(i)   S(B) from the seed DEFINITION (Ba+b | aB^n+b for all n>=1, tested for
      n up to a bound that covers preperiod+period of B^n mod (Ba+b)) equals
      the closed form (1.1), for B <= 60.
(ii)  For every B <= 3000 and every k in [1,B-2]:
        Step 1: D*E < B-1
        Step 2: floor((B-1)/L_k) < 2*sqrt(D*E)   (strict, checked exactly)
        Step 3: nonzero summand => k is the unique CRT solution in [1,DE)
      and the final inequality S(B) < 2 tau(B) tau(B-1) sqrt(B).
"""
from math import gcd, isqrt


def S_closed(B):
    s = 0
    for k in range(1, B - 1):
        s += ((B - 1) * gcd(k + 1, B) * gcd(k, B - 1)) // (B * k)
    return s


def S_def(B):
    cnt = 0
    for a in range(1, B):
        for b in range(1, B):
            m = B * a + b
            # B^n mod m is eventually periodic with preperiod <= log2(m) and
            # period <= m, so n in [1, 2m+64] covers every residue of B^n.
            ok = True
            p = B % m
            for n in range(1, 2 * m + 64):
                if (a * p + b) % m:
                    ok = False
                    break
                p = p * B % m
            cnt += ok
    return cnt


def tau(n):
    return sum(1 for d in range(1, n + 1) if n % d == 0)


def lt_2sqrt(val, X):
    """exact test val < 2*sqrt(X) for integers val>=0, X>=1"""
    return val * val < 4 * X


bad = []
for B in range(2, 61):
    if S_def(B) != S_closed(B):
        bad.append(B)
print("(i) definition vs closed form, B<=60, mismatches:", bad)

worst41 = 0.0
worst_step2 = 0.0
viol = []
for B in range(3, 3001):
    S = 0
    seen_pairs = set()
    for k in range(1, B - 1):
        D, E = gcd(k + 1, B), gcd(k, B - 1)
        if not D * E < B - 1:
            viol.append(("step1", B, k))
        term = ((B - 1) * D * E) // (B * k)
        if not lt_2sqrt(term, D * E):
            viol.append(("step2", B, k))
        if term:
            worst_step2 = max(worst_step2, term / (2 * (D * E) ** 0.5))
            if not D * E > k:
                viol.append(("step3-DE>k", B, k))
            if (D, E) in seen_pairs:
                viol.append(("step3-inj", B, k))
            seen_pairs.add((D, E))
        S += term
    tb = tau(B) * tau(B - 1)
    if not lt_2sqrt(S, tb * tb * B):
        viol.append(("thm4.1", B))
    worst41 = max(worst41, S / (2 * tb * B ** 0.5))
print("(ii) violations of Thm 4.1 steps, 3<=B<=3000:", viol[:10], len(viol))
print("     max summand/(2 sqrt(DE)) =", round(worst_step2, 4))
print("     max S(B)/(2 tau(B) tau(B-1) sqrt B) =", round(worst41, 4))
print("S(16) =", S_closed(16), " S(36) =", S_closed(36))
