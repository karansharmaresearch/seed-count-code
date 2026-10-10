"""check3 / c3: Lemma 5.4 (a)-(d), including the n = 1 support question.

For each x, every coprime (u',v') with n=u'v' <= x+2 and every t <= 3x+5:
  A_brute = sum over k in K(u',v') (k>=1) with (k+1)(tk+n) <= xn of
            floor((k+1)/n) - [v'=1]
compared with Q(Q+1)/2 - eps*Q (Q from the definition), and (b),(c),(d)
checked. Also random large (x,n,t,rho) tests of (b),(c) with exact integers.
"""
from math import gcd, isqrt
import random, sys

def crt_rho(up, vp):
    n = up * vp
    # k = -1 mod u', k = 0 mod v'; rho = (k+1) mod n
    for k in range(1, n + 1):
        if (k + 1) % up == 0 and k % vp == 0:
            return (k + 1) % n
    raise AssertionError


def Q_of(x, n, t, rho):
    lo, hi = 0, isqrt(x) + 5
    # largest q>=1 with (nq+rho)(t(nq+rho-1)+n) <= xn ; 0 if none
    Q = 0
    q = 1
    # monotone -> binary search
    lo, hi = 1, 2 * isqrt(x) + 10
    if (n + rho) * (t * (n + rho - 1) + n) > x * n:
        return 0
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if (n * mid + rho) * (t * (n * mid + rho - 1) + n) <= x * n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def A_brute(x, up, vp, t):
    n = up * vp
    eps = 1 if vp == 1 else 0
    s = 0
    k = 1
    while (k + 1) * (t * k + n) <= x * n:
        if (k + 1) % up == 0 and k % vp == 0:
            s += (k + 1) // n - eps
        k += 1
    return s


XMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
errs = []
maxratio = 0.0
max_t_n1 = {}
for x in range(3, XMAX + 1):
    for n in range(1, x + 3):
        pairs = [(up, n // up) for up in range(1, n + 1)
                 if n % up == 0 and gcd(up, n // up) == 1]
        for up, vp in pairs:
            rho = crt_rho(up, vp)
            eps = 1 if rho == 0 else 0
            if eps != (1 if vp == 1 else 0):
                errs.append(("rho", up, vp))
            for t in range(1, 3 * x + 6):
                A = A_brute(x, up, vp, t)
                Q = Q_of(x, n, t, rho)
                if A != Q * (Q + 1) // 2 - eps * Q:
                    errs.append(("a", x, up, vp, t, A, Q))
                R = (x / (n * t)) ** 0.5
                if not (max(0.0, R - 3) - 1e-12 <= Q <= R + 1 + 1e-12):
                    errs.append(("b", x, up, vp, t, Q, R))
                main = x / (2 * n * t) if t * n <= x else 0.0
                dev = abs(A - main)
                if dev > 3.5 * R + 2.5 + 1e-9:
                    errs.append(("c", x, up, vp, t))
                maxratio = max(maxratio, dev / (3.5 * R + 2.5))
                if A != 0:
                    if not 2 * n - 1 <= x:
                        errs.append(("d1", x, n, t))
                    if n >= 2 and not t * (n - 1) < x:
                        errs.append(("d2", x, n, t))
                    if n == 1:
                        max_t_n1[x] = max(max_t_n1.get(x, 0), t)
                if (A != 0 or t * n <= x) and t > 2 * x:
                    errs.append(("support t<=2x", x, n, t))
print("exhaustive x<=%d errors:" % XMAX, errs[:10], len(errs))
print("max |A - main| / (3.5R+2.5) =", round(maxratio, 4))
print("n=1: largest t with A_{1,1}(t) != 0, vs x/2-1 :",
      [(x, max_t_n1.get(x), x / 2 - 1) for x in (3, 4, 10, 50, XMAX)])

# random large tests of (b) and (c) directly from Q (any rho in [0,n-1]) --
# (a) makes A a function of Q, eps only.
random.seed(1)
bad = []
worst = 0.0
for _ in range(200000):
    x = random.randint(3, 10 ** random.randint(1, 12))
    n = random.randint(1, max(1, x))
    t = random.randint(1, max(1, 10 ** random.randint(0, 12)))
    rho = random.randint(0, n - 1)
    eps = 1 if rho == 0 else 0
    Q = Q_of(x, n, t, rho)
    R2 = x / (n * t)
    R = R2 ** 0.5
    if not (max(0.0, R - 3) - 1e-9 <= Q <= R + 1 + 1e-9):
        bad.append(("b", x, n, t, rho, Q, R))
    A = Q * (Q + 1) // 2 - eps * Q
    main = R2 / 2 if t * n <= x else 0.0
    dev = abs(A - main)
    worst = max(worst, dev / (3.5 * R + 2.5))
    if dev > 3.5 * R + 2.5 + 1e-6 * (1 + R):
        bad.append(("c", x, n, t, rho))
print("random (b),(c) failures:", bad[:5], len(bad), " worst ratio", round(worst, 4))
