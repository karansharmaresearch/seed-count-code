"""check3 / c6: constant G, Lemma 5.6 local factors, Prop 5.8, Cor 5.9, and
the numerical table/fits of Section 6 (from data/seed_count/partial_sums_1e8.txt, whose
values up to 1e7 are independently reproduced by c5_direct.c).
"""
import numpy as np
from math import log
from fractions import Fraction as F

# --- Lemma 5.6: local factor of h, exactly, for a few primes
for p in (2, 3, 5, 7, 11):
    # f(p^a) for a=0..6, convolve with mu*mu local series (1-X)^2
    f = [F(1)] + [2 * (1 - F(1, p))] * 6
    h = [f[a] - 2 * (f[a - 1] if a >= 1 else 0) + (f[a - 2] if a >= 2 else 0) for a in range(7)]
    assert h[1] == F(-2, p) and h[2] == -F(p - 2, p) and all(v == 0 for v in h[3:]), (p, h)
    assert 1 + h[1] / p + h[2] / p ** 2 == (1 - F(1, p)) ** 2 * (1 + F(2, p))
print("Lemma 5.6 local factors: OK for p in 2..11")

# --- G by Euler product to 2e7 plus tail estimate
N = 20_000_000
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
P = np.nonzero(s)[0].astype(np.float64)
logG = np.sum(np.log1p(-3 / P ** 2 + 2 / P ** 3))
tail = -3 / (N * log(N))  # sum_{p>N} -3/p^2 ~ -3/(N log N)
G = float(np.exp(logG + tail))
print("G = %.9f   G/12 = %.9f   (paper: 0.2867474, 0.0238956)" % (G, G / 12))

# --- Prop 5.8: sum_{n<=y} f(n)/n  vs  G/2 log^2 y, with f(n)/n = 2^w phi / n^2
M = 10_000_000
phi = np.arange(M + 1, dtype=np.float64)
w = np.zeros(M + 1, dtype=np.int8)
for p in np.nonzero(s[:M + 1])[0]:
    phi[p::p] *= (1 - 1 / p)
    w[p::p] += 1
n = np.arange(M + 1, dtype=np.float64); n[0] = 1
fn_over_n = (2.0 ** w) * phi / n ** 2
fn_over_n[0] = 0
cs = np.cumsum(fn_over_n)
print("Prop 5.8: (sum_{n<=y} f(n)/n - G/2 log^2 y) / (1+log y):")
for y in (10, 100, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    print("   y=%-9d %.5f" % (y, (cs[y] - G / 2 * log(y) ** 2) / (1 + log(y))))
# Cor 5.9: M(x) = x/2 sum f(n)/n log(x/n)
print("Cor 5.9: (M(x) - G/12 x log^3 x) / (x log^2 x):")
for x in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    Mx = x / 2 * np.sum(fn_over_n[1:x + 1] * np.log(x / n[1:x + 1]))
    print("   x=%-9d %.5f" % (x, (Mx - G / 12 * x * log(x) ** 3) / (x * log(x) ** 2)))

# --- Section 6 table and fits
data = np.loadtxt("../../data/seed_count/partial_sums_1e8.txt", usecols=(0, 1))
X, T = data[:, 0], data[:, 1]
for x0 in (1e6, 1e7, 1e8):
    i = np.where(X == x0)[0][0]
    L = log(x0)
    print("x=%.0e T=%d  T/(x L^3)=%.4f  (T-G/12 x L^3)/(x L^2)=%.4f"
          % (x0, T[i], T[i] / (x0 * L ** 3), (T[i] - G / 12 * x0 * L ** 3) / (x0 * L ** 2)))
sel = X >= 8000
L = np.log(X[sel]); y = T[sel] / X[sel]
A3 = np.vstack([L ** 2, L, np.ones_like(L)]).T
coef3 = np.linalg.lstsq(A3, y - G / 12 * L ** 3, rcond=None)[0]
print("fixed-leading fit (beta, a1, a0) =", np.round(coef3, 3), " paper: 0.193, -0.14, -0.16")
A4 = np.vstack([L ** 3, L ** 2, L, np.ones_like(L)]).T
coef4 = np.linalg.lstsq(A4, y, rcond=None)[0]
print("free 4-param fit, leading coeff = %.4f (paper: 0.0239); all:" % coef4[0], np.round(coef4, 3))
sel2 = (X >= 1e6) & (X <= 1e7)
th = np.polyfit(np.log(np.log(X[sel2])), np.log(T[sel2] / X[sel2]), 1)[0]
print("local theta over [1e6,1e7] = %.3f (paper: ~2.7)" % th)
