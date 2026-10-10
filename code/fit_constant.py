"""Fit sum_{B<=x} S(B) against C*x*(log x)^3 + lower-order terms.

Usage: python3 fit_constant.py [data/seed_count/partial_sums_2e6.txt]
Input file: lines "x T(x)" produced by seed_sum.c.
"""
import os
import sys
import numpy as np

# Predicted constant C = (1/12) * prod_p (1-1/p)^2 (1+2/p)
N = 10**7
sieve = np.ones(N, bool)
sieve[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
p = np.nonzero(sieve)[0].astype(float)
G = np.exp(np.sum(2 * np.log1p(-1 / p) + np.log1p(2 / p)))
C = G / 12
print(f"Euler product G = {G:.10f}")
print(f"Predicted C = G/12 = {C:.10f}")

d = np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "data", "seed_count", "partial_sums_2e6.txt"))
x, T = d[:, 0], d[:, 1]
L = np.log(x)
y = T / x
m = x >= 8000

A = np.vstack([L**3, L**2, L, np.ones_like(L)]).T[m]
c, *_ = np.linalg.lstsq(A, y[m], rcond=None)
print("Free cubic fit  [L^3, L^2, L, 1]:", c, " max resid:", np.abs(A @ c - y[m]).max())

A2 = np.vstack([L**2, L, np.ones_like(L)]).T[m]
r = (y - C * L**3)[m]
c2, *_ = np.linalg.lstsq(A2, r, rcond=None)
print("C fixed, fit    [L^2, L, 1]:     ", c2, " max resid:", np.abs(A2 @ c2 - r).max())
