"""Analyse T(x) = sum S(B) and the unfloored envelope U(x) from run_par.py output.

Usage: python analyze_avg.py data/seed_count/partial_sums_1e8.txt
Input lines: "x T(x) U(x)".

Reports
  - fit_constant.py-style fits (free cubic in L = log x; C fixed + quadratic),
    on the full range x >= 8000 and on the upper range x >= 1e5;
  - sensitivity: max residual when C is fixed at other values;
  - implied second-order coefficient (T/x - C L^3)/L^2 at each checkpoint;
  - rho(x) = T/U (the conjecture C = G/12 with c_T = G/4 needs rho -> 1/3);
  - U/(x L^3) against G/4.
"""
import sys
import numpy as np


def euler_G(N=10**7):
    s = np.ones(N, bool)
    s[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    p = np.nonzero(s)[0].astype(float)
    return float(np.exp(np.sum(2 * np.log1p(-1 / p) + np.log1p(2 / p))))


def fits(x, y, L, lo, C):
    m = x >= lo
    A = np.vstack([L**3, L**2, L, np.ones_like(L)]).T[m]
    c, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    A2 = np.vstack([L**2, L, np.ones_like(L)]).T[m]
    r = (y - C * L**3)[m]
    c2, *_ = np.linalg.lstsq(A2, r, rcond=None)
    return c, np.abs(A @ c - y[m]).max(), c2, np.abs(A2 @ c2 - r).max(), int(m.sum())


def main():
    d = np.loadtxt(sys.argv[1])
    x, T, U = d[:, 0], d[:, 1], d[:, 2]
    L, y = np.log(x), d[:, 1] / d[:, 0]
    G = euler_G()
    C = G / 12
    print(f"G = {G:.10f}   C = G/12 = {C:.10f}   G/4 = {G/4:.6f}\n")

    for lo in (8000, 1e5):
        c, r1, c2, r2, n = fits(x, y, L, lo, C)
        print(f"x >= {lo:g} ({n} points)")
        print(f"  free cubic  [L^3,L^2,L,1] = {np.round(c, 5)}   max resid {r1:.4f}")
        print(f"  C fixed     [L^2,L,1]     = {np.round(c2, 4)}   max resid {r2:.4f}")
        print("  sensitivity (C fixed -> max resid):", end="")
        for Cf in (0.018, 0.020, 0.022, C, 0.026, 0.028, 0.030):
            m = x >= lo
            A2 = np.vstack([L**2, L, np.ones_like(L)]).T[m]
            r = (y - Cf * L**3)[m]
            cc, *_ = np.linalg.lstsq(A2, r, rcond=None)
            print(f"  {Cf:.4f}:{np.abs(A2 @ cc - r).max():.4f}", end="")
        print("\n")

    print(f"{'x':>11} {'T/(xL^3)':>9} {'(T/x-CL^3)/L^2':>15} {'rho=T/U':>8} {'U/(xL^3)':>9}")
    for xi, Ti, Ui, Li in zip(x, T, U, L):
        print(f"{xi:11.0f} {Ti/xi/Li**3:9.5f} {(Ti/xi - C*Li**3)/Li**2:15.4f} {Ti/Ui:8.5f} {Ui/xi/Li**3:9.5f}")


if __name__ == "__main__":
    main()
