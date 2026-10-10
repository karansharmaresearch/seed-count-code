"""Seeds with a polynomial base P in F_q[x] (q prime here).

Digits are nonzero polynomials of degree < m = deg P; (a, b) is a seed if
aP + b divides aP^n + b for all n >= 1, equivalently aP + b | (P - 1) b.
The q^m - 1 pairs (a, -a) are always seeds ("trivial").

Checks, by brute force over all monic P of small degree:
  - nontrivial seeds exist  <=>  P and P - 1 are both reducible;
  - the count equals sum over h != 0, deg h < m, deg L_h < m of
    q^(m - deg L_h) - 1, with L_h = hP / (gcd(h+1, P) gcd(h, P-1)).
"""
import itertools


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return tuple(a)


def deg(a):
    return len(a) - 1


def add(a, b, q):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % q
                 for i in range(n)])


def mul(a, b, q):
    if not a or not b:
        return ()
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % q
    return trim(r)


def mod(a, b, q):
    a, inv = list(a), pow(b[-1], q - 2, q)
    while len(a) >= len(b) and a:
        c, s = a[-1] * inv % q, len(a) - len(b)
        for i, y in enumerate(b):
            a[s + i] = (a[s + i] - c * y) % q
        a = list(trim(a))
    return tuple(a)


def gcd(a, b, q):
    while b:
        a, b = b, mod(a, b, q)
    return a


def nonzero_polys(q, maxdeg):
    for n in range(maxdeg + 1):
        for co in itertools.product(range(q), repeat=n):
            for lead in range(1, q):
                yield co + (lead,)


def monics(q, n):
    for co in itertools.product(range(q), repeat=n):
        yield co + (1,)


def irreducible(P, q):
    return all(mod(P, f, q) for d in range(1, deg(P) // 2 + 1) for f in monics(q, d))


def seeds(P, q):
    m, Pm1 = deg(P), add(P, (q - 1,), q)
    out = []
    for a in nonzero_polys(q, m - 1):
        for b in nonzero_polys(q, m - 1):
            if not mod(mul(Pm1, b, q), add(mul(a, P, q), b, q), q):
                out.append((a, b))
    return out


def formula(P, q):
    m, Pm1, tot = deg(P), add(P, (q - 1,), q), 0
    for h in nonzero_polys(q, m - 1):
        dL = deg(h) + m - deg(gcd(add(h, (1,), q), P, q)) - deg(gcd(h, Pm1, q))
        if dL < m:
            tot += q ** (m - dL) - 1
    return tot


if __name__ == "__main__":
    for q, maxdeg in ((2, 5), (3, 4), (5, 3), (7, 2)):
        for m in range(1, maxdeg + 1):
            for P in monics(q, m):
                s = seeds(P, q)
                nontrivial = [ab for ab in s if add(*ab, q)]
                assert len(s) - len(nontrivial) == q ** m - 1
                Pm1 = add(P, (q - 1,), q)
                both_reducible = not irreducible(P, q) and not irreducible(Pm1, q)
                assert bool(nontrivial) == both_reducible, (q, P)
                assert len(s) == formula(P, q), (q, P)
        print(f"q={q}, deg P <= {maxdeg}: criterion and count formula hold")
    # paper's examples over F_2
    assert [ab for ab in seeds((0, 1, 1), 2) if add(*ab, 2)] == []   # x^2 + x
    assert ((1,), (0, 1)) in seeds((0, 0, 1), 2)                       # (1, x) in base x^2
    print("examples: x^2+x has only trivial seeds; (1, x) is a seed in base x^2")
