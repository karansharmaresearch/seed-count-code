"""check3 / c7: Proposition 3.3 lower bound S(pq) >= 2 gcd(p-1,q-1), pq <= 30000,
and the conditional family m+1, 2m+1 prime: S/sqrt(B) along it."""
from math import gcd, isqrt
def S(B): return sum(((B-1)*gcd(k+1,B)*gcd(k,B-1))//(B*k) for k in range(1,B-1))
def isp(n): return n>1 and all(n%d for d in range(2,isqrt(n)+1))
ps=[p for p in range(2,200) if isp(p)]
bad=[(p,q) for p in ps for q in ps if p<q and p*q<=30000 and S(p*q)<2*gcd(p-1,q-1)]
print("Prop 3.3 failures (pq<=30000):", bad)
fam=[(m,S((m+1)*(2*m+1))/((m+1)*(2*m+1))**0.5) for m in range(1,120) if isp(m+1) and isp(2*m+1)]
print("S(B)/sqrt(B) along m+1,2m+1 prime:", [(m,round(r,3)) for m,r in fam])
