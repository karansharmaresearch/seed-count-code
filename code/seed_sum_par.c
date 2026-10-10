/* Parallel-friendly exact computation of T(x) = sum_{B<=x} S(B) over a range
 * of bases, plus the unfloored envelope U(x) = sum_{B<=x} sum_{k=1}^{B-2} (B-1)/L_k.
 *
 * Usage: seed_sum_par B0 B1      (processes B0 <= B <= B1, B0 >= 2)
 * Output (partial sums start at B0):
 *   C c S U   at every checkpoint c in [B0,B1]
 *   E B1 S U  at the end of the range
 *   M r Bmax  largest S(B)/B over B != 36 in the range
 *
 * Differences from seed_sum.c (same S(B), validated against it):
 *  - B and B-1 are factored with a segmented sieve (small memory per process);
 *  - 128-bit products, so (B-1)*u*v cannot overflow for B up to 1e8 and beyond;
 *  - pairs (u,v) with u*v >= B-1 are skipped: for an exact pair
 *    u = gcd(k+1,B), v = gcd(k,B-1) one always has u*v < B-1 (paper,
 *    Theorem 11, first step), so such pairs contribute nothing.
 * U(B) uses gcd(n,m) = sum_{d | gcd} phi(d):
 *   U(B) = (B-1)/B * sum_{d|B, e|B-1} phi(d) phi(e) sum_{k<=B-2, k=-1 (d), k=0 (e)} 1/k,
 * the inner harmonic sum over an arithmetic progression via the digamma function.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef long long ll;
typedef __int128 i128;

static ll gcdll(ll a, ll b) { while (b) { ll t = a % b; a = b; b = t; } return a; }

static ll inv_mod(ll a, ll m) {           /* a^{-1} mod m, gcd(a,m)=1, m>=2 */
  ll t = 0, nt = 1, r = m, nr = a % m;
  while (nr) {
    ll q = r / nr, tmp = t - q * nt; t = nt; nt = tmp;
    tmp = r - q * nr; r = nr; nr = tmp;
  }
  if (t < 0) t += m;
  return t;
}

static double digamma(double z) {        /* z > 0 */
  double r = 0.0;
  while (z < 8.0) { r -= 1.0 / z; z += 1.0; }
  double z2 = 1.0 / (z * z);
  return r + log(z) - 0.5 / z
         - z2 * (1.0 / 12 - z2 * (1.0 / 120 - z2 * (1.0 / 252 - z2 * (1.0 / 240 - z2 / 132))));
}

/* sum_{j=0}^{J} 1/(r + j m) */
static double harm_ap(ll r, ll m, ll J) {
  if (J < 16) { double s = 0; for (ll j = J; j >= 0; j--) s += 1.0 / (double)(r + j * m); return s; }
  double a = (double)r / (double)m;
  return (digamma(a + (double)J + 1.0) - digamma(a)) / (double)m;
}

static int is_checkpoint(ll c) {
  static const ll dec[] = {1000000LL, 2000000LL, 5000000LL, 10000000LL, 20000000LL,
                           50000000LL, 100000000LL, 200000000LL, 0};
  for (int i = 0; dec[i]; i++) if (c == dec[i]) return 1;
  if (c < 1000 || c % 1000) return 0;
  ll q = c / 1000;
  return (q & (q - 1)) == 0;              /* c = 1000 * 2^j */
}

#define MAXF 12
#define MAXD 8192
#define BS 32768

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s B0 B1\n", argv[0]); return 1; }
  ll B0 = atoll(argv[1]), B1 = atoll(argv[2]);
  if (B0 < 2) B0 = 2;
  ll P = (ll)sqrt((double)B1) + 2;
  char *comp = calloc(P + 1, 1);
  ll *pr = malloc(sizeof(ll) * (P + 1)); int np = 0;
  for (ll i = 2; i <= P; i++) if (!comp[i]) { pr[np++] = i; for (ll j = i * i; j <= P; j += i) comp[j] = 1; }

  /* factorization buffers for n in [lo-1, hi] */
  ll *rem = malloc(sizeof(ll) * (BS + 1));
  int *nf = malloc(sizeof(int) * (BS + 1));
  ll (*fp)[MAXF] = malloc(sizeof(ll[MAXF]) * (BS + 1));
  int (*fe)[MAXF] = malloc(sizeof(int[MAXF]) * (BS + 1));
  static ll dB[MAXD], pB[MAXD], dC[MAXD], pC[MAXD];

  ll S = 0; double U = 0.0, Uc = 0.0; /* Kahan compensation for U */
  double bestr = 0; ll bestB = 0;
#ifdef RECORDS
  ll recS = -1;
#endif

  for (ll lo = B0; lo <= B1; lo += BS) {
    ll hi = lo + BS - 1; if (hi > B1) hi = B1;
    ll base = lo - 1, len = hi - base + 1;
    for (ll i = 0; i < len; i++) { rem[i] = base + i; nf[i] = 0; }
    for (int t = 0; t < np; t++) {
      ll p = pr[t]; if (p * p > hi) break;
      ll st = ((base + p - 1) / p) * p;
      for (ll n = st; n <= hi; n += p) {
        ll i = n - base; int e = 0;
        while (rem[i] % p == 0) { rem[i] /= p; e++; }
        fp[i][nf[i]] = p; fe[i][nf[i]] = e; nf[i]++;
      }
    }
    for (ll i = 0; i < len; i++) if (rem[i] > 1) { fp[i][nf[i]] = rem[i]; fe[i][nf[i]] = 1; nf[i]++; }

    for (ll B = lo; B <= hi; B++) {
      ll iB = B - base, iC = B - 1 - base;
      int nb = 1, nc = 1; dB[0] = 1; pB[0] = 1; dC[0] = 1; pC[0] = 1;
      for (int f = 0; f < nf[iB]; f++) {
        ll p = fp[iB][f]; int e = fe[iB][f]; int m = nb; ll pw = 1, ph = p - 1;
        for (int t = 1; t <= e; t++) { pw *= p; for (int i = 0; i < m; i++) { dB[nb] = dB[i] * pw; pB[nb] = pB[i] * ph; nb++; } ph *= p; }
      }
      for (int f = 0; f < nf[iC]; f++) {
        ll p = fp[iC][f]; int e = fe[iC][f]; int m = nc; ll pw = 1, ph = p - 1;
        for (int t = 1; t <= e; t++) { pw *= p; for (int i = 0; i < m; i++) { dC[nc] = dC[i] * pw; pC[nc] = pC[i] * ph; nc++; } ph *= p; }
      }
      ll SB = 0; double UB = 0.0;
      for (int i = 0; i < nb; i++) {
        ll u = dB[i];
        for (int j = 0; j < nc; j++) {
          ll v = dC[j];
          ll M = u * v;                               /* <= B(B-1): fits in 64 bits */
          ll r = (u == 1) ? v : v * ((u - (inv_mod(v % u, u))) % u);  /* r=-1 (u), r=0 (v), 1<=r<=M */
          if (r == 0) r = M;
          /* unfloored envelope term: all k <= B-2 in the progression r mod M */
          if (r <= B - 2) UB += (double)pB[i] * (double)pC[j] * harm_ap(r, M, (B - 2 - r) / M);
          /* floored count: exact pairs only, and those have u*v < B-1 */
          if (M >= B - 1) continue;
          for (ll k = r; k <= B - 2 && (i128)k * B <= (i128)(B - 1) * M; k += M) {
            if (gcdll(k + 1, B) != u || gcdll(k, B - 1) != v) continue;
            SB += (ll)(((i128)(B - 1) * M) / ((i128)B * k));
          }
        }
      }
      UB *= (double)(B - 1) / (double)B;
      S += SB;
      double y = UB - Uc, t = U + y; Uc = (t - U) - y; U = t;
      if (B != 36 && (double)SB / (double)B > bestr) { bestr = (double)SB / (double)B; bestB = B; }
#ifdef RECORDS
      /* chunk-local records of S(B); records_S.py merges them across chunks */
      if (SB > recS) { recS = SB; printf("R %lld %lld\n", B, SB); }
#endif
      if (is_checkpoint(B) && B != B1) printf("C %lld %lld %.6f\n", B, S, U);
    }
  }
  printf("E %lld %lld %.6f\n", B1, S, U);
  printf("M %.10f %lld\n", bestr, bestB);
  return 0;
}
