/* Brute-force check of the largest-seed theorem (Paper A, Theorem 12):
 * for composite B with smallest prime factor p, the largest seed is
 * B(B-1)/p, attained only by (a,b) = (B/p - 1, B - B/p).
 *
 * Uses only Theorem 1 ((a,b) is a seed iff Ba+b | (B-1)b) and tests every
 * pair 1 <= a,b <= B-1, independently of the parametrization of Lemma 3.
 * Also prints S(B) for comparison with the b-file of A394974.
 * Output for 2 <= B <= 1e4: data/seed_count/largest_seed_check_1e4.txt.
 *
 * Usage: largest_seed_check B0 B1     (prints "B S(B) maxseed a b ok")
 * Build: python -m ziglang cc -O2 -o largest_seed_check largest_seed_check.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s B0 B1\n", argv[0]); return 1; }
  uint32_t B0 = (uint32_t)atol(argv[1]), B1 = (uint32_t)atol(argv[2]);
  if (B0 < 2) B0 = 2;
  long bad = 0;
  for (uint32_t B = B0; B <= B1; B++) {
    uint32_t p = 0;
    for (uint32_t d = 2; (uint64_t)d * d <= B; d++) if (B % d == 0) { p = d; break; }
    uint64_t S = 0, best = 0; uint32_t ba = 0, bb = 0, nbest = 0;
    for (uint32_t a = 1; a < B; a++)
      for (uint32_t b = 1; b < B; b++) {
        uint64_t N1 = (uint64_t)B * a + b;
        if (((uint64_t)(B - 1) * b) % N1 == 0) {
          S++;
          if (N1 > best) { best = N1; ba = a; bb = b; nbest = 1; }
          else if (N1 == best) nbest++;
        }
      }
    int ok;
    if (p == 0) ok = (S == 0);                       /* prime: no seeds */
    else ok = (best == (uint64_t)B * (B - 1) / p && ba == B / p - 1 && bb == B - B / p && nbest == 1);
    if (!ok) bad++;
    printf("%u %llu %llu %u %u %d\n", B, (unsigned long long)S, (unsigned long long)best, ba, bb, ok);
  }
  fprintf(stderr, "range %u..%u: %ld failures\n", B0, B1, bad);
  return bad != 0;
}
