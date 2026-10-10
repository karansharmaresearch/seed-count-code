/* check3 / c5: independent direct computation of S(B) for all B <= X
 * (divisor-pair enumeration of formula (1.1), NOT the bijection of Sec. 5).
 *
 * A summand k of (1.1) is nonzero only if D*E > k, where D=gcd(k+1,B),
 * E=gcd(k,B-1) (Thm 4.1 Step 3 -- itself re-verified exhaustively for
 * B <= 3000 in c1). So for each divisor pair (D | B, E | B-1) take the CRT
 * solution k in [1, DE), keep it only if its TRUE gcds are (D,E) (so each k is
 * counted once), and add floor((B-1) D E / (B k)).
 *
 * Checks for every B: S(B) < 2 tau(B) tau(B-1) sqrt(B) (exactly, in integers).
 * Reports T(x) at checkpoints, S(16), S(36), all B with S(B)/B >= 4/9,
 * the largest S(B)/B with B != 36, and max S(B)/(tau tau sqrt B).
 * Usage: ./c5 X
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef long long ll;
static ll g(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t;} return a; }
static ll inv(ll a, ll m){ /* a^{-1} mod m, gcd=1, m>=1 */
    ll t=0,nt=1,r=m,nr=a%m; if(m==1) return 0;
    while(nr){ ll q=r/nr, tmp=t-q*nt; t=nt; nt=tmp; tmp=r-q*nr; r=nr; nr=tmp; }
    if(t<0) t+=m; return t;
}
int main(int argc, char**argv){
    ll X = atoll(argv[1]);
    int *cnt = calloc(X+2, sizeof(int));
    for(ll d=1; d<=X; d++) for(ll m=d; m<=X; m+=d) cnt[m]++;
    ll *start = malloc((X+3)*sizeof(ll)); start[0]=0;
    for(ll m=0;m<=X+1;m++) start[m+1]=start[m]+(m<=X?cnt[m]:0);
    int *dv = malloc(start[X+1]*sizeof(int));
    int *fill = calloc(X+2,sizeof(int));
    for(ll d=1; d<=X; d++) for(ll m=d; m<=X; m+=d) dv[start[m]+fill[m]++]=(int)d;
    ll T=0; ll next=1000;
    double bestr=0; ll bestB=0; double best2=0; ll best2B=0; double maxc=0; ll maxcB=0;
    ll viol=0;
    for(ll B=3; B<=X; B++){
        ll S=0;
        for(ll i=start[B]; i<start[B]+cnt[B]; i++){
            ll D=dv[i];
            for(ll j=start[B-1]; j<start[B-1]+cnt[B-1]; j++){
                ll E=dv[j];
                ll DE=D*E; if(DE>=B-1) continue;
                /* k = 0 mod E, k = -1 mod D  -> k = E*s, E s = -1 mod D */
                ll s = (D==1)?0:((D - inv(E%D, D)) % D);
                ll k = E*s; if(k==0) k=DE; /* representative in [1,DE] */
                if(k>=DE || k>B-2) continue;
                if(g(k+1,B)!=D || g(k,B-1)!=E) continue;
                S += ((B-1)*DE)/(B*k);
            }
        }
        ll tt = (ll)cnt[B]*cnt[B-1];
        /* S < 2 tt sqrt(B)  <=>  S^2 < 4 tt^2 B */
        if(!((__int128)S*S < (__int128)4*tt*tt*B)) { viol++; printf("VIOLATION B=%lld\n",B); }
        double r=(double)S/B;
        if(9*S>=4*B) printf("S(B)/B>=4/9: B=%lld S=%lld\n",B,S);
        if(r>bestr){bestr=r;bestB=B;}
        if(B!=36 && r>best2){best2=r;best2B=B;}
        double c=S/(tt*sqrt((double)B)); if(c>maxc){maxc=c;maxcB=B;}
        if(B==16||B==36) printf("S(%lld)=%lld\n",B,S);
        T+=S;
        if(B==next || B==X){ printf("T(%lld)=%lld\n",B,T); fflush(stdout); if(B==next) next*= (next%3==0? 1:1), next = (next==1000?10000:next*10); }
    }
    printf("violations of Thm 4.1: %lld\n",viol);
    printf("max S/B = %.6f at B=%lld; max over B!=36: %.6f at B=%lld\n",bestr,bestB,best2,best2B);
    printf("max S/(tau(B)tau(B-1)sqrt B) = %.5f at B=%lld (theorem: < 2)\n",maxc,maxcB);
    return 0;
}
