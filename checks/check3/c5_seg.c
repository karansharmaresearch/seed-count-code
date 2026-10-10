/* check3 / c5_seg: segmented version of c5_direct.c for large X.
 *
 * Same algorithm (divisor-pair enumeration of formula (1.1), independent of
 * Section 5), but divisors are produced block by block with a segmented
 * sieve (d <= sqrt(hi), cofactor m/d), so memory is O(block size) and
 * disjoint ranges can run as separate processes.
 *
 * Usage: ./c5_seg LO HI     (handles 3 <= LO <= B <= HI)
 * Output:
 *   "P <B> <sum of S over LO..B>"  for every B that is a multiple of 1000,
 *   "R B S" for every B with 9 S >= 4 B,  "V B" for any violation of
 *   S(B) < 2 tau(B) tau(B-1) sqrt(B),  and at the end the best S/B with
 *   B != 36, the max of S/(tau tau sqrt B), and the range total.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef long long ll;
static ll g(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t;} return a; }
static ll inv(ll a, ll m){
    ll t=0,nt=1,r=m,nr=a%m; if(m==1) return 0;
    while(nr){ ll q=r/nr, tmp=t-q*nt; t=nt; nt=tmp; tmp=r-q*nr; r=nr; nr=tmp; }
    if(t<0) t+=m; return t;
}
#define W 1000000LL

int main(int argc, char**argv){
    ll LO=atoll(argv[1]), HI=atoll(argv[2]);
    ll sq=(ll)sqrt((double)HI)+2;
    int *cnt=malloc((W+2)*sizeof(int)); ll *st=malloc((W+3)*sizeof(ll));
    ll cap=0; ll *dv=NULL;
    ll run=0, total=0, viol=0; double best2=0; ll best2B=0; double maxc=0; ll maxcB=0;
    for(ll a=LO-1; a<=HI; a+=W){           /* block covers m in [a, a+W] */
        ll b=a+W; if(b>HI) b=HI; ll len=b-a+1;
        for(ll i=0;i<len;i++) cnt[i]=0;
        for(ll d=1; d<=sq; d++){
            ll m0=((a+d-1)/d)*d; if(m0<d*d) m0=d*d;
            for(ll m=m0; m<=b; m+=d){ cnt[m-a]+= (m/d==d)?1:2; }
        }
        st[0]=0; for(ll i=0;i<len;i++) st[i+1]=st[i]+cnt[i];
        if(st[len]>cap){ cap=st[len]; dv=realloc(dv,cap*sizeof(ll)); }
        for(ll i=0;i<len;i++) cnt[i]=0;
        for(ll d=1; d<=sq; d++){
            ll m0=((a+d-1)/d)*d; if(m0<d*d) m0=d*d;
            for(ll m=m0; m<=b; m+=d){
                ll i=m-a; dv[st[i]+cnt[i]++]=d;
                if(m/d!=d) dv[st[i]+cnt[i]++]=m/d;
            }
        }
        for(ll B=(a+1>LO?a+1:LO); B<=b; B++){
            ll iB=B-a, iC=B-1-a; ll S=0;
            for(ll x=st[iB]; x<st[iB+1]; x++){
                ll D=dv[x];
                for(ll y=st[iC]; y<st[iC+1]; y++){
                    ll E=dv[y]; ll DE=D*E; if(DE>=B-1) continue;
                    ll s=(D==1)?0:((D-inv(E%D,D))%D);
                    ll k=E*s; if(k==0) k=DE;
                    if(k>=DE||k>B-2) continue;
                    if(g(k+1,B)!=D||g(k,B-1)!=E) continue;
                    S+=((B-1)*DE)/(B*k);
                }
            }
            ll tt=(st[iB+1]-st[iB])*(st[iC+1]-st[iC]);
            if(!((__int128)S*S<(__int128)4*tt*tt*B)){ viol++; printf("V %lld\n",B); }
            if(9*S>=4*B) printf("R %lld %lld\n",B,S);
            double r=(double)S/B; if(B!=36&&r>best2){best2=r;best2B=B;}
            double c=S/(tt*sqrt((double)B)); if(c>maxc){maxc=c;maxcB=B;}
            run+=S;
            if(B%1000==0) printf("P %lld %lld\n",B,run);
        }
        if(b==HI) break;
        a=b-W; /* next block starts at b (so B-1=b is included) */
        fflush(stdout);
    }
    printf("END total=%lld viol=%lld best_nonc36=%.6f@%lld maxc=%.5f@%lld\n",
           run,viol,best2,best2B,maxc,maxcB);
    return 0;
}
