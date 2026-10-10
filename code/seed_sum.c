#include <stdio.h>
#include <stdlib.h>
typedef long long ll;
static ll g(ll a, ll b){while(b){ll t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){
  ll X=atoll(argv[1]);
  // smallest prime factor sieve
  int *spf=calloc(X+2,sizeof(int));
  for(ll i=2;i<=X+1;i++) if(!spf[i]) for(ll j=i;j<=X+1;j+=i) if(!spf[j]) spf[j]=i;
  ll divB[4096], divC[4096];
  double T=0; ll next=1000;
  for(ll B=2;B<=X;B++){
    // divisors of B and B-1
    int nb=1,nc=1; divB[0]=1; divC[0]=1;
    ll n=B; while(n>1){ll p=spf[n];int e=0;while(n%p==0){n/=p;e++;}int m=nb;ll pw=1;for(int t=1;t<=e;t++){pw*=p;for(int i=0;i<m;i++)divB[nb++]=divB[i]*pw;}}
    n=B-1; while(n>1){ll p=spf[n];int e=0;while(n%p==0){n/=p;e++;}int m=nc;ll pw=1;for(int t=1;t<=e;t++){pw*=p;for(int i=0;i<m;i++)divC[nc++]=divC[i]*pw;}}
    ll S=0;
    for(int i=0;i<nb;i++){ ll u=divB[i];
      for(int j=0;j<nc;j++){ ll v=divC[j];
        ll M=u*v; // k ≡ -1 mod u, k ≡ 0 mod v, k <= (B-1)uv/B < uv
        // find k0 = v*t, v*t ≡ -1 mod u
        ll k0;
        if(u==1) k0=v; else { // inverse of v mod u
          ll a=v%u, m=u, x0=0,x1=1; ll aa=a,mm=m; while(aa>1){ll q=aa/mm,t=mm;mm=aa%mm;aa=t;t=x0;x0=x1-q*x0;x1=t;} if(x1<0)x1+=u;
          ll t=((u-1)%u)*x1%u; k0=v*t; if(k0==0)k0=M; }
        for(ll k=k0;k<=B-2 && k*B<=(B-1)*M; k+=M){
          if(g(k+1,B)!=u || g(k,B-1)!=v) continue;
          S += ((B-1)*M)/(B*k);
        }
      }
    }
    T+=S;
    if(B==next||B==X){printf("%lld %.0f\n",B,T); next*=2;}
  }
}
