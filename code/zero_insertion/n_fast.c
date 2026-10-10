// Usage: n_fast <X> <outfile>   (compile: python -m ziglang cc -O2 -o n_fast n_fast.c)
// N(b) for 3 <= b <= X via N(b) = sum_{g|b-1} sum_{s in S_b, s<n} [t=floor(n/s), r=n-s t,
//   gcd(t,b)=1, gcd(r,n)=1] floor((s-1)/r),  n=(b-1)/g.
// Output: "T x T(x)" at x = 1,2,5 x 10^k (T(x) = sum_{b<=x} N(b)), and "B b N(b) N(b)/b"
// for every b where N(b)/b sets a record or exceeds 1, and "R b N(b)" whenever N(b) sets a record.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint32_t *spf;
static uint64_t gcdu(uint64_t a, uint64_t b){ while(b){uint64_t t=a%b;a=b;b=t;} return a; }

int main(int argc, char**argv){
    uint32_t X = (uint32_t)atol(argv[1]);
    const char *out = argv[2];
    spf = calloc(X+2, sizeof(uint32_t));
    for (uint32_t i=2;i<=X+1;i++) if(!spf[i]) for(uint64_t j=i;j<=X+1;j+=i) if(!spf[j]) spf[j]=i;
    FILE *f = fopen(out,"w");
    static uint64_t S[1<<20]; static uint64_t D[1<<16];
    uint64_t T=0; uint64_t nextcp=1000;
    double best=0; uint64_t maxN=0;
    for (uint64_t b=3;b<=X;b++){
        // primes of b
        uint32_t P[16]; int k=0; uint64_t m=b, radb=1;
        while(m>1){ uint32_t p=spf[m]; P[k++]=p; radb*=p; while(m%p==0) m/=p; }
        // S_b up to b-1
        int ns=1; S[0]=1;
        for(int i=0;i<k;i++){ int cur=ns; for(int j=0;j<cur;j++){ uint64_t v=S[j]*P[i]; while(v<b){ S[ns++]=v; v*=P[i]; } } }
        // divisors of b-1
        int nd=1; D[0]=1; m=b-1;
        while(m>1){ uint32_t p=spf[m]; int e=0; while(m%p==0){m/=p;e++;} int cur=nd; uint64_t pp=1;
            for(int a=1;a<=e;a++){ pp*=p; for(int j=0;j<cur;j++) D[nd++]=D[j]*pp; } }
        uint64_t N=0;
        for(int i=0;i<nd;i++){ uint64_t n=(b-1)/D[i]; if(n<3) continue;
            for(int j=0;j<ns;j++){ uint64_t s=S[j]; if(s<2||s>=n) continue;
                uint64_t t=n/s, r=n-s*t;
                if(gcdu(t,radb)!=1) continue;
                if(gcdu(r,n)!=1) continue;
                uint64_t c=(s-1)/r; N+=c; } }
        T+=N;
        if(N>maxN){ maxN=N; fprintf(f,"R %llu %llu\n",(unsigned long long)b,(unsigned long long)N); }
        double ratio=(double)N/b;
        if(ratio>best || ratio>1.0){ if(ratio>best) best=ratio; fprintf(f,"B %llu %llu %.4f\n",(unsigned long long)b,(unsigned long long)N,ratio); }
        if(b==nextcp){ fprintf(f,"T %llu %llu\n",(unsigned long long)b,(unsigned long long)T); fflush(f);
            // checkpoints at 1,2,5 x 10^k
            uint64_t d=1; while(d*10<=b) d*=10; uint64_t lead=b/d;
            nextcp = (lead==1)?2*d:(lead==2)?5*d:10*d; }
    }
    fclose(f); return 0;
}
