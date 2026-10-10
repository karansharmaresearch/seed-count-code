// Statistics of N(B)/B for 3 <= B <= X, with the same algorithm as n_fast.c
// (N(B) via the second exact formula, Proposition 7.2 of the merged paper).
// Usage: n_ratio_stats <X> <outfile>
// Build: python -m ziglang cc -O2 -o n_ratio_stats n_ratio_stats.c
// Output lines:
//   G B N ratio   B > 1026 with N(B) > B           (flag: 1 if B = 2^j+2)
//   M B N ratio   overall maximum of N(B)/B
//   O B N ratio   running maximum of N(B)/B over B NOT of the form 2^j+2
//   D d B N ratio maximum of N(B)/B over B in (10^d, 10^(d+1)], B != 2^j+2
//   F j B N ratio B = 3^j+3
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static uint32_t *spf;
static uint64_t gcdu(uint64_t a, uint64_t b){ while(b){uint64_t t=a%b;a=b;b=t;} return a; }
static uint64_t S[1<<20], D[1<<16];

static uint64_t Ncount(uint64_t b){
    uint32_t P[16]; int k=0; uint64_t m=b, radb=1;
    while(m>1){ uint32_t p=spf[m]; P[k++]=p; radb*=p; while(m%p==0) m/=p; }
    int ns=1; S[0]=1;
    for(int i=0;i<k;i++){ int cur=ns; for(int j=0;j<cur;j++){ uint64_t v=S[j]*P[i]; while(v<b){ S[ns++]=v; v*=P[i]; } } }
    int nd=1; D[0]=1; m=b-1;
    while(m>1){ uint32_t p=spf[m]; int e=0; while(m%p==0){m/=p;e++;} int cur=nd; uint64_t pp=1;
        for(int a=1;a<=e;a++){ pp*=p; for(int j=0;j<cur;j++) D[nd++]=D[j]*pp; } }
    uint64_t N=0;
    for(int i=0;i<nd;i++){ uint64_t n=(b-1)/D[i]; if(n<3) continue;
        for(int j=0;j<ns;j++){ uint64_t s=S[j]; if(s<2||s>=n) continue;
            uint64_t t=n/s, r=n-s*t;
            if(gcdu(t,radb)!=1) continue;
            if(gcdu(r,n)!=1) continue;
            N+=(s-1)/r; } }
    return N;
}

static int is2j2(uint64_t b){ uint64_t m=b-2; return b>=3 && m>0 && (m&(m-1))==0; }

int main(int argc, char**argv){
    uint64_t X = (uint64_t)atoll(argv[1]);
    FILE *f = fopen(argv[2],"w");
    spf = calloc(X+2, sizeof(uint32_t));
    for (uint64_t i=2;i<=X+1;i++) if(!spf[i]) for(uint64_t j=i;j<=X+1;j+=i) if(!spf[j]) spf[j]=(uint32_t)i;
    double best=0, boff=0, dmax[12]={0}; uint64_t bB=0,bN=0, dB[12]={0}, dN[12]={0};
    uint64_t p3=1; int j3=0;
    for (uint64_t b=3;b<=X;b++){
        uint64_t N=Ncount(b); double r=(double)N/b;
        if(b>1026 && N>b) fprintf(f,"G %llu %llu %.6f %d\n",(unsigned long long)b,(unsigned long long)N,r,is2j2(b));
        if(r>best){best=r;bB=b;bN=N;}
        if(!is2j2(b)){
            if(r>boff){boff=r; fprintf(f,"O %llu %llu %.6f\n",(unsigned long long)b,(unsigned long long)N,r);}
            int d=0; uint64_t t=b-1; while(t>=10){t/=10;d++;}   // b in (10^d, 10^(d+1)]
            if(r>dmax[d]){dmax[d]=r;dB[d]=b;dN[d]=N;}
        }
        while(p3+3<b){p3*=3;j3++;}
        if(p3+3==b) fprintf(f,"F %d %llu %llu %.6f\n",j3,(unsigned long long)b,(unsigned long long)N,r);
    }
    fprintf(f,"M %llu %llu %.6f\n",(unsigned long long)bB,(unsigned long long)bN,best);
    for(int d=0;d<12;d++) if(dB[d]) fprintf(f,"D %d %llu %llu %.6f\n",d,(unsigned long long)dB[d],(unsigned long long)dN[d],dmax[d]);
    fclose(f); return 0;
}
