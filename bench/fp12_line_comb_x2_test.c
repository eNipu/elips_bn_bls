/* Correctness test for src/arith/fp12_line_comb_x2_x86_64.S. Build:
 *   gcc -O3 -std=c11 -Iinclude -Isrc -DELIPS_CURVE_BLS12_381 \
 *       bench/fp12_line_comb_x2_test.c src/arith/fp12_line_comb_x2_x86_64.S \
 *       -o linetest build/libelips_arith_BLS12_381.a -lm && ./linetest */
/* elips_fp12_line_comb_x2 must match the C combination limb for limb. */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <x86intrin.h>
#include "elips/fp.h"
#include "elips/fpx.h"
#include "elips/random.h"
typedef limb_t fpw_t[12];
void elips_fp12_line_comb_x2_6_x86_64(limb_t out[12][12],
                                      const limb_t prod[26][12],
                                      const limb_t p[6]);

static void addm(fpw_t r, const fpw_t a, const fpw_t b){
    unsigned char c=0,br=0; limb_t hi[6],cand[6];
    for(int i=0;i<6;i++) c=_addcarry_u64(c,a[i],b[i],(unsigned long long*)&r[i]);
    for(int i=0;i<6;i++) c=_addcarry_u64(c,a[i+6],b[i+6],(unsigned long long*)&hi[i]);
    for(int i=0;i<6;i++) br=_subborrow_u64(br,hi[i],FP_MODULUS[i],(unsigned long long*)&cand[i]);
    limb_t keep=(limb_t)0-(limb_t)(br&(c^1));
    for(int i=0;i<6;i++) r[i+6]=(hi[i]&keep)|(cand[i]&~keep);
}
static void subm(fpw_t r, const fpw_t a, const fpw_t b){
    unsigned char br=0,c=0; limb_t hi[6],cand[6];
    for(int i=0;i<6;i++) br=_subborrow_u64(br,a[i],b[i],(unsigned long long*)&r[i]);
    for(int i=0;i<6;i++) br=_subborrow_u64(br,a[i+6],b[i+6],(unsigned long long*)&hi[i]);
    for(int i=0;i<6;i++) c=_addcarry_u64(c,hi[i],FP_MODULUS[i],(unsigned long long*)&cand[i]);
    limb_t take=(limb_t)0-(limb_t)br;
    for(int i=0;i<6;i++) r[i+6]=(cand[i]&take)|(hi[i]&~take);
}
/* c0 = P1 + xi(P2), c1 = P3 + xi(P4), c2 = P5 - P1 - P4, on one five-product
 * block starting at slot base. */
static void block(fpw_t *out, const limb_t p[26][12], int base){
    fpw_t w0,w1;
    subm(w0,p[base+2],p[base+3]); addm(w1,p[base+2],p[base+3]);
    addm(out[0],p[base+0],w0); addm(out[1],p[base+1],w1);
    subm(w0,p[base+6],p[base+7]); addm(w1,p[base+6],p[base+7]);
    addm(out[2],p[base+4],w0); addm(out[3],p[base+5],w1);
    subm(w0,p[base+8],p[base+0]); subm(out[4],w0,p[base+6]);
    subm(w0,p[base+9],p[base+1]); subm(out[5],w0,p[base+7]);
}
static void ref(limb_t out[12][12], const limb_t p[26][12]){
    fpw_t t0[6], s[6], w0, w1;
    block(t0,p,0); block(s,p,10);
    for(int i=0;i<6;i++){ subm(w0,s[i],t0[i]); subm(out[6+i],w0,p[20+i]); }
    subm(w0,p[24],p[25]); addm(w1,p[24],p[25]);
    addm(out[0],t0[0],w0); addm(out[1],t0[1],w1);
    addm(out[2],t0[2],p[20]); addm(out[3],t0[3],p[21]);
    addm(out[4],t0[4],p[22]); addm(out[5],t0[5],p[23]);
}
/* random wide value in [0, p*R) */
static void rw(limb_t w[12]){
    limb_t k[FP_LIMBS];
    for(int i=0;i<6;i++){ elips_random_scalar(k); w[i]=k[0]; }
    fp_t h; elips_random_scalar(k); fp_from_limbs(h,k);
    for(int i=0;i<6;i++) w[i+6]=h[i];     /* high half < p, so w < p*R */
}
int main(int argc,char**argv){
    long n=(argc>1)?atol(argv[1]):200000; long bad=0;
    limb_t prod[26][12], r[12][12], g[12][12];
    for(long i=0;i<n;i++){
        for(int j=0;j<26;j++) rw(prod[j]);
        ref(r,(const limb_t(*)[12])prod);
        elips_fp12_line_comb_x2_6_x86_64(g,(const limb_t(*)[12])prod,FP_MODULUS);
        if(memcmp(r,g,sizeof r)){ bad++; if(bad<3){ for(int j=0;j<12;j++) if(memcmp(r[j],g[j],96)) printf("FAIL at %ld, out[%d]\n",i,j);} }
        if((i+1)%50000==0){printf("%8ld ok\n",i+1);fflush(stdout);}
    }
    printf("%s  %ld rounds, %ld mismatches\n",bad?"FAILED":"PASS",n,bad);
    return bad!=0;
}
