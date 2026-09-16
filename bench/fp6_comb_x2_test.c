/* Correctness test for src/arith/fp6_comb_x2_x86_64.S. Build:
 *   gcc -O3 -std=c11 -Iinclude -DELIPS_CURVE_BLS12_381 \
 *       bench/fp6_comb_x2_test.c src/arith/fp6_comb_x2_x86_64.S \
 *       -o combtest build/libelips_arith_BLS12_381.a -lm && ./combtest */
/* elips_fp6_comb_x2 must match the C combination limb for limb. */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <x86intrin.h>
#include "elips/fp.h"
#include "elips/fpx.h"
#include "elips/random.h"
typedef limb_t fpw_t[12];
void elips_fp6_comb_x2_6_x86_64(limb_t out[6][12], const limb_t prod[12][12], const limb_t p[6]);

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
/* prod index: v0=0,1  v1=2,3  v2=4,5  w01=6,7  w02=8,9  w12=10,11 */
static void ref(limb_t out[6][12], const limb_t prod[12][12]){
    fpw_t t0,t1,u0,u1;
    subm(t0,prod[10],prod[2]); subm(t0,t0,prod[4]);
    subm(t1,prod[11],prod[3]); subm(t1,t1,prod[5]);
    subm(u0,t0,t1); addm(u1,t0,t1);
    addm(out[0],u0,prod[0]); addm(out[1],u1,prod[1]);
    subm(t0,prod[6],prod[0]); subm(t0,t0,prod[2]);
    subm(t1,prod[7],prod[1]); subm(t1,t1,prod[3]);
    subm(u0,prod[4],prod[5]); addm(u1,prod[4],prod[5]);
    addm(out[2],t0,u0); addm(out[3],t1,u1);
    subm(t0,prod[8],prod[0]); subm(t0,t0,prod[4]);
    subm(t1,prod[9],prod[1]); subm(t1,t1,prod[5]);
    addm(out[4],t0,prod[2]); addm(out[5],t1,prod[3]);
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
    limb_t prod[12][12], r[6][12], g[6][12];
    for(long i=0;i<n;i++){
        for(int j=0;j<12;j++) rw(prod[j]);
        ref(r,prod);
        elips_fp6_comb_x2_6_x86_64(g,(const limb_t(*)[12])prod,FP_MODULUS);
        if(memcmp(r,g,sizeof r)){ bad++; if(bad<3){ for(int j=0;j<6;j++) if(memcmp(r[j],g[j],96)) printf("FAIL at %ld, out[%d]\n",i,j);} }
        if((i+1)%50000==0){printf("%8ld ok\n",i+1);fflush(stdout);}
    }
    printf("%s  %ld rounds, %ld mismatches\n",bad?"FAILED":"PASS",n,bad);
    return bad!=0;
}
