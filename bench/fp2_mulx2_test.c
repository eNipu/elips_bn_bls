/* Correctness test for bench/fp2_mulx2_x86_64.S. Build:
 *   gcc -O3 -std=c11 -Iinclude -DELIPS_CURVE_BLS12_381 \\
 *       bench/fp2_mulx2_test.c bench/fp2_mulx2_x86_64.S bench/lazy_kernels.S \\
 *       -o x2test build/libelips_arith_BLS12_381.a -lm && ./x2test
 *
 * elips_fp2_mulx2 must agree with the C Karatsuba it replaces, limb for limb. */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "elips/fp.h"
#include "elips/fpx.h"
#include "elips/random.h"
void elips_fp_mulw_6_x86_64(limb_t w[12], const limb_t a[6], const limb_t b[6]);
void elips_fp2_mulx2_6_x86_64(limb_t c0[12], limb_t c1[12],
                              const limb_t a[2][6], const limb_t b[2][6], const limb_t p[6]);
static void subm(limb_t r[12], const limb_t a[12], const limb_t b[12])
{
    limb_t br=0,c=0,hi[6],cand[6];
    for(int i=0;i<6;i++){unsigned __int128 s=(unsigned __int128)a[i]-b[i]-br;r[i]=(limb_t)s;br=(limb_t)((s>>64)&1);}
    for(int i=0;i<6;i++){unsigned __int128 s=(unsigned __int128)a[i+6]-b[i+6]-br;hi[i]=(limb_t)s;br=(limb_t)((s>>64)&1);}
    for(int i=0;i<6;i++){unsigned __int128 s=(unsigned __int128)hi[i]+FP_MODULUS[i]+c;cand[i]=(limb_t)s;c=(limb_t)(s>>64);}
    limb_t t=(limb_t)0-br;
    for(int i=0;i<6;i++) r[i+6]=(cand[i]&t)|(hi[i]&~t);
}
static void ref(limb_t c0[12], limb_t c1[12], const fp2_t a, const fp2_t b)
{
    limb_t t0[12],t1[12],s[12]; fp_t as,bs;
    elips_fp_mulw_6_x86_64(t0,a[0],b[0]);
    elips_fp_mulw_6_x86_64(t1,a[1],b[1]);
    fp_add(as,a[0],a[1]); fp_add(bs,b[0],b[1]);
    elips_fp_mulw_6_x86_64(s,as,bs);
    subm(c0,t0,t1); subm(s,s,t0); subm(c1,s,t1);
}
static void rf(fp_t x){limb_t k[FP_LIMBS];elips_random_scalar(k);fp_from_limbs(x,k);}
int main(int argc,char**argv){
    long n=(argc>1)?atol(argv[1]):200000; long bad=0;
    fp2_t a,b; limb_t r0[12],r1[12],g0[12],g1[12];
    /* edge: zero and one */
    fp2_set_zero(a); fp2_set_one(b);
    ref(r0,r1,a,b); elips_fp2_mulx2_6_x86_64(g0,g1,(const limb_t(*)[6])a,(const limb_t(*)[6])b,FP_MODULUS);
    if(memcmp(r0,g0,96)||memcmp(r1,g1,96)){printf("FAIL edge 0*1\n");bad++;}
    fp2_set_one(a);
    ref(r0,r1,a,b); elips_fp2_mulx2_6_x86_64(g0,g1,(const limb_t(*)[6])a,(const limb_t(*)[6])b,FP_MODULUS);
    if(memcmp(r0,g0,96)||memcmp(r1,g1,96)){printf("FAIL edge 1*1\n");bad++;}
    for(long i=0;i<n;i++){
        rf(a[0]);rf(a[1]);rf(b[0]);rf(b[1]);
        ref(r0,r1,a,b);
        elips_fp2_mulx2_6_x86_64(g0,g1,(const limb_t(*)[6])a,(const limb_t(*)[6])b,FP_MODULUS);
        if(memcmp(r0,g0,96)||memcmp(r1,g1,96)){ bad++; if(bad<3) printf("FAIL at %ld\n",i); }
        if((i+1)%50000==0){printf("%8ld ok\n",i+1);fflush(stdout);}
    }
    printf("%s  %ld rounds, %ld mismatches\n",bad?"FAILED":"PASS",n,bad);
    return bad!=0;
}
