/* Research-only prototype. No production files are modified.
 * The 8M2 kernel computes 2*f*L, NOT the exact raw Miller value.
 * Its extra Fp scalar is killed only by final exponentiation.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pairing/miller.c"

typedef struct { fp2_t B, C; } norm_line;
typedef struct { norm_line l[ELIPS_MILLER_LINES]; } norm_table;
static unsigned long long state = 0x137481;
static unsigned long long rng(void)
{
    state ^= state << 13; state ^= state >> 7; state ^= state << 17;
    return state;
}
static void rf(fp_t r)
{
    limb_t a[FP_LIMBS] = {0};
    /* Full-width canonical test values, not cryptographic randomness. */
    int less;
    do {
        for (int i=0;i<FP_LIMBS;i++) a[i]=rng();
        a[FP_LIMBS-1] &= (((limb_t)1 << (FP_BITS%64))-1);
        less=0;
        for (int i=FP_LIMBS-1;i>=0;i--) {
            if (a[i]!=FP_MODULUS[i]) { less=a[i]<FP_MODULUS[i]; break; }
        }
    } while (!less);
    fp_from_limbs(r, a);
}
static void r2(fp2_t a) { rf(a[0]); rf(a[1]); }
static void r12(fp12_t a)
{
    for (int i=0;i<2;i++) for (int j=0;j<3;j++) r2(a[i][j]);
}
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

/* 2*A*(B+C*v), v^3=xi, with four Fp2 products.
 * Polynomial evaluation at 0,+1,-1,infinity. Keeping the factor 2
 * eliminates interpolation divisions without requiring fp_half.
 */
static void linear4(fp6_t r, const fp6_t a, const fp2_t B,
                    const fp2_t C, const fp2_t plus, const fp2_t minus)
{
    fp2_t t0,t3,tp,tm,s,u,q0,q1,q2;
    fp2_mul(t0,a[0],B);
    fp2_mul(t3,a[2],C);
    fp2_add(s,a[0],a[2]);
    fp2_add(u,s,a[1]); fp2_mul(tp,u,plus);
    fp2_sub(u,s,a[1]); fp2_mul(tm,u,minus);
    fp2_mul_xi(u,t3); fp2_add(q0,t0,u); fp2_add(q0,q0,q0);
    fp2_add(u,t3,t3); fp2_sub(q1,tp,tm); fp2_sub(q1,q1,u);
    fp2_add(u,t0,t0); fp2_add(q2,tp,tm); fp2_sub(q2,q2,u);
    fp2_copy(r[0],q0); fp2_copy(r[1],q1); fp2_copy(r[2],q2);
}

/* L=1+B*w^3+C*w^5; f=A+D*w.
 * 2*f*L = (2*A+v^2*linear4(D)) + (2*D+v*linear4(A))*w.
 */
static void normalized8(fp12_t f, const fp2_t B, const fp2_t C)
{
    fp2_t plus,minus;
    fp6_t hA,hD,A,D;
    fp2_add(plus,B,C); fp2_sub(minus,B,C);
    linear4(hA,f[0],B,C,plus,minus);
    linear4(hD,f[1],B,C,plus,minus);
    fp6_add(A,f[0],f[0]); fp6_add(D,f[1],f[1]);
    fp6_mul_v(hD,hD); fp6_mul_v(hD,hD);
    fp6_mul_v(hA,hA);
    fp6_add(f[0],A,hD); fp6_add(f[1],D,hA);
}

/* Ordinary 10M2 normalized kernel, for a stronger control than the
 * repository's 15M2 general sparse kernel. Computes the EXACT product. */
/* Takes its two fp6 halves separately rather than an fp12_t.
 *
 * Not cosmetic. With an fp12_t parameter, GCC narrows what it believes f[1] to
 * be as soon as the body indexes f[1][1], then reports every later whole-fp6
 * access as a -Wstringop-overflow: "accessing 288 bytes in a region of size
 * 96". The object really is 288 bytes and the code was correct, but the repo
 * builds with -Werror on gcc, so this file did not compile there at all and
 * the experiment could not be reproduced off macOS.
 *
 * The library hit exactly this in fp12_mul_sparse035 and took the same way
 * out; see the note in PROGRESS.md. No runtime cost. */
static void normalized10(fp6_t f0, fp6_t f1, const fp2_t B, const fp2_t C)
{
    fp6_t t1,s,u;
    fp2_t b1,c2,x,y,z,q0,q1,q2,p;
    fp2_mul(b1,f1[1],B); fp2_mul(c2,f1[2],C);
    fp2_add(x,f1[1],f1[2]); fp2_add(y,B,C);
    fp2_mul(x,x,y); fp2_sub(x,x,b1); fp2_sub(x,x,c2);
    fp2_mul_xi(t1[0],x);
    fp2_mul(x,f1[0],B); fp2_mul_xi(y,c2); fp2_add(t1[1],x,y);
    fp2_mul(x,f1[0],C); fp2_add(t1[2],x,b1);
    fp6_add(s,f0,f1);
    fp2_mul(b1,s[1],B); fp2_mul(c2,s[2],C);
    fp2_add(x,s[1],s[2]); fp2_add(y,B,C);
    fp2_mul(x,x,y); fp2_sub(x,x,b1); fp2_sub(x,x,c2);
    fp2_mul_xi(x,x); fp2_add(q0,s[0],x);
    fp2_set_one(p);
    fp2_add(x,s[0],s[1]); fp2_add(y,p,B); fp2_mul(x,x,y);
    fp2_sub(x,x,s[0]); fp2_sub(x,x,b1); fp2_mul_xi(z,c2); fp2_add(q1,x,z);
    fp2_add(x,s[0],s[2]); fp2_add(y,p,C); fp2_mul(x,x,y);
    fp2_sub(x,x,s[0]); fp2_sub(x,x,c2); fp2_add(q2,x,b1);
    fp2_copy(s[0],q0); fp2_copy(s[1],q1); fp2_copy(s[2],q2);
    fp6_sub(s,s,f0); fp6_sub(s,s,t1);
    fp6_mul_v(u,t1); fp6_add(f0,f0,u); fp6_copy(f1,s);
}

static void normalize(norm_table *nt, const ep2_prec_t *pc)
{
    fp2_t prefix[ELIPS_MILLER_LINES], inv, ai;
    CHECK(!fp2_is_zero(pc->l[0].a));
    fp2_copy(prefix[0],pc->l[0].a);
    for (int i=1;i<ELIPS_MILLER_LINES;i++) {
        CHECK(!fp2_is_zero(pc->l[i].a));
        fp2_mul(prefix[i],prefix[i-1],pc->l[i].a);
    }
    fp2_inv(inv,prefix[ELIPS_MILLER_LINES-1]);
    for (int i=ELIPS_MILLER_LINES-1;i>=0;i--) {
        if (i) fp2_mul(ai,inv,prefix[i-1]); else fp2_copy(ai,inv);
        fp2_mul(nt->l[i].B,pc->l[i].c,ai);
        fp2_mul(nt->l[i].C,pc->l[i].b,ai);
        if (i) fp2_mul(inv,inv,pc->l[i].a);
    }
}

static void norm_apply(fp12_t f, const norm_line *l, const fp_t yi,
                       const fp_t xy, int mode)
{
    fp2_t B,C;
    fp2_mul_fp(B,l->B,yi); fp2_mul_fp(C,l->C,xy);
    if (mode==8) normalized8(f,B,C); else normalized10(f[0],f[1],B,C);
}
static void norm_miller(fp12_t f, const norm_table *pc,
                        const fp_t *px, const fp_t *py, int n, int mode)
{
    fp_t yi[8],xy[8];
    CHECK(n<=8);
    for (int j=0;j<n;j++) { fp_inv(yi[j],py[j]); fp_mul(xy[j],px[j],yi[j]); }
    int k=0;
    fp12_set_one(f);
    for (int i=ELIPS_LOOP_TOP-1;i>=0;i--) {
        fp12_sqr(f,f);
        for (int j=0;j<n;j++) {
            norm_apply(f,&pc[j].l[k],yi[j],xy[j],mode);
            if (ELIPS_LOOP[i]) norm_apply(f,&pc[j].l[k+1],yi[j],xy[j],mode);
        }
        k+=1+(ELIPS_LOOP[i]!=0);
    }
    CHECK(k==ELIPS_MILLER_LINES);
}

static void tangent_checks(void)
{
    ep2_t Q,T,U;
    ep2_generator(&Q);
    for (int i=0;i<32;i++) {
        limb_t k[1]={(rng()%100000)+1};
        ep2_mul(&T,&Q,k,64);
        fp2_t z,xx,yy,zz,yz,b3,t;
        r2(z);
        fp2_mul(T.x,T.x,z); fp2_mul(T.y,T.y,z); fp2_mul(T.z,T.z,z);
        ep2_copy(&U,&T);
        ep2_line_t old,neu;
        fp2_sqr(xx,T.x); fp2_sqr(yy,T.y); fp2_sqr(zz,T.z);
        fp2_mul(yz,T.y,T.z);
        ep2_curve_b(b3); fp2_add(t,b3,b3); fp2_add(b3,b3,t);
        fp2_mul(t,b3,zz); fp2_sub(neu.c,yy,t);
        fp2_add(t,xx,xx); fp2_add(t,t,xx); fp2_neg(neu.b,t);
        fp2_add(t,yz,yz); fp2_mul_xi(neu.a,t);
        dbl_line(&old,&U);
        fp2_mul(t,neu.a,T.z); CHECK(fp2_eq(t,old.a));
        fp2_mul(t,neu.b,T.z); CHECK(fp2_eq(t,old.b));
        fp2_mul(t,neu.c,T.z); CHECK(fp2_eq(t,old.c));
    }
    puts("PASS 32 projectively rescaled tangent-line checks");
}

static double clockus(void)
{
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec*1e6+ts.tv_nsec/1e3;
}
static int cmp(const void *a,const void *b)
{ double x=*(const double*)a,y=*(const double*)b; return (x>y)-(x<y); }
static volatile limb_t sink;
static void bench_paths(const ep2_prec_t *pc,const norm_table *nt,
                         const fp_t *px,const fp_t *py,int n,int full)
{
    enum {REPS=15,INNER=30};
    double times[3][REPS];
    fp12_t f,out;
    for (int rep=-1;rep<REPS;rep++) for (int j=0;j<3;j++) {
        /* Rotate the order to avoid always penalizing the same method. */
        int mode=(j+rep+REPS)%3;
        double start=clockus();
        for (int k=0;k<INNER;k++) {
            if (!mode) miller_prec_chunk(f,pc,px,py,(size_t)n);
            else norm_miller(f,nt,px,py,n,mode==1?10:8);
            if (full) { pairing_final_exp_fast(out,f); sink^=out[0][0][0][0]; }
            else sink^=f[0][0][0][0];
        }
        if (rep>=0) times[mode][rep]=(clockus()-start)/INNER;
    }
    const char *names[]={"repository15","normalized10","scaled8"};
    for (int mode=0;mode<3;mode++) {
        qsort(times[mode],REPS,sizeof(double),cmp);
        printf("BENCH n=%d %s %s median=%.3f us min=%.3f max=%.3f\n",
               n,full?"miller+final":"miller",names[mode],times[mode][REPS/2],
               times[mode][0],times[mode][REPS-1]);
    }
}

int main(int argc,char **argv)
{
    (void)argv;
    fp12_t f,a,b,L,t;
    fp2_t B,C,one;
    fp2_set_one(one);
    for (int i=0;i<2000;i++) {
        r12(f); r2(B); r2(C);
        if (i%11==0) fp12_set_one(f);
        if (i%13==0) fp12_set_zero(f);
        if (i%17==0) fp2_copy(C,B);
        if (i%19==0) fp2_neg(C,B);
        if (i%5==0) fp2_set_zero(B);
        if (i%7==0) fp2_set_zero(C);
        fp12_set_one(L); fp2_copy(L[1][1],B); fp2_copy(L[1][2],C);
        fp12_mul(t,f,L);
        fp12_copy(a,f); normalized10(a[0],a[1],B,C); CHECK(fp12_eq(a,t));
        fp12_copy(a,f); normalized8(a,B,C);
        fp12_add(t,t,t); CHECK(fp12_eq(a,t));
    }
    puts("PASS 2000 kernels against dense multiplication, including zero coefficients");
    tangent_checks();
    ep_t G,P; ep2_t Q,R;
    ep_generator(&G); ep2_generator(&Q);
    ep2_prec_t pcs[8]; norm_table nts[8];
    fp_t px[8],py[8]; fp2_t qx[8],qy[8];
    for (int trial=0;trial<4;trial++) {
        for (int j=0;j<8;j++) {
            limb_t kp[1]={(rng()%100000)+1},kq[1]={(rng()%100000)+1};
            ep_mul(&P,&G,kp,64); ep2_mul(&R,&Q,kq,64);
            CHECK(ep_to_affine(px[j],py[j],&P));
            CHECK(ep2_to_affine(qx[j],qy[j],&R));
            CHECK(ep2_precompute(&pcs[j],&R)); normalize(&nts[j],&pcs[j]);
        }
        for (int n=0;n<=8;n++) {
            pairing_miller_multi(f,qx,qy,px,py,(size_t)n);
            pairing_final_exp_fast(a,f);
            norm_miller(f,nts,px,py,n,8);
            pairing_final_exp_fast(b,f); CHECK(fp12_eq(a,b));
            norm_miller(f,nts,px,py,n,10);
            pairing_final_exp_fast(b,f); CHECK(fp12_eq(a,b));
        }
    }
    puts("PASS 36 full multi-pairing comparisons, n=0..8, two candidate kernels");
    /* Check the exact exponent as well, not just the library's cubed pairing. */
    miller_prec_chunk(f,pcs,px,py,2); pairing_final_exp_plain(a,f);
    norm_miller(f,nts,px,py,2,8); pairing_final_exp_plain(b,f);
    CHECK(fp12_eq(a,b)); puts("PASS exact final-exponent comparison");
    printf("TABLE bytes old=%zu normalized=%zu\n",sizeof(ep2_prec_t),sizeof(norm_table));
    if (argc>1) {
        bench_paths(pcs,nts,px,py,1,0); bench_paths(pcs,nts,px,py,1,1);
        bench_paths(pcs,nts,px,py,2,0); bench_paths(pcs,nts,px,py,2,1);
        bench_paths(pcs,nts,px,py,8,0); bench_paths(pcs,nts,px,py,8,1);
    }
    return 0;
}
