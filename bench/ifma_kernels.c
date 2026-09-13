/*
 * IFMA kernels for bench/ifma_probe.c. Separate translation unit on purpose.
 *
 * These are measured, not shipped. They live apart from the driver because
 * that is the only thing that reliably stopped gcc from measuring nothing:
 * inside one file it infers that a pure function of unchanged inputs need not
 * be called twice, and hoists the whole thing out of the timing loop. That
 * produced 6 cycles for a Montgomery multiply, and a 37-cycle reading for a
 * loop that really takes 181. noinline did not stop it; a memory clobber did
 * not stop it; a call across a translation unit boundary does.
 */
#include "ifma_kernels.h"
#include <string.h>
#include <immintrin.h>
#define M52 ((uint64_t)0xFFFFFFFFFFFFFull)
static const uint64_t P52[8]={0x0effffffffaaabull,0x0feb153ffffb9full,0x06b0f6241eabffull,0x012bf6730d2a0full,0x0764774b84f385ull,0x01ba7b6434bacdull,0x01ea397fe69a4bull,0x0000000001a011ull};
static const uint64_t K0=0x3fffcfffcfffdull;
static void to52(uint64_t o[8],const uint64_t in[6],unsigned s){uint64_t t[7];
 if(s){t[0]=in[0]<<s;for(int i=1;i<6;i++)t[i]=(in[i]<<s)|(in[i-1]>>(64-s));t[6]=in[5]>>(64-s);}else{memcpy(t,in,48);t[6]=0;}
 for(int j=0;j<8;j++){unsigned b=52u*(unsigned)j,w=b>>6,q=b&63;uint64_t v=t[w]>>q;if(q)v|=t[w+1]<<(64-q);o[j]=v&M52;}}
static void from52(uint64_t o[6],const uint64_t in[8]){uint64_t t[7]={0};
 for(int j=0;j<8;j++){unsigned b=52u*(unsigned)j,w=b>>6,q=b&63;t[w]|=in[j]<<q;if(q)t[w+1]|=in[j]>>(64-q);}memcpy(o,t,48);}
void part_to52(uint64_t o[8],const uint64_t in[6],unsigned s){to52(o,in,s);}
uint64_t part_convert(const uint64_t in[6]){uint64_t ca[8],cb[8],back[6];
 to52(ca,in,32); to52(cb,in,0); from52(back,ca); return back[0]+ca[7]+cb[7];}
__attribute__((target("avx512f,avx512vl,avx512ifma")))
uint64_t part_loop(const uint64_t a52[8],const uint64_t b52[8]){
 const __m512i A=_mm512_loadu_si512(a52),P=_mm512_loadu_si512(P52);
 __m512i acc=_mm512_setzero_si512();
 for(int i=0;i<8;i++){const __m512i bi=_mm512_set1_epi64((long long)b52[i]);
  acc=_mm512_madd52lo_epu64(acc,A,bi);
  uint64_t l0=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc));
  uint64_t m=(l0*K0)&M52; const __m512i mv=_mm512_set1_epi64((long long)m);
  acc=_mm512_madd52lo_epu64(acc,P,mv);
  uint64_t c=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc))>>52;
  acc=_mm512_alignr_epi64(_mm512_setzero_si512(),acc,1);
  acc=_mm512_mask_add_epi64(acc,0x01,acc,_mm512_set1_epi64((long long)c));
  acc=_mm512_madd52hi_epu64(acc,A,bi); acc=_mm512_madd52hi_epu64(acc,P,mv);}
 uint64_t d[8]; _mm512_storeu_si512(d,acc); return d[0]+d[7];}
__attribute__((target("avx512f,avx512vl,avx512ifma")))
uint64_t part_loop_nodep(const uint64_t a52[8],const uint64_t b52[8]){
 const __m512i A=_mm512_loadu_si512(a52),P=_mm512_loadu_si512(P52);
 __m512i acc=_mm512_setzero_si512(); const __m512i mv=_mm512_set1_epi64(0x123456789abull);
 for(int i=0;i<8;i++){const __m512i bi=_mm512_set1_epi64((long long)b52[i]);
  acc=_mm512_madd52lo_epu64(acc,A,bi); acc=_mm512_madd52lo_epu64(acc,P,mv);
  acc=_mm512_alignr_epi64(_mm512_setzero_si512(),acc,1);
  acc=_mm512_madd52hi_epu64(acc,A,bi); acc=_mm512_madd52hi_epu64(acc,P,mv);}
 uint64_t d[8]; _mm512_storeu_si512(d,acc); return d[0]+d[7];}

/* K independent Montgomery multiplications interleaved. Same work per product
 * as part_loop; the only difference is that K accumulator chains are in flight,
 * so one product's latency is another's throughput. This is the measurement
 * that decides whether batching is the route. */
__attribute__((target("avx512f,avx512vl,avx512ifma")))
uint64_t part_loop_x2(const uint64_t a52[8],const uint64_t b52[8]){
 const __m512i A=_mm512_loadu_si512(a52),P=_mm512_loadu_si512(P52);
 __m512i acc[2]; for(int k=0;k<2;k++) acc[k]=_mm512_setzero_si512();
 for(int i=0;i<8;i++){
  const __m512i bi=_mm512_set1_epi64((long long)b52[i]);
  __m512i mv[2];
  for(int k=0;k<2;k++) acc[k]=_mm512_madd52lo_epu64(acc[k],A,bi);
  for(int k=0;k<2;k++){
   uint64_t l0=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k]));
   mv[k]=_mm512_set1_epi64((long long)((l0*K0)&M52)); }
  for(int k=0;k<2;k++) acc[k]=_mm512_madd52lo_epu64(acc[k],P,mv[k]);
  for(int k=0;k<2;k++){
   uint64_t c=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k]))>>52;
   acc[k]=_mm512_alignr_epi64(_mm512_setzero_si512(),acc[k],1);
   acc[k]=_mm512_mask_add_epi64(acc[k],0x01,acc[k],_mm512_set1_epi64((long long)c)); }
  for(int k=0;k<2;k++) acc[k]=_mm512_madd52hi_epu64(acc[k],A,bi);
  for(int k=0;k<2;k++) acc[k]=_mm512_madd52hi_epu64(acc[k],P,mv[k]); }
 uint64_t r=0,d[8];
 for(int k=0;k<2;k++){ _mm512_storeu_si512(d,acc[k]); r+=d[0]+d[7]; }
 return r;}

/* K independent Montgomery multiplications interleaved. Same work per product
 * as part_loop; the only difference is that K accumulator chains are in flight,
 * so one product's latency is another's throughput. This is the measurement
 * that decides whether batching is the route. */
__attribute__((target("avx512f,avx512vl,avx512ifma")))
uint64_t part_loop_x4(const uint64_t a52[8],const uint64_t b52[8]){
 const __m512i A=_mm512_loadu_si512(a52),P=_mm512_loadu_si512(P52);
 __m512i acc[4]; for(int k=0;k<4;k++) acc[k]=_mm512_setzero_si512();
 for(int i=0;i<8;i++){
  const __m512i bi=_mm512_set1_epi64((long long)b52[i]);
  __m512i mv[4];
  for(int k=0;k<4;k++) acc[k]=_mm512_madd52lo_epu64(acc[k],A,bi);
  for(int k=0;k<4;k++){
   uint64_t l0=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k]));
   mv[k]=_mm512_set1_epi64((long long)((l0*K0)&M52)); }
  for(int k=0;k<4;k++) acc[k]=_mm512_madd52lo_epu64(acc[k],P,mv[k]);
  for(int k=0;k<4;k++){
   uint64_t c=(uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k]))>>52;
   acc[k]=_mm512_alignr_epi64(_mm512_setzero_si512(),acc[k],1);
   acc[k]=_mm512_mask_add_epi64(acc[k],0x01,acc[k],_mm512_set1_epi64((long long)c)); }
  for(int k=0;k<4;k++) acc[k]=_mm512_madd52hi_epu64(acc[k],A,bi);
  for(int k=0;k<4;k++) acc[k]=_mm512_madd52hi_epu64(acc[k],P,mv[k]); }
 uint64_t r=0,d[8];
 for(int k=0;k<4;k++){ _mm512_storeu_si512(d,acc[k]); r+=d[0]+d[7]; }
 return r;}

/* A COMPLETE Montgomery multiply, four in flight, from 52-bit inputs to
 * reduced 52-bit outputs: the loop plus carry normalisation plus the
 * conditional subtraction. No radix conversion -- that is the part step 3
 * removes by keeping the tower in 52-bit form. This is the number that is
 * comparable with fp_mul, which also includes its own final reduction. */
__attribute__((target("avx512f,avx512vl,avx512ifma")))
uint64_t part_full_x4(const uint64_t a52[8], const uint64_t b52[8])
{
    const __m512i A = _mm512_loadu_si512(a52), P = _mm512_loadu_si512(P52);
    __m512i acc[4];
    for (int k = 0; k < 4; k++) acc[k] = _mm512_setzero_si512();

    for (int i = 0; i < 8; i++) {
        const __m512i bi = _mm512_set1_epi64((long long)b52[i]);
        __m512i mv[4];
        for (int k = 0; k < 4; k++) acc[k] = _mm512_madd52lo_epu64(acc[k], A, bi);
        for (int k = 0; k < 4; k++) {
            uint64_t l0 = (uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k]));
            mv[k] = _mm512_set1_epi64((long long)((l0 * K0) & M52));
        }
        for (int k = 0; k < 4; k++) acc[k] = _mm512_madd52lo_epu64(acc[k], P, mv[k]);
        for (int k = 0; k < 4; k++) {
            uint64_t c = (uint64_t)_mm_cvtsi128_si64(_mm512_castsi512_si128(acc[k])) >> 52;
            acc[k] = _mm512_alignr_epi64(_mm512_setzero_si512(), acc[k], 1);
            acc[k] = _mm512_mask_add_epi64(acc[k], 0x01, acc[k], _mm512_set1_epi64((long long)c));
        }
        for (int k = 0; k < 4; k++) acc[k] = _mm512_madd52hi_epu64(acc[k], A, bi);
        for (int k = 0; k < 4; k++) acc[k] = _mm512_madd52hi_epu64(acc[k], P, mv[k]);
    }

    uint64_t out = 0;
    for (int k = 0; k < 4; k++) {
        uint64_t d[8], c = 0;
        _mm512_storeu_si512(d, acc[k]);
        for (int j = 0; j < 8; j++) { uint64_t v = d[j] + c; d[j] = v & M52; c = v >> 52; }
        /* conditional subtraction of p, by mask, in 52-bit digits */
        uint64_t s[8], borrow = 0;
        for (int j = 0; j < 8; j++) {
            uint64_t x = d[j] - P52[j] - borrow;
            borrow = (x >> 63) & 1;
            s[j] = x & M52;
        }
        uint64_t keep = (uint64_t)0 - (borrow ^ 1);
        for (int j = 0; j < 8; j++) d[j] = (s[j] & keep) | (d[j] & ~keep);
        out += d[0] + d[7];
    }
    return out;
}
void part_from52(uint64_t o[6],const uint64_t in[8]){from52(o,in);}
