//
//  fp4.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/22/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef fp4_h
#define fp4_h

#include <ELiPS_bn_bls/bn_fp2.h>

extern void Fp4_init(Fp4 *A);
extern void Fp4_set(Fp4 *ANS, Fp4 *A);
extern void Fp4_set_ui(Fp4 *A,signed long int B);
extern void Fp4_set_random(Fp4 *A, gmp_randstate_t state);
extern void Fp4_clear(Fp4 *A);
extern void Fp4_printf(Fp4 *A);
extern void Fp4_add(Fp4 *ANS, Fp4 *A, Fp4 *B);
extern void Fp4_add_ui (Fp4 *ANS, Fp4 *A,unsigned long int B);
extern void Fp4_sub(Fp4 *ANS, Fp4 *A, Fp4 *B);
extern void Fp4_sqr(Fp4 *ANS, Fp4 *A);
extern void Fp4_sqr_complex(Fp4 *ANS, Fp4 *A);
extern void Fp4_mul(Fp4 *ANS, Fp4 *A, Fp4 *B);
extern void Fp4_mul_basis(Fp4 *ANS, Fp4 *A);
extern void Fp4_mul_ui(Fp4 *ANS, Fp4 *A,unsigned long int B);
extern void Fp4_mul_mpz(Fp4 *ANS, Fp4 *A,mpz_t B);
extern void Fp4_mul_Fp (Fp4 *ANS, Fp4 *A, Fp *B);
extern void Fp4_neg(Fp4 *ANS, Fp4 *A);
extern void Fp4_invert(Fp4 *ANS, Fp4 *A);
extern void Fp4_div(Fp4 *ANS, Fp4 *A, Fp4 *B);
extern void Fp4_pow(Fp4 *ANS, Fp4 *A,mpz_t B);
extern void Fp4_sqrt(Fp4 *ANS, Fp4 *A);//x^2=a mod p
extern int  Fp4_legendre(Fp4 *A);
extern int  Fp4_cmp(Fp4 *A, Fp4 *B);
extern int  Fp4_cmp_mpz(Fp4 *A,mpz_t B);
extern void Fp4_neg(Fp4 *ANS, Fp4 *A);
extern void Fp4_frobenius_map(Fp4 *ANS, Fp4 *A);
extern void Fp4_mul_betainv(Fp4 *ANS);
#endif /* fp4_h */
