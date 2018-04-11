//
//  bn_fp12.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_fp12_h
#define bn_fp12_h

#include <ELiPS_bn_bls/bn_fp6.h>

extern void Fp12_init(Fp12 *A);
extern void Fp12_clear(Fp12 *A);
extern void Fp12_printf(Fp12 *A,char *str);
extern void Fp12_set(Fp12 *ANS,Fp12 *A);
extern void Fp12_set_ui(Fp12 *ANS,unsigned long int UI);
extern void Fp12_set_mpz(Fp12 *ANS,mpz_t A);
extern void Fp12_set_neg(Fp12 *ANS,Fp12 *A);
extern void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state);
extern void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);
extern void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
extern void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
extern void Fp12_sqr(Fp12 *ANS,Fp12 *A);
extern void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B);
extern void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
extern void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
extern void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B);
extern void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI);
extern void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B);
extern void Fp12_inv(Fp12 *ANS,Fp12 *A);
extern int  Fp12_legendre(Fp12 *A);
extern int  Fp12_isCNR(Fp12 *A);
extern void Fp12_sqrt(Fp12 *ANS,Fp12 *A);
extern void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar);
extern int  Fp12_cmp(Fp12 *A,Fp12 *B);
extern int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI);
extern int  Fp12_cmp_mpz(Fp12 *A,mpz_t B);
extern int  Fp12_cmp_zero(Fp12 *A);
extern int  Fp12_cmp_one(Fp12 *A);

#endif /* bn_fp12_h */
