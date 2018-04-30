//
//  bn_fp2.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/29/18.

/*
 * ELiPS is an Efficient Library for Pairing-based Systems
 * Copyright (C) 2008-2018 ELiPS Authors. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * This file is part of ELiPS Project. ELiPS is legal property of its
 * developers, whose names are not listed here. 
 *
 * ELiPS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * ELiPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ELiPS. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef bn_fp2_h
#define bn_fp2_h


#include <ELiPS_bn_bls/bn_fp.h>
#include <ELiPS_bn_bls/field_dtype.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>

//Fp2
extern void Fp2_init(Fp2 *A);
extern void Fp2_clear(Fp2 *A);
extern void Fp2_printf(Fp2 *A,char *str);
extern void Fp2_set(Fp2 *ANS,Fp2 *A);
extern void Fp2_set_ui(Fp2 *ANS,unsigned long int UI);
extern void Fp2_set_mpz(Fp2 *ANS,mpz_t A);
extern void Fp2_set_neg(Fp2 *ANS,Fp2 *A);
extern void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state);
extern void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_mul_Fp(Fp2 *ANS, Fp2 *A, Fp *B);
extern void Fp2_mul_basis(Fp2 *ANS,Fp2 *A);
extern void Fp2_inv_basis(Fp2 *ANS,Fp2 *A);
extern void Fp2_sqr(Fp2 *ANS,Fp2 *A);
extern void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_inv(Fp2 *ANS,Fp2 *A);
extern int  Fp2_legendre(Fp2 *A);
extern int  Fp2_isCNR(Fp2 *A);
extern void Fp2_sqrt(Fp2 *ANS,Fp2 *A);
extern void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar);
extern int  Fp2_cmp(Fp2 *A,Fp2 *B);
extern int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI);
extern int  Fp2_cmp_mpz(Fp2 *A,mpz_t B);
extern int  Fp2_cmp_zero(Fp2 *A);
extern int  Fp2_cmp_one(Fp2 *A);
extern void Fp2_neg(Fp2 *ans,Fp2 *a);

extern void Fp2_mul_KSS16(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_mul_basis_KSS16(Fp2 *ANS, Fp2 *A);
extern void Fp2_invert_kss16(Fp2 *ANS,Fp2 *A);


// #pragma mark Fp2 methods
//void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_sqr_kss16(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_sqr_complex(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
//void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
//void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B);
//void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);
extern void Fp2_inv(Fp2 *ANS, Fp2 *A);
//void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
//void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
//void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);//x^2=a mod p
//int  Fp2_legendre(struct Fp2 *A);
//int  Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
//int  Fp2_cmp_mpz(struct Fp2 *A,mpz_t B);
//void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);

#endif /* bn_fp2_h */
