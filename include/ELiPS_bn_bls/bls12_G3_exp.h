//
//  bls12_G3_exp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_G3_exp_h
#define bls12_G3_exp_h

#include <ELiPS_bn_bls/bn_utils.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bls12_frobenius.h>

//G3
extern void BLS12_plain_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar);
extern void BLS12_2split_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar);
extern void BLS12_4split_G3_exp( Fp12 *ANS,Fp12 *A,mpz_t scalar);

#endif /* bls12_G3_exp_h */
