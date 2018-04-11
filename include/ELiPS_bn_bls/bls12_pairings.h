//
//  bls12_pairings.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_pairings_h
#define bls12_pairings_h

#include <ELiPS_bn_bls/bls12_miller_ate.h>
#include <ELiPS_bn_bls/bls12_miller_tate.h>
#include <ELiPS_bn_bls/bls12_miller_optate.h>
#include <ELiPS_bn_bls/bls12_finalexp.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bn_utils.h>

extern void BLS12_tate(Fp12 *ANS,EFp12 *P, EFp12 *Q);
extern void BLS12_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BLS12_opt_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q);

#endif /* bls12_pairings_h */
