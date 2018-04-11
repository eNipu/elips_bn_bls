//
//  bn_pairings.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_pairings_h
#define bn_pairings_h

#include <ELiPS_bn_bls/bn_miller_ate.h>
#include <ELiPS_bn_bls/bn_miller_tate.h>
#include <ELiPS_bn_bls/bn_miller_xate.h>
#include <ELiPS_bn_bls/bn_miller_optate.h>
#include <ELiPS_bn_bls/bn_final_exp.h>
#include <ELiPS_bn_bls/bn_utils.h>

extern void BN12_tate(Fp12 *ANS,EFp12 *P, EFp12 *Q);
extern void BN12_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);
extern void BN12_opt_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q);
extern void BN12_x_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q);

#endif /* bn_pairings_h */
