//
//  bn_miller_xate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_miller_xate_h
#define bn_miller_xate_h

#include <ELiPS_bn_bls/bn_line_ate.h>
#include <ELiPS_bn_bls/bn_frobenius.h>
#include <ELiPS_bn_bls/bn_twist.h>
#include <ELiPS_bn_bls/bn_skew_frobenius.h>

//miller algo for x-ate
extern void Miller_algo_for_x_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);

#endif /* bn_miller_xate_h */
