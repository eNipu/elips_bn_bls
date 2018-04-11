//
//  bn_miller_ate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_miller_ate_h
#define bn_miller_ate_h

#include <ELiPS_bn_bls/bn_twist.h>
#include <ELiPS_bn_bls/bn_p8sparse.h>
#include <ELiPS_bn_bls/bn_line_ate.h>

//miller algo for plain ate
extern void Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);

#endif /* bn_miller_ate_h */
