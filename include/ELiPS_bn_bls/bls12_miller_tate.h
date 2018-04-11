//
//  bls12_miller_tate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_miller_tate_h
#define bls12_miller_tate_h

#include <ELiPS_bn_bls/bls12_line_tate.h>
#include <ELiPS_bn_bls/bls12_twist.h>

//miller algo for tate
extern void BLS12_Miller_algo_for_tate(Fp12 *ANS,EFp12 *Q,EFp12 *P);
#endif /* bls12_miller_tate_h */
