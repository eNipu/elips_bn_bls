//
//  bn_final_exp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_final_exp_h
#define bn_final_exp_h

#include <ELiPS_bn_bls/bn_frobenius.h>

extern void Final_exp_plain(Fp12 *ANS,Fp12 *A);
extern void Final_exp_optimal(Fp12 *ANS,Fp12 *A);
extern void Fp12_pow_X(Fp12 *ANS,Fp12 *A);

#endif /* bn_final_exp_h */
