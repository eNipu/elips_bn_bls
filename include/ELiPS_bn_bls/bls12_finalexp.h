//
//  bls12_finalexp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_finalexp_h
#define bls12_finalexp_h

#include <ELiPS_bn_bls/bls12_frobenius.h>

//final exp
extern void BLS12_Final_exp_plain(Fp12 *ANS,Fp12 *A);
extern void BLS12_Final_exp_optimal(Fp12 *ANS,Fp12 *A);

#endif /* bls12_finalexp_h */
