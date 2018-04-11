//
//  bls12_p8sparse.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_p8sparse_h
#define bls12_p8sparse_h

#include <ELiPS_bn_bls/bn_efp2.h>

extern void BLS12_Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);
extern void BLS12_Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);

#endif /* bls12_p8sparse_h */
