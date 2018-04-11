//
//  bn_p8sparse.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_p8sparse_h
#define bn_p8sparse_h

#include <ELiPS_bn_bls/bn_efp2.h>

extern void Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);
extern void Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);

#endif /* bn_p8sparse_h */
