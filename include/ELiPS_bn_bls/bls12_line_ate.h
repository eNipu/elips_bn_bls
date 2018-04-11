//
//  bls12_line_ate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_line_ate_h
#define bls12_line_ate_h

#include <ELiPS_bn_bls/bls12_p8sparse.h>

extern void BLS12_ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
extern void BLS12_f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);

#endif /* bls12_line_ate_h */
