//
//  bls12_line_tate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_line_tate_h
#define bls12_line_tate_h

#include <ELiPS_bn_bls/bn_efp12.h>

extern void BLS12_EFp_ECD_return_lambda(EFp *ANS,Fp *lambda,EFp *P);
extern void BLS12_EFp_ECA_return_lambda(EFp *ANS,Fp *lambda,EFp *P1, EFp *P2);
extern void BLS12_ff_ltt_vtt_for_tate(Fp12 *f,EFp *T,EFp12 *Q);
extern void BLS12_f_ltp_vtp_for_tate(Fp12 *f,EFp *T,EFp *P,EFp12 *Q);

#endif /* bls12_line_tate_h */
