//
//  bn_line_tate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_line_tate_h
#define bn_line_tate_h

#include <ELiPS_bn_bls/bn_efp12.h>


extern void EFp_ECD_return_lambda(EFp *ANS,Fp *lambda,EFp *P);
extern void EFp_ECA_return_lambda(EFp *ANS,Fp *lambda,EFp *P1, EFp *P2);
extern void ff_ltt_vtt_for_tate(Fp12 *f,EFp *T,EFp12 *Q);
extern void f_ltp_vtp_for_tate(Fp12 *f,EFp *T,EFp *P,EFp12 *Q);

#endif /* bn_line_tate_h */
