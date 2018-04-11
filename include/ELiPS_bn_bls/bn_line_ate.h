//
//  bn_line_ate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_line_ate_h
#define bn_line_ate_h

#include <ELiPS_bn_bls/bn_p8sparse.h>

extern void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
extern void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);

#endif /* bn_line_ate_h */
