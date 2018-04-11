//
//  bn_efp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_efp_h
#define bn_efp_h

#include <ELiPS_bn_bls/curve_dtypes.h>

extern void EFp_init(EFp *P);
extern void EFp_clear(EFp *P);
extern void EFp_printf(EFp *P,char *str);
extern void EFp_set(EFp *ANS,EFp *P);
extern void EFp_set_ui(EFp *ANS,unsigned long int UI);
extern void EFp_set_mpz(EFp *ANS,mpz_t A);
extern void EFp_set_neg(EFp *ANS,EFp *P);
extern void EFp_rational_point_BN(EFp *P);
extern void EFp_rational_point_BLS12(EFp *P);
extern void EFp_ECD(EFp *ANS,EFp *P);
extern void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
extern void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);

#endif /* bn_efp_h */
