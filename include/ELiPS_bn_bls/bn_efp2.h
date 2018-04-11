//
//  bn_efp2.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_efp2_h
#define bn_efp2_h

#include <ELiPS_bn_bls/bn_efp.h>

extern void EFp2_init(EFp2 *P);
extern void EFp2_clear(EFp2 *P);
extern void EFp2_printf(EFp2 *P,char *str);
extern void EFp2_set(EFp2 *ANS,EFp2 *P);
extern void EFp2_set_ui(EFp2 *ANS,unsigned long int UI);
extern void EFp2_set_mpz(EFp2 *ANS,mpz_t A);
extern void EFp2_set_neg(EFp2 *ANS,EFp2 *P);
extern void EFp2_rational_point(EFp2 *P);
extern void EFp2_ECD(EFp2 *ANS,EFp2 *P);
extern void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);
extern void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);

#endif /* bn_efp2_h */
