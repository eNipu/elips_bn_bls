//
//  bn_efp12.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_efp12_h
#define bn_efp12_h

#include <ELiPS_bn_bls/bn_efp6.h>

extern void EFp12_init(EFp12 *P);
extern void EFp12_clear(EFp12 *P);
extern void EFp12_printf(EFp12 *P,char *str);
extern void EFp12_set(EFp12 *ANS,EFp12 *P);
extern void EFp12_set_ui(EFp12 *ANS,unsigned long int UI);
extern void EFp12_set_mpz(EFp12 *ANS,mpz_t A);
extern void EFp12_set_neg(EFp12 *ANS,EFp12 *P);
extern void EFp12_rational_point_BN(EFp12 *P);
extern void EFp12_rational_point_BLS12(EFp12 *P);
extern void EFp12_ECD(EFp12 *ANS,EFp12 *P);
extern void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2);
extern void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar);

#endif /* bn_efp12_h */
