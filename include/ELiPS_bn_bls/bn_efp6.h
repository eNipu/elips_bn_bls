//
//  bn_efp6.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_efp6_h
#define bn_efp6_h

#include <ELiPS_bn_bls/bn_efp2.h>


extern void EFp6_init(EFp6 *P);
extern void EFp6_clear(EFp6 *P);
extern void EFp6_printf(EFp6 *P,char *str);
extern void EFp6_set(EFp6 *ANS,EFp6 *P);
extern void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);
extern void EFp6_set_mpz(EFp6 *ANS,mpz_t A);
extern void EFp6_set_neg(EFp6 *ANS,EFp6 *P);
extern void EFp6_rational_point(EFp6 *P);
extern void EFp6_ECD(EFp6 *ANS,EFp6 *P);
extern void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);
extern void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar);

#endif /* bn_efp6_h */
