//
//  bls12_twist.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_twist_h
#define bls12_twist_h

#include <ELiPS_bn_bls/bn_efp12.h>
//sextic twist
extern void BLS12_EFp12_to_EFp2(EFp2 *ANS,EFp12 *P);
extern void BLS12_EFp2_to_EFp12(EFp12 *ANS,EFp2 *P);
//primary twist
extern void BLS12_EFp12_to_EFp(EFp *ANS,EFp12 *P);
extern void BLS12_EFp_to_EFp12(EFp12 *ANS,EFp *P);

#endif /* bls12_twist_h */
