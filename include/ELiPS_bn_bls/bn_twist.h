//
//  bn_twist.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_twist_h
#define bn_twist_h

#include <ELiPS_bn_bls/bn_efp12.h>

extern void EFp12_to_EFp2(EFp2 *ANS,EFp12 *P);
extern void EFp2_to_EFp12(EFp12 *ANS,EFp2 *P);

#endif /* bn_twist_h */
