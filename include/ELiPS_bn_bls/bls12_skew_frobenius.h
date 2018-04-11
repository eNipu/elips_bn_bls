//
//  bls12_skew_frobenius.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_skew_frobenius_h
#define bls12_skew_frobenius_h

#include <ELiPS_bn_bls/bn_bls12_precoms.h>
#include <ELiPS_bn_bls/bn_efp12.h>


//skew frobenius
extern void BLS12_EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);
extern void BLS12_EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);
extern void BLS12_EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);
extern void BLS12_EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);
extern void BLS12_EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);


#endif /* bls12_skew_frobenius_h */
