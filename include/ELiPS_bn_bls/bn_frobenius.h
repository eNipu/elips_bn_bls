//
//  bn_frobenius.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_frobenius_h
#define bn_frobenius_h

#include <ELiPS_bn_bls/bn_efp12.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>

//frobenius
extern void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A);
extern void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A);

#endif /* bn_frobenius_h */
