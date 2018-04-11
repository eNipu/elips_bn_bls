//
//  bn_generate_points.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_generate_points_h
#define bn_generate_points_h

#include <ELiPS_bn_bls/bn_frobenius.h>


extern EFp12 G1_P;
extern EFp12 G2_Q;
extern EFp12 Random_R;

//rational points
extern void BN12_generate_G1_point(EFp12 *P);
extern void BN12_generate_G2_point(EFp12 *Q);
extern void BN12_generate_random_point(EFp12 *R);

#endif /* bn_generate_points_h */
