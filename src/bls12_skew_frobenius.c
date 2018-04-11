//
//  bls12_skew_frobenius.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_skew_frobenius.h>

void BLS12_EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A){
    Fp_mul_mpz(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
}

void BLS12_EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&d12_skew_frobenius_constant[f_p1][0]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&d12_skew_frobenius_constant[f_p1][1]);
}

void BLS12_EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&d12_skew_frobenius_constant[f_p2][0]);
    //y
    Fp2_mul(&ANS->y,&A->y,&d12_skew_frobenius_constant[f_p2][1]);
}

void BLS12_EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&d12_skew_frobenius_constant[f_p3][0]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&d12_skew_frobenius_constant[f_p3][1]);
}

void BLS12_EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&d12_skew_frobenius_constant[f_p10][0]);
    //y
    Fp2_mul(&ANS->y,&A->y,&d12_skew_frobenius_constant[f_p10][1]);
}
