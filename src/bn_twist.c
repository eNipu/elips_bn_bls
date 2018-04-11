//
//  bn_twist.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_twist.h>

void EFp12_to_EFp2(EFp2 *ANS,EFp12 *P){
    Fp2_set_ui(&ANS->x,0);
    Fp2_set(&ANS->x,&P->x.x0.x2);
    Fp2_mul_basis(&ANS->x,&ANS->x);
    Fp2_set_ui(&ANS->y,0);
    Fp2_set(&ANS->y,&P->y.x1.x1);
    Fp2_mul_basis(&ANS->y,&ANS->y);
    ANS->infinity=P->infinity;
}

void EFp2_to_EFp12(EFp12 *ANS,EFp2 *P){
    Fp12_set_ui(&ANS->x,0);
    Fp2_set(&ANS->x.x0.x2,&P->x);
    Fp2_inv_basis(&ANS->x.x0.x2,&ANS->x.x0.x2);
    Fp12_set_ui(&ANS->y,0);
    Fp2_set(&ANS->y.x1.x1,&P->y);
    Fp2_inv_basis(&ANS->y.x1.x1,&ANS->y.x1.x1);
    ANS->infinity=P->infinity;
}
