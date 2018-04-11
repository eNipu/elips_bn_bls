//
//  bn_generate_points.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_generate_points.h>

EFp12 G1_P;
EFp12 G2_Q;
EFp12 Random_R;

void BN12_generate_G1_point(EFp12 *P){
    EFp12_init(P);
    EFp Tmp;
    EFp_init(&Tmp);
    
    EFp_rational_point_BN(&Tmp);
    EFp12_set_ui(P,0);
    Fp_set(&P->x.x0.x0.x0,&Tmp.x);
    Fp_set(&P->y.x0.x0.x0,&Tmp.y);
    P->infinity=Tmp.infinity;
    
    EFp_clear(&Tmp);
}

//TODO
void BN12_generate_G2_point(EFp12 *Q){
    EFp12_init(Q);
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point_BN(&random_P);
    mpz_pow_ui(exp,curve_parameters.order,2);
    mpz_tdiv_q(exp,curve_parameters.EFpd_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void BN12_generate_random_point(EFp12 *R){
    EFp12_init(R);
    EFp12_rational_point_BN(R);
}
